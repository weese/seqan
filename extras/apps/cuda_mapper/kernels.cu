// ==========================================================================
//                                cuda_mapper
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================

// ============================================================================
// Prerequisites
// ============================================================================

#include "kernels.h"
#include <thrust/sort.h>

using namespace seqan;

// ============================================================================
// Kernels
// ============================================================================

// --------------------------------------------------------------------------
// Kernel _hashReadsKernel()
// --------------------------------------------------------------------------

template <typename TReadSeqsView, typename THashes, typename TIdx>
__global__ void
_hashReadsKernel(TReadSeqsView readSeqs, THashes hashes, TIdx idxs)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    // Return silently if there is no job left.
    if (idx >= length(readSeqs)) return;

    // Compute the hash of a read.
    hashes[idx] = hashRead(infix(readSeqs[idx], 0, 20));

    // Fill idxs with the identity permutation.
    idxs[idx] = idx;
}

// --------------------------------------------------------------------------
// Kernel _mapReadsKernel()
// --------------------------------------------------------------------------

template <typename TIndexView, typename TReadSeqsView, typename TIdxView>
__global__ void
_mapReadsKernel(TIndexView index, TReadSeqsView readSeqs, TIdxView idxs)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    // Return silently if there is no job left.
    if (idx >= length(readSeqs)) return;

    // Pick the read id to map.
    idx = idxs[idx];

    // Map the read.
    mapRead(index, readSeqs[idx], (unsigned)idx);
}

// ============================================================================
// Functions
// ============================================================================

template <typename TReadSeq>
SEQAN_FUNC __uint32
hashRead(TReadSeq const & readSeq)
{
    __uint32 result = 0;

    for (int i = 0; i < _min(16,length(readSeq)); ++i)
    {
        __uint32 c = ordValue(readSeq[i]);
        result |= c << i*2;
    }

    return result;
}

// --------------------------------------------------------------------------
// Function mapReads()                                                  [GPU]
// --------------------------------------------------------------------------

void _mapReads(TGenomeIndex & index, TReadSeqs & readSeqs, GPU const & /* tag */)
{
    typedef typename Device<TGenomeIndex>::Type         TDeviceIndex;
    typedef typename Device<TReadSeqs>::Type            TDeviceReadSeqs;
    typedef typename View<TDeviceIndex>::Type           TDeviceIndexView;
    typedef typename View<TDeviceReadSeqs>::Type        TDeviceReadSeqsView;

    // Copy index to device.
    TDeviceIndex deviceIndex;
    assign(deviceIndex, index);

    // Copy read seqs to device.
    TDeviceReadSeqs deviceReadSeqs;
    assign(deviceReadSeqs, readSeqs);

    cudaPrintFreeMemory();

    // Instantiate views of device objects.
    TDeviceIndexView deviceIndexView = view(deviceIndex);
    TDeviceReadSeqsView deviceReadSeqsView = view(deviceReadSeqs);

    // Setup kernel parameters.
    unsigned threadsPerBlock = 256;
    unsigned blocksPerGrid = (length(readSeqs) + threadsPerBlock - 1) / threadsPerBlock;

    thrust::device_vector<__uint32> hashes(length(readSeqs));
    thrust::device_vector<__uint32> idx(length(readSeqs));

    // Launch kernel 10 times!
    double start, finish;
    start = sysTime();

//    for (unsigned i = 0; i < 10; i++)
//    {
        // Sort reads.
        _hashReadsKernel<<<blocksPerGrid, threadsPerBlock>>>(deviceReadSeqsView, view(hashes), view(idx));
        cudaDeviceSynchronize();
        thrust::sort_by_key(hashes.begin(), hashes.end(), idx.begin());

        // Map reads.
        _mapReadsKernel<<<blocksPerGrid, threadsPerBlock>>>(deviceIndexView, deviceReadSeqsView, view(idx));
//    }

    cudaDeviceSynchronize();

    finish = sysTime();
    std::cout << (finish - start)/10 << " sec" << std::endl;
}

void mapReads(TGenomeIndex & index, TReadSeqs & readSeqs, GPU const & /* tag */)
{
    _mapReads(index, readSeqs, GPU());
    cudaDeviceReset();
}
