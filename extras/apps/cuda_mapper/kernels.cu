// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2013 NVIDIA Corporation
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
//     * Neither the name of NVIDIA Corporation nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL NVIDIA CORPORATION BE LIABLE
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
#include <thrust/reduce.h>

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
    typedef typename Value<TReadSeqsView>::Type         TReadSeq;
    typedef typename Value<TReadSeq>::Type              TAlphabet;
    typedef Shape<TAlphabet, UngappedShape<16> >        TShape;

    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    // Return silently if there is no job left.
    if (idx >= length(readSeqs)) return;

    // Compute the hash of a read.
    TShape shape;
    hashes[idx] = hash(shape, begin(readSeqs[idx], Standard()));

    // Fill idxs with the identity permutation.
    idxs[idx] = idx;
}

// --------------------------------------------------------------------------
// Kernel _mapReadsKernel()
// --------------------------------------------------------------------------

template <typename TIndexView, typename TReadSeqsView, typename TIdxView, typename TOccView>
__global__ void
_mapReadsKernel(TIndexView index, TReadSeqsView readSeqs, TIdxView idxs, TOccView occs)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    // Return silently if there is no job left.
    if (idx >= length(readSeqs)) return;

    // Pick the read id to map.
    idx = idxs[idx];

    // Map the read.
    occs[idx] = mapRead(index, readSeqs[idx], (unsigned)idx);
}

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function mapReads()                                                  [GPU]
// --------------------------------------------------------------------------

void mapReads(TGenomeIndex & index, TReadSeqs & readSeqs, GPU const & /* tag */)
{
    typedef typename Device<TGenomeIndex>::Type         TDeviceIndex;
    typedef typename Device<TReadSeqs>::Type            TDeviceReadSeqs;
    typedef typename View<TDeviceIndex>::Type           TDeviceIndexView;
    typedef typename View<TDeviceReadSeqs>::Type        TDeviceReadSeqsView;

    typedef typename Value<TReadSeqs>::Type             TReadSeq;
    typedef typename Value<TReadSeq>::Type              TAlphabet;
    typedef Shape<TAlphabet, UngappedShape<16> >        TShape;
    typedef typename Value<TShape>::Type                THash;

    // Copy index to device.
    TDeviceIndex deviceIndex;
    assign(deviceIndex, index);

    // Copy read seqs to device.
    TDeviceReadSeqs deviceReadSeqs;
    assign(deviceReadSeqs, readSeqs);

    cudaDeviceSynchronize();

    // Instantiate views of device objects.
    TDeviceIndexView deviceIndexView = view(deviceIndex);
    TDeviceReadSeqsView deviceReadSeqsView = view(deviceReadSeqs);

    // Setup kernel parameters.
    unsigned threadsPerBlock = 256;
    unsigned blocksPerGrid = (length(readSeqs) + threadsPerBlock - 1) / threadsPerBlock;

    thrust::device_vector<THash> hashes(length(readSeqs));
    thrust::device_vector<__uint32> idx(length(readSeqs));
    thrust::device_vector<__uint64> occs(length(readSeqs));

    cudaDeviceSynchronize();
    cudaPrintFreeMemory();

    // Launch kernel 10 times!
    double start, finish;
    start = sysTime();

    unsigned n_tests = 10;
    for (unsigned i = 0; i < n_tests; i++)
    {
        // Sort reads.
        _hashReadsKernel<<<blocksPerGrid, threadsPerBlock>>>(deviceReadSeqsView, view(hashes), view(idx));
        cudaDeviceSynchronize();
        thrust::sort_by_key(hashes.begin(), hashes.end(), idx.begin());
        cudaDeviceSynchronize();

        // Map reads.
        _mapReadsKernel<<<blocksPerGrid, threadsPerBlock>>>(deviceIndexView, deviceReadSeqsView, view(idx), view(occs));
        cudaDeviceSynchronize();
    }

    finish = sysTime();
    std::cout << (finish - start)/n_tests << " sec" << std::endl;

    thrust::host_vector<__uint64> h_occs(occs);
    const __uint64 n_occ = thrust::reduce( h_occs.begin(), h_occs.end() );
    std::cout << n_occ << std::endl;

    cudaPrintFreeMemory();
}
