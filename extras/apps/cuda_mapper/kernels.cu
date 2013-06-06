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

using namespace seqan;

// --------------------------------------------------------------------------
// Function mapReadsGPU()
// --------------------------------------------------------------------------

template <typename TIndexView, typename TReadSeqsView>
__global__ void
mapReadsGPU(TIndexView index, TReadSeqsView readSeqs)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

//    mapRead(index, readSeqs[idx]);

    typename Iterator<TIndexView, TopDown<> >::Type it(index);

    unsigned occurrences = goDown(it, readSeqs[idx]) ? countOccurrences(it) : 0;

    printf("index=%i, occurrences=%d\n", idx, occurrences);
}

// --------------------------------------------------------------------------
// Function mapReads()                                                  [GPU]
// --------------------------------------------------------------------------

void
mapReads(Index<StringSet<String<Dna>, Owner<ConcatDirect<> > >, FMIndex<> > & index,
         StringSet<String<Dna>, Owner<ConcatDirect<> > > & readSeqs,
         GPU const & /* tag */)
{
    typedef Index<StringSet<String<Dna>, Owner<ConcatDirect<> > >, FMIndex<> >  TIndex;
    typedef StringSet<String<Dna>, Owner<ConcatDirect<> > >                     TReadSeqs;

    typedef typename Device<TIndex>::Type               TDeviceIndex;
    typedef typename Device<TReadSeqs>::Type            TDeviceReadSeqs;

    // Copy index to device.
    TDeviceIndex deviceIndex;
    assign(deviceIndex, index);

    // Copy read seqs to device.
    TDeviceReadSeqs deviceReadSeqs;
    assign(deviceReadSeqs, readSeqs);

    // Launch kernel.
    mapReadsGPU<<<10,100>>>(view(deviceIndex), view(deviceReadSeqs));
    cudaDeviceSynchronize();
}
