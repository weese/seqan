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

#include "cuda_demo.h"

using namespace seqan;

// ==========================================================================
// Kernels
// ==========================================================================

// --------------------------------------------------------------------------
// Kernel _findKernel()
// --------------------------------------------------------------------------

template <typename TIndex, typename TString>
__global__ void
_findKernel(TIndex index, TString pattern)
{
    typedef typename Iterator<TIndex, TopDown<> >::Type TIterator;
    typedef typename EdgeLabel<TIterator>::Type         TEdgeLabel;
    typedef typename Fibre<TIndex, FibreText>::Type     TTextFibre;
    typedef typename Infix<TTextFibre const>::Type      TRepresentative;

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    printf("index=%i\n", idx);

    // Print the Compressed SA values.
    printf("lengthSA=%ld\n", length(indexSA(index)));
    for (unsigned i = 0; i < length(indexSA(index)); ++i)
        printf("<%ld,%ld>\n", indexSA(index)[i].i1, indexSA(index)[i].i2);
//        printf("%ld\n", indexSA(index)[i]);

    // Instantiate a virtual suffix tree iterator.
    TIterator it(index);

    // At root.
    printf("isRoot()=%d\n", isRoot(it));
    printf("repLength()=%ld\n", repLength(it));
    printf("countOccurrences()=%ld\n", countOccurrences(it));

    // Visit the leftmost children of the root.
    if (goDown(it))
    {
        // Visit all the siblings at depth one.
        do
        {
            printf("repLength()=%ld\n", repLength(it));
            printf("parentEdgeLabel()=%c\n", static_cast<char>(parentEdgeLabel(it)));
            printf("countOccurrences()=%ld\n", countOccurrences(it));
            printf("isLeaf()=%d\n", isLeaf(it));
        }
        while (goRight(it));
    }

    // Restart from root.
    goRoot(it);
    printf("goRoot()\n");

    // Search the pattern.
    printf("goDown(pattern)=%d\n", goDown(it, pattern));
    printf("repLength()=%ld\n", repLength(it));
    printf("countOccurrences()=%ld\n", countOccurrences(it));
    printf("isLeaf()=%d\n", isLeaf(it));
}

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function findKernel()
// --------------------------------------------------------------------------

void findKernel(Index<StringSet<DnaString, Owner<ConcatDirect<> > >, FMIndex<> > & index, DnaString & pattern)
{
    typedef DnaString                                   TString;
    typedef StringSet<TString, Owner<ConcatDirect<> > > TStringSet;
    typedef Index<TStringSet, FMIndex<> >               TIndex;

//template <typename TIndex, typename TString>
//void findKernel(TIndex & index, TString & pattern)
//{

    typedef typename Device<TIndex>::Type               TDeviceIndex;
    typedef typename Device<TString>::Type              TDeviceString;

    // Copy index to device.
    TDeviceIndex deviceIndex;
    assign(deviceIndex, index);

    // Copy pattern to device.
    TDeviceString devicePattern;
    assign(devicePattern, pattern);

    _findKernel<<< 1,1 >>>(view(deviceIndex), view(devicePattern));
    cudaDeviceSynchronize();
}

//void findKernel(Index<StringSet<DnaString, Owner<ConcatDirect<> > >, FMIndex<> > & index, DnaString & pattern);
