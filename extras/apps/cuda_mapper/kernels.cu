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
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class HitsCounter
// ----------------------------------------------------------------------------

template <>
struct HitsCounter<GPU> : HitsCounter<>
{
    template <typename TFinder>
    SEQAN_FUNC void
    operator() (TFinder const & /* finder */)
    {
//        SEQAN_OMP_PRAGMA(atomic)
//        hits += countOccurrences(finder._pool[omp_get_thread_num()]._textIt);
    }
};

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function mapReads()                                                  [GPU]
// --------------------------------------------------------------------------

void mapReads(TGenomeIndex & index, TReadSeqs & readSeqs, GPU const & /* tag */)
{
    typedef typename Device<TGenomeIndex>::Type                             TDeviceIndex;
    typedef typename Device<TReadSeqs>::Type                                TDeviceReadSeqs;
    typedef Finder2<TDeviceIndex, TDeviceReadSeqs, Multiple<FinderSTree> >  TFinder;

    // Copy index to device.
    TDeviceIndex deviceIndex;
    assign(deviceIndex, index);

    // Copy read seqs to device.
    TDeviceReadSeqs deviceReadSeqs;
    assign(deviceReadSeqs, readSeqs);

    // Wait for the copy to finish.
    cudaDeviceSynchronize();

    // Instantiate a multiple finder.
    TFinder finder(deviceIndex);

    // Count hits.
    HitsCounter<GPU> counter;
    find(finder, deviceReadSeqs, counter);
    std::cout << "Hits count:\t\t\t" << counter.hits << std::endl;
}
