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

#ifndef SEQAN_EXTRAS_CUDAMAPPER_MAPPER_H_
#define SEQAN_EXTRAS_CUDAMAPPER_MAPPER_H_

// ============================================================================
// Prerequisites
// ============================================================================

#include "types.h"

using namespace seqan;


// ============================================================================
// Tags
// ============================================================================

// ----------------------------------------------------------------------------
// Tag Ranges_
// ----------------------------------------------------------------------------

struct Ranges_;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class Hits
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec = void>
struct Hits
{
    typename Member<Hits, Ranges_>::Type    ranges;

    template <typename TFinder>
    SEQAN_FUNC void
    operator() (TFinder const & finder)
    {
        appendRange(*this, finder);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Member Ranges_
// ----------------------------------------------------------------------------

struct Ranges_;

namespace seqan {
template <typename TIndex, typename TSpec>
struct Member<Hits<TIndex, TSpec>, Ranges_>
{
    typedef Pair<typename Size<TIndex>::Type>   TRange_;
    typedef String<TRange_>                     Type;
};

#ifdef PLATFORM_CUDA
template <typename TIndex, typename TSpec>
struct Member<Hits<TIndex, Device<TSpec> >, Ranges_>
{
    typedef Pair<typename Size<TIndex>::Type>   TRange_;
    typedef thrust::device_vector<TRange_>      Type;
};
#endif

template <typename TIndex, typename TSpec>
struct Member<Hits<TIndex, View<TSpec> >, Ranges_>
{
    typedef typename Member<Hits<TIndex, TSpec>, Ranges_>::Type TRanges_;
    typedef ContainerView<TRanges_, Resizable<TSpec> >          Type;
};
}

// ----------------------------------------------------------------------------
// Metafunction View
// ----------------------------------------------------------------------------

namespace seqan {
template <typename TIndex, typename TSpec>
struct View<Hits<TIndex, TSpec> >
{
    typedef Hits<TIndex, View<TSpec> >  Type;
};
}

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function reserve()
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec, typename TSize>
inline void
reserve(Hits<TIndex, TSpec> & hits, TSize newCapacity)
{
    reserve(hits.ranges, newCapacity, Exact());
}

// ----------------------------------------------------------------------------
// Function view()
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec>
typename View<Hits<TIndex, TSpec> >::Type
view(Hits<TIndex, TSpec> & hits)
{
    typename View<Hits<TIndex, TSpec> >::Type hitsView;

    hitsView.ranges = view(hits.ranges);

    return hitsView;
}

// ----------------------------------------------------------------------------
// Function appendRange()
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec, typename TFinder>
inline void
appendRange(Hits<TIndex, TSpec> & hits, TFinder const & finder)
{
    SEQAN_OMP_PRAGMA(critical(Hits_appendRange))
    appendValue(hits.ranges, range(textIterator(finder)));
}

#ifdef PLATFORM_CUDA
template <typename TIndex, typename TSpec, typename TFinder>
SEQAN_FUNC void
appendRange(Hits<TIndex, View<Device<TSpec> > > & hits, TFinder const & finder)
{
    // TODO(esiragusa): Global lock.
    appendValue(hits.ranges, range(textIterator(finder)));
}
#endif

// ----------------------------------------------------------------------------
// Function _mapReads()
// ----------------------------------------------------------------------------

template <typename TIndex, typename TReadSeqs, typename TSpace>
inline void
_mapReads(TIndex & index, TReadSeqs & readSeqs)
{
    typedef Multiple<Backtracking<HammingDistance> >    TFinderSpec;
//    typedef Multiple<FinderSTree>                       TFinderSpec;
    typedef Finder2<TIndex, TReadSeqs, TFinderSpec>     TFinder;

    // Instantiate a multiple finder.
    TFinder finder(index);

    // Instantiate an object to save the hits.
    Hits<TIndex, TSpace> hits;

    // Reserve space for hits.
    reserve(hits, length(readSeqs));

    // Find hits.
    find(finder, readSeqs, hits);

    std::cout << "Ranges count:\t\t\t" << length(hits.ranges) << std::endl;
}

// ----------------------------------------------------------------------------
// Function mapReads()
// ----------------------------------------------------------------------------

template <typename TIndex, typename TReadSeqs>
inline void
mapReads(TIndex & index, TReadSeqs & readSeqs, CPU const & /* tag */)
{
    // Map reads.
    _mapReads<TIndex, TReadSeqs, void>(index, readSeqs);
}

#endif  // #ifndef SEQAN_EXTRAS_CUDAMAPPER_MAPPER_H_
