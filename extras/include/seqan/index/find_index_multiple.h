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

#ifndef SEQAN_EXTRAS_INDEX_FIND_INDEX_MULTIPLE_H
#define SEQAN_EXTRAS_INDEX_FIND_INDEX_MULTIPLE_H

namespace seqan {

// ============================================================================
// Tags
// ============================================================================

// ----------------------------------------------------------------------------
// Tag Multiple
// ----------------------------------------------------------------------------

template <typename TSpec>
struct Multiple;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Member FinderText_
// ----------------------------------------------------------------------------

struct FinderText_;

template <typename TText, typename TPattern, typename TSpec>
struct Member<Finder2<TText, TPattern, TSpec>, FinderText_>
{
    typedef Holder<TText>   Type;
};

template <typename TText, typename TViewSpec, typename TIndexSpec, typename TPattern, typename TSpec>
struct Member<Finder2<Index<ContainerView<TText, TViewSpec>, TIndexSpec>, TPattern, TSpec>, FinderText_>
{
    typedef typename View<Index<TText, TIndexSpec> >::Type  Type;
};

template <typename TText, typename TViewSpec, typename TSSetSpec, typename TIndexSpec, typename TPattern, typename TSpec>
struct Member<Finder2<Index<StringSet<ContainerView<TText, TViewSpec>, TSSetSpec>, TIndexSpec>, TPattern, TSpec>, FinderText_>
{
    typedef typename View<Index<StringSet<ContainerView<TText, TViewSpec>, TSSetSpec>, TIndexSpec> >::Type  Type;
};

// ----------------------------------------------------------------------------
// Member FinderHistory_
// ----------------------------------------------------------------------------

struct FinderHistory_;

template <typename TText, typename TPattern, typename TSpec>
struct Member<Finder2<TText, TPattern, Multiple<TSpec> >, FinderHistory_>
{
    typedef typename TextIterator_<TText, TSpec>::Type      TTextIterator_;
    typedef typename HistoryStack_<TTextIterator_>::Type    Type;
};

template <typename TText, typename TViewSpec, typename TIndexSpec, typename TPattern, typename TSpec>
struct Member<Finder2<Index<ContainerView<TText, TViewSpec>, TIndexSpec>, TPattern, Multiple<TSpec> >, FinderHistory_>
{
    typedef Index<TText, TIndexSpec>                        TIndex_;
    typedef Finder2<TIndex_, TPattern, Multiple<TSpec> >    TFinder_;
    typedef typename Member<TFinder_, FinderHistory_>::Type TFinderHistory_;
    typedef typename View<TFinderHistory_>::Type            Type;
};

template <typename TText, typename TViewSpec, typename TSSetSpec, typename TIndexSpec, typename TPattern, typename TSpec>
struct Member<Finder2<Index<StringSet<ContainerView<TText, TViewSpec>, TSSetSpec>, TIndexSpec>, TPattern, Multiple<TSpec> >, FinderHistory_>
{
    typedef Index<StringSet<TText, TSSetSpec>, TIndexSpec>  TIndex_;
    typedef Finder2<TIndex_, TPattern, Multiple<TSpec> >    TFinder_;
    typedef typename Member<TFinder_, FinderHistory_>::Type TFinderHistory_;
    typedef typename View<TFinderHistory_>::Type            Type;
};

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class Finder
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec>
struct Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<TSpec> >
{
    typedef Index<TText, TIndexSpec>                TIndex;

    typename Member<Finder2, FinderText_>::Type     _index;

    Finder2() :
        _index()
    {}

    Finder2(TIndex & index) :
        _index(index)
    {}
};

template <typename TText, typename TIndexSpec, typename TPattern, typename TDistance, typename TSpec>
struct Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<Backtracking<TDistance, TSpec> > >
{
    typedef Index<TText, TIndexSpec>                TIndex;

    typename Member<Finder2, FinderText_>::Type     _index;
    typename Member<Finder2, FinderHistory_>::Type  _history;
    typename Size<TPattern>::Type                   _historyLength;
    // TODO(esiragusa): change type to typename Size<typename Value<TPattern>::Type>::Type.

    Finder2() :
        _index()
    {}

    Finder2(TIndex & index) :
        _index(index)
    {}
};

// ----------------------------------------------------------------------------
// Class Proxy                                                         [Finder]
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec>
class Proxy<Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<TSpec> > >
{
public:
    typedef Index<TText, TIndexSpec>                    TIndex;
    typedef typename TextIterator_<TIndex, TSpec>::Type TTextIterator;
    typedef typename Iterator<TPattern, Standard>::Type TPatternIterator;

    TTextIterator const & _textIt;
    unsigned _patternIt;

    template <typename TFinder>
    SEQAN_FUNC
    Proxy(TFinder const & finder) :
        _textIt(finder._textIt)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Delegated
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec>
struct Delegated<Finder2<TText, TPattern, Multiple<TSpec> > >
{
    typedef Proxy<Finder2<TText, TPattern, Multiple<TSpec> > > Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function textIterator()
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec>
SEQAN_FUNC typename TextIterator_<TText, Multiple<TSpec> >::Type &
textIterator(Proxy<Finder2<TText, TPattern, Multiple<TSpec> > > & finder)
{
    return finder._textIt;
}

template <typename TText, typename TPattern, typename TSpec>
SEQAN_FUNC typename TextIterator_<TText, Multiple<TSpec> >::Type const &
textIterator(Proxy<Finder2<TText, TPattern, Multiple<TSpec> > > const & finder)
{
    return finder._textIt;
}

// ----------------------------------------------------------------------------
// Function find()
// ----------------------------------------------------------------------------

//template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec, typename TDelegate>
//inline void
//find(Finder2<Index<TText, TIndexSpec>, TPattern,  Multiple<TSpec> > & finder, TPattern & pattern, TDelegate & delegate)
//{
//    typedef Index<TText, TIndexSpec>                                        TIndex;
//    typedef Finder2<TIndex, TPattern,  Multiple<TSpec> >                    TFinder;
//    typedef Delegator<TFinder, TDelegate>                                   TDelegator;
//    typedef typename Iterator<TPattern, Standard>::Type                     TPatternIter;
//
//    // Use a delegator object to delegate this finder instead of the pooled finders.
//    TDelegator delegator(finder, delegate);
//
//    // Initialize the pool.
//    _initPool(finder, omp_get_max_threads());
//
//    // Find all patterns in parallel.
//    TPatternIter patternBegin = begin(pattern, Standard());
//    TPatternIter patternEnd = end(pattern, Standard());
//    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic))
//    for (TPatternIter patternIt = patternBegin; patternIt != patternEnd; ++patternIt)
//    {
////        finder._patternIt[omp_get_thread_num()] = patternIt;
//        clear(finder._pool[omp_get_thread_num()]);
//        find(finder._pool[omp_get_thread_num()], value(patternIt), delegator);
//    }
//}

// ----------------------------------------------------------------------------
// Function view()
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec>
inline typename View<Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<TSpec> > >::Type
view(Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<TSpec> > & finder)
{
    typename View<Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<TSpec> > >::Type finderView;

    finderView._index = view(value(finder._index));

    return finderView;
}

template <typename TText, typename TIndexSpec, typename TPattern, typename TDistance, typename TSpec>
inline typename View<Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<Backtracking<TDistance, TSpec> > > >::Type
view(Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<Backtracking<TDistance, TSpec> > > & finder)
{
    typename View<Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<Backtracking<TDistance, TSpec> > > >::Type finderView;

    finderView._index = view(value(finder._index));
    finderView._history = view(finder._history);
    finderView._historyLength = finder._historyLength;
//    finderView._idxs = view(finder._idxs);
//    finderView._hashes = view(finder._hashes);

    return finderView;
}

}

#endif  // #ifndef SEQAN_EXTRAS_INDEX_FIND_INDEX_MULTIPLE_H
