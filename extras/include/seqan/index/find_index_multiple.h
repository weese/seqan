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

// ----------------------------------------------------------------------------
// Tag Factory_
// ----------------------------------------------------------------------------

struct Factory_;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class Finder
// ----------------------------------------------------------------------------

template <typename TText, typename TPatterns, typename TSpec>
struct Finder2<TText, TPatterns, Multiple<TSpec> >
{
    typename Member<Finder2, Factory_>::Type    _factory;

    SEQAN_HOST_DEVICE
    Finder2() {}

    Finder2(TText & text) :
        _factory(text)
    {}
};

template <typename TText, typename TIndexSpec, typename TPatterns, typename TSpec>
struct Finder2<Index<TText, TIndexSpec>, TPatterns, Multiple<TSpec> >
{
    typedef Index<TText, TIndexSpec>            TIndex;

    typename Member<Finder2, Factory_>::Type    _factory;

    SEQAN_HOST_DEVICE
    Finder2() {}

    Finder2(TIndex & index) :
        _factory(index)
    {}
};

// ----------------------------------------------------------------------------
// Class Proxy                                                         [Finder]
// ----------------------------------------------------------------------------

template <typename TText, typename TPatterns, typename TSpec>
class Proxy<Finder2<TText, TPatterns, Multiple<TSpec> > >
{
public:
    typedef typename TextIterator_<TText, TSpec>::Type      TTextIterator;
    typedef typename Iterator<TPatterns, Standard>::Type    TPatternIterator;

    TTextIterator const &   _textIt;
    unsigned                _patternIt;

    template <typename TFinder>
    SEQAN_HOST_DEVICE
    Proxy(TFinder const & finder) :
        _textIt(textIterator(finder))
    {}
};

// ============================================================================
// Members
// ============================================================================

// ----------------------------------------------------------------------------
// Member Factory_
// ----------------------------------------------------------------------------

template <typename TText, typename TPatterns, typename TSpec>
struct Member<Finder2<TText, TPatterns, Multiple<TSpec> >, Factory_>
{
    typedef Factory<typename TextIterator_<TText, TSpec>::Type>  Type;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Delegated
// ----------------------------------------------------------------------------

template <typename TText, typename TPatterns, typename TSpec>
struct Delegated<Finder2<TText, TPatterns, Multiple<TSpec> > >
{
    typedef Proxy<Finder2<TText, TPatterns, Multiple<TSpec> > > Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function textIterator()
// ----------------------------------------------------------------------------

template <typename TText, typename TPatterns, typename TSpec>
inline SEQAN_HOST_DEVICE typename TextIterator_<TText, Multiple<TSpec> >::Type &
textIterator(Proxy<Finder2<TText, TPatterns, Multiple<TSpec> > > & finder)
{
    return finder._textIt;
}

template <typename TText, typename TPatterns, typename TSpec>
inline SEQAN_HOST_DEVICE typename TextIterator_<TText, Multiple<TSpec> >::Type const &
textIterator(Proxy<Finder2<TText, TPatterns, Multiple<TSpec> > > const & finder)
{
    return finder._textIt;
}

// ----------------------------------------------------------------------------
// Function _find(); ExecHost
// ----------------------------------------------------------------------------

template <typename TText, typename TPatterns, typename TSpec, typename TDelegate>
inline void
_find(Finder2<TText, TPatterns, Multiple<TSpec> > & finder,
      TPatterns & patterns,
      TDelegate & delegate,
      ExecHost const & /* tag */)
{
    typedef Finder2<TText, TPatterns,  Multiple<TSpec> >    TFinder;
    typedef typename Value<TPatterns>::Type                 TPattern;
    typedef Finder2<TText, TPattern, TSpec>                 TFinderSimple;
    typedef typename View<TFinder>::Type                    TFinderView;
    typedef typename View<TFinderSimple>::Type              TFinderSimpleView;
    typedef Proxy<TFinderView>                              TFinderProxy;
    typedef Delegator<TFinderProxy, TDelegate>              TDelegator;
    typedef typename Iterator<TPatterns, Standard>::Type    TPatternsIter;

    // Initialize the iterator factory.
    setMaxHistoryLength(finder._factory, length(back(patterns)));
    setMaxObjects(finder._factory, omp_get_max_threads());
    build(finder._factory);

    // Use a finder view.
    TFinderView finderView = view(finder);

    // Instantiate a finder.
    TFinderSimpleView simpleFinder = getObject(finderView._factory, omp_get_thread_num());

    // Instantiate a finder proxy to be delegated.
    TFinderProxy finderProxy(simpleFinder);

    // Instantiate a delegator object to delegate the finder proxy instead of the simple finder.
    TDelegator delegator(finderProxy, delegate);

    // Find all patterns in parallel.
    TPatternsIter patternsBegin = begin(patterns, Standard());
    TPatternsIter patternsEnd = end(patterns, Standard());
//    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic))
    for (TPatternsIter patternsIt = patternsBegin; patternsIt != patternsEnd; ++patternsIt)
    {
        clear(simpleFinder);
        find(simpleFinder, value(patternsIt), delegator);
    }
}

// ----------------------------------------------------------------------------
// Function find()
// ----------------------------------------------------------------------------

template <typename TText, typename TPatterns, typename TSpec, typename TDelegate>
inline void
find(Finder2<TText, TPatterns, Multiple<TSpec> > & finder, TPatterns & patterns, TDelegate & delegate)
{
    typedef Finder2<TText, TPatterns,  Multiple<TSpec> > TFinder;

    _find(finder, patterns, delegate, typename ExecSpace<TFinder>::Type());
}

// ----------------------------------------------------------------------------
// Function view()
// ----------------------------------------------------------------------------

template <typename TText, typename TPatterns, typename TSpec>
inline typename View<Finder2<TText, TPatterns, Multiple<TSpec> > >::Type
view(Finder2<TText, TPatterns, Multiple<TSpec> > & finder)
{
    typename View<Finder2<TText, TPatterns, Multiple<TSpec> > >::Type finderView;

    finderView._factory = view(finder._factory);
//    finderView._idxs = view(finder._idxs);
//    finderView._hashes = view(finder._hashes);

    return finderView;
}

}

#endif  // #ifndef SEQAN_EXTRAS_INDEX_FIND_INDEX_MULTIPLE_H
