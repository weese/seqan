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

template <typename TText, typename TPattern, typename TSpec>
struct Finder2<TText, TPattern, Multiple<TSpec> >
{
    typename Member<Finder2, Factory_>::Type    _factory;

    Finder2() {}

    Finder2(TText & text) :
        _factory(text)
    {}
};

template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec>
struct Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<TSpec> >
{
    typedef Index<TText, TIndexSpec>            TIndex;

    typename Member<Finder2, Factory_>::Type    _factory;

    Finder2() {}

    Finder2(TIndex & index) :
        _factory(index)
    {}
};

// ----------------------------------------------------------------------------
// Class Proxy                                                         [Finder]
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec>
class Proxy<Finder2<TText, TPattern, Multiple<TSpec> > >
{
public:
    typedef typename TextIterator_<TText, TSpec>::Type  TTextIterator;
    typedef typename Iterator<TPattern, Standard>::Type TPatternIterator;

    TTextIterator const &   _textIt;
    unsigned                _patternIt;

    template <typename TFinder>
    SEQAN_FUNC
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

template <typename TText, typename TPattern, typename TSpec>
struct Member<Finder2<TText, TPattern, Multiple<TSpec> >, Factory_>
{
    typedef Factory<typename TextIterator_<TText, TSpec>::Type>  Type;
};

template <typename TText, typename TIndexSpec, typename TPattern, typename TDistance, typename TSpec>
struct Member<Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<Backtracking<TDistance, TSpec> > >, Factory_>
{
    typedef Factory<typename TextIterator_<Index<TText, TIndexSpec>, Backtracking<TDistance, View<TSpec> > >::Type>  Type;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction FinderSerial_
// ----------------------------------------------------------------------------

template <typename TFinder>
struct FinderSerial_;

template <typename TText, typename TPattern, typename TSpec>
struct FinderSerial_<Finder2<TText, TPattern, Multiple<TSpec> > >
{
    typedef Finder2<TText, typename Value<TPattern>::Type, TSpec>    Type;
};

template <typename TText, typename TIndexSpec, typename TPattern, typename TDistance, typename TSpec>
struct FinderSerial_<Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<Backtracking<TDistance, TSpec> > > >
{
    typedef Finder2<Index<TText, TIndexSpec>, typename Value<TPattern>::Type, Backtracking<TDistance, View<TSpec> > >    Type;
};

// ----------------------------------------------------------------------------
// Metafunction TextIterator_
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TDistance, typename TSpec>
struct TextIterator_<Index<TText, TIndexSpec>, Backtracking<TDistance, View<TSpec> > >
{
    typedef Index<TText, TIndexSpec>                            TIndex_;
//    typedef Backtracking<TDistance, TSpec>                      TIterSpec_;
//    typedef typename TextIterator_<TIndex_, TIterSpec_>::Type   TIter_;

//    typedef typename View<TIter_>::Type                         Type;

    typedef Iter<TIndex_, VSTree<TopDown<ParentLinks<VSTreeIteratorTraits<Preorder_, True, True> > > > >  Type;
};

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
// Function _find(); ExecHost
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec, typename TDelegate>
inline void
_find(Finder2<TText, TPattern, Multiple<TSpec> > & finder,
      TPattern & pattern,
      TDelegate & delegate,
      ExecHost const & /* tag */)
{
    typedef Finder2<TText, TPattern,  Multiple<TSpec> >     TFinder;
    typedef Proxy<TFinder>                                  TFinderProxy;
    typedef Delegator<TFinderProxy, TDelegate>              TDelegator;
    typedef typename FinderSerial_<TFinder>::Type           TSerialFinder;
    typedef typename Iterator<TPattern, Standard>::Type     TPatternIter;

    // Initialize the iterator factory.
    setMaxHistoryLength(finder._factory, length(back(pattern)));
    setMaxObjects(finder._factory, omp_get_max_threads());
    build(finder._factory);

//    static_cast<Nothing>(finder._factory);
//    static_cast<Nothing>(finder._factory._history);
//    static_cast<Nothing>(getObject(finder._factory, 0));

    // Instantiate a serial finder.
    TSerialFinder serialFinder = getObject(finder._factory, omp_get_thread_num());

    // Instantiate a finder proxy to be delegated.
    TFinderProxy finderProxy(serialFinder);

    // Instantiate a delegator object to delegate the finder proxy instead of the serial finder.
    TDelegator delegator(finderProxy, delegate);

    // Find all patterns in parallel.
    TPatternIter patternBegin = begin(pattern, Standard());
    TPatternIter patternEnd = end(pattern, Standard());
//    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic))
    for (TPatternIter patternIt = patternBegin; patternIt != patternEnd; ++patternIt)
    {
        clear(serialFinder);
        find(serialFinder, value(patternIt), delegator);
    }
}

// ----------------------------------------------------------------------------
// Function find()
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec, typename TDelegate>
inline void
find(Finder2<TText, TPattern,  Multiple<TSpec> > & finder, TPattern & pattern, TDelegate & delegate)
{
    typedef Finder2<TText, TPattern,  Multiple<TSpec> > TFinder;

    _find(finder, pattern, delegate, typename ExecSpace<TFinder>::Type());
}

// ----------------------------------------------------------------------------
// Function view()
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec>
inline typename View<Finder2<TText, TPattern, Multiple<TSpec> > >::Type
view(Finder2<TText, TPattern, Multiple<TSpec> > & finder)
{
    typename View<Finder2<TText, TPattern, Multiple<TSpec> > >::Type finderView;

    finderView._factory = view(finder._factory);
//    finderView._idxs = view(finder._idxs);
//    finderView._hashes = view(finder._hashes);

    return finderView;
}

}

#endif  // #ifndef SEQAN_EXTRAS_INDEX_FIND_INDEX_MULTIPLE_H
