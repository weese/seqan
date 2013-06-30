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

#ifndef SEQAN_EXTRAS_INDEX_FIND_INDEX_H_
#define SEQAN_EXTRAS_INDEX_FIND_INDEX_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// TODO(esiragusa): Remove this when there will be a base finder class.
template <typename TText, typename TPattern, typename TSpec = void>
struct Finder2;

// TODO(esiragusa): Remove this.
template <typename TDistance, typename TSpec>
struct Backtracking;

// TODO(esiragusa): Remove this.
template <typename TSpec>
struct View;

// TODO(esiragusa): Remove this.
template <typename TSpec>
struct RemoveView;

// TODO(esiragusa): Remove this.
template <typename TSpec>
struct IsView;

// TODO(esiragusa): Remove this.
template <typename TSpec>
struct IsDevice;

// TODO(esiragusa): Remove this.
template <typename TObject, typename T1, typename T2>
struct IfView;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction View
// ----------------------------------------------------------------------------
// TODO(esiragusa): move this function in the base finder class.

template <typename TText, typename TPattern, typename TSpec>
struct View<Finder2<TText, TPattern, TSpec> >
{
    typedef Finder2<typename View<TText>::Type, typename View<TPattern>::Type, TSpec>   Type;
};

// ----------------------------------------------------------------------------
// Metafunction RemoveView
// ----------------------------------------------------------------------------
// TODO(esiragusa): move this function in the base finder class.

template <typename TText, typename TPattern, typename TSpec>
struct RemoveView<Finder2<TText, TPattern, TSpec> >
{
    typedef Finder2<typename RemoveView<TText>::Type, typename RemoveView<TPattern>::Type, TSpec>   Type;
};

// ----------------------------------------------------------------------------
// Metafunction IsView
// ----------------------------------------------------------------------------
// TODO(esiragusa): move this function in the base finder class.

template <typename TText, typename TPattern, typename TSpec>
struct IsView<Finder2<TText, TPattern, TSpec> > : IsView<TText> {};

// ----------------------------------------------------------------------------
// Metafunction IsDevice
// ----------------------------------------------------------------------------
// TODO(esiragusa): move this function in the base finder class.

template <typename TText, typename TPattern, typename TSpec>
struct IsDevice<Finder2<TText, TPattern, TSpec> > : IsDevice<TText> {};

// ----------------------------------------------------------------------------
// Metafunction TextIterator_
// ----------------------------------------------------------------------------

// TODO(esiragusa): move this function in the base finder class.
template <typename TText, typename TSpec>
struct TextIterator_
{
    typedef typename Iterator<TText>::Type  Type;
};

// ----------------------------------------------------------------------------
// Metafunction TextIterator_
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TSpec>
struct TextIterator_<Index<TText, TIndexSpec>, TSpec>
{
    typedef typename Iterator<Index<TText, TIndexSpec>, TopDown<> >::Type  Type;
};

template <typename TText, typename TIndexSpec, typename TDistance, typename TSpec>
struct TextIterator_<Index<TText, TIndexSpec>, Backtracking<TDistance, TSpec> >
{
    typedef typename Iterator<Index<TText, TIndexSpec>, TopDown<ParentLinks<> > >::Type  Type;
};

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class Finder
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec>
struct Finder2<Index<TText, TIndexSpec>, TPattern, TSpec>
{
    typedef Index<TText, TIndexSpec>                    TIndex;
    typedef typename TextIterator_<TIndex, TSpec>::Type TTextIterator;

    TTextIterator _textIt;

    SEQAN_FUNC
    Finder2() {}

    SEQAN_FUNC
    Finder2(TIndex const & index) :
        _textIt(index)
    {}

    SEQAN_FUNC
    Finder2(TTextIterator const & textIt) :
        _textIt(textIt)
    {}
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function textIterator()
// ----------------------------------------------------------------------------
// TODO(esiragusa): move this function in the base finder class.

template <typename TText, typename TPattern, typename TSpec>
SEQAN_FUNC typename TextIterator_<TText, TSpec>::Type &
textIterator(Finder2<TText, TPattern, TSpec> & finder)
{
    return finder._textIt;
}

template <typename TText, typename TPattern, typename TSpec>
SEQAN_FUNC typename TextIterator_<TText, TSpec>::Type const &
textIterator(Finder2<TText, TPattern, TSpec> const & finder)
{
    return finder._textIt;
}

// ----------------------------------------------------------------------------
// Function preprocess()
// ----------------------------------------------------------------------------
// TODO(esiragusa): move this function in the base finder class.

template <typename TText, typename TPattern, typename TSpec>
SEQAN_FUNC void
preprocess(Finder2<TText, TPattern, TSpec> & /* finder */, TPattern const & /* pattern */)
{}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec>
SEQAN_FUNC void
clear(Finder2<Index<TText, TIndexSpec>, TPattern, TSpec> & finder)
{
    goRoot(textIterator(finder));
}

// ----------------------------------------------------------------------------
// Function find()
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPattern, typename TDelegate>
SEQAN_FUNC void
find(Finder2<Index<TText, TIndexSpec>, TPattern, FinderSTree> & finder,
     TPattern const & pattern,
     TDelegate & delegate)
{
    if (goDown(textIterator(finder), pattern)) delegate(finder);
}

template <typename TText, typename TIndexSpec, typename TPattern, typename TDistance, typename TSpec, typename TDelegate>
SEQAN_FUNC void
find(Finder2<Index<TText, TIndexSpec>, TPattern, Backtracking<TDistance, TSpec> > & finder,
     TPattern const & pattern,
     TDelegate & delegate)
{
    if (goDown(textIterator(finder), pattern)) delegate(finder);
}

}

#endif  // #ifndef SEQAN_EXTRAS_INDEX_FIND_INDEX_H_
