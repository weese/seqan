// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
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

#ifndef EXTRAS_TESTS_INDEX_VIEW_H_
#define EXTRAS_TESTS_INDEX_VIEW_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index_extras.h>

#include <seqan/misc/misc_view.h>
#include <seqan/index/index_view.h>

using namespace seqan;

// ============================================================================
// Functions
// ============================================================================

// ============================================================================
// Tests
// ============================================================================

// ----------------------------------------------------------------------------
// Test test_index_view_basic
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_index_view_basic)
{
    typedef CharString                                  TText;
    typedef View<CharString>                            TTextView;
    typedef Index<TText, IndexSa<> >                    TIndex;
    typedef Index<TTextView, IndexSa<> >                TIndexView;
    typedef typename Fibre<TIndex, FibreText>::Type     TIndexTextFibre;
    typedef typename Fibre<TIndexView, FibreText>::Type TIndexViewTextFibre;
    typedef typename Fibre<TIndex, FibreSA>::Type       TIndexSAFibre;
    typedef typename Fibre<TIndexView, FibreSA>::Type   TIndexViewSAFibre;

    // ----------------------------------------------------------------------

    TText text("text");
    TTextView textView(text);

    SEQAN_ASSERT_EQ(length(text), length(textView));
    for (unsigned i = 0; i < length(text); ++i)
        SEQAN_ASSERT_EQ(text[i], textView[i]);

    // ----------------------------------------------------------------------

    TIndex index(text);
    indexCreate(index, FibreSA());

    TIndexView indexView;
    indexText(indexView) = TIndexViewTextFibre(indexText(index));
    indexSA(indexView) = TIndexViewSAFibre(indexSA(index));

//    TIndexView indexView2 = toView(index);

    // NOTE(esiragusa): Index view constructor must not take the text view, but the original index!
//    TIndexView indexView(textView);

    // NOTE(esiragusa): Calling the index view constructor with the original index calls the generic TText constructor.
//    TIndexView indexView3(index);

    // NOTE(esiragusa): indexCreate() doesn't work on views, it should never be called.
//    indexCreate(indexView, FibreSA());

    // ----------------------------------------------------------------------

    TIndexTextFibre & textFibre = indexText(index);
    TIndexViewTextFibre & textFibreView = indexText(indexView);

    SEQAN_ASSERT_EQ(length(textFibre), length(textFibreView));
    for (unsigned i = 0; i < length(textFibre); ++i)
        SEQAN_ASSERT_EQ(textFibre[i], textFibreView[i]);

    // ----------------------------------------------------------------------

    TIndexSAFibre & saFibre = indexSA(index);
    TIndexViewSAFibre & saFibreView = indexSA(indexView);

    SEQAN_ASSERT_EQ(length(saFibre), length(saFibreView));
    for (unsigned i = 0; i < length(saFibre); ++i)
        SEQAN_ASSERT_EQ(saFibre[i], saFibreView[i]);
}

#endif  // EXTRAS_TESTS_INDEX_VIEW_H_
