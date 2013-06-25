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
// Forwards
// ==========================================================================

extern void findKernel(Index<StringSet<DnaString, Owner<ConcatDirect<> > >, FMIndex<> > & index, DnaString & pattern);

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function testInfix()
// --------------------------------------------------------------------------

//template <typename TAlphabet>
//void testInfix()
//{
//    typedef String<TAlphabet>                           TString;
//    typedef typename View<TString>::Type                TStringView;
//    typedef typename Infix<TString>::Type               TStringInfix;
//    typedef typename Infix<TStringView>::Type           TStringViewInfix;
//
//    TString s = "AAACCCGGGTTT";
//    TStringInfix sInfix = infix(s, 3, 6);
//
//    TStringView sView = view(s);
//    TStringViewInfix sViewInfix = infix(sView, 3, 6);
//
//    SEQAN_ASSERT(isEqual(sInfix, sViewInfix));
//}

// --------------------------------------------------------------------------
// Function testStringSet()
// --------------------------------------------------------------------------

//template <typename TAlphabet>
//void testStringSet()
//{
//    typedef String<TAlphabet>                           TString;
//    typedef StringSet<TString, Owner<ConcatDirect<> > > TStringSet;
//    typedef typename Device<TStringSet>::Type           TDeviceStringSet;
//    typedef typename View<TString>::Type                TStringView;
//    typedef typename View<TStringSet>::Type             TStringSetView;
//    typedef typename View<TDeviceStringSet>::Type       TDeviceStringSetView;
//
//    TStringSet ss;
//    appendValue(ss, "AAAAAAAA");
//    appendValue(ss, "CCCCCCC");
//    appendValue(ss, "GGGGGGGGGGGGGG");
//    appendValue(ss, "T");
//
//    TStringSetView ssView = view(ss);
//
//    SEQAN_ASSERT_EQ(length(ss), length(ssView));
//    for (unsigned i = 0; i < length(ss); ++i)
//        SEQAN_ASSERT(isEqual(ss[i], ssView[i]));
//
//    TDeviceStringSet deviceSs;
//    assign(deviceSs, ss);
//    TDeviceStringSetView deviceSsView = view(deviceSs);
//}

// --------------------------------------------------------------------------
// Function testIndex()
// --------------------------------------------------------------------------

template <typename TAlphabet, typename TIndexSpec>
void testIndex()
{
    typedef String<TAlphabet>                           TString;
    typedef StringSet<TString, Owner<ConcatDirect<> > > TStringSet;
    typedef Index<TStringSet, TIndexSpec>               TIndex;

    // Create a text.
    TStringSet text;
    appendValue(text, "ATAAAAAAAAA");
    appendValue(text, "CCCCTACCC");

    // Create an index on the reversed text.
    TIndex index(text);
    reverse(text);
    indexCreate(index);
    reverse(text);

    // Create a pattern.
    TString pattern("TA");

    // Launch find kernel.
    findKernel(index, pattern);
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

int main(int /* argc */, char const ** /* argv */)
{
//    testInfix<Dna>();
//    testStringSet<Dna>();
    testIndex<Dna, FMIndex<> >();
};
