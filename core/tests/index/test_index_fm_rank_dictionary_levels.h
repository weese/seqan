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

#ifndef TEST_INDEX_FM_RANK_DICTIONARY_LEVELS_H_
#define TEST_INDEX_FM_RANK_DICTIONARY_LEVELS_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/random.h>
#include <seqan/index.h>

using namespace seqan;

// ----------------------------------------------------------------------------
// Metafunction Size                                           [RankDictionary]
// ----------------------------------------------------------------------------

namespace seqan {
template <typename TValue>
struct Size<RankDictionary<TwoLevels<TValue, unsigned> > >
{
    typedef unsigned    Type;
};
}

// ----------------------------------------------------------------------------
// Function testRankDictionaryGetValue()
// ----------------------------------------------------------------------------

template <typename TAlphabet, typename TStringSpec, typename TSpec>
void testRankDictionaryGetValue(String<TAlphabet, TStringSpec> const & text, TSpec const & /* tag */)
{
    typedef String<TAlphabet, TStringSpec> const        TText;
    typedef typename Iterator<TText>::Type              TTextIterator;
    typedef TwoLevels<TAlphabet, unsigned>              TRankDictionarySpec;
    typedef RankDictionary<TRankDictionarySpec>         TRankDictionary;

    // Build the RankDictionary object to test.
    TRankDictionary dict(text);

    // Scan the text.
    TTextIterator textBegin = begin(text, Standard());
    TTextIterator textEnd = end(text, Standard());
    for (TTextIterator textIt = textBegin; textIt != textEnd; ++textIt)
        SEQAN_ASSERT_EQ(value(textIt), getValue(dict, textIt - textBegin));
}

// ----------------------------------------------------------------------------
// Function testRankDictionaryGetRank()
// ----------------------------------------------------------------------------

template <typename TAlphabet, typename TStringSpec, typename TSpec>
void testRankDictionaryGetRank(String<TAlphabet, TStringSpec> const & text, TSpec const & /* tag */)
{
    typedef String<TAlphabet, TStringSpec> const        TText;
    typedef typename Iterator<TText>::Type              TTextIterator;
    typedef typename Size<TText>::Type                  TTextSize;
    typedef typename ValueSize<TAlphabet>::Type         TAlphabetSize;
    typedef String<TTextSize>                           TRankMap;
    typedef TwoLevels<TAlphabet, unsigned>              TRankDictionarySpec;
    typedef RankDictionary<TRankDictionarySpec>         TRankDictionary;

    TAlphabetSize alphabetSize = ValueSize<TAlphabet>::VALUE;

    // Build the RankDictionary object to test.
    TRankDictionary dict(text);

    // The naive rank dictionary is built while scanning the text.
    TRankMap rank;
    resize(rank, alphabetSize, 0);

    // Scan the text.
    TTextIterator textBegin = begin(text, Standard());
    TTextIterator textEnd = end(text, Standard());
    for (TTextIterator textIt = textBegin; textIt != textEnd; ++textIt)
    {
        // Update the naive rank.
        rank[ordValue(value(textIt))]++;

        // Check the rank for all alphabet symbols.
        for (TAlphabetSize c = 0; c < alphabetSize; ++c)
            SEQAN_ASSERT_EQ(rank[c], getRank(dict, textIt - textBegin, c));
    }
}

// ----------------------------------------------------------------------------
// Test test_rss_getvalue                                      [RankDictionary]
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_rss_getvalue)
{
    {
        String<char> text = "testestestestestestestestestestestestestestestestestest";
        testRankDictionaryGetValue(text, TwoLevels<char, unsigned>());
    }
    {
        String<Dna> text = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
        testRankDictionaryGetValue(text, TwoLevels<Dna, unsigned>());
    }
    {
        String<bool> text;
        for (unsigned i = 0; i < 10; i++)
        {
            appendValue(text, 1);
            appendValue(text, 0);
            appendValue(text, 1);
            appendValue(text, 0);
            appendValue(text, 1);
            appendValue(text, 1);
            appendValue(text, 1);
        }
        testRankDictionaryGetValue(text, TwoLevels<bool, unsigned>());
    }
}

// ----------------------------------------------------------------------------
// Test test_rss_getrank                                       [RankDictionary]
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_rss_getrank)
{
    {
        String<char> text = "testestestestestestestestestestestestestestestestestest";
        testRankDictionaryGetRank(text, TwoLevels<char, unsigned>());
    }
    {
        String<Dna> text = "ACGTACGTACGTACGTACGTACGTACGTACGTCCCCCCCCCCCCCCCCCCGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
        testRankDictionaryGetRank(text, TwoLevels<Dna, unsigned>());
    }
    {
        String<bool> text;
        for (unsigned i = 0; i < 10; i++)
        {
            appendValue(text, 1);
            appendValue(text, 0);
            appendValue(text, 1);
            appendValue(text, 0);
            appendValue(text, 1);
            appendValue(text, 1);
            appendValue(text, 1);
        }
        for (unsigned i = 0; i < 150; i++)
            appendValue(text, 1);
        testRankDictionaryGetRank(text, TwoLevels<bool, unsigned>());
    }
}

// ----------------------------------------------------------------------------
// Test ...                                                    [RankDictionary]
// ----------------------------------------------------------------------------

//SEQAN_DEFINE_TEST(test_rss_sizeof)
//{
//    typedef Dna                                         TAlphabet;
//    typedef Alloc<unsigned>                             TTextSpec;
//    typedef String<TAlphabet, TTextSpec>                TText;
//
//    typedef TwoLevels<TAlphabet, unsigned>              TRankDictionarySpec;
//    typedef RankDictionary<TRankDictionarySpec>         TRankDictionary;
//
//    TRankSupport rs;
//
//    std::cout << "sizeof(Block): " << sizeof(rs.block) << std::endl;
//
//    std::cout << "bits(Block): " << BitsPerValue<TRankSupport::TBlock>::VALUE << std::endl;
//    std::cout << "length(Block): " << length(rs.block) << std::endl;
//    std::cout << "capacity(Block): " << capacity(rs.block) << std::endl;
//    std::cout << std::endl;
//
//    std::cout << "sizeof(SuperBlock): " << sizeof(rs.sblock) << std::endl;
//    std::cout << "bits(SuperBlock): " << BitsPerValue<TRankSupport::TSuperBlock>::VALUE << std::endl;
//    std::cout << "length(SuperBlock): " << length(rs.sblock) << std::endl;
//    std::cout << "capacity(SuperBlock): " << capacity(rs.sblock) << std::endl;
//    std::cout << std::endl;
//
//    std::cout << "sizeof(RankSupport): " << sizeof(rs) << std::endl;
//    std::cout << "bits(RankSupport): " << BitsPerValue<TRankSupport>::VALUE << std::endl;
//    std::cout << std::endl;
//}

//SEQAN_DEFINE_TEST(test_rss_resize)
//{
//    typedef Dna                                         TAlphabet;
//    typedef Alloc<unsigned>                             TTextSpec;
//    typedef String<TAlphabet, TTextSpec>                TText;
//
//    typedef TwoLevels<TAlphabet, unsigned>              TRankDictionarySpec;
//    typedef RankDictionary<TRankDictionarySpec>         TRankDictionary;
//
////    TText text = "ACGTNACGTNACGTNACGTNA";
//    TText text = "ACGTACGTACGTACGTACGTACGTACGTACGT";
////    TText text = "ACGTACGTACGTACGTACGTACGTACGTACGTCCCCCCCCCCCCCCC";
//
//    TRankDictionary dict(text);
////    createRankDictionary(dict, text);
//
//    std::cout << "Text: " << text << std::endl;
////    std::cout << "Block: " << rs.block << std::endl;
//
//    for (unsigned i = 0; i < 10; i++)
//        for (unsigned char c = 0; c < 4; c++)
//            std::cout << "getRank(" << Dna(c) << ", " << i << "): " << getRank(dict, i, Dna(c)) << std::endl;
//
//    std::cout << std::endl;
//}

#endif  // TEST_INDEX_FM_RANK_DICTIONARY_LEVELS_H_

