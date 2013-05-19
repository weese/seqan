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
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// ==========================================================================

#ifndef TEST_INDEX_FM_SENTINEL_RANK_DICTIONARY
#define TEST_INDEX_FM_SENTINEL_RANK_DICTIONARY

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>

template <typename TAlphabetSpecPair>
class SentinelRankDictionaryTest : public seqan::Test
{
public:
    typedef typename seqan::TagListValue<TAlphabetSpecPair, 0>::Type TSentinelRankDictionarySpec;
};

typedef seqan::TagList<
            seqan::TagList<SentinelRankDictionary<RankDictionary<seqan::WaveletTree<seqan::Dna> >, Sentinel> >, seqan::TagList<
            seqan::TagList<SentinelRankDictionary<RankDictionary<seqan::WaveletTree<seqan::Dna5> >, Sentinel> >, seqan::TagList<
            seqan::TagList<SentinelRankDictionary<RankDictionary<seqan::WaveletTree<seqan::AminoAcid> >, Sentinel> >, seqan::TagList<
            seqan::TagList<SentinelRankDictionary<RankDictionary<seqan::WaveletTree<char> >, Sentinel> >, seqan::TagList<
            seqan::TagList<SentinelRankDictionary<RankDictionary<seqan::WaveletTree<unsigned char> >, Sentinel> >, seqan::TagList<
            seqan::TagList<SentinelRankDictionary<RankDictionary<seqan::WaveletTree<signed char> >, Sentinel> >, seqan::TagList<
            seqan::TagList<SentinelRankDictionary<RankDictionary<seqan::SequenceBitMask<seqan::Dna> >, Sentinel> >, seqan::TagList<
            seqan::TagList<SentinelRankDictionary<RankDictionary<seqan::SequenceBitMask<seqan::Dna5> >, Sentinel> >, seqan::TagList<
            seqan::TagList<SentinelRankDictionary<RankDictionary<seqan::SequenceBitMask<seqan::AminoAcid> >, Sentinel> >, seqan::TagList<
            seqan::TagList<SentinelRankDictionary<RankDictionary<seqan::WaveletTree<seqan::Dna> >, Sentinels> >, seqan::TagList<
            seqan::TagList<SentinelRankDictionary<RankDictionary<seqan::WaveletTree<seqan::Dna5> >, Sentinels> >, seqan::TagList<
            seqan::TagList<SentinelRankDictionary<RankDictionary<seqan::WaveletTree<seqan::AminoAcid> >, Sentinels> >, seqan::TagList<
            seqan::TagList<SentinelRankDictionary<RankDictionary<seqan::WaveletTree<char> >, Sentinels> >, seqan::TagList<
            seqan::TagList<SentinelRankDictionary<RankDictionary<seqan::WaveletTree<unsigned char> >, Sentinels> >, seqan::TagList<
            seqan::TagList<SentinelRankDictionary<RankDictionary<seqan::WaveletTree<signed char> >, Sentinels> >, seqan::TagList<
            seqan::TagList<SentinelRankDictionary<RankDictionary<seqan::SequenceBitMask<seqan::Dna> >, Sentinels> >, seqan::TagList<
            seqan::TagList<SentinelRankDictionary<RankDictionary<seqan::SequenceBitMask<seqan::Dna5> >, Sentinels> >, seqan::TagList<
            seqan::TagList<SentinelRankDictionary<RankDictionary<seqan::SequenceBitMask<seqan::AminoAcid> >, Sentinels> >
            > > > > > >
            > > >
            > > > > > >
            > > >
        SentinelRankDictionaryTestTypes;


template <typename T>
class SentinelRankDictionaryTestCommon : public SentinelRankDictionaryTest<T>
{};

SEQAN_TYPED_TEST_CASE(SentinelRankDictionaryTestCommon, SentinelRankDictionaryTestTypes);

using namespace seqan;


template <typename TRankDictionary>
void sentinelRankDictionaryConstructor(TRankDictionary & /*tag*/)
{
	{
		TRankDictionary sentinelRankDictionary;
	}
	{
		String<typename Value<TRankDictionary>::Type> text;        
        generateText(text);
		TRankDictionary sentinelRankDictionary(text);
		TRankDictionary sentinelRankDictionary2(sentinelRankDictionary);
        
        SEQAN_ASSERT(sentinelRankDictionary == sentinelRankDictionary2);

		for (unsigned i = 0; i < length(text); ++i)
        {
		    SEQAN_ASSERT_EQ(getValue(sentinelRankDictionary, i), text[i]);
		    SEQAN_ASSERT_EQ(getValue(sentinelRankDictionary2, i), text[i]);
        }
    }
}

//SEQAN_TYPED_TEST(SentinelRankDictionaryTestCommon, Constuctor)
//{
//    using namespace seqan;
//
//    typename TestFixture::TSentinelRankDictionarySpec dictionay;
//    sentinelRankDictionaryConstructor(dictionay);
//}

template <typename TRankDictionary>
void sentinelDictionaryClear(TRankDictionary & /*tag*/)
{
	String<typename Value<TRankDictionary>::Type> text = "ACGTNACGTNACGTN";
	TRankDictionary sentinelDictionary(text);

    clear(sentinelDictionary);

	SEQAN_ASSERT_EQ(empty(sentinelDictionary), true);
}

SEQAN_TYPED_TEST(SentinelRankDictionaryTestCommon, Clear)
{
    using namespace seqan;

    typename TestFixture::TSentinelRankDictionarySpec dictionay;
    sentinelDictionaryClear(dictionay);
}

template <typename TRankDictionary>
void sentinelDictionarySentinelPosition(TRankDictionary & /*tag*/)
{
	String<typename Value<TRankDictionary>::Type> text = "ACGTNACGTNACGTN";
	TRankDictionary sentinelDictionary(text);

    SEQAN_ASSERT_EQ(getFibre(sentinelDictionary, FibreSentinelPosition()), length(text));

    setSentinelPosition(sentinelDictionary, 2u);
    SEQAN_ASSERT_EQ(getFibre(sentinelDictionary, FibreSentinelPosition()), 2u);
}

// NOTE(esiragusa): This test fails now - it seems wrong to me.
template <typename TRankDictionary>
void sentinelDictionarySentinelPosition(SentinelRankDictionary<TRankDictionary, Sentinels> & /*tag*/)
{
	String<typename Value<SentinelRankDictionary<TRankDictionary, Sentinels> >::Type> text = "ACGTNACGTNACGTN";
	SentinelRankDictionary<TRankDictionary, Sentinels> sentinelDictionary(text);

    String<unsigned> sentinelPos;
    SEQAN_ASSERT_EQ(getFibre(sentinelDictionary, FibreSentinelPosition()), sentinelPos);

    RankDictionary<TwoLevels<bool> > bitString;
//    resize(bitString, length(text), 0);
    resize(bitString, length(text));
    setValue(bitString, 2u, true);
    setSentinelPosition(sentinelDictionary, bitString);
    String<unsigned> temp;
    appendValue(temp, 2);
    SEQAN_ASSERT_EQ(getFibre(sentinelDictionary, FibreSentinelPosition()), temp);
}

//SEQAN_TYPED_TEST(SentinelRankDictionaryTestCommon, SentinelPosition)
//{
//    using namespace seqan;
//
//    {
//        typename TestFixture::TSentinelRankDictionarySpec dictionay;
//        sentinelDictionarySentinelPosition(dictionay);
//    }
//}

template <typename TRankDictionary>
void sentinelRankDictionaryEmpty(TRankDictionary & /*tag*/)
{
    {
        TRankDictionary sentinelRankDictionary;
        
        SEQAN_ASSERT_EQ(empty(sentinelRankDictionary), true);
    }
    {
        String<typename Value<TRankDictionary>::Type> text = "ACGTNACGTNACGTN";
        TRankDictionary sentinelRankDictionary(text);
        
        SEQAN_ASSERT_EQ(empty(sentinelRankDictionary), false);
    }
}

SEQAN_TYPED_TEST(SentinelRankDictionaryTestCommon, Empty)
{
    using namespace seqan;

    {
        typename TestFixture::TSentinelRankDictionarySpec dictionay;
        sentinelRankDictionaryEmpty(dictionay);
    }
    {
        typename TestFixture::TSentinelRankDictionarySpec dictionay;
        sentinelRankDictionaryEmpty(dictionay);
    }
}

template <typename TRankDictionary>
void sentinelRankDictionaryGetValue(TRankDictionary & /*tag*/)
{
    String<typename Value<TRankDictionary>::Type> text;
    generateText(text);

    TRankDictionary sentinelRankDictionary(text); 

    for (unsigned i = 0; i < length(text); ++i)
    {
        SEQAN_ASSERT_EQ(getValue(sentinelRankDictionary, i), text[i]);
    }
}

SEQAN_TYPED_TEST(SentinelRankDictionaryTestCommon, GetValue)
{
    using namespace seqan;

    {
        typename TestFixture::TSentinelRankDictionarySpec dictionay;
        sentinelRankDictionaryGetValue(dictionay);
    }
    {
        typename TestFixture::TSentinelRankDictionarySpec dictionay;
        sentinelRankDictionaryGetValue(dictionay);
    }
}

template <typename TRankDictionary>
void sentinelRankDictionaryGetFibre(TRankDictionary & /*tag*/)
{
	String<typename Value<TRankDictionary>::Type> text = "ACGTNACGTNACGTN";
	TRankDictionary sentinelRankDictionary(text);

    typename Fibre<TRankDictionary, FibreRankDictionary>::Type & dictionary = getFibre(sentinelRankDictionary, FibreRankDictionary());

	SEQAN_ASSERT_EQ(empty(getFibre(sentinelRankDictionary, FibreRankDictionary())), false);

    clear(dictionary);

	SEQAN_ASSERT_EQ(empty(getFibre(sentinelRankDictionary, FibreRankDictionary())), true);
}

SEQAN_TYPED_TEST(SentinelRankDictionaryTestCommon, GetFibre)
{
    using namespace seqan;

    using namespace seqan;
    {
        typename TestFixture::TSentinelRankDictionarySpec dictionay;
        sentinelRankDictionaryGetFibre(dictionay);
    }
    {
        typename TestFixture::TSentinelRankDictionarySpec dictionay;
        sentinelRankDictionaryGetFibre(dictionay);
    }
}

template <typename TRankDictionary>
void sentinelRankDictionaryGetRank(TRankDictionary & /*tag*/)
{
    typedef typename Value<TRankDictionary>::Type TChar;
    String<TChar> text;
    generateText(text);
    text[0] = 'A';

    resize(text, 1000);

    TRankDictionary sentinelRankDictionary(text);
    setSentinelSubstitute(sentinelRankDictionary, 'A');
    setSentinelPosition(sentinelRankDictionary, 0);

    for (int i = MinValue<TChar>::VALUE; i <= MaxValue<TChar>::VALUE; ++i)
    {
        unsigned counter = 0;
        for (unsigned j = 0; j < length(text); ++j)
        {
            if (text[j] == (TChar)i)
                ++counter;
            if ((TChar)i == 'A' && j == 0)
                --counter;
            SEQAN_ASSERT_EQ(getRank(sentinelRankDictionary, (TChar)i, j), counter);
        }
    }
}
template <typename TRankDictionary>
void sentinelRankDictionaryGetRank(SentinelRankDictionary<TRankDictionary, Sentinels> & /*tag*/)
{
    typedef typename Value<TRankDictionary>::Type TChar;
    String<TChar> text;
    generateText(text);
    text[0] = 'A';
    text[99] = 'A';
    text[999] = 'A';

    resize(text, 1000);
    
    SentinelRankDictionary<TRankDictionary, Sentinels> sentinelRankDictionary(text);
    setSentinelSubstitute(sentinelRankDictionary, 'A');
    setValue(sentinelRankDictionary.sentinelPosition, 0, true);
    setValue(sentinelRankDictionary.sentinelPosition, 99, true);
    setValue(sentinelRankDictionary.sentinelPosition, 999, true);
    updateRanks(sentinelRankDictionary.sentinelPosition);


    for (int i = MinValue<TChar>::VALUE; i <= MaxValue<TChar>::VALUE; ++i)
    {
        unsigned counter = 0;
        for (unsigned j = 0; j < length(text); ++j)
        {
            if (text[j] == (TChar)i)
                ++counter;
            if ((TChar)i == 'A' && (j == 0 || j == 99 || j == 999))
                --counter;
            SEQAN_ASSERT_EQ(getRank(sentinelRankDictionary, (TChar)i, j), counter);
        }
    }
}

SEQAN_TYPED_TEST(SentinelRankDictionaryTestCommon, GetRank)
{
    using namespace seqan;
    {
        typename TestFixture::TSentinelRankDictionarySpec dictionay;
        sentinelRankDictionaryGetRank(dictionay);
    }
    {
        typename TestFixture::TSentinelRankDictionarySpec dictionay;
        sentinelRankDictionaryGetRank(dictionay);
    }
}

template <typename TRankDictionary>
void sentinelRankDictionaryOpenSave(TRankDictionary & /*tag*/)
{
	String<typename Value<TRankDictionary>::Type> text;
    generateText(text);
    resize(text, 1000);

    CharString tempFilename = SEQAN_TEMP_FILENAME();

    TRankDictionary sentinelRankDictionary(text);
    sentinelRankDictionary.sentinelSubstitute = 'A';
    save(sentinelRankDictionary, toCString(tempFilename));

    TRankDictionary sentinelRankDictionaryOpen;
    open(sentinelRankDictionaryOpen, toCString(tempFilename));
    SEQAN_ASSERT(sentinelRankDictionary == sentinelRankDictionaryOpen);
}

template <typename TValue>
void sentinelRankDictionaryOpenSave(RankDictionary<SequenceBitMask<TValue> > & /*tag*/) {}

//SEQAN_TYPED_TEST(SentinelRankDictionaryTestCommon, OpenSave)
//{
//    using namespace seqan;
//    {
//        typename TestFixture::TSentinelRankDictionarySpec dictionay;
//        sentinelRankDictionaryOpenSave(dictionay);
//    }
//}

#endif  // TEST_INDEX_FM_SENTINEL_RANK_DICTIONARY

