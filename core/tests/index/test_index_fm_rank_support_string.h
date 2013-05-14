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

#ifndef TEST_INDEX_FM_RANK_SUPPORT_STRING_H_
#define TEST_INDEX_FM_RANK_SUPPORT_STRING_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/random.h>
#include <seqan/index.h>

#include <seqan/index/index_fm_rank_support_string.h>

namespace seqan {

template <typename TValue>
struct Size<String<TValue, Alloc<Nothing> > >
{
    typedef unsigned Type;
};

}

using namespace seqan;

SEQAN_DEFINE_TEST(test_rss_resize)
{
    typedef Dna                                         TAlphabet;
    typedef Alloc<Nothing>                              TTextSpec;
    typedef String<TAlphabet, TTextSpec>                TText;
    typedef typename Iterator<TText, Standard>::Type    TTextIterator;
    typedef RankSupport<TText>                          TRankSupport;
    typedef String<TRankSupport>                        TRankSupportString;

    std::cout << std::endl;

// ==========================================================================

    std::cout << "sizeof(text): " << sizeof(Size<TText>::Type) << std::endl;
    std::cout << "sizeof(alphabet): " << sizeof(Size<TAlphabet>::Type) << std::endl;

    std::cout << "bits(alphabet): " << static_cast<unsigned>(BitsPerValue<typename Value<TText>::Type>::VALUE) << std::endl;
    std::cout << "bits(word): " << BitsPerValue<unsigned long>::VALUE << std::endl;

    std::cout << std::endl;

// ==========================================================================

    TRankSupport rs;

    std::cout << "sizeof(Block): " << sizeof(rs.block) << std::endl;
    std::cout << "bits(Block): " << BitsPerValue<TRankSupport::TBlock>::VALUE << std::endl;
    std::cout << "length(Block): " << length(rs.block) << std::endl;
    std::cout << "capacity(Block): " << capacity(rs.block) << std::endl;
    std::cout << std::endl;

    std::cout << "sizeof(SuperBlock): " << sizeof(rs.sblock) << std::endl;
    std::cout << "bits(SuperBlock): " << BitsPerValue<TRankSupport::TSuperBlock>::VALUE << std::endl;
    std::cout << "length(SuperBlock): " << length(rs.sblock) << std::endl;
    std::cout << "capacity(SuperBlock): " << capacity(rs.sblock) << std::endl;
    std::cout << std::endl;

    std::cout << "sizeof(RankSupport): " << sizeof(rs) << std::endl;
    std::cout << "bits(RankSupport): " << BitsPerValue<TRankSupport>::VALUE << std::endl;
    std::cout << std::endl;

// ==========================================================================

    TText text = "ACGTACGTACGTACGTACGTACGTACGTACGT";
//    TText text = "ACGTNACGTNACGTNACGTNA";

    for (TTextIterator it = begin(text, Standard()); it != end(text, Standard()); ++it)
        assignValue(rs.block, position(it, text), value(it));

    std::cout << "Text: " << text << std::endl;
    std::cout << "Block: " << rs.block << std::endl;

    std::cout << std::endl;

// ==========================================================================

    clear(rs.sblock);
    for (TTextIterator it = begin(text, Standard()); it != end(text, Standard()); ++it)
        rs.sblock[ordValue(value(it))]++;

    std::cout << "SuperBlock: " << rs.sblock << std::endl;
    std::cout << std::endl;

// ==========================================================================

    TRankSupportString rss;
    reserve(rss, 1);

    std::cout << "length(RankSupportString): " << length(rss) << std::endl;
    std::cout << "capacity(RankSupportString): " << capacity(rss) << std::endl;

    appendValue(rss, rs);

    std::cout << "length(RankSupportString): " << length(rss) << std::endl;
    std::cout << "capacity(RankSupportString): " << capacity(rss) << std::endl;

}

#endif  // TEST_INDEX_FM_RANK_SUPPORT_STRING_H_

