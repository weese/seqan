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
struct Size<String<TValue, Alloc<unsigned> > >
{
    typedef unsigned Type;
};

}

using namespace seqan;

SEQAN_DEFINE_TEST(test_rss_size)
{
    typedef Dna                                         TAlphabet;
    typedef Alloc<Nothing>                              TTextSpec;
    typedef String<TAlphabet, TTextSpec>                TText;
    typedef RankSupport<TText>                          TRankSupport;
    typedef String<TRankSupport>                        TRankSupportString;

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
}

SEQAN_DEFINE_TEST(test_rss_resize)
{
    typedef Dna                                         TAlphabet;
    typedef Alloc<Nothing>                              TTextSpec;
    typedef String<TAlphabet, TTextSpec>                TText;
    typedef RankSupport<TText>                          TRankSupport;
    typedef String<TRankSupport>                        TRankSupportString;

//    TText text = "ACGTNACGTNACGTNACGTNA";
//    TText text = "ACGTACGTACGTACGTACGTACGTACGTACGT";
    TText text = "ACGTACGTACGTACGTACGTACGTACGTACGTCCCCCCCCCCCCCCC";

    TRankSupportString rss;
    fillRankSupportString(rss, text);

    std::cout << "Text: " << text << std::endl;
//    std::cout << "Block: " << rs.block << std::endl;

    for (unsigned i = 0; i < length(rss); i++)
        std::cout << rss[i].block << std::endl;

    for (unsigned i = 0; i < length(rss); i++)
        std::cout << rss[i].sblock << std::endl;

    std::cout << std::endl;

    std::cout << "getRank(A, 5): " << getRank(rss, 5u, Dna('A')) << std::endl;
    std::cout << "getRank(A, 31): " << getRank(rss, 31u, Dna('A')) << std::endl;
    std::cout << "getRank(G, 9): " << getRank(rss, 9u, Dna('G')) << std::endl;
    std::cout << "getRank(C, 40): " << getRank(rss, 40u, Dna('C')) << std::endl;

    std::cout << std::endl;
}

#endif  // TEST_INDEX_FM_RANK_SUPPORT_STRING_H_

