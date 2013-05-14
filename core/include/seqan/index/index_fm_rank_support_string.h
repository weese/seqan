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

#ifndef INDEX_FM_RANK_SUPPORT_STRING_H_
#define INDEX_FM_RANK_SUPPORT_STRING_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

template <typename T>
struct BlockSize;

template <typename TValue>
struct BlockSize
{
    static const unsigned VALUE = BitsPerValue<unsigned long>::VALUE / BitsPerValue<TValue>::VALUE;
};

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Struct RankSupport
// ----------------------------------------------------------------------------

template <typename TText>
struct RankSupport
{
    typedef typename Value<TText>::Type                             TValue;
    typedef typename Size<TText>::Type                              TSize;
    typedef Tuple<TValue, BlockSize<TValue>::VALUE, BitPacked<> >   TBlock;
    typedef Tuple<TSize, ValueSize<TValue>::VALUE>                  TSuperBlock;

    TBlock      block;      // A bit-compressed snippet of TText long BlockSize symbols.
    TSuperBlock sblock;     // A summary of counts for each TText symbol.
};

// ============================================================================
// Functions
// ============================================================================

template <typename TText, typename TSpec, typename TPos, typename TChar>
inline typename Size<TText>::Type
_getRankInSuperBlock(String<RankSupport<TText>, TSpec> const & me, TPos sblockPos, TChar c)
{
    return me[sblockPos].sblock[ordValue(c)];
}

template <typename TText, typename TSpec, typename TPos, typename TChar>
inline typename Size<TText>::Type
_getRankInBlock(String<RankSupport<TText>, TSpec> const & me, TPos sblockPos, TPos blockPos, TChar c)
{
    typename Size<TText>::Type rankInBlock = 0;

    // TODO(esiragusa): Use popCount().
    for (unsigned i = 0; i < blockPos; ++i)
        rankInBlock += isEqual(me[sblockPos].block[i], c);

    return rankInBlock;
}

template <typename TText, typename TSpec, typename TPos, typename TChar>
inline typename Size<TText>::Type
getRank(String<RankSupport<TText>, TSpec> const & me, TPos pos, TChar c)
{
    typedef typename Value<TText>::Type     TValue;
    typedef typename Size<TText>::Type      TSize;

    // TODO(esiragusa): Use bit shifts to derive positions.
    TSize sblockPos = pos / BlockSize<TValue>::VALUE;
    TSize blockPos = pos % BlockSize<TValue>::VALUE;

    return _getRankInSuperBlock(me, sblockPos, c) + _getRankInBlock(me, sblockPos + 1, blockPos, c);
}

}


#endif  // INDEX_FM_RANK_SUPPORT_STRING_H_
