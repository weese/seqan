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

template <typename TSpec = void>
struct TwoLevels;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction BlockSize
// ----------------------------------------------------------------------------

template <typename TValue>
struct BlockSize
{
    static const unsigned VALUE = BitsPerValue<unsigned long>::VALUE / BitsPerValue<TValue>::VALUE;
};

// ----------------------------------------------------------------------------
// Metafunction RankSupportBitMask_
// ----------------------------------------------------------------------------

template <typename TValue>
struct RankSupportBitMask_;

template <>
struct RankSupportBitMask_<unsigned short>
{
    static const unsigned short VALUE = 0x5555;
};

template <>
struct RankSupportBitMask_<unsigned int>
{
    static const unsigned int VALUE = 0x55555555;
};

template <>
struct RankSupportBitMask_<__uint64>
{
    static const __uint64 VALUE = 0x5555555555555555;
};

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Struct TwoLevelRankDictionaryValue_
// ----------------------------------------------------------------------------

template <typename TText>
struct TwoLevelRankDictionaryValue_
{
    typedef typename Value<TText>::Type                             TValue;
    typedef typename Size<TText>::Type                              TSize;
    typedef Tuple<TSize, ValueSize<TValue>::VALUE>                  TBlock;
    typedef Tuple<TValue, BlockSize<TValue>::VALUE, BitPacked<> >   TBits;

    TBlock  block;     // A summary of counts for each TText symbol.
    TBits   bits;      // A bit-compressed snippet of TText long BlockSize symbols.
};

// ----------------------------------------------------------------------------
// Class RankDictionary_
// ----------------------------------------------------------------------------

template <typename TText, typename = TwoLevels<> >
struct RankDictionary_ {};

// ----------------------------------------------------------------------------
// Class TwoLevels RankDictionary_
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec>
struct RankDictionary_<TText, TwoLevels<TSpec> >
{
    typedef TwoLevelRankDictionaryValue_<TText>         TRankDictionaryValue;
    typedef String<TRankDictionaryValue>                TRanksString;

    TRanksString ranks;

    RankDictionary_() {};

    RankDictionary_(TText const & text)
    {
        createRankDictionary(*this, text);
    };
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _clear()                                               [RankEntry_]
// ----------------------------------------------------------------------------

template <typename TText>
inline void _clear(TwoLevelRankDictionaryValue_<TText> & dict)
{
    clear(dict.block);
    clear(dict.bits);
}

// ----------------------------------------------------------------------------
// Function _assignBits()                                          [RankEntry_]
// ----------------------------------------------------------------------------

template <typename TBits, typename TTextIterator>
inline void _assignBits(TBits & bits, TTextIterator const & blockBegin, TTextIterator const & blockEnd)
{
    // Assign the text character by character.
    for (TTextIterator blockIt = blockBegin; blockIt != blockEnd; ++blockIt)
        assignValue(bits, blockIt - blockBegin, value(blockIt));
}

// ----------------------------------------------------------------------------
// Function _updateBlockRank()                                     [RankEntry_]
// ----------------------------------------------------------------------------

template <typename TBlock, typename TTextIterator>
inline void _updateBlockRank(TBlock & block, TTextIterator const & blockBegin, TTextIterator const & blockEnd)
{
    // Update in place the sum character by character.
    for (TTextIterator blockIt = blockBegin; blockIt != blockEnd; ++blockIt)
        block[ordValue(value(blockIt))]++;
}

// ----------------------------------------------------------------------------
// Function bitsAt()                                          [RankDictionary_]
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TPos>
inline typename TwoLevelRankDictionaryValue_<TText>::TBits const &
bitsAt(RankDictionary_<TText, TwoLevels<TSpec> > const & dict, TPos pos)
{
    return dict.ranks[pos].bits;
}

// ----------------------------------------------------------------------------
// Function blockAt()                                         [RankDictionary_]
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TPos>
inline typename TwoLevelRankDictionaryValue_<TText>::TBlock const &
blockAt(RankDictionary_<TText, TwoLevels<TSpec> > const & dict, TPos pos)
{
    return dict.ranks[pos].block;
}

// ----------------------------------------------------------------------------
// Function createRankDictionary()
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec>
inline void
createRankDictionary(RankDictionary_<TText, TwoLevels<TSpec> > & dict, TText const & text)
{
    typedef TwoLevelRankDictionaryValue_<TText>                     TRankEntry;
    typedef typename Value<TText>::Type                             TValue;
    typedef typename Iterator<TText const, Standard>::Type          TTextIterator;

    // Reserve space in the RankSupport String.
//    reserve(dict, length(text), Exact());

    // Get an empty RankSupport entry.
    TRankEntry entry;
    _clear(entry);

    // Last bits might be smaller than BlockSize.
    TTextIterator textEnd = end(text, Standard());
    TTextIterator lastBlockBegin = textEnd - length(text) % BlockSize<TValue>::VALUE;
    TTextIterator blockBegin = begin(text, Standard());
    TTextIterator blockEnd = blockBegin + BlockSize<TValue>::VALUE;

    // Scan the text blockwise.
    while (blockBegin != lastBlockBegin)
    {
        // Assign the text bits to the RankSupport entry.
        _assignBits(entry.bits, blockBegin, blockEnd);

        // Append entry to the RankSupport String.
        appendValue(dict.ranks, entry);

        // Update the ranks.
        _updateBlockRank(entry.block, blockBegin, blockEnd);

        blockBegin = blockEnd;
        blockEnd += BlockSize<TValue>::VALUE;
    }

    // Scan last text bits.
    if (blockBegin != textEnd)
    {
        // Assign the text bits to the RankSupport entry.
        _assignBits(entry.bits, blockBegin, textEnd);

        // Append entry to the RankSupport String.
        appendValue(dict.ranks, entry);
    }
}

// ----------------------------------------------------------------------------
// Function _getBlockRank()
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TPos, typename TChar>
inline typename Size<TText>::Type
_getBlockRank(RankDictionary_<TText, TwoLevels<TSpec> > const & dict, TPos blockPos, TChar c)
{
    return blockAt(dict, blockPos)[ordValue(c)];
}

// ----------------------------------------------------------------------------
// Function _getBitsRank()
// ----------------------------------------------------------------------------

//template <typename TText, typename TSpec, typename TPos, typename TChar>
//inline typename Size<TText>::Type
//_getBitsRank(RankDictionary_<TText, TwoLevels<TSpec> > const & dict, TPos blockPos, TPos bitsPos, TChar c)
//{
//    typename Size<TText>::Type bitsRank = 0;
//
//    // NOTE(esiragusa): This version is generic but absymally slow.
//    for (unsigned i = 0; i < bitsPos; ++i)
//        rankInBlock += isEqual(value(bitsAt(dict, blockPos), i), c);
//
//    return bitsRank;
//}

template <typename TText, typename TSpec, typename TPos, typename TChar>
inline typename Size<TText>::Type
_getBitsRank(RankDictionary_<TText, TwoLevels<TSpec> > const & dict, TPos blockPos, TPos bitsPos, TChar c)
{
    typedef TwoLevelRankDictionaryValue_<TText>         TRankEntry;
    typedef typename TRankEntry::TBits                  TBits;
    typedef typename TBits::TBitVector                  TBlockBitVector;
    typedef typename Value<TText>::Type                 TValue;
    typedef typename Size<TText>::Type                  TSize;

    TBits bits = bitsAt(dict, blockPos);

    // Clear the last blockPos positions.
    // TODO(esiragusa): Change blockPos << 1 to blockPos * BitsPerValue in the generic case.
    TBlockBitVector word = bits.i & ~(MaxValue<TBlockBitVector>::VALUE >> (bitsPos << 1));

    // And matches when c == G|T.
    TBlockBitVector odd  = ((ordValue(c) & ordValue(Dna('G'))) ? word : ~word) >> 1;

    // And matches when c == C|T.
    TBlockBitVector even = ((ordValue(c) & ordValue(Dna('C'))) ? word : ~word);

    TBlockBitVector mask = odd & even & RankSupportBitMask_<TBlockBitVector>::VALUE;

    // The rank is the sum of bits on.
    TSize bitsRank = popCount(mask);

    // If c == A then masked character positions must be subtracted from the count.
    if (c == Dna('A')) bitsRank -= BlockSize<TValue>::VALUE - bitsPos;

    return bitsRank;
}

// ----------------------------------------------------------------------------
// Function getRank()
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TPos, typename TChar>
inline typename Size<TText>::Type
getRank(RankDictionary_<TText, TwoLevels<TSpec> > const & dict, TPos pos, TChar c)
{
    typedef typename Value<TText>::Type     TValue;
    typedef typename Size<TText>::Type      TSize;

    // TODO(esiragusa): Use bit shifts to derive positions.
    TSize blockPos = pos / BlockSize<TValue>::VALUE;
    TSize bitsPos = pos % BlockSize<TValue>::VALUE;

    return _getBlockRank(dict, blockPos, c) + _getBitsRank(dict, blockPos, bitsPos, c);
}

}


#endif  // INDEX_FM_RANK_SUPPORT_STRING_H_
