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

#ifndef INDEX_FM_RANK_DICTIONARY_LEVELS_H_
#define INDEX_FM_RANK_DICTIONARY_LEVELS_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ----------------------------------------------------------------------------
// Struct RankDictionaryValue_
// ----------------------------------------------------------------------------

template <typename TSpec>
struct RankDictionaryValue_;

// ----------------------------------------------------------------------------
// Metafunction RankDictionaryBlock_
// ----------------------------------------------------------------------------

template <typename TSpec>
struct RankDictionaryBlock_;

// ----------------------------------------------------------------------------
// Metafunction RankDictionaryBits_
// ----------------------------------------------------------------------------

template <typename TSpec>
struct RankDictionaryBits_;

// ----------------------------------------------------------------------------
// Metafunction RankDictionaryBitMask_
// ----------------------------------------------------------------------------

template <typename TSpec>
struct RankDictionaryBitMask_;

// ============================================================================
// Tags
// ============================================================================

// ----------------------------------------------------------------------------
// Tag FibreRanks
// ----------------------------------------------------------------------------

struct FibreRanks_;

typedef Tag<FibreRanks_>
const FibreRanks;

// ----------------------------------------------------------------------------
// Tag TwoLevels
// ----------------------------------------------------------------------------

template <typename TValue = bool, typename TSpec = void>
struct TwoLevels {};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Size                                           [RankDictionary]
// ----------------------------------------------------------------------------
 
template <typename TValue, typename TSpec>
struct Size<RankDictionary<TwoLevels<TValue, TSpec> > >
{
    typedef typename Size<String<TValue, TSpec> >::Type Type;
};

template <typename TValue, typename TSpec>
struct Size<RankDictionary<TwoLevels<TValue, TSpec> > const> :
    public Size<RankDictionary<TwoLevels<TValue, TSpec> > > {};

// ----------------------------------------------------------------------------
// Metafunction BlockSize
// ----------------------------------------------------------------------------

template <typename TValue>
struct BlockSize
{
    static const unsigned VALUE = BitsPerValue<unsigned long>::VALUE / BitsPerValue<TValue>::VALUE;
};

// ----------------------------------------------------------------------------
// Metafunction RankDictionaryBlock_                                [TwoLevels]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct RankDictionaryBlock_<TwoLevels<TValue, TSpec> >
{
    typedef RankDictionary<TwoLevels<TValue, TSpec> >               TRankDictionary_;
    typedef typename Size<TRankDictionary_>::Type                   TSize_;
    
    typedef Tuple<TSize_, ValueSize<TValue>::VALUE>                 Type;
};

template <typename TSpec>
struct RankDictionaryBlock_<TwoLevels<bool, TSpec> >
{
    typedef RankDictionary<TwoLevels<bool, TSpec> >                 TRankDictionary_;
    typedef typename Size<TRankDictionary_>::Type                   Type;
};

// ----------------------------------------------------------------------------
// Metafunction RankDictionaryBits_                                 [TwoLevels]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct RankDictionaryBits_<TwoLevels<TValue, TSpec> >
{
    typedef Tuple<TValue, BlockSize<TValue>::VALUE, BitPacked<> >   Type;
};

// ----------------------------------------------------------------------------
// Metafunction RankDictionaryBitMask_
// ----------------------------------------------------------------------------

template <>
struct RankDictionaryBitMask_<unsigned short>
{
    static const unsigned short VALUE = 0x5555;
};

template <>
struct RankDictionaryBitMask_<unsigned int>
{
    static const unsigned int VALUE = 0x55555555;
};

template <>
struct RankDictionaryBitMask_<__uint64>
{
    static const __uint64 VALUE = 0x5555555555555555;
};

// ----------------------------------------------------------------------------
// Metafunction Fibre
// ----------------------------------------------------------------------------

template <typename TSpec>
struct Fibre<RankDictionary<TSpec>, FibreRanks>
{
    typedef RankDictionary<TSpec>                                           TRankDictionary_;
    typedef typename RankDictionaryFibreSpec<TRankDictionary_>::Type        TRankDictionaryFibreSpec_;

    typedef String<RankDictionaryValue_<TSpec>, TRankDictionaryFibreSpec_>  Type;
};

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Struct RankDictionaryValue_
// ----------------------------------------------------------------------------

template <typename TSpec = TwoLevels<> >
struct RankDictionaryValue_ {};

// ----------------------------------------------------------------------------
// Struct TwoLevels RankDictionaryValue_
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct RankDictionaryValue_<TwoLevels<TValue, TSpec> >
{
    typedef TwoLevels<TValue, TSpec>    TRankDictionarySpec;

    // A summary of counts for each block of TValue symbols.
    typename RankDictionaryBlock_<TRankDictionarySpec>::Type    block;

   // A bit-compressed block of TValue symbols.
    typename RankDictionaryBits_<TRankDictionarySpec>::Type     bits;
};

// ----------------------------------------------------------------------------
// Class TwoLevels RankDictionary
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct RankDictionary<TwoLevels<TValue, TSpec> >
{
    typedef RankDictionaryValue_<TwoLevels<TValue, TSpec> >     TRankDictionaryValue;
    typedef typename Fibre<RankDictionary, FibreRanks>::Type    TRankDictionaryFibre;

    TRankDictionaryFibre ranks;

    RankDictionary() {};

    template <typename TText>
    RankDictionary(TText const & text)
    {
        createRankDictionary(*this, text);
    };
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _clearBlockAt()                                    [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TPos>
inline void _clearBlockAt(RankDictionary<TwoLevels<TValue, TSpec> > & dict, TPos pos)
{
    clear(blockAt(dict, pos));
}

// ----------------------------------------------------------------------------
// Function _clearBlockAt()                                    [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TPos>
inline void _clearBlockAt(RankDictionary<TwoLevels<bool, TSpec> > & dict, TPos pos)
{
    blockAt(dict, pos) = 0u;
}

// ----------------------------------------------------------------------------
// Function _assignBits()                                [RankDictionaryValue_]
// ----------------------------------------------------------------------------

template <typename TBits, typename TTextIterator>
inline void _assignBits(TBits & bits,
                        TTextIterator const & blockBegin,
                        TTextIterator const & blockEnd)
{
    // Assign the text character by character.
    for (TTextIterator blockIt = blockBegin; blockIt != blockEnd; ++blockIt)
        assignValue(bits, blockIt - blockBegin, value(blockIt));
}

// ----------------------------------------------------------------------------
// Function bitsAt()                                           [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TPos>
inline typename RankDictionaryBits_<TwoLevels<TValue, TSpec> >::Type &
bitsAt(RankDictionary<TwoLevels<TValue, TSpec> > & dict, TPos pos)
{
    return dict.ranks[pos].bits;
}

template <typename TValue, typename TSpec, typename TPos>
inline typename RankDictionaryBits_<TwoLevels<TValue, TSpec> >::Type const &
bitsAt(RankDictionary<TwoLevels<TValue, TSpec> > const & dict, TPos pos)
{
    return dict.ranks[pos].bits;
}

// ----------------------------------------------------------------------------
// Function blockAt()                                          [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TPos>
inline typename RankDictionaryBlock_<TwoLevels<TValue, TSpec> >::Type &
blockAt(RankDictionary<TwoLevels<TValue, TSpec> > & dict, TPos pos)
{
    return dict.ranks[pos].block;
}

template <typename TValue, typename TSpec, typename TPos>
inline typename RankDictionaryBlock_<TwoLevels<TValue, TSpec> >::Type const &
blockAt(RankDictionary<TwoLevels<TValue, TSpec> > const & dict, TPos pos)
{
    return dict.ranks[pos].block;
}

// ----------------------------------------------------------------------------
// Function reserve()                                          [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TSize, typename TExpand>
inline void
reserve(RankDictionary<TwoLevels<TValue, TSpec> > & dict, TSize size, Tag<TExpand> const tag)
{
    reserve(dict.ranks, std::ceil(size / static_cast<double>(BlockSize<TValue>::VALUE)), tag);
}

// ----------------------------------------------------------------------------
// Function resize()                                           [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TSize, typename TExpand>
inline void
resize(RankDictionary<TwoLevels<TValue, TSpec> > & dict, TSize size, Tag<TExpand> const tag)
{
    resize(dict.ranks, std::ceil(size / static_cast<double>(BlockSize<TValue>::VALUE)), tag);
}

// ----------------------------------------------------------------------------
// Function _getBlockRank()                                    [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TPos, typename TChar>
inline typename Size<RankDictionary<TwoLevels<TValue, TSpec> > const>::Type
_getBlockRank(RankDictionary<TwoLevels<TValue, TSpec> > const & dict, TPos blockPos, TChar c)
{
    return blockAt(dict, blockPos)[ordValue(c)];
}

// ----------------------------------------------------------------------------
// Function _getBlockRank(bool)                                [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TPos>
inline typename Size<RankDictionary<TwoLevels<bool, TSpec> > const>::Type
_getBlockRank(RankDictionary<TwoLevels<bool, TSpec> > const & dict, TPos blockPos, bool c)
{
    // Not blockPos but starting position of the block.
    return c ? blockAt(dict, blockPos) : blockPos - blockAt(dict, blockPos);
}

// ----------------------------------------------------------------------------
// Function _getBitsRank()                                     [RankDictionary]
// ----------------------------------------------------------------------------
// NOTE(esiragusa): This version is generic but absymally slow.

template <typename TValue, typename TSpec, typename TBlockPos, typename TBitsPos>
inline typename Size<RankDictionary<TwoLevels<TValue, TSpec> > const>::Type
_getBitsRank(RankDictionary<TwoLevels<TValue, TSpec> > const & dict, TBlockPos blockPos, TBitsPos bitsPos, TValue c)
{
    typedef TwoLevels<Dna, TSpec>                                   TRankDictionarySpec;
    typedef RankDictionary<TRankDictionarySpec>                    TRankDictionary;
    typedef typename Size<TRankDictionary>::Type                    TSize;

    TSize bitsRank = 0;

    for (TSize i = 0; i < bitsPos; ++i)
        bitsRank += isEqual(bitsAt(dict, blockPos)[i], c);

    return bitsRank;
}

// ----------------------------------------------------------------------------
// Function _getBitsRank(Dna)                                  [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TBlockPos, typename TBitsPos>
inline typename Size<RankDictionary<TwoLevels<Dna, TSpec> > const>::Type
_getBitsRank(RankDictionary<TwoLevels<Dna, TSpec> > const & dict, TBlockPos blockPos, TBitsPos bitsPos, Dna c)
{
    typedef TwoLevels<Dna, TSpec>                                   TRankDictionarySpec;
    typedef RankDictionary<TRankDictionarySpec>                    TRankDictionary;
    typedef typename Size<TRankDictionary>::Type                    TSize;
    typedef typename RankDictionaryBits_<TRankDictionarySpec>::Type TBits;
    typedef typename TBits::TBitVector                              TBitVector;

    TBits bits = bitsAt(dict, blockPos);

    // Clear the last 2 * bitsPos positions.
    TBitVector word = bits.i & ~(MaxValue<TBitVector>::VALUE >> (bitsPos << 1));

    // And matches when c == G|T.
    TBitVector odd  = ((ordValue(c) & ordValue(Dna('G'))) ? word : ~word) >> 1;

    // And matches when c == C|T.
    TBitVector even = ((ordValue(c) & ordValue(Dna('C'))) ? word : ~word);

    // Apply the interleaved mask.
    TBitVector mask = odd & even & RankDictionaryBitMask_<TBitVector>::VALUE;

    // The rank is the sum of bits on.
    TSize bitsRank = popCount(mask);

    // If c == A then masked character positions must be subtracted from the count.
    if (c == Dna('A')) bitsRank -= BlockSize<Dna>::VALUE - bitsPos;

    return bitsRank;
}

// ----------------------------------------------------------------------------
// Function _getBitsRank(bool)                                 [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TBlockPos, typename TBitsPos>
inline typename Size<RankDictionary<TwoLevels<bool, TSpec> > const>::Type
_getBitsRank(RankDictionary<TwoLevels<bool, TSpec> > const & dict, TBlockPos blockPos, TBitsPos bitsPos, bool c)
{
    typedef TwoLevels<bool, TSpec>                                  TRankDictionarySpec;
    typedef RankDictionary<TRankDictionarySpec>                    TRankDictionary;
    typedef typename Size<TRankDictionary>::Type                    TSize;
    typedef typename RankDictionaryBits_<TRankDictionarySpec>::Type TBits;
    typedef typename TBits::TBitVector                              TBitVector;

    TBits bits = bitsAt(dict, blockPos);

    // Clear the last bitsPos positions.
    TBitVector word = bits.i & ~(MaxValue<TBitVector>::VALUE >> bitsPos);

    // Get the sum of the bits on.
    TSize bitsRank = popCount(word);

    // Return either the rank for c == true or its complement for c == false.
    return c ? bitsRank : bitsPos - bitsRank;
}

// ----------------------------------------------------------------------------
// Function _getBitsRank()                                     [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TPos>
inline typename Size<RankDictionary<TwoLevels<TValue, TSpec> > const>::Type
_getBitsRank(RankDictionary<TwoLevels<TValue, TSpec> > const & dict, TPos blockPos, TValue c)
{
    typedef TwoLevels<TValue, TSpec>                                TRankDictionarySpec;
    typedef typename RankDictionaryBits_<TRankDictionarySpec>::Type TBits;
    typedef typename TBits::TBitVector                              TBitVector;

    return _getBitsRank(dict, blockPos, MaxValue<TBitVector>::VALUE, c);
}

// ----------------------------------------------------------------------------
// Function _getBitsRanks()                                    [RankDictionary]
// ----------------------------------------------------------------------------
// TODO(esiragusa): Specialize _getBitsRanks() for Dna.

template <typename TValue, typename TSpec, typename TBlockPos, typename TBitsPos>
inline typename RankDictionaryBlock_<TwoLevels<TValue, TSpec> >::Type
_getBitsRanks(RankDictionary<TwoLevels<TValue, TSpec> > const & dict, TBlockPos blockPos, TBitsPos bitsPos)
{
    typedef TwoLevels<TValue, TSpec>                                    TRankDictionarySpec;
    typedef typename RankDictionaryBlock_<TRankDictionarySpec>::Type    TBlock;
    typedef typename ValueSize<TValue>::Type                            TValueSize;

    TBlock blockRank;

    for (TValueSize c = 0; c < ValueSize<TValue>::VALUE; ++c)
        assignValue(blockRank, c, _getBitsRank(dict, blockPos, bitsPos, TValue(c)));

    return blockRank;
}

// ----------------------------------------------------------------------------
// Function _getBitsRanks(bool)                                [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TBlockPos, typename TBitsPos>
inline typename RankDictionaryBlock_<TwoLevels<bool, TSpec> >::Type
_getBitsRanks(RankDictionary<TwoLevels<bool, TSpec> > const & dict, TBlockPos blockPos, TBitsPos bitsPos)
{
    return _getBitsRank(dict, blockPos, bitsPos, true);
}

// ----------------------------------------------------------------------------
// Function _getBitsRanks()                                    [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TPos>
inline typename RankDictionaryBlock_<TwoLevels<TValue, TSpec> >::Type
_getBitsRanks(RankDictionary<TwoLevels<TValue, TSpec> > const & dict, TPos blockPos)
{
    typedef TwoLevels<TValue, TSpec>                                TRankDictionarySpec;
    typedef typename RankDictionaryBits_<TRankDictionarySpec>::Type TBits;
    typedef typename TBits::TBitVector                              TBitVector;

    return _getBitsRanks(dict, blockPos, MaxValue<TBitVector>::VALUE);
}

// ----------------------------------------------------------------------------
// Function getRank()                                          [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TPos, typename TChar>
inline typename Size<RankDictionary<TwoLevels<TValue, TSpec> > const>::Type
getRank(RankDictionary<TwoLevels<TValue, TSpec> > const & dict, TPos pos, TChar c)
{
    typedef TwoLevels<TValue, TSpec>                                TRankDictionarySpec;
    typedef RankDictionary<TRankDictionarySpec>                    TRankDictionary;
    typedef typename Size<TRankDictionary>::Type                    TSize;

    // TODO(esiragusa): Use bit shifts to derive positions.
    TSize blockPos = pos / BlockSize<TValue>::VALUE;
    TSize bitsPos = pos % BlockSize<TValue>::VALUE;

    return _getBlockRank(dict, blockPos, convert<TValue>(c)) +
           _getBitsRank(dict, blockPos, bitsPos, convert<TValue>(c));
}

// ----------------------------------------------------------------------------
// Function getValue()                                         [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TPos>
inline TValue getValue(RankDictionary<TwoLevels<TValue, TSpec> > const & dict, TPos pos)
{
    typedef TwoLevels<TValue, TSpec>                                TRankDictionarySpec;
    typedef RankDictionary<TRankDictionarySpec>                    TRankDictionary;
    typedef typename Size<TRankDictionary>::Type                    TSize;

    // TODO(esiragusa): Use bit shifts to derive positions.
    TSize blockPos = pos / BlockSize<TValue>::VALUE;
    TSize bitsPos = pos % BlockSize<TValue>::VALUE;

    return bitsAt(dict, blockPos)[bitsPos];
}

// ----------------------------------------------------------------------------
// Function setValue()                                         [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TPos, typename TChar>
inline void setValue(RankDictionary<TwoLevels<TValue, TSpec> > & dict, TPos pos, TChar c)
{
    typedef TwoLevels<TValue, TSpec>                                TRankDictionarySpec;
    typedef RankDictionary<TRankDictionarySpec>                    TRankDictionary;
    typedef typename Size<TRankDictionary>::Type                    TSize;

    // TODO(esiragusa): Use bit shifts to derive positions.
    TSize blockPos = pos / BlockSize<TValue>::VALUE;
    TSize bitsPos = pos % BlockSize<TValue>::VALUE;

    bitsAt(dict, blockPos)[bitsPos] = convert<TValue>(c);
}

// ----------------------------------------------------------------------------
// Function appendValue()                                      [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TChar, typename TExpand>
inline void appendValue(RankDictionary<TwoLevels<TValue, TSpec> > & dict, TChar c, Tag<TExpand> const tag)
{
    resize(dict, length(dict) + 1, tag);
    setValue(dict, length(dict), c);
}

// ----------------------------------------------------------------------------
// Function updateRanks()                                      [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TPos>
inline void updateRanks(RankDictionary<TwoLevels<TValue, TSpec> > & dict, TPos pos)
{
    typedef TwoLevels<TValue, TSpec>                                TRankDictionarySpec;
    typedef RankDictionary<TRankDictionarySpec>                    TRankDictionary;
    typedef typename Size<TRankDictionary>::Type                    TSize;
    typedef typename Fibre<TRankDictionary, FibreRanks>::Type       TFibreRanks;
    typedef typename Iterator<TFibreRanks, Standard>::Type          TFibreRanksIter;

    TFibreRanksIter ranksBegin = begin(dict.ranks, Standard());
    TFibreRanksIter ranksEnd = end(dict.ranks, Standard());

    _clearBlockAt(dict, 0u);

    for (TFibreRanksIter ranksIt = ranksBegin + pos / BlockSize<TValue>::VALUE; ranksIt != ranksEnd - 1; ++ranksIt)
    {
        TSize blockPos = ranksIt - ranksBegin;

        blockAt(dict, blockPos + 1) = blockAt(dict, blockPos) + _getBitsRanks(dict, blockPos);
    }
}

// ----------------------------------------------------------------------------
// Function updateRanks()                                      [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline void updateRanks(RankDictionary<TwoLevels<TValue, TSpec> > & dict)
{
    updateRanks(dict, 0u);
}

// ----------------------------------------------------------------------------
// Function createRankDictionary()                             [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TText>
inline void
createRankDictionary(RankDictionary<TwoLevels<TValue, TSpec> > & dict, TText const & text)
{
    typedef TwoLevels<TValue, TSpec>                                TRankDictionarySpec;
    typedef RankDictionary<TRankDictionarySpec>                    TRankDictionary;
    typedef typename Size<TRankDictionary>::Type                    TSize;
    typedef typename Iterator<TText const, Standard>::Type          TTextIterator;

    // Resize the RankDictionary.
    resize(dict, length(text), Exact());

    // Last bits might be smaller than BlockSize.
    TTextIterator textEnd = end(text, Standard());
    TTextIterator lastBlockBegin = textEnd - length(text) % BlockSize<TValue>::VALUE;
    TTextIterator blockBegin = begin(text, Standard());
    TTextIterator blockEnd = blockBegin + BlockSize<TValue>::VALUE;

    TSize blockPos = 0;

    // Scan the text blockwise.
    while (blockBegin != lastBlockBegin)
    {
        // Assign the text bits to the RankSupport entry.
        _assignBits(bitsAt(dict, blockPos), blockBegin, blockEnd);

        blockBegin = blockEnd;
        blockEnd += BlockSize<TValue>::VALUE;
        blockPos++;
    }

    // Assign the last text bits to the RankSupport entry.
    if (blockBegin != textEnd)
        _assignBits(bitsAt(dict, blockPos), blockBegin, blockEnd);

    // Update all ranks.
    updateRanks(dict);
}

}

#endif  // INDEX_FM_RANK_DICTIONARY_LEVELS_H_
