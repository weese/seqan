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
// TODO(esiragusa): Move Base RankDictionary stuff to a separete header file.

// ----------------------------------------------------------------------------
// Class RankDictionary_
// ----------------------------------------------------------------------------

template <typename TSpec>
struct RankDictionary_ {};

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

// ----------------------------------------------------------------------------
// Metafunction RankDictionaryFibreSpec
// ----------------------------------------------------------------------------

template <typename TRankDictionary>
struct RankDictionaryFibreSpec
{
    typedef Alloc<> Type;
};

// ============================================================================
// Tags
// ============================================================================

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
struct Fibre<RankDictionary_<TSpec>, FibreValueString>
{
    typedef RankDictionary_<TSpec>                                          TRankDictionary_;
    typedef typename RankDictionaryFibreSpec<TRankDictionary_>::Type        TRankDictionaryFibreSpec_;

    typedef String<RankDictionaryValue_<TSpec>, TRankDictionaryFibreSpec_>  Type;
};

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Struct RankDictionaryValue_
// ----------------------------------------------------------------------------
// TODO(esiragusa): Move Base RankDictionaryValue_ to a separete header file.

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
// Class TwoLevels RankDictionary_
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct RankDictionary_<TwoLevels<TValue, TSpec> >
{
    typedef RankDictionaryValue_<TwoLevels<TValue, TSpec> >         TRankDictionaryValue;
    typedef typename Fibre<RankDictionary_, FibreValueString>::Type TRankDictionaryFibre;

    TRankDictionaryFibre ranks;

    RankDictionary_() {};

    template <typename TText>
    RankDictionary_(TText const & text)
    {
        createRankDictionary(*this, text);
    };
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _clear()                                     [RankDictionaryValue_]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline void _clear(RankDictionaryValue_<TwoLevels<TValue, TSpec> > & entry)
{
    clear(entry.block);
    clear(entry.bits);
}

template <typename TSpec>
inline void _clear(RankDictionaryValue_<TwoLevels<bool, TSpec> > & entry)
{
    entry.block = 0;
    clear(entry.bits);
}

// ----------------------------------------------------------------------------
// Function _assignBits()                                [RankDictionaryValue_]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TTextIterator>
inline void _assignBits(RankDictionaryValue_<TwoLevels<TValue, TSpec> > & rank,
                        TTextIterator const & blockBegin,
                        TTextIterator const & blockEnd)
{
    // Assign the text character by character.
    for (TTextIterator blockIt = blockBegin; blockIt != blockEnd; ++blockIt)
        assignValue(rank.bits, blockIt - blockBegin, value(blockIt));
}

// ----------------------------------------------------------------------------
// Function _updateBlockRank()                           [RankDictionaryValue_]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TTextIterator>
inline void _updateBlockRank(RankDictionaryValue_<TwoLevels<TValue, TSpec> > & rank,
                             TTextIterator const & blockBegin,
                             TTextIterator const & blockEnd)
{
    // Update in place the sum character by character.
    for (TTextIterator blockIt = blockBegin; blockIt != blockEnd; ++blockIt)
        rank.block[ordValue(value(blockIt))]++;
}

template <typename TSpec, typename TTextIterator>
inline void _updateBlockRank(RankDictionaryValue_<TwoLevels<bool, TSpec> > & rank,
                             TTextIterator const & blockBegin,
                             TTextIterator const & blockEnd)
{
    // Update in place the sum character by character.
    for (TTextIterator blockIt = blockBegin; blockIt != blockEnd; ++blockIt)
        if (value(blockIt)) rank.block++;
}

// ----------------------------------------------------------------------------
// Function bitsAt()                                           [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TPos>
inline typename RankDictionaryBits_<TwoLevels<TValue, TSpec> >::Type const &
bitsAt(RankDictionary_<TwoLevels<TValue, TSpec> > const & dict, TPos pos)
{
    return dict.ranks[pos].bits;
}

// ----------------------------------------------------------------------------
// Function blockAt()                                          [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TPos>
inline typename RankDictionaryBlock_<TwoLevels<TValue, TSpec> >::Type const &
blockAt(RankDictionary_<TwoLevels<TValue, TSpec> > const & dict, TPos pos)
{
    return dict.ranks[pos].block;
}

// ----------------------------------------------------------------------------
// Function createRankDictionary()                             [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TSize, typename TExpand>
inline void
reserve(RankDictionary_<TwoLevels<TValue, TSpec> > & dict, TSize size, Tag<TExpand> const tag)
{
    reserve(dict.ranks, size / BlockSize<TValue>::VALUE, tag);
}

// ----------------------------------------------------------------------------
// Function createRankDictionary()                             [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TText>
inline void
createRankDictionary(RankDictionary_<TwoLevels<TValue, TSpec> > & dict, TText const & text)
{
    typedef TwoLevels<TValue, TSpec>                            TRankDictionarySpec;
    typedef RankDictionaryValue_<TRankDictionarySpec>           TRankDictionaryValue;
    typedef typename Iterator<TText const, Standard>::Type      TTextIterator;

    // Reserve space in the RankDictionary.
    reserve(dict, length(text), Exact());

    // Get an empty RankSupport entry.
    TRankDictionaryValue rank;
    _clear(rank);

    // Last bits might be smaller than BlockSize.
    TTextIterator textEnd = end(text, Standard());
    TTextIterator lastBlockBegin = textEnd - length(text) % BlockSize<TValue>::VALUE;
    TTextIterator blockBegin = begin(text, Standard());
    TTextIterator blockEnd = blockBegin + BlockSize<TValue>::VALUE;

    // Scan the text blockwise.
    while (blockBegin != lastBlockBegin)
    {
        // Assign the text bits to the RankSupport entry.
        _assignBits(rank, blockBegin, blockEnd);

        // Append entry to the RankSupport String.
        appendValue(dict.ranks, rank);

        // Update the ranks.
        _updateBlockRank(rank, blockBegin, blockEnd);

        blockBegin = blockEnd;
        blockEnd += BlockSize<TValue>::VALUE;
    }

    // Scan last text bits.
    if (blockBegin != textEnd)
    {
        // Assign the text bits to the RankSupport entry.
        _assignBits(rank, blockBegin, textEnd);

        // Append entry to the RankSupport String.
        appendValue(dict.ranks, rank);
    }
}

// ----------------------------------------------------------------------------
// Function _getBlockRank()                                    [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TPos, typename TChar>
inline typename Size<RankDictionary_<TwoLevels<TValue, TSpec> > const>::Type
_getBlockRank(RankDictionary_<TwoLevels<TValue, TSpec> > const & dict, TPos blockPos, TChar c)
{
    return blockAt(dict, blockPos)[ordValue(c)];
}

template <typename TSpec, typename TPos>
inline typename Size<RankDictionary_<TwoLevels<bool, TSpec> > const>::Type
_getBlockRank(RankDictionary_<TwoLevels<bool, TSpec> > const & dict, TPos blockPos, bool c)
{
    // Not blockPos but starting position of the block.
    return c ? blockAt(dict, blockPos) : blockPos - blockAt(dict, blockPos);
}

// ----------------------------------------------------------------------------
// Function _getBitsRank()                                     [RankDictionary]
// ----------------------------------------------------------------------------
// NOTE(esiragusa): This version is generic but absymally slow.

template <typename TValue, typename TSpec, typename TPos>
inline typename Size<RankDictionary_<TwoLevels<TValue, TSpec> > const>::Type
_getBitsRank(RankDictionary_<TwoLevels<TValue, TSpec> > const & dict, TPos blockPos, TPos bitsPos, TValue c)
{
    typedef TwoLevels<Dna, TSpec>                                   TRankDictionarySpec;
    typedef RankDictionary_<TRankDictionarySpec>                    TRankDictionary;
    typedef typename Size<TRankDictionary>::Type                    TSize;

    TSize bitsRank = 0;

    for (TSize i = 0; i < bitsPos; ++i)
        bitsRank += isEqual(bitsAt(dict, blockPos)[i], c);

    return bitsRank;
}

// ----------------------------------------------------------------------------
// Function _getBitsRank()                                     [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TPos>
inline typename Size<RankDictionary_<TwoLevels<Dna, TSpec> > const>::Type
_getBitsRank(RankDictionary_<TwoLevels<Dna, TSpec> > const & dict, TPos blockPos, TPos bitsPos, Dna c)
{
    typedef TwoLevels<Dna, TSpec>                                   TRankDictionarySpec;
    typedef RankDictionary_<TRankDictionarySpec>                    TRankDictionary;
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
// Function _getBitsRank()                                     [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TPos>
inline typename Size<RankDictionary_<TwoLevels<bool, TSpec> > const>::Type
_getBitsRank(RankDictionary_<TwoLevels<bool, TSpec> > const & dict, TPos blockPos, TPos bitsPos, bool c)
{
    typedef TwoLevels<bool, TSpec>                                  TRankDictionarySpec;
    typedef RankDictionary_<TRankDictionarySpec>                    TRankDictionary;
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
// Function getRank()                                          [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TPos, typename TChar>
inline typename Size<RankDictionary_<TwoLevels<TValue, TSpec> > const>::Type
getRank(RankDictionary_<TwoLevels<TValue, TSpec> > const & dict, TPos pos, TChar c)
{
    typedef TwoLevels<TValue, TSpec>                                TRankDictionarySpec;
    typedef RankDictionary_<TRankDictionarySpec>                    TRankDictionary;
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
inline TValue getValue(RankDictionary_<TwoLevels<TValue, TSpec> > const & dict, TPos pos)
{
    typedef TwoLevels<TValue, TSpec>                                TRankDictionarySpec;
    typedef RankDictionary_<TRankDictionarySpec>                    TRankDictionary;
    typedef typename Size<TRankDictionary>::Type                    TSize;

    // TODO(esiragusa): Use bit shifts to derive positions.
    TSize blockPos = pos / BlockSize<TValue>::VALUE;
    TSize bitsPos = pos % BlockSize<TValue>::VALUE;

    return bitsAt(dict, blockPos)[bitsPos];
}

}


#endif  // INDEX_FM_RANK_SUPPORT_STRING_H_
