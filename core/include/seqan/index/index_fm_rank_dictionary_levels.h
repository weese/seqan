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
    typedef TwoLevels<TValue, TSpec>                            TRankDictionarySpec;
    typedef RankDictionaryValue_<TRankDictionarySpec>           TRankDictionaryValue;
    typedef typename Fibre<RankDictionary, FibreRanks>::Type    TRankDictionaryFibre;
    typedef typename Size<RankDictionary>::Type                 TSize;

    TRankDictionaryFibre ranks;
//    TSize _length;

    // ------------------------------------------------------------------------
    // Constructors
    // ------------------------------------------------------------------------

    RankDictionary() /* : _length(0) */ {};

    template <typename TText>
    RankDictionary(TText const & text)
    {
        createRankDictionary(*this, text);
    };

    // ------------------------------------------------------------------------
    // Operator==()
    // ------------------------------------------------------------------------

    inline bool operator==(RankDictionary const & other) const
    {
        return ranks == other.ranks;// &&
//               _length == other.length;
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _toBitPos()                                        [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TPos>
inline typename Size<RankDictionary<TwoLevels<TValue, TSpec> > >::Type
_toBitPos(RankDictionary<TwoLevels<TValue, TSpec> > const & /* dict */, TPos pos)
{
    // TODO(esiragusa): Use bit shifts.
    return pos % BlockSize<TValue>::VALUE;
}

// ----------------------------------------------------------------------------
// Function _toBlockPos()                                      [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TPos>
inline typename Size<RankDictionary<TwoLevels<TValue, TSpec> > >::Type
_toBlockPos(RankDictionary<TwoLevels<TValue, TSpec> > const & /* dict */, TPos pos)
{
    // TODO(esiragusa): Use bit shifts.
    return pos / BlockSize<TValue>::VALUE;
}

// ----------------------------------------------------------------------------
// Function _toPos()                                           [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TBlockPos>
inline typename Size<RankDictionary<TwoLevels<TValue, TSpec> > >::Type
_toPos(RankDictionary<TwoLevels<TValue, TSpec> > const & /* dict */, TBlockPos blockPos)
{
    // TODO(esiragusa): Use bit shifts.
    return blockPos * BlockSize<TValue>::VALUE;
}

// ----------------------------------------------------------------------------
// Function _bitsAt()                                          [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TPos>
inline typename RankDictionaryBits_<TwoLevels<TValue, TSpec> >::Type &
_bitsAt(RankDictionary<TwoLevels<TValue, TSpec> > & dict, TPos pos)
{
    return dict.ranks[_toBlockPos(dict, pos)].bits;
}

template <typename TValue, typename TSpec, typename TPos>
inline typename RankDictionaryBits_<TwoLevels<TValue, TSpec> >::Type const &
_bitsAt(RankDictionary<TwoLevels<TValue, TSpec> > const & dict, TPos pos)
{
    return dict.ranks[_toBlockPos(dict, pos)].bits;
}

// ----------------------------------------------------------------------------
// Function _blockAt()                                         [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TPos>
inline typename RankDictionaryBlock_<TwoLevels<TValue, TSpec> >::Type &
_blockAt(RankDictionary<TwoLevels<TValue, TSpec> > & dict, TPos pos)
{
    return dict.ranks[_toBlockPos(dict, pos)].block;
}

template <typename TValue, typename TSpec, typename TPos>
inline typename RankDictionaryBlock_<TwoLevels<TValue, TSpec> >::Type const &
_blockAt(RankDictionary<TwoLevels<TValue, TSpec> > const & dict, TPos pos)
{
    return dict.ranks[_toBlockPos(dict, pos)].block;
}

// ----------------------------------------------------------------------------
// Function _clearBlockAt()                                    [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TPos>
inline void _clearBlockAt(RankDictionary<TwoLevels<TValue, TSpec> > & dict, TPos pos)
{
    clear(_blockAt(dict, pos));
}

// ----------------------------------------------------------------------------
// Function _clearBlockAt(bool)                                [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TPos>
inline void _clearBlockAt(RankDictionary<TwoLevels<bool, TSpec> > & dict, TPos pos)
{
    _blockAt(dict, pos) = 0u;
}

// ----------------------------------------------------------------------------
// Function _getBlockRank()                                    [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TPos, typename TChar>
inline typename Size<RankDictionary<TwoLevels<TValue, TSpec> > const>::Type
_getBlockRank(RankDictionary<TwoLevels<TValue, TSpec> > const & dict, TPos pos, TChar c)
{
    return _blockAt(dict, pos)[ordValue(c)];
}

// ----------------------------------------------------------------------------
// Function _getBlockRank(bool)                                [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TPos>
inline typename Size<RankDictionary<TwoLevels<bool, TSpec> > const>::Type
_getBlockRank(RankDictionary<TwoLevels<bool, TSpec> > const & dict, TPos pos, bool c)
{
    // If c == false then return the complementary rank.
    return c ? _blockAt(dict, pos) : pos - _toBitPos(dict, pos) - _blockAt(dict, pos);
}

// ----------------------------------------------------------------------------
// Function _getBitsRank()                                     [RankDictionary]
// ----------------------------------------------------------------------------
// NOTE(esiragusa): This version is generic but absymally slow.

template <typename TValue, typename TSpec, typename TPos>
inline typename Size<RankDictionary<TwoLevels<TValue, TSpec> > const>::Type
_getBitsRank(RankDictionary<TwoLevels<TValue, TSpec> > const & dict, TPos pos, TValue c)
{
    typedef TwoLevels<Dna, TSpec>                                   TRankDictionarySpec;
    typedef RankDictionary<TRankDictionarySpec>                     TRankDictionary;
    typedef typename Size<TRankDictionary>::Type                    TSize;

    TSize bitsRank = 0;

    TSize bitPos = _toBitPos(dict, pos);
    for (TSize i = 0; i <= bitPos; ++i)
        bitsRank += isEqual(_bitsAt(dict, pos)[i], c);

    return bitsRank;
}

// ----------------------------------------------------------------------------
// Function _getBitsRank(Dna)                                  [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TPos>
inline typename Size<RankDictionary<TwoLevels<Dna, TSpec> > const>::Type
_getBitsRank(RankDictionary<TwoLevels<Dna, TSpec> > const & dict, TPos pos, Dna c)
{
    typedef TwoLevels<Dna, TSpec>                                   TRankDictionarySpec;
    typedef RankDictionary<TRankDictionarySpec>                     TRankDictionary;
    typedef typename Size<TRankDictionary>::Type                    TSize;
    typedef typename RankDictionaryBits_<TRankDictionarySpec>::Type TBits;
    typedef typename TBits::TBitVector                              TBitVector;

    TSize bitPos = _toBitPos(dict, pos);
    TBits bits = _bitsAt(dict, pos);

    // Clear the last 2 * bitsPos positions.
    TBitVector word = bits.i & ~(MaxValue<TBitVector>::VALUE >> (bitPos << 1));

    // And matches when c == G|T.
    TBitVector odd  = ((ordValue(c) & ordValue(Dna('G'))) ? word : ~word) >> 1;

    // And matches when c == C|T.
    TBitVector even = ((ordValue(c) & ordValue(Dna('C'))) ? word : ~word);

    // Apply the interleaved mask.
    TBitVector mask = odd & even & RankDictionaryBitMask_<TBitVector>::VALUE;

    // The rank is the sum of bits on.
    TSize bitsRank = popCount(mask);

    // If c == A then masked character positions must be subtracted from the count.
    if (c == Dna('A')) bitsRank -= BlockSize<Dna>::VALUE - bitPos;

    return bitsRank;
}

// ----------------------------------------------------------------------------
// Function _getBitsRank(bool)                                 [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TPos>
inline typename Size<RankDictionary<TwoLevels<bool, TSpec> > const>::Type
_getBitsRank(RankDictionary<TwoLevels<bool, TSpec> > const & dict, TPos pos, bool c)
{
    typedef TwoLevels<bool, TSpec>                                  TRankDictionarySpec;
    typedef RankDictionary<TRankDictionarySpec>                     TRankDictionary;
    typedef typename Size<TRankDictionary>::Type                    TSize;
    typedef typename RankDictionaryBits_<TRankDictionarySpec>::Type TBits;
    typedef typename TBits::TBitVector                              TBitVector;

    TSize bitPos = _toBitPos(dict, pos);
    TBits bits = _bitsAt(dict, pos);

    // Clear the last bitsPos positions.
    TBitVector word = bits.i & ~(MaxValue<TBitVector>::VALUE >> bitPos);

    // Get the sum of the bits on.
    TSize bitsRank = popCount(word);

    // Return either the rank for c == true or its complement for c == false.
    return c ? bitsRank : bitPos - bitsRank;
}

// ----------------------------------------------------------------------------
// Function _getBitsRanks()                                    [RankDictionary]
// ----------------------------------------------------------------------------
// TODO(esiragusa): Specialize _getBitsRanks() for Dna.

template <typename TValue, typename TSpec, typename TPos>
inline typename RankDictionaryBlock_<TwoLevels<TValue, TSpec> >::Type
_getBitsRanks(RankDictionary<TwoLevels<TValue, TSpec> > const & dict, TPos pos)
{
    typedef TwoLevels<TValue, TSpec>                                    TRankDictionarySpec;
    typedef typename RankDictionaryBlock_<TRankDictionarySpec>::Type    TBlock;
    typedef typename ValueSize<TValue>::Type                            TValueSize;

    TBlock blockRank;

    for (TValueSize c = 0; c < ValueSize<TValue>::VALUE; ++c)
        assignValue(blockRank, c, _getBitsRank(dict, pos, TValue(c)));

    return blockRank;
}

// ----------------------------------------------------------------------------
// Function _getBitsRanks(bool)                                [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TPos>
inline typename RankDictionaryBlock_<TwoLevels<bool, TSpec> >::Type
_getBitsRanks(RankDictionary<TwoLevels<bool, TSpec> > const & dict, TPos pos)
{
    return _getBitsRank(dict, pos, true);
}

// ----------------------------------------------------------------------------
// Function getRank()                                          [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TPos, typename TChar>
inline typename Size<RankDictionary<TwoLevels<TValue, TSpec> > const>::Type
getRank(RankDictionary<TwoLevels<TValue, TSpec> > const & dict, TPos pos, TChar c)
{
//    SEQAN_ASSERT_LT(pos, length(dict));
    return _getBlockRank(dict, pos, convert<TValue>(c)) +
           _getBitsRank(dict, pos, convert<TValue>(c));
}

// ----------------------------------------------------------------------------
// Function getRank(bool)                                      [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TPos>
inline typename Size<RankDictionary<TwoLevels<bool, TSpec> > const>::Type
getRank(RankDictionary<TwoLevels<bool, TSpec> > const & dict, TPos pos)
{
//    SEQAN_ASSERT_LT(pos, length(dict));
    return getRank(dict, pos, true);
}

// ----------------------------------------------------------------------------
// Function getValue()                                         [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TPos>
inline TValue getValue(RankDictionary<TwoLevels<TValue, TSpec> > const & dict, TPos pos)
{
//    SEQAN_ASSERT_LT(pos, length(dict));
    return _bitsAt(dict, pos)[_toBitPos(dict, pos)];
}

// ----------------------------------------------------------------------------
// Function setValue()                                         [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TPos, typename TChar>
inline void setValue(RankDictionary<TwoLevels<TValue, TSpec> > & dict, TPos pos, TChar c)
{
//    SEQAN_ASSERT_LT(pos, length(dict));
    assignValue(_bitsAt(dict, pos), _toBitPos(dict, pos), convert<TValue>(c));
}

// ----------------------------------------------------------------------------
// Function appendValue()                                      [RankDictionary]
// ----------------------------------------------------------------------------
//
//template <typename TValue, typename TSpec, typename TChar, typename TExpand>
//inline void appendValue(RankDictionary<TwoLevels<TValue, TSpec> > & dict, TChar c, Tag<TExpand> const tag)
//{
//    resize(dict, length(dict) + 1, tag);
//    setValue(dict, length(dict) - 1, c);
//}

// ----------------------------------------------------------------------------
// Function updateRanks()                                      [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TPos>
inline void updateRanks(RankDictionary<TwoLevels<TValue, TSpec> > & dict, TPos /* pos */)
{
    typedef TwoLevels<TValue, TSpec>                                TRankDictionarySpec;
    typedef RankDictionary<TRankDictionarySpec>                     TRankDictionary;
    typedef typename Size<TRankDictionary>::Type                    TSize;
    typedef typename Fibre<TRankDictionary, FibreRanks>::Type       TFibreRanks;
    typedef typename Iterator<TFibreRanks, Standard>::Type          TFibreRanksIter;

//    SEQAN_ASSERT_LT(pos, length(dict));

    TFibreRanksIter ranksBegin = begin(dict.ranks, Standard());
    TFibreRanksIter ranksEnd = end(dict.ranks, Standard());

    // Insures the first block ranks start from zero.
    _clearBlockAt(dict, 0u);

    // Iterate through the blocks.
    for (TFibreRanksIter ranksIt = ranksBegin /* + _toBlockPos(dict, pos) */; ranksIt != ranksEnd - 1; ++ranksIt)
    {
        TSize blockPos = ranksIt - ranksBegin;
        TSize curr = _toPos(dict, blockPos);
        TSize next = _toPos(dict, blockPos + 1);

        _blockAt(dict, next) = _blockAt(dict, curr) + _getBitsRanks(dict, next - 1);

        std::cout << "Block " << blockPos + 1 << ": [" << curr << "," << next - 1 << "]" << std::endl;
        std::cout << _blockAt(dict, next) << std::endl;
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
    typedef RankDictionary<TRankDictionarySpec>                     TRankDictionary;
    typedef typename Size<TRankDictionary>::Type                    TSize;
    typedef typename Iterator<TText const, Standard>::Type          TTextIterator;

    // Resize the RankDictionary.
    resize(dict, length(text), Exact());

    // Assign the text value by value.
    TTextIterator textBegin = begin(text, Standard());
    TTextIterator textEnd = end(text, Standard());
    for (TTextIterator textIt = textBegin; textIt != textEnd; ++textIt)
        setValue(dict, textIt - textBegin, value(textIt));

    // Update all ranks.
    updateRanks(dict);
}

// ----------------------------------------------------------------------------
// Function clear()                                            [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline void clear(RankDictionary<TwoLevels<TValue, TSpec> > & dict)
{
    clear(dict.ranks);
}

// ----------------------------------------------------------------------------
// Function empty()                                            [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline bool empty(RankDictionary<TwoLevels<TValue, TSpec> > const & dict)
{
    return empty(dict.ranks);
}

// ----------------------------------------------------------------------------
// Function length()                                           [RankDictionary]
// ----------------------------------------------------------------------------

//template <typename TValue, typename TSpec>
//inline typename Size<RankDictionary<TwoLevels<TValue, TSpec> > >::Type
//length(RankDictionary<TwoLevels<TValue, TSpec> > const & dict)
//{
//    return dict._length;
//}

// ----------------------------------------------------------------------------
// Function reserve()                                          [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TSize, typename TExpand>
inline typename Size<RankDictionary<TwoLevels<TValue, TSpec> > >::Type
reserve(RankDictionary<TwoLevels<TValue, TSpec> > & dict, TSize newCapacity, Tag<TExpand> const tag)
{
   return reserve(dict.ranks, std::ceil(newCapacity / static_cast<double>(BlockSize<TValue>::VALUE)), tag);
}

// ----------------------------------------------------------------------------
// Function resize()                                           [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TSize, typename TExpand>
inline typename Size<RankDictionary<TwoLevels<TValue, TSpec> > >::Type
resize(RankDictionary<TwoLevels<TValue, TSpec> > & dict, TSize newLength, Tag<TExpand> const tag)
{
//    dict._length = newLength;
    return resize(dict.ranks, std::ceil(newLength / static_cast<double>(BlockSize<TValue>::VALUE)), tag);
}

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline bool open(RankDictionary<TwoLevels<TValue, TSpec> > & dict, const char * fileName, int openMode)
{
    // TODO(esiragusa): Open _length.
    return open(dict.ranks, fileName, openMode);
}

template <typename TValue, typename TSpec>
inline bool open(RankDictionary<TwoLevels<TValue, TSpec> > & dict, const char * fileName)
{
    return open(dict, fileName, DefaultOpenMode<RankDictionary<TwoLevels<TValue, TSpec> > >::VALUE);
}

// ----------------------------------------------------------------------------
// Function save()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline bool save(RankDictionary<TwoLevels<TValue, TSpec> > const & dict, const char * fileName, int openMode)
{
    return save(dict.ranks, fileName, openMode);
}

template <typename TValue, typename TSpec>
inline bool save(RankDictionary<TwoLevels<TValue, TSpec> > const & dict, const char * fileName)
{
    // TODO(esiragusa): Save _length.
    return save(dict, fileName, DefaultOpenMode<RankDictionary<TwoLevels<TValue, TSpec> > >::VALUE);
}

}

#endif  // INDEX_FM_RANK_DICTIONARY_LEVELS_H_
