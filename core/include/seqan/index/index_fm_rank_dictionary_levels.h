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
// Struct RankDictionaryEntry_
// ----------------------------------------------------------------------------

template <typename TSpec>
struct RankDictionaryEntry_;

// ----------------------------------------------------------------------------
// Metafunction RankDictionaryValuesPerBlock_
// ----------------------------------------------------------------------------

template <typename TSpec>
struct RankDictionaryValuesPerBlock_;

// ----------------------------------------------------------------------------
// Metafunction RankDictionaryBlock_
// ----------------------------------------------------------------------------

template <typename TSpec>
struct RankDictionaryBlock_;

// ----------------------------------------------------------------------------
// Metafunction RankDictionaryValues_
// ----------------------------------------------------------------------------

template <typename TSpec>
struct RankDictionaryValues_;

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
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct Value<RankDictionary<TwoLevels<TValue, TSpec> > >
{
    typedef TValue  Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size                                           [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct Size<RankDictionary<TwoLevels<TValue, TSpec> > >
{
    // TODO(esiragusa): Choose a better RankDictinonary size type.
    typedef typename Size<String<TValue, TSpec> >::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction RankDictionaryValuesPerBlock_
// ----------------------------------------------------------------------------
// TODO(esiragusa): Specialize it for CUDA.

template <typename TValue, typename TSpec>
struct RankDictionaryValuesPerBlock_<TwoLevels<TValue, TSpec> >
{
//    static const unsigned VALUE = BitsPerValue<unsigned long>::VALUE / BitsPerValue<TValue>::VALUE;
    static const unsigned VALUE = BitsPerValue<unsigned>::VALUE / BitsPerValue<TValue>::VALUE;
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
// Metafunction RankDictionaryValues_                               [TwoLevels]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct RankDictionaryValues_<TwoLevels<TValue, TSpec> >
{
    typedef Tuple<TValue, RankDictionaryValuesPerBlock_<TwoLevels<TValue, TSpec> >::VALUE, BitPacked<> >   Type;
};

// ----------------------------------------------------------------------------
// Metafunction RankDictionaryBitMask_
// ----------------------------------------------------------------------------

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
// Metafunction Fibre                                          [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct Fibre<RankDictionary<TwoLevels<TValue, TSpec> >, FibreRanks>
{
    typedef TwoLevels<TValue, TSpec>                                        TRankDictionarySpec_;
    typedef RankDictionary<TRankDictionarySpec_>                            TRankDictionary_;
    typedef typename RankDictionaryFibreSpec<TRankDictionary_>::Type        TRankDictionaryFibreSpec_;

    typedef String<RankDictionaryEntry_<TRankDictionarySpec_>, TRankDictionaryFibreSpec_>   Type;
};

template <typename TValue, typename TSpec>
struct Fibre<RankDictionary<TwoLevels<TValue, TSpec> > const, FibreRanks>
{
    typedef TwoLevels<TValue, TSpec>                                        TRankDictionarySpec_;
    typedef RankDictionary<TRankDictionarySpec_> const                      TRankDictionary_;
    typedef typename RankDictionaryFibreSpec<TRankDictionary_>::Type        TRankDictionaryFibreSpec_;

    typedef String<RankDictionaryEntry_<TRankDictionarySpec_>, TRankDictionaryFibreSpec_> const Type;
};

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Struct RankDictionaryEntry_
// ----------------------------------------------------------------------------

template <typename TSpec = TwoLevels<> >
struct RankDictionaryEntry_ {};

// ----------------------------------------------------------------------------
// Struct TwoLevels RankDictionaryEntry_
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct RankDictionaryEntry_<TwoLevels<TValue, TSpec> >
{
    typedef TwoLevels<TValue, TSpec>    TRankDictionarySpec;

    // A summary of counts for each block of TValue symbols.
    typename RankDictionaryBlock_<TRankDictionarySpec>::Type    block;

    // A bit-compressed block of TValue symbols.
    typename RankDictionaryValues_<TRankDictionarySpec>::Type   values;
};

// ----------------------------------------------------------------------------
// Class TwoLevels RankDictionary
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct RankDictionary<TwoLevels<TValue, TSpec> >
{
    typedef TwoLevels<TValue, TSpec>                            TRankDictionarySpec;
    typedef RankDictionaryEntry_<TRankDictionarySpec>           TRankDictionaryEntry;
    typedef typename Fibre<RankDictionary, FibreRanks>::Type    TRankDictionaryFibre;
    typedef typename Size<RankDictionary>::Type                 TSize;

    TRankDictionaryFibre ranks;
    TSize _length;

    // ------------------------------------------------------------------------
    // Constructors
    // ------------------------------------------------------------------------

    // NOTE(esiragusa): NVCC cyclic SEQAN_FUNC problem.
    RankDictionary() {};

//    RankDictionary() : _length(0) {};

    template <typename TText>
    RankDictionary(TText const & text)
    {
        createRankDictionary(*this, text);
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function getFibre()                                         [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
SEQAN_FUNC typename Fibre<RankDictionary<TwoLevels<TValue, TSpec> >, FibreRanks>::Type &
getFibre(RankDictionary<TwoLevels<TValue, TSpec> > & dict, FibreRanks)
{
    return dict.ranks;
}

template <typename TValue, typename TSpec>
SEQAN_FUNC typename Fibre<RankDictionary<TwoLevels<TValue, TSpec> >, FibreRanks>::Type const &
getFibre(RankDictionary<TwoLevels<TValue, TSpec> > const & dict, FibreRanks)
{
    return dict.ranks;
}

// ----------------------------------------------------------------------------
// Function _toValuePos()                                      [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TPos>
SEQAN_FUNC typename Size<RankDictionary<TwoLevels<TValue, TSpec> > >::Type
_toValuePos(RankDictionary<TwoLevels<TValue, TSpec> > const & /* dict */, TPos pos)
{
    return pos % RankDictionaryValuesPerBlock_<TwoLevels<TValue, TSpec> >::VALUE;
}

// ----------------------------------------------------------------------------
// Function _toBlockPos()                                      [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TPos>
SEQAN_FUNC typename Size<RankDictionary<TwoLevels<TValue, TSpec> > >::Type
_toBlockPos(RankDictionary<TwoLevels<TValue, TSpec> > const & /* dict */, TPos pos)
{
    return pos / RankDictionaryValuesPerBlock_<TwoLevels<TValue, TSpec> >::VALUE;
}

// ----------------------------------------------------------------------------
// Function _toPos()                                           [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TBlockPos>
inline typename Size<RankDictionary<TwoLevels<TValue, TSpec> > >::Type
_toPos(RankDictionary<TwoLevels<TValue, TSpec> > const & /* dict */, TBlockPos blockPos)
{
    return blockPos * RankDictionaryValuesPerBlock_<TwoLevels<TValue, TSpec> >::VALUE;
}

// ----------------------------------------------------------------------------
// Function _valueAt()                                         [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TPos>
SEQAN_FUNC typename RankDictionaryValues_<TwoLevels<TValue, TSpec> >::Type &
_valueAt(RankDictionary<TwoLevels<TValue, TSpec> > & dict, TPos pos)
{
    return dict.ranks[_toBlockPos(dict, pos)].values;
}

template <typename TValue, typename TSpec, typename TPos>
SEQAN_FUNC typename RankDictionaryValues_<TwoLevels<TValue, TSpec> >::Type const &
_valueAt(RankDictionary<TwoLevels<TValue, TSpec> > const & dict, TPos pos)
{
    return dict.ranks[_toBlockPos(dict, pos)].values;
}

// ----------------------------------------------------------------------------
// Function _blockAt()                                         [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TPos>
SEQAN_FUNC typename RankDictionaryBlock_<TwoLevels<TValue, TSpec> >::Type &
_blockAt(RankDictionary<TwoLevels<TValue, TSpec> > & dict, TPos pos)
{
    return dict.ranks[_toBlockPos(dict, pos)].block;
}

template <typename TValue, typename TSpec, typename TPos>
SEQAN_FUNC typename RankDictionaryBlock_<TwoLevels<TValue, TSpec> >::Type const &
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
SEQAN_FUNC typename Size<RankDictionary<TwoLevels<TValue, TSpec> > const>::Type
_getBlockRank(RankDictionary<TwoLevels<TValue, TSpec> > const & dict, TPos pos, TChar c)
{
    return _blockAt(dict, pos)[ordValue(c)];
}

// ----------------------------------------------------------------------------
// Function _getBlockRank(bool)                                [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TPos>
SEQAN_FUNC typename Size<RankDictionary<TwoLevels<bool, TSpec> > const>::Type
_getBlockRank(RankDictionary<TwoLevels<bool, TSpec> > const & dict, TPos pos, bool c)
{
    // If c == false then return the complementary rank.
    return c ? _blockAt(dict, pos) : pos - _toValuePos(dict, pos) - _blockAt(dict, pos);
}

// ----------------------------------------------------------------------------
// Function _getValueRank()                                    [RankDictionary]
// ----------------------------------------------------------------------------
// NOTE(esiragusa): This version is generic but absymally slow.

template <typename TValue, typename TSpec, typename TPos>
inline typename Size<RankDictionary<TwoLevels<TValue, TSpec> > const>::Type
_getValueRank(RankDictionary<TwoLevels<TValue, TSpec> > const & dict, TPos pos, TValue c)
{
    typedef TwoLevels<TValue, TSpec>                                TRankDictionarySpec;
    typedef RankDictionary<TRankDictionarySpec>                     TRankDictionary;
    typedef typename Size<TRankDictionary>::Type                    TSize;

    TSize valueRank = 0;

    TSize bitPos = _toValuePos(dict, pos);
    for (TSize i = 0; i <= bitPos; ++i)
        valueRank += isEqual(_valueAt(dict, pos)[i], c);

    return valueRank;
}

// ----------------------------------------------------------------------------
// Function _getValueRank(Dna)                                 [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TPos>
SEQAN_FUNC typename Size<RankDictionary<TwoLevels<Dna, TSpec> > const>::Type
_getValueRank(RankDictionary<TwoLevels<Dna, TSpec> > const & dict, TPos pos, Dna c)
{
    typedef TwoLevels<Dna, TSpec>                                       TRankDictionarySpec;
    typedef RankDictionary<TRankDictionarySpec>                         TRankDictionary;
    typedef typename Size<TRankDictionary>::Type                        TSize;
    typedef typename RankDictionaryValues_<TRankDictionarySpec>::Type   TValues;
    typedef typename TValues::TBitVector                                TBitVector;

    TSize bitPos = _toValuePos(dict, pos);
    TValues values = _valueAt(dict, pos);

    // Clear the last bitsPos positions.
    TBitVector word = values.i & ~((TBitVector(1) << (BitsPerValue<TBitVector>::VALUE - 2 - bitPos*2)) - TBitVector(1));

    // And matches when c == G|T.
    TBitVector odd  = ((ordValue(c) & ordValue(Dna('G'))) ? word : ~word) >> 1;

    // And matches when c == C|T.
    TBitVector even = ((ordValue(c) & ordValue(Dna('C'))) ? word : ~word);

    // Apply the interleaved mask.
    TBitVector mask = odd & even & RankDictionaryBitMask_<TBitVector>::VALUE;

    // The rank is the sum of values on.
    TSize valueRank = popCount(mask);

    // If c == A then masked character positions must be subtracted from the count.
    if (c == Dna('A')) valueRank -= RankDictionaryValuesPerBlock_<TwoLevels<Dna, TSpec> >::VALUE - 1 - bitPos;

    return valueRank;
}

// ----------------------------------------------------------------------------
// Function _getValueRank(bool)                                [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TPos>
SEQAN_FUNC typename Size<RankDictionary<TwoLevels<bool, TSpec> > const>::Type
_getValueRank(RankDictionary<TwoLevels<bool, TSpec> > const & dict, TPos pos, bool c)
{
    typedef TwoLevels<bool, TSpec>                                  TRankDictionarySpec;
    typedef RankDictionary<TRankDictionarySpec>                     TRankDictionary;
    typedef typename Size<TRankDictionary>::Type                    TSize;
    typedef typename RankDictionaryValues_<TRankDictionarySpec>::Type TValues;
    typedef typename TValues::TBitVector                              TBitVector;

    TSize bitPos = _toValuePos(dict, pos);
    TValues values = _valueAt(dict, pos);

    // Negate the values to compute the rank of zero.
    TBitVector word = c ? values.i : ~values.i;

    // Clear the last bitsPos positions.
    TBitVector mask = word & ~((TBitVector(1) << (BitsPerValue<TBitVector>::VALUE - 1 - bitPos)) - TBitVector(1));

    // Get the sum of the values on.
    return popCount(mask);
}

// ----------------------------------------------------------------------------
// Function _getValuesRanks()                                  [RankDictionary]
// ----------------------------------------------------------------------------
// TODO(esiragusa): Specialize _getValuesRanks() for Dna.

template <typename TValue, typename TSpec, typename TPos>
inline typename RankDictionaryBlock_<TwoLevels<TValue, TSpec> >::Type
_getValuesRanks(RankDictionary<TwoLevels<TValue, TSpec> > const & dict, TPos pos)
{
    typedef TwoLevels<TValue, TSpec>                                    TRankDictionarySpec;
    typedef typename RankDictionaryBlock_<TRankDictionarySpec>::Type    TBlock;
    typedef typename ValueSize<TValue>::Type                            TValueSize;

    TBlock blockRank;

    for (TValueSize c = 0; c < ValueSize<TValue>::VALUE; ++c)
        assignValue(blockRank, c, _getValueRank(dict, pos, TValue(c)));

    return blockRank;
}

// ----------------------------------------------------------------------------
// Function _getValuesRanks(bool)                              [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TPos>
inline typename RankDictionaryBlock_<TwoLevels<bool, TSpec> >::Type
_getValuesRanks(RankDictionary<TwoLevels<bool, TSpec> > const & dict, TPos pos)
{
    return _getValueRank(dict, pos, true);
}

// ----------------------------------------------------------------------------
// Function getRank()                                          [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TPos, typename TChar>
SEQAN_FUNC typename Size<RankDictionary<TwoLevels<TValue, TSpec> > const>::Type
getRank(RankDictionary<TwoLevels<TValue, TSpec> > const & dict, TPos pos, TChar c)
{
//    SEQAN_ASSERT_LT(pos, length(dict));
    return _getBlockRank(dict, pos, static_cast<TValue>(c)) + _getValueRank(dict, pos, static_cast<TValue>(c));
}

// ----------------------------------------------------------------------------
// Function getRank(bool)                                      [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TPos>
SEQAN_FUNC typename Size<RankDictionary<TwoLevels<bool, TSpec> > const>::Type
getRank(RankDictionary<TwoLevels<bool, TSpec> > const & dict, TPos pos)
{
//    SEQAN_ASSERT_LT(pos, length(dict));
    return getRank(dict, pos, true);
}

// ----------------------------------------------------------------------------
// Function getValue()                                         [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TPos>
SEQAN_FUNC TValue getValue(RankDictionary<TwoLevels<TValue, TSpec> > & dict, TPos pos)
{
//    SEQAN_ASSERT_LT(pos, length(dict));
    return _valueAt(dict, pos)[_toValuePos(dict, pos)];
}

template <typename TValue, typename TSpec, typename TPos>
SEQAN_FUNC TValue getValue(RankDictionary<TwoLevels<TValue, TSpec> > const & dict, TPos pos)
{
//    SEQAN_ASSERT_LT(pos, length(dict));
    return _valueAt(dict, pos)[_toValuePos(dict, pos)];
}

// ----------------------------------------------------------------------------
// Function setValue()                                         [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TPos, typename TChar>
inline void setValue(RankDictionary<TwoLevels<TValue, TSpec> > & dict, TPos pos, TChar c)
{
//    SEQAN_ASSERT_LT(pos, length(dict));
    assignValue(_valueAt(dict, pos), _toValuePos(dict, pos), static_cast<TValue>(c));
}

// ----------------------------------------------------------------------------
// Function appendValue()                                      [RankDictionary]
// ----------------------------------------------------------------------------
// NOTE(esiragusa): Better not to have appendValue() - it is not efficient - and thus neither length().

template <typename TValue, typename TSpec, typename TChar, typename TExpand>
inline void appendValue(RankDictionary<TwoLevels<TValue, TSpec> > & dict, TChar c, Tag<TExpand> const tag)
{
    resize(dict, length(dict) + 1, tag);
    setValue(dict, length(dict) - 1, c);
}

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

    if (empty(dict)) return;

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

        _blockAt(dict, next) = _blockAt(dict, curr) + _getValuesRanks(dict, next - 1);
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
SEQAN_FUNC bool empty(RankDictionary<TwoLevels<TValue, TSpec> > const & dict)
{
    return empty(dict.ranks);
}

// ----------------------------------------------------------------------------
// Function length()                                           [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline typename Size<RankDictionary<TwoLevels<TValue, TSpec> > >::Type
length(RankDictionary<TwoLevels<TValue, TSpec> > const & dict)
{
    return dict._length;
}

// ----------------------------------------------------------------------------
// Function reserve()                                          [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TSize, typename TExpand>
inline typename Size<RankDictionary<TwoLevels<TValue, TSpec> > >::Type
reserve(RankDictionary<TwoLevels<TValue, TSpec> > & dict, TSize newCapacity, Tag<TExpand> const tag)
{
   return reserve(dict.ranks, std::ceil(newCapacity / static_cast<double>(RankDictionaryValuesPerBlock_<TwoLevels<TValue, TSpec> >::VALUE)), tag);
}

// ----------------------------------------------------------------------------
// Function resize()                                           [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TSize, typename TExpand>
inline typename Size<RankDictionary<TwoLevels<TValue, TSpec> > >::Type
resize(RankDictionary<TwoLevels<TValue, TSpec> > & dict, TSize newLength, Tag<TExpand> const tag)
{
    dict._length = newLength;
    return resize(dict.ranks, std::ceil(newLength / static_cast<double>(RankDictionaryValuesPerBlock_<TwoLevels<TValue, TSpec> >::VALUE)), tag);
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
