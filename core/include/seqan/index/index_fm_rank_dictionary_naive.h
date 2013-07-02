// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2013 NVIDIA Corporation
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
//     * Neither the name of NVIDIA Corporation nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL NVIDIA CORPORATION BE LIABLE
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

#ifndef INDEX_FM_RANK_DICTIONARY_NAIVE_H_
#define INDEX_FM_RANK_DICTIONARY_NAIVE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

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
// Tag Naive
// ----------------------------------------------------------------------------

template <typename TValue = bool, typename TSpec = void>
struct Naive {};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Value                                          [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct Value<RankDictionary<Naive<TValue, TSpec> > >
{
    typedef TValue  Type;
};

template <typename TValue, typename TSpec>
struct Value<RankDictionary<Naive<TValue, TSpec> > const> :
    Value<RankDictionary<Naive<TValue, TSpec> > > {};

// ----------------------------------------------------------------------------
// Metafunction Fibre                                          [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct Fibre<RankDictionary<Naive<TValue, TSpec> >, FibreRanks>
{
    typedef Naive<TValue, TSpec>                                        TRankDictionarySpec_;
    typedef RankDictionary<TRankDictionarySpec_>                        TRankDictionary_;
    typedef typename Size<TRankDictionary_>::Type                       TSize_;
    typedef typename RankDictionaryFibreSpec<TRankDictionary_>::Type    TRankDictionaryFibreSpec_;

    typedef String<TSize_, TRankDictionaryFibreSpec_>   Type;
};

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class Naive RankDictionary
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct RankDictionary<Naive<TValue, TSpec> >
{
    // ------------------------------------------------------------------------
    // Fibres
    // ------------------------------------------------------------------------

    typename Fibre<RankDictionary, FibreRanks>::Type    ranks;

    // ------------------------------------------------------------------------
    // Constructors
    // ------------------------------------------------------------------------

    RankDictionary() {};

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
SEQAN_FUNC typename Fibre<RankDictionary<Naive<TValue, TSpec> >, FibreRanks>::Type &
getFibre(RankDictionary<Naive<TValue, TSpec> > & dict, FibreRanks)
{
    return dict.ranks;
}

template <typename TValue, typename TSpec>
SEQAN_FUNC typename Fibre<RankDictionary<Naive<TValue, TSpec> >, FibreRanks>::Type const &
getFibre(RankDictionary<Naive<TValue, TSpec> > const & dict, FibreRanks)
{
    return dict.ranks;
}

// ----------------------------------------------------------------------------
// Function getRank()                                          [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TPos>
SEQAN_FUNC typename Size<RankDictionary<Naive<bool, TSpec> > const>::Type
getRank(RankDictionary<Naive<bool, TSpec> > const & dict, TPos pos, bool c)
{
    typedef Naive<bool, TSpec>                                              TRankDictionarySpec;
    typedef RankDictionary<TRankDictionarySpec> const                       TRankDictionary;
    typedef typename Fibre<TRankDictionary, FibreRanks>::Type               TFibreRanks;
    typedef typename Iterator<TFibreRanks, Standard>::Type                  TFibreRanksIterator;

    TFibreRanksIterator ranksBegin = begin(dict.ranks, Standard());
    TFibreRanksIterator ranksEnd   = end(dict.ranks, Standard());
    TFibreRanksIterator ranksIt    = ranksBegin;

    for (; ranksIt != ranksEnd && value(ranksIt) < pos; ++ranksIt);

    return c ? ranksIt - ranksBegin : pos - (ranksIt - ranksBegin);
}

// ----------------------------------------------------------------------------
// Function getRank(bool)                                      [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TPos>
SEQAN_FUNC typename Size<RankDictionary<Naive<bool, TSpec> > const>::Type
getRank(RankDictionary<Naive<bool, TSpec> > const & dict, TPos pos)
{
    return getRank(dict, pos, true);
}

// ----------------------------------------------------------------------------
// Function getValue()                                         [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TPos>
SEQAN_FUNC typename Value<RankDictionary<Naive<bool, TSpec> > >::Type
getValue(RankDictionary<Naive<bool, TSpec> > & dict, TPos pos)
{
    typedef Naive<bool, TSpec>                                              TRankDictionarySpec;
    typedef RankDictionary<TRankDictionarySpec> const                       TRankDictionary;
    typedef typename Fibre<TRankDictionary, FibreRanks>::Type               TFibreRanks;
    typedef typename Iterator<TFibreRanks, Standard>::Type                  TFibreRanksIterator;

    TFibreRanksIterator ranksBegin = begin(dict.ranks, Standard());
    TFibreRanksIterator ranksEnd   = end(dict.ranks, Standard());
    TFibreRanksIterator ranksIt    = ranksBegin;

    for (; ranksIt != ranksEnd && value(ranksIt) < pos; ++ranksIt);

    return ranksIt != ranksEnd && value(ranksIt) == pos;
}

template <typename TSpec, typename TPos>
SEQAN_FUNC typename Value<RankDictionary<Naive<bool, TSpec> > const>::Type
getValue(RankDictionary<Naive<bool, TSpec> > const & dict, TPos pos)
{
    typedef Naive<bool, TSpec>                                              TRankDictionarySpec;
    typedef RankDictionary<TRankDictionarySpec> const                       TRankDictionary;
    typedef typename Fibre<TRankDictionary, FibreRanks>::Type               TFibreRanks;
    typedef typename Iterator<TFibreRanks, Standard>::Type                  TFibreRanksIterator;

    TFibreRanksIterator ranksBegin = begin(dict.ranks, Standard());
    TFibreRanksIterator ranksEnd   = end(dict.ranks, Standard());
    TFibreRanksIterator ranksIt    = ranksBegin;

    for (; ranksIt != ranksEnd && value(ranksIt) < pos; ++ranksIt);

    return ranksIt != ranksEnd && value(ranksIt) == pos;
}

// ----------------------------------------------------------------------------
// Function setValue()                                         [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TPos, typename TChar>
inline void setValue(RankDictionary<Naive<TValue, TSpec> > & dict, TPos pos, TChar c)
{
    SEQAN_ASSERT_GT(pos, back(dict.ranks));

    if (c == false) return;

    appendValue(dict.ranks, pos);
}

// ----------------------------------------------------------------------------
// Function appendValue()                                      [RankDictionary]
// ----------------------------------------------------------------------------
// NOTE(esiragusa): Better not to have appendValue() - it is not efficient - and thus neither length().

template <typename TSpec, typename TChar, typename TExpand>
inline void appendValue(RankDictionary<Naive<bool, TSpec> > & dict, TChar c, Tag<TExpand> const tag)
{
    if (c == false) return;

// NOTE(esiragusa): RankDictionary's resize() is desabled.
//    resize(dict, length(dict) + 1, tag);
    resize(dict.ranks, length(dict) + 1, tag);
    setValue(dict, length(dict) - 1, c);
}

// ----------------------------------------------------------------------------
// Function updateRanks()                                      [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline void updateRanks(RankDictionary<Naive<TValue, TSpec> > & /* dict */)
{}

// ----------------------------------------------------------------------------
// Function createRankDictionary()                             [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TText>
inline void
createRankDictionary(RankDictionary<Naive<bool, TSpec> > & dict, TText const & text)
{
    typedef Naive<bool, TSpec>                                      TRankDictionarySpec;
    typedef RankDictionary<TRankDictionarySpec>                     TRankDictionary;
    typedef typename Iterator<TText const, Standard>::Type          TTextIterator;

    // Assign the text value by value.
    TTextIterator textBegin = begin(text, Standard());
    TTextIterator textEnd = end(text, Standard());
    for (TTextIterator textIt = textBegin; textIt != textEnd; ++textIt)
        appendValue(dict, textIt - textBegin, value(textIt));
}

// ----------------------------------------------------------------------------
// Function clear()                                            [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline void clear(RankDictionary<Naive<TValue, TSpec> > & dict)
{
    clear(dict.ranks);
}

template <typename TValue, typename TSpec>
SEQAN_FUNC bool empty(RankDictionary<Naive<TValue, TSpec> > const & dict)
{
    return empty(dict.ranks);
}

// ----------------------------------------------------------------------------
// Function length()                                           [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline typename Size<RankDictionary<Naive<TValue, TSpec> > >::Type
length(RankDictionary<Naive<TValue, TSpec> > const & dict)
{
    return length(dict.ranks);
}

// ----------------------------------------------------------------------------
// Function reserve()                                          [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TSize, typename TExpand>
inline typename Size<RankDictionary<Naive<TValue, TSpec> > >::Type
reserve(RankDictionary<Naive<TValue, TSpec> > & dict, TSize newCapacity, Tag<TExpand> const tag)
{
   return reserve(dict.ranks, newCapacity, tag);
}

// ----------------------------------------------------------------------------
// Function resize()                                           [RankDictionary]
// ----------------------------------------------------------------------------
// NOTE(esiragusa): disabled because LfTable::_createBwt() resizes the rank dict to the bwt length.

template <typename TValue, typename TSpec, typename TSize, typename TExpand>
inline typename Size<RankDictionary<Naive<TValue, TSpec> > >::Type
resize(RankDictionary<Naive<TValue, TSpec> > & dict, TSize /* newLength */, Tag<TExpand> const /* tag */)
{
    return length(dict);
//    return resize(dict.ranks, newLength, tag);
}

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline bool open(RankDictionary<Naive<TValue, TSpec> > & dict, const char * fileName, int openMode)
{
    return open(dict.ranks, fileName, openMode);
}

template <typename TValue, typename TSpec>
inline bool open(RankDictionary<Naive<TValue, TSpec> > & dict, const char * fileName)
{
    return open(dict, fileName, DefaultOpenMode<RankDictionary<Naive<TValue, TSpec> > >::VALUE);
}

// ----------------------------------------------------------------------------
// Function save()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline bool save(RankDictionary<Naive<TValue, TSpec> > const & dict, const char * fileName, int openMode)
{
    return save(dict.ranks, fileName, openMode);
}

template <typename TValue, typename TSpec>
inline bool save(RankDictionary<Naive<TValue, TSpec> > const & dict, const char * fileName)
{
    return save(dict, fileName, DefaultOpenMode<RankDictionary<Naive<TValue, TSpec> > >::VALUE);
}

}

#endif  // INDEX_FM_RANK_DICTIONARY_NAIVE_H_
