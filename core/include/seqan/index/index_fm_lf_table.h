// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// LfTable is an object storing all necessary information for the LF-mapping.
// To be more precise, the occurrence-table data structure as well as the
// prefix-sum table are stored.
// ============================================================================

#ifndef INDEX_FM_LF_TABLE_H_
#define INDEX_FM_LF_TABLE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename TText, typename TSpec>
struct LfTable;

// ============================================================================
// Tags
// ============================================================================

/**
.Tag.LF Table Fibres
..summary:Tag to select a specific fibre of a @Spec.FMIndex@.
..remarks:These tags can be used to get @Metafunction.Fibre.Fibres@ of a FM index.
..cat:Index

..tag.FibreOccTable:The occurrence table of the lf table.
..tag.FibrePrefixSumTable:The prefix sum table of the lf table.

..see:Metafunction.Fibre
..see:Function.getFibre
..include:seqan/index_fm.h
*/

struct FibrePrefixSum_;
struct FibreValues_;
struct FibreSentinels_;
struct FibreTempBwt_;

typedef Tag<FibrePrefixSum_>    const FibrePrefixSum;
typedef Tag<FibreValues_>       const FibreValues;
typedef Tag<FibreSentinels_>    const FibreSentinels;
typedef Tag<FibreTempBwt_>      const FibreTempBwt;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec>
struct Value<LfTable<TText, TSpec> >
{
    typedef typename Value<TText>::Type     Type;
};

template <typename TText, typename TSpec>
struct Value<LfTable<TText, TSpec> const> :
    public Value<LfTable<TText, TSpec> > {};

template <typename TText, typename TSSetSpec, typename TSpec>
struct Value<LfTable<StringSet<TText, TSSetSpec>, TSpec> >
{
    typedef typename Value<TText>::Type     Type;
};

template <typename TText, typename TSSetSpec, typename TSpec>
struct Value<LfTable<StringSet<TText, TSSetSpec>, TSpec> const > :
    public Value<LfTable<StringSet<TText, TSSetSpec>, TSpec> > {};

// ----------------------------------------------------------------------------
// Metafunction Fibre
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec>
struct Fibre<LfTable<TText, TSpec>, FibrePrefixSum>
{
    typedef typename Value<LfTable<TText, TSpec> >::Type    TValue_;
    typedef typename MakeUnsigned<TValue_>::Type            TUValue_;
    typedef PrefixSumTable<TUValue_, TSpec>                 Type;
};

template <typename TText, typename TSpec>
struct Fibre<LfTable<TText, TSpec>, FibreValues>
{
    typedef typename Value<LfTable<TText, TSpec> >::Type    TValue_;
    typedef RankDictionary<TwoLevels<TValue_, TSpec> >      Type;
//    typedef RankDictionary<WaveletTree<TValue_> >           Type;
};

template <typename TText, typename TSpec>
struct Fibre<LfTable<TText, TSpec>, FibreSentinels>
{
    typedef typename Size<TText>::Type   Type;
};

template <typename TText, typename TSSetSpec, typename TSpec>
struct Fibre<LfTable<StringSet<TText, TSSetSpec>, TSpec>, FibreSentinels>
{
//    typedef RankDictionary<TwoLevels<bool, TSpec> >         Type;
    typedef RankDictionary<Naive<bool, TSpec> >         Type;
};

template <typename TText, typename TSpec>
struct Fibre<LfTable<TText, TSpec>, FibreTempBwt>
{
    typedef typename Value<LfTable<TText, TSpec> >::Type        TValue_;

    typedef String<TValue_, External<ExternalConfigLarge<> > >  Type;
};

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class LfTable
// ----------------------------------------------------------------------------

/**
.Class.LfTable:
..cat:Index
..summary:LfTable is an object storing all necessary information for the LF-mapping.
..signature:LfTable<TText, TSpec>
..param.TOccTable:The occurrence table data structure.
...type:Spec.WaveletTree
..param.TPrefixSumTable:The specialisation tag.
...default:String
..include:seqan/Index.h
*/
template <typename TText, typename TSpec>
struct LfTable
{
    typename Fibre<LfTable, FibrePrefixSum>::Type   prefixSum;
    typename Fibre<LfTable, FibreValues>::Type      values;
    typename Fibre<LfTable, FibreSentinels>::Type   sentinels;
    typename Value<LfTable>::Type                   sentinelSubstitute;

    // NOTE(esiragusa): NVCC cyclic SEQAN_FUNC problem.
    LfTable() {}

//    LfTable() :
//        sentinelSubstitute(0)
//    {}

//    LfTable(TText const & text)
//    {
//        createLFTable(text);
//    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function getFibre()
// ----------------------------------------------------------------------------

/**
.Function.LfTable#getFibre:
..summary:Returns a specific fibre of a container.
..signature:getFibre(container, fibreTag)
..class:Class.LfTable
..cat:Index
..param.container:The container holding the fibre.
...type:Class.LfTable
..param.fibreTag:A tag that identifies the @Metafunction.Fibre@.
...type:Tag.LF Table Fibres
..returns:A reference to the @Metafunction.Fibre@ object.
..include:seqan/index.h
*/

template <typename TText, typename TSpec>
SEQAN_FUNC typename Fibre<LfTable<TText, TSpec>, FibrePrefixSum>::Type &
getFibre(LfTable<TText, TSpec> & lfTable, FibrePrefixSum)
{
    return lfTable.prefixSum;
}

template <typename TText, typename TSpec>
SEQAN_FUNC typename Fibre<LfTable<TText, TSpec>, FibrePrefixSum>::Type const &
getFibre(LfTable<TText, TSpec> const & lfTable, FibrePrefixSum)
{
    return lfTable.prefixSum;
}

template <typename TText, typename TSpec>
SEQAN_FUNC typename Fibre<LfTable<TText, TSpec>, FibreValues>::Type &
getFibre(LfTable<TText, TSpec> & lfTable, FibreValues)
{
    return lfTable.values;
}

template <typename TText, typename TSpec>
SEQAN_FUNC typename Fibre<LfTable<TText, TSpec>, FibreValues>::Type const &
getFibre(LfTable<TText, TSpec> const & lfTable, FibreValues)
{
    return lfTable.values;
}

template <typename TText, typename TSpec>
SEQAN_FUNC typename Fibre<LfTable<TText, TSpec>, FibreSentinels>::Type &
getFibre(LfTable<TText, TSpec> & lfTable, FibreSentinels)
{
    return lfTable.sentinels;
}

template <typename TText, typename TSpec>
SEQAN_FUNC typename Fibre<LfTable<TText, TSpec>, FibreSentinels>::Type const &
getFibre(LfTable<TText, TSpec> const & lfTable, FibreSentinels)
{
    return lfTable.sentinels;
}

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

/**
.Function.LfTable#empty
..class:Class.LfTable
..summary:Clears the LF table.
..signature:empty(lfTable)
..param.lfTable:The LF table to be cleared.
...type:Class.LfTable
..returns:$true$ if the LF table is empty, $false$ otherwise.
...type:nolink:$bool$
..include:seqan/index.h
*/

template <typename TText, typename TSpec>
SEQAN_FUNC bool empty(LfTable<TText, TSpec> const & lfTable)
{
    return empty(lfTable.values) &&
           empty(lfTable.prefixSum);
}

template <typename TText, typename TSSetSpec, typename TSpec>
SEQAN_FUNC bool empty(LfTable<StringSet<TText, TSSetSpec>, TSpec> const & lfTable)
{
    return empty(lfTable.values) &&
           empty(lfTable.sentinels) &&
           empty(lfTable.prefixSum);
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

/**
.Function.LfTable#clear
..class:Class.LfTable
..summary:Clears the LF table.
..signature:clear(lfTable)
..param.lfTable:The LF table to be cleared.
...type:Class.LfTable
..include:seqan/index.h
*/

template <typename TText, typename TSpec>
inline void clear(LfTable<TText, TSpec> & lfTable)
{
    clear(lfTable.values);
    _clearSentinels(lfTable);
    clear(lfTable.prefixSum);
    lfTable.sentinelSubstitute = 0;
}

// ----------------------------------------------------------------------------
// Function _clearSentinels()
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec>
inline void _clearSentinels(LfTable<TText, TSpec> & lfTable)
{
    getFibre(lfTable, FibreSentinels()) = 0;
}

template <typename TText, typename TSSetSpec, typename TSpec>
inline void _clearSentinels(LfTable<StringSet<TText, TSSetSpec>, TSpec> & lfTable)
{
    clear(getFibre(lfTable, FibreSentinels()));
}

// ----------------------------------------------------------------------------
// Function sentinelAt()
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TPos>
SEQAN_FUNC bool sentinelAt(LfTable<TText, TSpec> const & lfTable, TPos pos)
{
    return getFibre(lfTable, FibreSentinels()) == pos;
}

template <typename TText, typename TSSetSpec, typename TSpec, typename TPos>
SEQAN_FUNC bool sentinelAt(LfTable<StringSet<TText, TSSetSpec>, TSpec> const & lfTable, TPos pos)
{
    return getValue(getFibre(lfTable, FibreSentinels()), pos);
}

// ----------------------------------------------------------------------------
// Function _getRank()
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TPos, typename TValue>
SEQAN_FUNC typename Size<typename Fibre<LfTable<TText, TSpec>, FibreValues>::Type>::Type
_getRank(LfTable<TText, TSpec> const & lfTable, TPos pos, TValue character)
{
    typedef LfTable<TText, TSpec> const                         TLfTable;
    typedef typename Fibre<TLfTable, FibreValues>::Type         TValues;
    typedef typename Size<TValues>::Type                        TSize;

    TValues const & values = getFibre(lfTable, FibreValues());
    TSize rank = getRank(values, pos, character);

    // TODO(esiragusa): Make single sentinel case transparent.
    if (ordEqual(lfTable.sentinelSubstitute, character) && pos >= getFibre(lfTable, FibreSentinels()))
         --rank;

    return rank;
}

template <typename TText, typename TSSetSpec, typename TSpec, typename TPos, typename TValue>
SEQAN_FUNC typename Size<typename Fibre<LfTable<StringSet<TText, TSSetSpec>, TSpec>, FibreValues>::Type>::Type
_getRank(LfTable<StringSet<TText, TSSetSpec>, TSpec> const & lfTable, TPos pos, TValue character)
{
    typedef LfTable<StringSet<TText, TSSetSpec>, TSpec> const   TLfTable;
    typedef typename Fibre<TLfTable, FibreValues>::Type         TValues;
    typedef typename Size<TValues>::Type                        TSize;

    TValues const & values = getFibre(lfTable, FibreValues());
    TSize rank = getRank(values, pos, character);

    if (ordEqual(lfTable.sentinelSubstitute, character))
        return rank - getRank(getFibre(lfTable, FibreSentinels()), pos);

    return rank;
}

// ----------------------------------------------------------------------------
// Function lfMapping()
// ----------------------------------------------------------------------------

/**
.Function.LfTable#lfMapping:
..summary:Returns the position of an character at a specified position of L in F. L corresponds to the last column of 
the sorted cyclic rotations of the original text, while F correspond to the first column.
..cat:Index
..signature:lfMapping(lfTable, pos)
..param.lfTable:The @Class.LfTable@ holding the occurrence and prefix sum table.
...type:Class.LfTable
..param.pos:The position in L
..returns:Returns the position of the character L[c] in F. The returned position is of the same type as pos.
..include:seqan/index.h
*/

// TODO(esiragusa): rename lfMapping() as getValue() or getFrontPos()?
template <typename TLfTable, typename TPos>
SEQAN_FUNC TPos
lfMapping(TLfTable const & lfTable, TPos pos)
{
    typedef typename Fibre<TLfTable const, FibreValues>::Type       TValues;
    typedef typename Fibre<TLfTable const, FibreSentinels>::Type    TSentinels;
    typedef typename Fibre<TLfTable const, FibrePrefixSum>::Type    TPrefixSum;
    typedef typename Value<TLfTable const>::Type                    TValue;

    TValues const & values = getFibre(lfTable, FibreValues());
    TPrefixSum const & prefixSum = getFibre(lfTable, FibrePrefixSum());

    TValue c = getValue(values, pos);

    return _getRank(lfTable, pos, c) + getPrefixSum(prefixSum, ordValue(c)) - 1;
}

// ----------------------------------------------------------------------------
// Function _setSentinelSubstitute()
// ----------------------------------------------------------------------------
// This function determines the '$' substitute.
// The character with the smallest number of occurrences greater 0 is chosen.

template <typename TText, typename TSpec>
inline void _setSentinelSubstitute(LfTable<TText, TSpec> & lfTable)
{
    typedef LfTable<TText, TSpec>                                   TLfTable;
    typedef typename Fibre<TLfTable, FibrePrefixSum>::Type          TPrefixSum;
    typedef typename Size<TPrefixSum>::Type                         TSize;
    typedef typename Value<TPrefixSum>::Type                        TValue;

    TPrefixSum & prefixSum = getFibre(lfTable, FibrePrefixSum());

    TValue minRank = MaxValue<TValue>::VALUE;
    TSize ordVal = length(prefixSum) - 1;

    // NOTE(esiragusa): This doesn't work for maps!!
    for (TSize i = 0; i < length(prefixSum) - 1; ++i)
    {
        TSize diff = prefixSum[i + 1] - prefixSum[i];
        if (diff != 0 && diff < minRank)
        {
            minRank = diff;
            ordVal = i;
        }
    }

    lfTable.sentinelSubstitute = ordVal;
}

// ----------------------------------------------------------------------------
// Function _computeBwtLength()
// ----------------------------------------------------------------------------

// This function computes the length of the bwt string.
template <typename TText>
inline typename Size<TText>::Type
_computeBwtLength(TText const & text)
{
    return length(text) + 1;
}

// This function computes the length of the bwt string.
template <typename TText, typename TSSetSpec>
inline typename Size<StringSet<TText, TSSetSpec> >::Type
_computeBwtLength(StringSet<TText, TSSetSpec> const & text)
{
    return lengthSum(text) + countSequences(text);
}

// ----------------------------------------------------------------------------
// Function _createBwt()
// ----------------------------------------------------------------------------
// This function computes the BWT of a text. Note that the sentinel sign is substituted and its position stored.

template <typename TText, typename TSpec, typename TBwt, typename TOtherText, typename TSA>
inline void
_createBwt(LfTable<TText, TSpec> & lfTable, TBwt & bwt, TOtherText const & text, TSA const & sa)
{
    typedef typename GetValue<TSA>::Type                    TSAValue;
    typedef typename Iterator<TSA const, Standard>::Type    TSAIter;
    typedef typename Iterator<TBwt, Standard>::Type         TBwtIter;

    TSAIter saIt = begin(sa, Standard());
    TSAIter saItEnd = end(sa, Standard());
    TBwtIter bwtIt = begin(bwt, Standard());

	assignValue(bwtIt, back(text));
    ++bwtIt;

    for (; saIt != saItEnd; ++saIt, ++bwtIt)
    {
        TSAValue pos = getValue(saIt);

        if (pos != 0)
        {
            assignValue(bwtIt, getValue(text, pos - 1));
        }
        else
        {
            assignValue(bwtIt, lfTable.sentinelSubstitute);
            lfTable.sentinels = bwtIt - begin(bwt, Standard());
        }
    }
}

// ----------------------------------------------------------------------------
// Function _createBwt()
// ----------------------------------------------------------------------------
// This function computes the BWT of a text. Note that the sentinel sign is substituted and its position stored.

template <typename TText, typename TSSetSpec, typename TSpec, typename TBwt, typename TOtherText, typename TSA>
inline void
_createBwt(LfTable<StringSet<TText, TSSetSpec>, TSpec > & lfTable, TBwt & bwt, TOtherText const & text, TSA const & sa)
{
    typedef typename Value<TSA>::Type                       TSAValue;
    typedef typename Size<TSA>::Type                        TSize;
    typedef typename Iterator<TSA const, Standard>::Type    TSAIter;
    typedef typename Iterator<TBwt, Standard>::Type         TBwtIter;

    TSize seqNum = countSequences(text);
    TSize totalLen = lengthSum(text);

    resize(lfTable.sentinels, seqNum + totalLen, Exact());

    TSAIter saIt = begin(sa, Standard());
    TSAIter saItEnd = end(sa, Standard());
    TBwtIter bwtItBeg = begin(bwt, Standard());
    TBwtIter bwtIt = bwtItBeg;

    // Fill the sentinel positions (they are all at the beginning of the bwt).
    for (TSize i = 1; i <= seqNum; ++i, ++bwtIt)
        assignValue(bwtIt, back(text[seqNum - i]));

    // Compute the rest of the bwt.
    for (; saIt != saItEnd; ++saIt, ++bwtIt)
    {
        TSAValue pos;    // = SA[i];
        posLocalize(pos, getValue(saIt), stringSetLimits(text));
        
        if (getSeqOffset(pos) != 0)
        {
            assignValue(bwtIt, getValue(getValue(text, getSeqNo(pos)), getSeqOffset(pos) - 1));
            setValue(lfTable.sentinels, bwtIt - bwtItBeg, false);
        }
        else
        {
            assignValue(bwtIt, lfTable.sentinelSubstitute);
            setValue(lfTable.sentinels, bwtIt - bwtItBeg, true);
        }
    }

    // Update the auxiliary RankDictionary of sentinel positions.
    updateRanks(lfTable.sentinels);
}

// ----------------------------------------------------------------------------
// Function _createRankDictionary()
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TBwt, typename TPrefixSum>
inline void
_createRankDictionary(LfTable<TText, TSpec> & lfTable,
                      TBwt const & bwt,
                      TPrefixSum & /* prefixSum */)
{
    createRankDictionary(getFibre(lfTable, FibreValues()), bwt);
}

// TODO(esiragusa): Specialize _createRankDictionary() for WaveletTree.
//template <typename TText, typename TSpec, typename TBwt, typename TPrefixSum>
//inline void
//_createRankDictionary(LfTable<TText, TSpec> & lfTable,
//                      TBwt const & bwt,
//                      TPrefixSum & prefixSum)
//{
//    createRankDictionary(getFibre(lfTable, FibreValues()), bwt, prefixSum);
//}

// ----------------------------------------------------------------------------
// Function _insertSentinels()
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TNumSentinel>
void _insertSentinels(LfTable<TText, TSpec> & lfTable, TNumSentinel numSentinel)
{
    typedef LfTable<TText, TSpec>                                   TLfTable;
    typedef typename Fibre<TLfTable, FibrePrefixSum>::Type          TPrefixSum;
    typedef typename Size<TPrefixSum>::Type                         TSize;

    TPrefixSum & prefixSum_ = getFibre(lfTable, FibrePrefixSum());

    // NOTE(esiragusa): This doesn't work for maps!!
    for (TSize i = 0; i < length(prefixSum_); ++i)
        prefixSum(prefixSum_, i) = getPrefixSum(prefixSum_, i) + numSentinel;
}

// ----------------------------------------------------------------------------
// Function createLFTable()
// ----------------------------------------------------------------------------

// This function creates all table of the lf table given a text and a suffix array.
template <typename TText, typename TSpec, typename TOtherText, typename TSA>
inline void createLFTable(LfTable<TText, TSpec> & lfTable, TOtherText const & text, TSA const & sa)
{
    typedef LfTable<TText, TSpec>                                   TLfTable;
    typedef typename Fibre<TLfTable, FibrePrefixSum>::Type          TPrefixSum;
    typedef typename Fibre<TLfTable, FibreValues>::Type             TValues;
    typedef typename Value<TLfTable>::Type                          TValue;
    typedef typename Fibre<TLfTable, FibreTempBwt>::Type            TBwt;

    TPrefixSum & prefixSum = getFibre(lfTable, FibrePrefixSum());

    createPrefixSumTable(prefixSum, text);

    _setSentinelSubstitute(lfTable);

    // TODO(esiragusa): Rename this to fillValues() / fillSentinels().
    TBwt bwt;
    resize(bwt, _computeBwtLength(text), Exact());
    _createBwt(lfTable, bwt, text, sa);

    _createRankDictionary(lfTable, bwt, prefixSum);

    _insertSentinels(lfTable, countSequences(text));
}

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

/**
.Function.LfTable#open
..class:Class.LfTable
..summary:This functions loads a LF table from disk.
..signature:open(lfTable, fileName [, openMode])
..param.lfTable:The lfTable.
...type:Class.LfTable
..param.fileName:C-style character string containing the file name.
..param.openMode:The combination of flags defining how the file should be opened.
...remarks:To open a file read-only, write-only or to read and write use $OPEN_RDONLY$, $OPEN_WRONLY$, or $OPEN_RDWR$.
...remarks:To create or overwrite a file add $OPEN_CREATE$.
...remarks:To append a file if existing add $OPEN_APPEND$.
...remarks:To circumvent problems, files are always opened in binary mode.
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A nolink:$bool$ which is $true$ on success.
..include:seqan/index.h
*/

template <typename TText, typename TSpec>
inline bool open(LfTable<TText, TSpec> & lfTable, const char * fileName, int openMode)
{
    String<char> name;

    name = fileName;
    append(name, ".pst");
    if (!open(getFibre(lfTable, FibrePrefixSum()), toCString(name), openMode)) return false;

    name = fileName;
    append(name, ".drv");
    if (!open(getFibre(lfTable, FibreValues()), toCString(name), openMode)) return false;

    name = fileName;
    append(name, ".drp");
    if (!open(getFibre(lfTable, FibreSentinels()), toCString(name), openMode)) return false;

    name = fileName;
    append(name, ".drs");
    if (!open(lfTable.sentinelSubstitute, toCString(name), openMode)) return false;

    return true;
}

template <typename TText, typename TSpec>
inline bool open(LfTable<TText, TSpec> & lfTable, const char * fileName)
{
    return open(lfTable, fileName, DefaultOpenMode<LfTable<TText, TSpec> >::VALUE);
}

// ----------------------------------------------------------------------------
// Function save()
// ----------------------------------------------------------------------------

/**
.Function.LfTable#save
..class:Class.LfTable
..summary:This functions saves a LF table to disk.
..signature:save(lfTable, fileName [, openMode])
..param.lfTable:The dictionary.
...type:Class.LfTable
..param.fileName:C-style character string containing the file name.
..param.openMode:The combination of flags defining how the file should be opened.
...remarks:To open a file read-only, write-only or to read and write use $OPEN_RDONLY$, $OPEN_WRONLY$, or $OPEN_RDWR$.
...remarks:To create or overwrite a file add $OPEN_CREATE$.
...remarks:To append a file if existing add $OPEN_APPEND$.
...remarks:To circumvent problems, files are always opened in binary mode.
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A nolink:$bool$ which is $true$ on success.
..include:seqan/index.h
*/

template <typename TText, typename TSpec>
inline bool save(LfTable<TText, TSpec> const & lfTable, const char * fileName, int openMode)
{
    String<char> name;

    name = fileName;
    append(name, ".pst");
    if (!save(getFibre(lfTable, FibrePrefixSum()), toCString(name), openMode)) return false;

    name = fileName;
    append(name, ".drv");
    if (!save(getFibre(lfTable, FibreValues()), toCString(name), openMode)) return false;

    name = fileName;
    append(name, ".drp");
    if (!save(getFibre(lfTable, FibreSentinels()), toCString(name), openMode)) return false;

    name = fileName;
    append(name, ".drs");
    if (!save(lfTable.sentinelSubstitute, toCString(name), openMode)) return false;

    return true;
}

template <typename TText, typename TSpec>
inline bool save(LfTable<TText, TSpec> const & lfTable, const char * fileName)
{
    return save(lfTable, fileName, DefaultOpenMode<LfTable<TText, TSpec> >::VALUE);
}

}
#endif // INDEX_FM_LF_TABLE_H_
