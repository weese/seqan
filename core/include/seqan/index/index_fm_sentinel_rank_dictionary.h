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
// ==========================================================================

#ifndef INDEX_FM_SENTINEL_RANK_DICTIONARY_H_
#define INDEX_FM_SENTINEL_RANK_DICTIONARY_H_

namespace seqan {

// ==========================================================================
// Forwards
// ==========================================================================

template <typename TValue>
struct RankDictionary;

template<typename TRankDictionarySpec, typename TSpec> 
struct SentinelRankDictionary;

// ==========================================================================
// Tags
// ==========================================================================

struct Sentinel_;
struct Sentinels_;
struct FibreRankDictionary_;
struct FibreSentinelPosition_;

typedef Tag<Sentinel_> const                Sentinel;
typedef Tag<Sentinels_> const               Sentinels;
typedef Tag<FibreRankDictionary_> const     FibreRankDictionary;
typedef Tag<FibreSentinelPosition_> const   FibreSentinelPosition;

// ----------------------------------------------------------------------------
// Spec SentinelRankDictionary Fibres
// ----------------------------------------------------------------------------

/**
.Tag.SentinelRankDictionary Fibres
..summary:Tag to select a specific fibre (e.g. table, object, ...) of a @Class.SentinelRankDictionary@.
..remarks:These tags can be used to get @Metafunction.Fibre.Fibres@ of a SentinelRankDictionary.
..cat:SentinelRankDictionary

..tag.FibreRankDictionary:The rank dictionary.

..tag.FibreSentinelPosition:The bit string encoding the position of the sentinel sign.

..see:Metafunction.Fibre
..see:Function.getFibre
..include:seqan/index.h
*/

// ==========================================================================
// Metafunctions
// ==========================================================================

// ----------------------------------------------------------------------------
// Metafunction Fibre
// ----------------------------------------------------------------------------

// TODO(DOC)
///.Metafunction.Fibre.param.TSpec.type:Tag.SentinelRankDictionary Fibres

template <typename TValue, typename TSpec>
struct Fibre<SentinelRankDictionary<RankDictionary<WaveletTree<TValue> >, TSpec>, FibreRankDictionary>
{
    typedef RankDictionary<WaveletTree<TValue> > Type;
};

template <typename TValue, typename TRankDictSpec, typename TSpec>
struct Fibre<SentinelRankDictionary<RankDictionary<TwoLevels<TValue, TRankDictSpec> >, TSpec>, FibreRankDictionary>
{
    typedef RankDictionary<TwoLevels<TValue, TRankDictSpec> >   Type;
};

template <typename TRankDictionary, typename TSpec>
struct Fibre<SentinelRankDictionary<TRankDictionary, TSpec> const, FibreRankDictionary>
{
    typedef typename Fibre<SentinelRankDictionary<TRankDictionary, TSpec>, FibreRankDictionary>::Type const Type;
};

template <typename TRankDictionary>
struct Fibre<SentinelRankDictionary<TRankDictionary, Sentinel>, FibreSentinelPosition>
{
    typedef typename Size<TRankDictionary>::Type Type;
};

template <typename TRankDictionary>
struct Fibre<SentinelRankDictionary<TRankDictionary, Sentinel> const, FibreSentinelPosition>
{
    typedef typename Fibre<SentinelRankDictionary<TRankDictionary, Sentinel>, FibreSentinelPosition>::Type const Type;
};

template <typename TRankDictionary>
struct Fibre<SentinelRankDictionary<TRankDictionary, Sentinels>, FibreSentinelPosition>
{
    typedef RankDictionary<TwoLevels<bool> >    Type;
};

template <typename TRankDictionary>
struct Fibre<SentinelRankDictionary<TRankDictionary, Sentinels> const, FibreSentinelPosition>
{
    typedef typename Fibre<SentinelRankDictionary<TRankDictionary, Sentinels>, FibreSentinelPosition>::Type const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TRankDictionary, typename TSpec>
struct Value<SentinelRankDictionary<TRankDictionary, TSpec> > :
    public Value<TRankDictionary> {};

template <typename TRankDictionary, typename TSpec>
struct Value<SentinelRankDictionary<TRankDictionary, TSpec> const> :    
    public Value<TRankDictionary const> {};

// ==========================================================================
// Classes
// ==========================================================================

// ----------------------------------------------------------------------------
// Class SentinelRankDictionary
// ----------------------------------------------------------------------------

/**
.Class.SentinelRankDictionary:
..cat:Index
..summary:A rank dictionary, additional storing sentinel character which are not accounted for in a rank querry.
..signature:SentinelRankDictionary<TRankDictionary, TSpec>
..param.TRankDictionary:The rank dictionary of a text.
...type:Class.RankDictionary
..param.TSpec:Specialisation
..include:seqan/index.h
*/

template <typename TRankDictionary, typename TSpec>
struct SentinelRankDictionary
{
    typedef typename Value<SentinelRankDictionary>::Type                            TValue;
    typedef typename Fibre<SentinelRankDictionary, FibreSentinelPosition>::Type     TSentinelPosition;
//    typedef typename Fibre<SentinelRankDictionary, FibreBwt>::Type                  TBwt;

//    TBwt                bwt;
    TRankDictionary     rankDictionary;
    TSentinelPosition   sentinelPosition;
    TValue              sentinelSubstitute;

    SentinelRankDictionary() :
        sentinelSubstitute()
    {}

// NOTE(esiragusa): createSentinelRankDictionary() needs also a PrefixSumTable object.
//    template <typename TText, typename TSA>
//    SentinelRankDictionary(TText const & text, TSA const & sa) :
//        sentinelSubstitute()
//    {
//        createSentinelRankDictionary(*this, text, sa);
//    }
};

// ==========================================================================
// Functions
// ==========================================================================

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

/**
.Function.SentinelRankDictionary#clear
..class:Class.SentinelRankDictionary
..summary:Clears the dictionary.
..signature:clear(dictionary)
..param.dictionary:The rank dictionary to be cleared.
...type:Class.SentinelRankDictionary
..include:seqan/index.h
*/

template <typename TRankDictionary>
inline void _clearSentinel(SentinelRankDictionary<TRankDictionary, Sentinel> &)
{}

template <typename TRankDictionary>
inline void _clearSentinel(SentinelRankDictionary<TRankDictionary, Sentinels> & dict)
{
    clear(getFibre(dict, FibreSentinelPosition()));
}

template <typename TRankDictionary, typename TSpec>
inline void clear(SentinelRankDictionary<TRankDictionary, TSpec> & dict)
{
    clear(getFibre(dict, FibreRankDictionary()));
    _clearSentinel(dict);
}

// ----------------------------------------------------------------------------
// Function isSentinelPosition()
// ----------------------------------------------------------------------------

/**
.Function.SentinelRankDictionary#sentinelPosition
..class:Class.SentinelRankDictionary
..summary:Returns whether a specified position is a sentinel position.
..signature:isSentinelPosition(dictionary, pos)
..param.dictionary:The dictionary.
...type:Class.SentinelRankDictionary
..param.pos:The position.
...type:Concept.UnsignedIntegerConcept
..include:seqan/index.h
*/

//.Function.sentinelPosition.param.type:Class.RankDictionary
template <typename TRankDictionary, typename TPos>
inline bool isSentinelPosition(SentinelRankDictionary<TRankDictionary, Sentinel> const & dict, TPos pos)
{
    return dict.sentinelPosition == pos;
}

template <typename TRankDictionary, typename TPos>
inline bool isSentinelPosition(SentinelRankDictionary<TRankDictionary, Sentinels> const & dict, TPos pos)
{
    return getValue(getFibre(dict, FibreSentinelPosition()), pos);
}

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

/**
.Function.SentinelRankDictionary#empty
..class:Class.SentinelRankDictionary
..summary:Returns whether or not the dictionary is empty.
..signature:empty(dictionary)
..param.dictionary:The rank dictionary to be checked.
...type:Class.SentinelRankDictionary
..include:seqan/index.h
*/
template <typename TRankDictionary, typename TSpec>
inline bool empty(SentinelRankDictionary<TRankDictionary, TSpec> const & dict)
{
    return empty(getFibre(dict, FibreRankDictionary()));
}

// ----------------------------------------------------------------------------
// Function getValue()
// ----------------------------------------------------------------------------

/**
.Function.SentinelRankDictionary#getValue
..class:Class.SentinelRankDictionary
..summary:Returns the character of a specified position.
..signature:getCharacter(dictionary, pos)
..param.dictionary:The rank dictionary.
...type:Class.RankDictionary
..param.pos:The position
..include:seqan/index.h
..example.code:
*/
template <typename TRankDictionary, typename TSpec, typename TPos>
inline typename Value<TRankDictionary>::Type
getValue(SentinelRankDictionary<TRankDictionary, TSpec > const & dict, TPos pos)
{
    SEQAN_ASSERT_NEQ(isSentinelPosition(dict, pos), true);
    return getValue(getFibre(dict, FibreRankDictionary()), pos);
}

template <typename TRankDictionary, typename TSpec, typename TPos>
inline typename Value<TRankDictionary>::Type
getValue(SentinelRankDictionary<TRankDictionary, TSpec > & dict, TPos pos)
{
    SEQAN_ASSERT_NEQ(isSentinelPosition(dict, pos), true);
    return getValue(getFibre(dict, FibreRankDictionary()), pos);
}
//
//// ----------------------------------------------------------------------------
//// Function getFibre()
//// ----------------------------------------------------------------------------
//
///**
//.Function.SentinelRankDictionary#getFibre:
//..class:Class.SentinelRankDictionary
//..summary:Returns a specific fibre of a dictionary.
//..signature:getFibre(dictionary, fibreTag)
//..class:Class.RankDictionary
//..cat:Index
//..param.dictionary:The dictionary holding the fibre.
//...type:Class.RankDictionary
//..param.fibreTag:A tag that identifies the @Metafunction.Fibre@.
//...type:Tag.SentinelRankDictionary Fibres
//..returns:A reference to the @Metafunction.Fibre@ object.
//..include:seqan/index.h
//*/
//template <typename TRankDictionary, typename TSpec>
//inline typename Fibre<SentinelRankDictionary<TRankDictionary, TSpec>, FibreRankDictionary>::Type &
//getFibre(SentinelRankDictionary<TRankDictionary, TSpec> & dict, FibreRankDictionary)
//{
//    return dict.rankDictionary;
//}
//
//template <typename TRankDictionary, typename TSpec>
//inline typename Fibre<SentinelRankDictionary<TRankDictionary, TSpec>, FibreRankDictionary>::Type const &
//getFibre(SentinelRankDictionary<TRankDictionary, TSpec> const & dict, FibreRankDictionary)
//{
//    return dict.rankDictionary;
//}
//
//template <typename TRankDictionary, typename TSpec>
//inline typename Fibre<SentinelRankDictionary<TRankDictionary, TSpec>, FibreSentinelPosition>::Type &
//getFibre(SentinelRankDictionary<TRankDictionary, TSpec> & dict, FibreSentinelPosition)
//{
//    return dict.sentinelPosition;
//}
//
//template <typename TRankDictionary, typename TSpec>
//inline typename Fibre<SentinelRankDictionary<TRankDictionary, TSpec>, FibreSentinelPosition>::Type const &
//getFibre(SentinelRankDictionary<TRankDictionary, TSpec> const & dict, FibreSentinelPosition)
//{
//    return dict.sentinelPosition;
//}
//
//// ----------------------------------------------------------------------------
//// Function getSentinelSubstitute()
//// ----------------------------------------------------------------------------
//
///**
//.Function.SentinelRankDictionary#getSentinelSubstitute
//..class:Class.SentinelRankDictionary
//..summary:Returns the character used to substitute the sentinel sign.
//..signature:getSentinelSubstitute(dictionary)
//..param.dictionary:The dictionary.
//...type:Class.RankDictionary
//..include:seqan/index.h
//*/
//template <typename TRankDictionary, typename TSpec>
//inline typename Value<TRankDictionary>::Type
//getSentinelSubstitute(SentinelRankDictionary<TRankDictionary, TSpec> const & dict /*tag*/)
//{
//    return dict.sentinelSubstitute;
//}
//
//// ----------------------------------------------------------------------------
//// Function setSentinelSubstitute()
//// ----------------------------------------------------------------------------
//
///**
//.Function.SentinelRankDictionary#setSentinelSubstitute
//..class:Class.SentinelRankDictionary
//..summary:Sets the character used to substitute the sentinel sign.
//..signature:setSentinelSubstitute(dictionary, character)
//..param.dictionary:The dictionary.
//...type:Class.RankDictionary
//..param.character:The sentinel substitute.
//..include:seqan/index.h
//*/
//
//template <typename TRankDictionary, typename TSpec, typename TChar>
//inline void setSentinelSubstitute(SentinelRankDictionary<TRankDictionary, TSpec> & dict,
//                                  TChar sentinelSubstitute)
//{
//    dict.sentinelSubstitute = sentinelSubstitute;
//}
//
//// ----------------------------------------------------------------------------
//// Function setSentinelPosition()
//// ----------------------------------------------------------------------------
//
///**
//.Function.SentinelRankDictionary#setSentinelPosition
//..class:Class.SentinelRankDictionary
//..summary:Sets the sentinel position..
//..signature:setSentinelPosition(dictionary, pos)
//..param.dictionary:The dictionary.
//...type:Class.RankDictionary
//..param.pos:The sentinel position.
//..include:seqan/index.h
//*/
//
//template <typename TRankDictionary, typename TSpec, typename TPos>
//inline void setSentinelPosition(SentinelRankDictionary<TRankDictionary, TSpec> & dict,
//                                TPos const & position)
//{
//    dict.sentinelPosition = position;
//}

// ----------------------------------------------------------------------------
// Function getRank()
// ----------------------------------------------------------------------------

/**
.Function.SentinelRankDictionary#countOccurrences
..class:Class.SentinelRankDictionary
..summary:Returns the number of occurrences of a specified character from the start
to a specified position.
..signature:getRank(dictionary, pos, character)
..param.dictionary:The dictionary.
...type:Class.RankDictionary
..param.character:The character.
..param.pos:The position (which is included in the counting).
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGT";
RankDictionary<String<Dna5> > dictionary(genome);

std::cout << getRank(dictionary, 3, 'a') << std::endl; // 1
std::cout << getRank(dictionary, 4, 'a') << std::endl; // 2
*/

template <typename TRankDictionary, typename TPos, typename TChar>
inline typename Size<SentinelRankDictionary<TRankDictionary, Sentinel> >::Type
getRank(SentinelRankDictionary<TRankDictionary, Sentinel> const & dict, TPos pos, TChar character)
{
    typedef SentinelRankDictionary<TRankDictionary, Sentinel>   TSentinelRankDictionary;
    typedef typename Size<TSentinelRankDictionary>::Type        TSize;

    TSize rank = getRank(getFibre(dict, FibreRankDictionary()), pos, character);

    if (ordEqual(getSentinelSubstitute(dict), character) && pos >= dict.sentinelPosition)
         --rank;

    return rank;
}

template <typename TRankDictionary, typename TPos, typename TChar>
inline  typename Size<SentinelRankDictionary<TRankDictionary, Sentinels> >::Type
getRank(SentinelRankDictionary<TRankDictionary, Sentinels> const & dict, TPos pos, TChar character)
{
    typedef SentinelRankDictionary<TRankDictionary, Sentinels>  TSentinelRankDictionary;
    typedef typename Size<TSentinelRankDictionary>::Type        TSize;

    TSize rank = getRank(getFibre(dict, FibreRankDictionary()), pos, character);

    if (ordEqual(getSentinelSubstitute(dict), character))
        return rank - getRank(getFibre(dict, FibreSentinelPosition()), pos);

    return rank;
}

// ----------------------------------------------------------------------------
// Function _computeBwtLength()
// ----------------------------------------------------------------------------

// This function computes the length of the bwt string.
template <typename TText>
typename Size<TText>::Type _computeBwtLength(TText const & text)
{
    return length(text) + 1;
}

// This function computes the length of the bwt string.
template <typename TText, typename TSetSpec>
typename Size<StringSet<TText, TSetSpec> >::Type _computeBwtLength(StringSet<TText, TSetSpec> const & text)
{
    return lengthSum(text) + countSequences(text);
}

// ----------------------------------------------------------------------------
// Function _createBwt()
// ----------------------------------------------------------------------------

// This function computes the BWT of a text. Note that the sentinel sign is substituted and its position stored.
template <typename TRankDictionary, typename TSpec, typename TBwt, typename TText, typename TSA>
inline void
_createBwtFromSA(SentinelRankDictionary<TRankDictionary, TSpec> & dict, TBwt & bwt, TText const & text, TSA const & sa)
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
            assignValue(bwtIt, dict.sentinelSubstitute);
            dict.sentinelPosition = bwtIt - begin(bwt, Standard());
        }
    }
}

// This function computes the BWT of a text. Note that the sentinel sign is substituted and its position stored.
template <typename TRankDictionary, typename TSpec, typename TBwt, typename TText, typename TSetSpec, typename TSA>
inline void
_createBwtFromSA(SentinelRankDictionary<TRankDictionary, TSpec> & dict, TBwt & bwt, StringSet<TText, TSetSpec> const & text, TSA const & sa)
{
    typedef typename Value<TSA>::Type                       TSAValue;
    typedef typename Size<TSA>::Type                        TSize;
    typedef typename Iterator<TSA const, Standard>::Type    TSAIter;
    typedef typename Iterator<TBwt, Standard>::Type         TBwtIter;

    TSize seqNum = countSequences(text);
    TSize totalLen = lengthSum(text);

    resize(dict.sentinelPosition, seqNum + totalLen, Exact());
    
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
            setValue(dict.sentinelPosition, bwtIt - bwtItBeg, false);
        }
        else
        {
            assignValue(bwtIt, dict.sentinelSubstitute);
            setValue(dict.sentinelPosition, bwtIt - bwtItBeg, true);
        }
    }

    // Update the auxiliary RankDictionary of sentinel positions.
    updateRanks(dict.sentinelPosition);
}

// ----------------------------------------------------------------------------
// Function _createRankDictionary()
// ----------------------------------------------------------------------------

template <typename TRankDictionary, typename TSpec, typename TBwt, typename TPrefixSumTable>
inline void
_createRankDictionary(SentinelRankDictionary<TRankDictionary, TSpec> & dict,
                      TBwt const & bwt,
                      TPrefixSumTable & /* prefixSumTable */)
{
    createRankDictionary(getFibre(dict, FibreRankDictionary()), bwt);
}

template <typename TValue, typename TSpec, typename TBwt, typename TPrefixSumTable>
inline void
_createRankDictionary(SentinelRankDictionary<RankDictionary<WaveletTree<TValue> >, TSpec> & dict,
                      TBwt const & bwt,
                      TPrefixSumTable & prefixSumTable)
{
    createRankDictionary(getFibre(dict, FibreRankDictionary()), bwt, prefixSumTable);
}

// ----------------------------------------------------------------------------
// Function createSentinelRankDictionary()
// ----------------------------------------------------------------------------

/**
.Function.SentinelRankDictionary#createSentinelRankDictionary
..class:Class.SentinelRankDictionary
..summary:This functions creates the dictionary structure.
..signature:void createSentinelRankDictionary(dictionary, text)
..param.dictionary:The dictionary.
...type:Class.RankDictionary.
..param.text:A text to be represented by the dictionary.
...type:Class.String
..include:seqan/index.h
*/

template <typename TRankDictionary, typename TSpec, typename TPrefixSumTable, typename TText, typename TSA>
inline void
createSentinelRankDictionary(SentinelRankDictionary<TRankDictionary, TSpec> & dict,
                             TPrefixSumTable & prefixSumTable,
                             TText const & text,
                             TSA const & sa)
{
    typedef SentinelRankDictionary<TRankDictionary, TSpec>          TSentinelRankDictionary;
    typedef typename Value<TSentinelRankDictionary>::Type           TValue;
    typedef String<TValue>                                          TBwt;

    createPrefixSumTable(prefixSumTable, text);

    dict.sentinelSubstitute = determineSentinelSubstitute(prefixSumTable);

    TBwt bwt;
    resize(bwt, _computeBwtLength(text), Exact());
    _createBwtFromSA(dict, bwt, text, sa);

    _createRankDictionary(dict, bwt, prefixSumTable);

    insertSentinels(prefixSumTable, countSequences(text));
}

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

/**
.Function.SentinelRankDictionary#open
..class:Class.SentinelRankDictionary
..summary:This functions loads a dictionary from disk.
..signature:open(dictionary, fileName [, openMode])
..param.dictionary:The dictionary.
...type:Class.SentinelRankDictionary
..param.fileName:C-style character string containing the file name.
..param.openMode:The combination of flags defining how the file should be opened.
...remarks:To open a file read-only, write-only or to read and write use $OPEN_RDONLY$, $OPEN_WRONLY$, or $OPEN_RDWR$.
...remarks:To create or overwrite a file add $OPEN_CREATE$.
...remarks:To append a file if existing add $OPEN_APPEND$.
...remarks:To circumvent problems, files are always opened in binary mode.
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A $bool$ which is $true$ on success.
..include:seqan/index.h
*/

template <typename TRankDictionary, typename TSpec>
inline bool
_open(SentinelRankDictionary<TRankDictionary, TSpec> & dict, const char * fileName, int openMode)
{
    String<char> name;

    name = fileName;
    append(name, ".drs");
    if (!open(dict.sentinelSubstitute, toCString(name), openMode)) return false;

    name = fileName;
    append(name, ".drp");
    if (!open(dict.sentinelPosition, toCString(name), openMode)) return false;

    return true;
}

template <typename TRankDictionary, typename TSpec>
inline bool open(SentinelRankDictionary<TRankDictionary, TSpec> & dict, const char * fileName, int openMode)
{
    if (!open(getFibre(dict, FibreRankDictionary()), fileName, openMode)) return false;
    if (!_open(dict, fileName, openMode)) return false;

    return true;
}

template <typename TRankDictionary, typename TSpec>
inline bool open(SentinelRankDictionary<TRankDictionary, TSpec> & dict, const char * fileName)
{
    return open(dict, fileName, DefaultOpenMode<SentinelRankDictionary<TRankDictionary, TSpec> >::VALUE);
}

// ----------------------------------------------------------------------------
// Function save()
// ----------------------------------------------------------------------------

/**
.Function.SentinelRankDictionary#save
..class:Class.SentinelRankDictionary
..summary:This functions saves a dictionary to disk.
..signature:open(dictionary, fileName [, openMode])
..param.dictionary:The dictionary.
...type:Class.SentinelRankDictionary
..param.fileName:C-style character string containing the file name.
..param.openMode:The combination of flags defining how the file should be opened.
...remarks:To open a file read-only, write-only or to read and write use $OPEN_RDONLY$, $OPEN_WRONLY$, or $OPEN_RDWR$.
...remarks:To create or overwrite a file add $OPEN_CREATE$.
...remarks:To append a file if existing add $OPEN_APPEND$.
...remarks:To circumvent problems, files are always opened in binary mode.
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A $bool$ which is $true$ on success.
..include:seqan/index.h
*/

template <typename TRankDictionary, typename TSpec>
inline bool
_save(SentinelRankDictionary<TRankDictionary, TSpec> const & dict, const char * fileName, int openMode)
{
    String<char> name;

    name = fileName;
    append(name, ".drs");
    if (!save(dict.sentinelSubstitute, toCString(name), openMode)) return false;

    name = fileName;
    append(name, ".drp");
    if (!save(dict.sentinelPosition, toCString(name), openMode)) return false;
    
    return true;
}

template <typename TRankDictionary, typename TSpec>
inline bool save(SentinelRankDictionary<TRankDictionary, TSpec> const & dict, const char * fileName, int openMode)
{
    if (!save(getFibre(dict, FibreRankDictionary()), fileName, openMode)) return false;
    if (!_save(dict, fileName, openMode)) return false;

    return true;
}

template <typename TRankDictionary, typename TSpec>
inline bool save(SentinelRankDictionary<TRankDictionary, TSpec> const & dict, const char * fileName)
{
    return save(dict, fileName, DefaultOpenMode<SentinelRankDictionary<TRankDictionary, TSpec> >::VALUE);
}

}
#endif  // INDEX_FM_SENTINEL_RANK_DICTIONARY_H_

