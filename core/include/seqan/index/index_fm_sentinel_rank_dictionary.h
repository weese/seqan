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
class SentinelRankDictionary;

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
class SentinelRankDictionary
{
public:
    typedef typename Value<SentinelRankDictionary>::Type TChar;
    typedef typename Fibre<SentinelRankDictionary<TRankDictionary, TSpec>, FibreSentinelPosition>::Type TSentinelPosition;

    TRankDictionary     rankDictionary;
    TSentinelPosition   sentinelPosition;
    TChar               sentinelSubstitute;

    SentinelRankDictionary() :
        sentinelSubstitute()
    {
        _initSentinelPosition(*this, 0u);
    }

    // TODO(singer): Use concept sequence when available
    template <typename TText>
    SentinelRankDictionary(TText const & text) :
        rankDictionary(text),
        sentinelSubstitute()
    {
        _initSentinelPosition(*this, length(text));
    }
};

// ==========================================================================
// Functions
// ==========================================================================

// ----------------------------------------------------------------------------
// Function _initSentinelPosition()
// ----------------------------------------------------------------------------

template <typename TRankDictionary, typename TSize>
inline void _initSentinelPosition(SentinelRankDictionary<TRankDictionary, Sentinel> & dictionary, TSize sentinels)
{
    setSentinelPosition(dictionary, sentinels);
}

template <typename TRankDictionary, typename TSize>
inline void _initSentinelPosition(SentinelRankDictionary<TRankDictionary, Sentinels> & dictionary, TSize sentinels)
{
//    resize(dictionary.sentinelPosition, sentinels, 0, Exact());
    resize(dictionary.sentinelPosition, sentinels, Exact());
}

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
inline void _clearSentinel(SentinelRankDictionary<TRankDictionary, Sentinels> & dictionary)
{
    clear(getFibre(dictionary, FibreSentinelPosition()));
}

template <typename TRankDictionary, typename TSpec>
inline void clear(SentinelRankDictionary<TRankDictionary, TSpec> & dictionary)
{
    clear(getFibre(dictionary, FibreRankDictionary()));
    _clearSentinel(dictionary);
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
inline bool isSentinelPosition(SentinelRankDictionary<TRankDictionary, Sentinel> const & dictionary, TPos pos)
{
    return dictionary.sentinelPosition == pos;
}

template <typename TRankDictionary, typename TPos>
inline bool isSentinelPosition(SentinelRankDictionary<TRankDictionary, Sentinels> const & dictionary, TPos pos)
{
    return getValue(getFibre(dictionary, FibreSentinelPosition()), pos);
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
inline bool empty(SentinelRankDictionary<TRankDictionary, TSpec> const & dictionary)
{
    return empty(getFibre(dictionary, FibreRankDictionary()));
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
getValue(SentinelRankDictionary<TRankDictionary, TSpec > const & dictionary, TPos pos)
{
    SEQAN_ASSERT_NEQ(isSentinelPosition(dictionary, pos), true);
    return getValue(getFibre(dictionary, FibreRankDictionary()), pos);
}

template <typename TRankDictionary, typename TSpec, typename TPos>
inline typename Value<TRankDictionary>::Type
getValue(SentinelRankDictionary<TRankDictionary, TSpec > & dictionary, TPos pos)
{
    SEQAN_ASSERT_NEQ(isSentinelPosition(dictionary, pos), true);
    return getValue(getFibre(dictionary, FibreRankDictionary()), pos);
}

// ----------------------------------------------------------------------------
// Function getFibre()
// ----------------------------------------------------------------------------

/**
.Function.SentinelRankDictionary#getFibre:
..class:Class.SentinelRankDictionary
..summary:Returns a specific fibre of a dictionary.
..signature:getFibre(dictionary, fibreTag)
..class:Class.RankDictionary
..cat:Index
..param.dictionary:The dictionary holding the fibre.
...type:Class.RankDictionary
..param.fibreTag:A tag that identifies the @Metafunction.Fibre@.
...type:Tag.SentinelRankDictionary Fibres
..returns:A reference to the @Metafunction.Fibre@ object.
..include:seqan/index.h
*/
template <typename TRankDictionary, typename TSpec>
inline typename Fibre<SentinelRankDictionary<TRankDictionary, TSpec>, FibreRankDictionary>::Type &
getFibre(SentinelRankDictionary<TRankDictionary, TSpec>& dictionary, FibreRankDictionary)
{
    return dictionary.rankDictionary;
}

template <typename TRankDictionary, typename TSpec>
inline typename Fibre<SentinelRankDictionary<TRankDictionary, TSpec>, FibreRankDictionary>::Type const &
getFibre(SentinelRankDictionary<TRankDictionary, TSpec> const & dictionary, FibreRankDictionary)
{
    return dictionary.rankDictionary;
}

template <typename TRankDictionary, typename TSpec>
inline typename Fibre<SentinelRankDictionary<TRankDictionary, TSpec>, FibreSentinelPosition>::Type &
getFibre(SentinelRankDictionary<TRankDictionary, TSpec> & dictionary, FibreSentinelPosition)
{
    return dictionary.sentinelPosition;
}

template <typename TRankDictionary, typename TSpec>
inline typename Fibre<SentinelRankDictionary<TRankDictionary, TSpec>, FibreSentinelPosition>::Type const &
getFibre(SentinelRankDictionary<TRankDictionary, TSpec> const & dictionary, FibreSentinelPosition)
{
    return dictionary.sentinelPosition;
}

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
getRank(SentinelRankDictionary<TRankDictionary, Sentinel> const & dictionary, TPos pos, TChar character)
{
    typedef SentinelRankDictionary<TRankDictionary, Sentinel>   TSentinelRankDictionary;
    typedef typename Size<TSentinelRankDictionary>::Type        TSize;

    TSize rank = getRank(getFibre(dictionary, FibreRankDictionary()), pos, character);

    if (ordEqual(getSentinelSubstitute(dictionary), character) && pos >= dictionary.sentinelPosition)
         --rank;

    return rank;
}

template <typename TRankDictionary, typename TPos, typename TChar>
inline  typename Size<SentinelRankDictionary<TRankDictionary, Sentinels> >::Type
getRank(SentinelRankDictionary<TRankDictionary, Sentinels> const & dictionary, TPos pos, TChar character)
{
    typedef SentinelRankDictionary<TRankDictionary, Sentinels>  TSentinelRankDictionary;
    typedef typename Size<TSentinelRankDictionary>::Type        TSize;

    TSize rank = getRank(getFibre(dictionary, FibreRankDictionary()), pos, character);

    if (ordEqual(getSentinelSubstitute(dictionary), character))
        return rank - getRank(getFibre(dictionary, FibreSentinelPosition()), pos);

    return rank;
}

// ----------------------------------------------------------------------------
// Function getSentinelSubstitute()
// ----------------------------------------------------------------------------

/**
.Function.SentinelRankDictionary#getSentinelSubstitute
..class:Class.SentinelRankDictionary
..summary:Returns the character used to substitute the sentinel sign.
..signature:getSentinelSubstitute(dictionary)
..param.dictionary:The dictionary.
...type:Class.RankDictionary
..include:seqan/index.h
*/
template <typename TRankDictionary, typename TSpec>
inline typename Value<TRankDictionary>::Type
getSentinelSubstitute(SentinelRankDictionary<TRankDictionary, TSpec> const & dictionary /*tag*/)
{
    return dictionary.sentinelSubstitute;
}

// ----------------------------------------------------------------------------
// Function setSentinelSubstitute()
// ----------------------------------------------------------------------------

/**
.Function.SentinelRankDictionary#setSentinelSubstitute
..class:Class.SentinelRankDictionary
..summary:Sets the character used to substitute the sentinel sign.
..signature:setSentinelSubstitute(dictionary, character)
..param.dictionary:The dictionary.
...type:Class.RankDictionary
..param.character:The sentinel substitute.
..include:seqan/index.h
*/

template <typename TRankDictionary, typename TSpec, typename TChar>
inline void setSentinelSubstitute(SentinelRankDictionary<TRankDictionary, TSpec> & dictionary,
                                  TChar sentinelSubstitute)
{
    dictionary.sentinelSubstitute = sentinelSubstitute;
}

// ----------------------------------------------------------------------------
// Function setSentinelPosition()
// ----------------------------------------------------------------------------

/**
.Function.SentinelRankDictionary#setSentinelPosition
..class:Class.SentinelRankDictionary
..summary:Sets the sentinel position..
..signature:setSentinelPosition(dictionary, pos)
..param.dictionary:The dictionary.
...type:Class.RankDictionary
..param.pos:The sentinel position.
..include:seqan/index.h
*/

template <typename TRankDictionary, typename TSpec, typename TPos>
inline void setSentinelPosition(SentinelRankDictionary<TRankDictionary, TSpec> & dictionary,
                                TPos const & position)
{
    dictionary.sentinelPosition = position;
}

// ----------------------------------------------------------------------------
// Function sentinelRankDictionaryCreate()
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

template <typename TRankDictionary, typename TSpec, typename TText, typename TSentinelSub, typename TSentinelPos>
inline void createSentinelRankDictionary(SentinelRankDictionary<TRankDictionary, TSpec> & dictionary,
                                         TText const & text,
                                         TSentinelSub const & sentinelSub,
                                         TSentinelPos const & sentinelPos)
{
    setSentinelSubstitute(dictionary, sentinelSub);
    setSentinelPosition(dictionary, sentinelPos);

    createRankDictionary(getFibre(dictionary, FibreRankDictionary()), text);
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

template <typename TRankDictionary>
inline bool _openSentinelInformation(SentinelRankDictionary<TRankDictionary, Sentinel> & dictionary,
                                     const char * fileName,
                                     int openMode)
{
    typedef typename Fibre<SentinelRankDictionary<TRankDictionary, Sentinel>, FibreSentinelPosition>::Type TSentinelString;
    typedef typename Value<SentinelRankDictionary<TRankDictionary, Sentinel> >::Type TChar;

    String<Pair<TChar, TSentinelString, Pack> > sentinelValues;

    String<char> name = fileName;
    append(name, ".dr");

    if (!open(sentinelValues, toCString(name), openMode) || empty(sentinelValues)) return false;

    dictionary.sentinelSubstitute = sentinelValues[0].i1;
    dictionary.sentinelPosition = sentinelValues[0].i2;
    
    return true;
}

template <typename TRankDictionary>
inline bool _openSentinelInformation(SentinelRankDictionary<TRankDictionary, Sentinels> & dictionary,
                                     const char * fileName,
                                     int openMode)
{
    typedef typename Value<SentinelRankDictionary<TRankDictionary, Sentinels> >::Type TChar;

    String<TChar> sentinelSub;

    String<char> name;

    name = fileName;
    append(name, ".drs");
    if (!open(sentinelSub, toCString(name), openMode)) return false;
    
    name = fileName;
    append(name, ".drp");
    if (!open(dictionary.sentinelPosition, toCString(name), openMode)) return false;
    
    if (empty(sentinelSub))
        return false;

    dictionary.sentinelSubstitute = sentinelSub[0];
    
    return true;
}

template <typename TRankDictionary, typename TSpec>
inline bool open(SentinelRankDictionary<TRankDictionary, TSpec> & dictionary, const char * fileName, int openMode)
{
    if (!open(getFibre(dictionary, FibreRankDictionary()), fileName, openMode)) return false;
    if (!_openSentinelInformation(dictionary, fileName, openMode)) return false;
    
    return true;
}


template <typename TRankDictionary, typename TSpec>
inline bool open(SentinelRankDictionary<TRankDictionary, TSpec> & dictionary, const char * fileName)
{
    return open(dictionary, fileName, DefaultOpenMode<SentinelRankDictionary<TRankDictionary, Sentinel> >::VALUE);
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

template <typename TRankDictionary>
inline bool _saveSentinelInformation(SentinelRankDictionary<TRankDictionary, Sentinel> const & dictionary,
                                     const char * fileName,
                                     int openMode)
{
    typedef typename Value<SentinelRankDictionary<TRankDictionary, Sentinel> >::Type TChar;
    typedef typename Fibre<SentinelRankDictionary<TRankDictionary, Sentinel>, FibreSentinelPosition>::Type TSentinelString;

    String<char> name;

    String<Pair<TChar, TSentinelString, Pack> > sentinelValues;
    appendValue(sentinelValues, Pair<TChar, TSentinelString>(dictionary.sentinelSubstitute, dictionary.sentinelPosition));
    
    name = fileName;
    append(name, ".dr");
    return save(sentinelValues, toCString(name), openMode);
}

template <typename TRankDictionary>
inline bool _saveSentinelInformation(SentinelRankDictionary<TRankDictionary, Sentinels> const & dictionary,
                                     const char * fileName,
                                     int openMode)
{
    typedef typename Value<SentinelRankDictionary<TRankDictionary, Sentinels> >::Type TChar;

    String<char> name;

    String<TChar> sentinelSub;
    appendValue(sentinelSub, dictionary.sentinelSubstitute);

    name = fileName;
    append(name, ".drs");
    if (!save(sentinelSub, toCString(name), openMode)) return false;

    name = fileName;
    append(name, ".drp");
    if (!save(dictionary.sentinelPosition, toCString(name), openMode)) return false;
    
    return true;
}

template <typename TRankDictionary, typename TSpec>
inline bool save(SentinelRankDictionary<TRankDictionary, TSpec> const & dictionary, const char * fileName, int openMode)
{
    if (!save(getFibre(dictionary, FibreRankDictionary()), fileName, openMode)) return false;
    if (!_saveSentinelInformation(dictionary, fileName, openMode)) return false;
    
    return true;
}

template <typename TRankDictionary, typename TSpec>
inline bool save(SentinelRankDictionary<TRankDictionary, TSpec> const & dictionary, const char * fileName)
{
    return save(dictionary, fileName, DefaultOpenMode<SentinelRankDictionary<TRankDictionary, Sentinel> >::VALUE);
}

template <typename TRankDictionary, typename TSpec>
inline bool save(SentinelRankDictionary<TRankDictionary, TSpec> & dictionary, const char * fileName)
{
    return save(dictionary, fileName, DefaultOpenMode<SentinelRankDictionary<TRankDictionary, Sentinel> >::VALUE);
}

}
#endif  // INDEX_FM_SENTINEL_RANK_DICTIONARY_H_

