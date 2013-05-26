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
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================

#ifndef INDEX_FM_H
#define INDEX_FM_H

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// WT = WaveletTree
/**
.Tag.WT
..summary:Tag that specifies the @Spec.FMIndex@ to use a wavelet tree as the occurrence table.
..cat:Index
*/

template <typename TSpec = void>
class WT;

template <typename TSpec = void>
class TL;

template <typename TOccSpec = WT<>, typename TSpec = void>
class FMIndex;

// ============================================================================
// Tags
// ============================================================================

// FM index fibres
/**
.Tag.FM Index Fibres
..summary:Tag to select a specific fibre of a @Spec.FMIndex@.
..remarks:These tags can be used to get @Metafunction.Fibre.Fibres@ of a FM index.
..cat:Index

..tag.FibrePrefixSum:The prefix sum table of the index.
..tag.FibreSA:The compressed suffix array of the text.
..tag.FibreText:The original text of the index.
..tag.FibreLF:The lf table.
..tag.FibreSALF:The lf table as well as the compressed suffix array.
...remarks:This tag can only be used with the functions @Function.indexRequire@ or @Function.indexSupplied@.

..see:Metafunction.Fibre
..see:Function.getFibre
..include:seqan/index_fm.h
*/

struct FibreText_;
struct FibreSA_;
struct FibreTempSA_;
struct FibreLF_;
struct FibreSALF_;

typedef Tag<FibreText_> const           FibreText;
typedef Tag<FibreSA_> const             FibreSA;
typedef Tag<FibreTempSA_> const         FibreTempSA;
typedef Tag<FibreLF_> const             FibreLF;
typedef Tag<FibreSALF_> const           FibreSALF;

// ----------------------------------------------------------------------------
// Tag FinderFMIndex
// ----------------------------------------------------------------------------

struct FinderFMIndex_;
typedef Tag<FinderFMIndex_> FinderFMIndex;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Fibre
// ----------------------------------------------------------------------------

/*
.Metafunction.Fibre:
..summary:Type of a specific FMIndex member (fibre).
..signature:Fibre<Index<TText, FMIndex<TSentinelRankDictionary, TSpec> >, TFibreSpec>::Type
..class:Spec.FMIndex
..cat:Index
..param.TText:The text type.
...type:Class.String
...type:Class.StringSet
..param.TSentinelRankDictionary:The type of the sentinel rank dictionary.
...type:Class.SentinelRankDictionary
..param.TSpec:Tag to specify a certain variant of the FM index.
...default;$void$
..param.TFibreSpec:Tag to specify the fibre.
...type:Tag.FM Index Fibres
..returns:Fibre type.
..remarks:Some containers, such as @Spec.FMIndex@, can be seen as a bundle consisting of various fibres. Because not 
every table is a fibre we did not call them tables, however, in many cases one can think of fibres as tables. The 
fibre interface was designed to unify the access to the members of the different fibres.
To get a reference or the type of a specific fibre use @Function.getFibre@ or @Metafunction.Fibre@.		
..include:seqan/index.h
*/

template <typename TText, typename TOccSpec, typename TSpec>
struct Fibre<Index<TText, FMIndex<TOccSpec, TSpec> >, FibreLF>
{
    typedef LfTable<TText, TSpec>   Type;
};

template <typename TText, typename TOccSpec, typename TSpec>
struct Fibre<Index<TText, FMIndex<TOccSpec, TSpec> >, FibreSA>
{
    typedef CompressedSA<TText, TSpec>  Type;
};

template <typename TText, typename TOccSpec, typename TSpec>
struct Fibre<Index<TText, FMIndex<TOccSpec, TSpec> >, FibreTempSA>
{
    typedef Index<TText, FMIndex<TOccSpec, TSpec> >                 TIndex_;
    typedef typename SAValue<TIndex_>::Type                         TSAValue_;

    // TODO(esiragusa): Revert this to External, now it is causing problems on device code.
//    typedef String<TSAValue_, External<ExternalConfigLarge<> > >    Type;
    typedef String<TSAValue_, typename DefaultIndexStringSpec<TText>::Type>     Type;
};

// ----------------------------------------------------------------------------
// Metafunction DefaultFinder
// ----------------------------------------------------------------------------

template <typename TText, typename TOccSpec, typename TSpec>
struct DefaultFinder<Index<TText, FMIndex<TOccSpec, TSpec> > >
{
    typedef FinderFMIndex Type;
};

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class FMIndex 
// ----------------------------------------------------------------------------

/**
.Spec.FMIndex:
..summary:An index based on the Burrows-Wheeler transform.
..cat:Index
..general:Class.Index
..signature:Index<TText, FMIndex<TOccSpec, TSpec> >
..param.TText:The text type.
...type:Class.String
...type:Class.StringSet
..param.TOccSpec:Occurrence table specialisation. 
...type:Tag.WT
...type:Tag.SBM
...remarks:The tags are really shortcuts for the different @Class.SentinelRankDictionary@s
...default:Tag.WT
..param.TSpec:FM index specialisation.
...default:void
..include:seqan/index.h
*/

template <typename TText, typename TOccSpec, typename TSpec>
class Index<TText, FMIndex<TOccSpec, TSpec> >
{
public:
    typename Member<Index, FibreText>::Type         text;
    typename Fibre<Index, FibreLF>::Type            lfTable;
    typename Fibre<Index, FibreSA>::Type            compressedSA;

//  TODO(esiragusa): Remove bwtLength and move compressionFactor into CSA.
    typename Size<TText>::Type                      bwtLength;
    unsigned                                        compressionFactor;

    // NOTE(esiragusa): NVCC cyclic SEQAN_FUNC problem.
    Index() {};

//    Index() :
//        bwtLength(0),
//        compressionFactor(10)
//    {}

    // TODO(esiragusa): Move compression factor inside CompressedSA.
    Index(TText & text, unsigned compressionFactor = 10) :
        text(text),
        bwtLength(_computeBwtLength(text)),
        compressionFactor(compressionFactor)
    {}
};

// ----------------------------------------------------------------------------
// Class FmIndexInfo_ 
// ----------------------------------------------------------------------------
//  TODO(esiragusa): Move compressionFactor into CSA and remove FmIndexInfo_.

// Stores the information about an FM index file bundle and is written to the .fma file.

#ifdef PLATFORM_WINDOWS
    #pragma pack(push,1)
#endif

struct FmIndexInfo_
{
    // The compression factor.
    __uint32 compressionFactor;
    // The sizeof(TSAEntry) values for suffix array entries.
    __uint32 sizeOfSAEntry;
    // The length of the BWT.
    __uint64 bwtLength;
}

#ifndef PLATFORM_WINDOWS
    __attribute__((packed))
#endif
    ;
#ifdef PLATFORM_WINDOWS
      #pragma pack(pop)
#endif

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

// Already documented
template <typename TText, typename TOccSpec, typename TSpec>
inline void clear(Index<TText, FMIndex<TOccSpec, TSpec> > & index)
{
    clear(getFibre(index, FibreLF()));
    clear(getFibre(index, FibreSA()));
    index.bwtLength = 0;
    index.compressionFactor = 0;
}

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

// This function checks whether the index is empty. Its already documented.
template <typename TText, typename TOccSpec, typename TSpec>
inline bool empty(Index<TText, FMIndex<TOccSpec, TSpec> > const & index)
{
    return empty(getFibre(index, FibreLF())) && empty(getFibre(index, FibreSA()));
}

// ----------------------------------------------------------------------------
// Function getFibre()
// ----------------------------------------------------------------------------

/**
.Function.FMIndex#getFibre:
..summary:Returns a specific fibre of a fm index.
..signature:getFibre(index, fibreTag)
..class:Spec.FMIndex
..cat:Index
..param.index:The index holding the fibre.
...type:Spec.FMIndex
..param.fibreTag:A tag that identifies the @Metafunction.Fibre@.
...type:Tag.FM Index Fibres
..returns:A reference to the @Metafunction.Fibre@ object.
..include:seqan/index.h
*/

template <typename TText, typename TOccSpec, typename TSpec>
SEQAN_FUNC typename Fibre<Index<TText, FMIndex<TOccSpec, TSpec> >, FibreLF >::Type &
getFibre(Index<TText, FMIndex<TOccSpec, TSpec> > & index, FibreLF /*tag*/)
{
    return index.lfTable;
}

template <typename TText, typename TOccSpec, typename TSpec>
SEQAN_FUNC typename Fibre<Index<TText, FMIndex<TOccSpec, TSpec> >, FibreLF >::Type const &
getFibre(Index<TText, FMIndex<TOccSpec, TSpec> > const & index, FibreLF /*tag*/)
{
    return index.lfTable;
}

template <typename TText, typename TOccSpec, typename TSpec>
SEQAN_FUNC typename Fibre<Index<TText, FMIndex<TOccSpec, TSpec> >, FibreSA >::Type &
getFibre(Index<TText, FMIndex<TOccSpec, TSpec> > & index, FibreSA /*tag*/)
{
    return index.compressedSA;
}

template <typename TText, typename TOccSpec, typename TSpec>
SEQAN_FUNC typename Fibre<Index<TText, FMIndex<TOccSpec, TSpec> >, FibreSA >::Type const &
getFibre(Index<TText, FMIndex<TOccSpec, TSpec> > const & index, FibreSA /*tag*/)
{
    return index.compressedSA;
}

// ----------------------------------------------------------------------------
// Function indexLF()
// ----------------------------------------------------------------------------

template <typename TText, typename TOccSpec, typename TSpec>
SEQAN_FUNC typename Fibre<Index<TText, FMIndex<TOccSpec, TSpec> >, FibreLF >::Type &
indexLF(Index<TText, FMIndex<TOccSpec, TSpec> > & index)
{
    return getFibre(index, FibreLF());
}

template <typename TText, typename TOccSpec, typename TSpec>
SEQAN_FUNC typename Fibre<Index<TText, FMIndex<TOccSpec, TSpec> >, FibreLF >::Type const &
indexLF(Index<TText, FMIndex<TOccSpec, TSpec> > const & index)
{
    return getFibre(index, FibreLF());
}

// ----------------------------------------------------------------------------
// Function toSuffixPosition()
// ----------------------------------------------------------------------------

/**
.Function.FMIndex#toSuffixPosition
..class:Spec.FMIndex
..summary:This function computes the position of a specified position in the suffix array (additionally containing 
entries for the sentinels. The returned position correspond to the suffix array of the original text without sentinels.
..signature:toSuffixPosition(fmIndex, pos, offset)
..param.fmIndex:The FM index.
...type:Spec.FMIndex
..param.pos:The position in the suffix array of the fm index (with sentinels).
...type:Concept.UnsignedIntegerConcept
..param.offset:The number of sequences in the original text.
...type:Concept.UnsignedIntegerConcept
..include:seqan/index.h
*/

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TPos, typename TSize>
inline typename SAValue<Index<TText, FMIndex<TOccSpec, TIndexSpec > > >::Type
toSuffixPosition(Index<TText, FMIndex<TOccSpec, TIndexSpec > > & index, TPos i, TSize offset)
{
    SEQAN_ASSERT_GEQ(suffixLength(i, index), offset);
    setSeqOffset(i, suffixLength(i, index) - offset);
    return i;
}

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TPos, typename TSize>
inline typename SAValue<Index<TText, FMIndex<TOccSpec, TIndexSpec > > const>::Type
toSuffixPosition(Index<TText, FMIndex<TOccSpec, TIndexSpec > > const & index, TPos i, TSize offset)
{
    SEQAN_ASSERT_GEQ(suffixLength(i, index), offset);
    setSeqOffset(i, suffixLength(i, index) - offset);
    return i;
}

// ----------------------------------------------------------------------------
// Function indexCreate()
// ----------------------------------------------------------------------------

/**
.Function.FMIndex#indexCreate
..summary:Creates a specific @Metafunction.Fibre@.
..signature:indexCreate(index, fibreTag)
..param.index:The index to be created.
...type:Spec.FMIndex
..param.fibreTag:The fibre of the index to be computed.
...type:Tag.FM Index Fibres.tag.FibreSALF
*/

template <typename TText, typename TIndexSpec, typename TSpec>
inline bool indexCreate(Index<TText, FMIndex<TIndexSpec, TSpec> > & index, FibreSALF const)
{
    typedef Index<TText, FMIndex<TIndexSpec, TSpec> >   TIndex;
    typedef typename Fibre<TIndex, FibreTempSA>::Type   TTempSA;
    typedef typename Size<TIndex>::Type                 TSize;
    
    TText const & text = getFibre(index, FibreText());

    if (empty(text))
        return false;

    TTempSA tempSA;

    // Create the full SA.
    resize(tempSA, length(text), Exact());
    createSuffixArray(tempSA, text, Skew7());

    // Create the LF table.
    createLFTable(getFibre(index, FibreLF()), text, tempSA);

    // Set the FMIndex LfTable as the CompressedSA LfTable.
    setLfTable(getFibre(index, FibreSA()), getFibre(index, FibreLF()));

    // Create the compressed SA.
    TSize numSentinel = countSequences(text);
    createCompressedSa(getFibre(index, FibreSA()), tempSA, index.compressionFactor, numSentinel);

    return true;
}

template <typename TText, typename TIndexSpec, typename TSpec>
inline bool indexCreate(Index<TText, FMIndex<TIndexSpec, TSpec> > & index)
{
    return indexCreate(index, FibreSALF());
}

// ----------------------------------------------------------------------------
// Function indexSupplied()
// ----------------------------------------------------------------------------

/**
.Function.FMIndex#indexSupplied:
..summary:Returns whether a specific @Metafunction.Fibre@ is present.
..param.fibreTag:
...type:Tag.FM Index Fibres
*/
template <typename TText, typename TIndexSpec, typename TSpec>
SEQAN_FUNC bool indexSupplied(Index<TText, FMIndex<TIndexSpec, TSpec > > & index, FibreSALF const)
{
    return !(empty(getFibre(index, FibreSA())) || empty(getFibre(index, FibreLF())));
}

template <typename TText, typename TIndexSpec, typename TSpec>
SEQAN_FUNC bool indexSupplied(Index<TText, FMIndex<TIndexSpec, TSpec > > const & index, FibreSALF const)
{
    return !(empty(getFibre(index, FibreSA())) || empty(getFibre(index, FibreLF())));
}

// ----------------------------------------------------------------------------
// Function _range()
// ----------------------------------------------------------------------------

// This function computes a range in the suffix array whose entries point to location
// in the text where the pattern occurs.

// TODO(esiragusa): Merge this with VSTree code.
template <typename TText, typename TOccSpec, typename TSpec, typename TPattern, typename TIter, typename TPairSpec>
inline void _range(Index<TText, FMIndex<TOccSpec, TSpec> > const & index, TPattern const & pattern,
                   Pair<TIter, TPairSpec> & range)
{
    typedef Index<TText, FMIndex<TOccSpec, TSpec> >                 TIndex;
    typedef typename Fibre<TIndex, FibreLF>::Type                   TLF;
    typedef typename Fibre<TLF, FibrePrefixSum>::Type               TPrefixSum;
    typedef typename Fibre<TLF, FibreValues>::Type                  TValues;
    typedef typename Value<TIndex>::Type                            TAlphabet;
    typedef typename ValueSize<TAlphabet>::Type                     TAlphabetSize;
    typedef typename Size<TIndex>::Type                             TSize;
    typedef typename Value<TPattern>::Type                      	TChar;

    TLF const & lfTable = getFibre(index, FibreLF());
    TPrefixSum const & prefixSumTable = getFibre(lfTable, FibrePrefixSum());

    if (empty(pattern))
    {
	    setPosition(range.i1, countSequences(index));
	    setPosition(range.i2, index.bwtLength);
    }

    TSize i = length(pattern) - 1;
    TChar letter = pattern[i];

    // Initilization.
    TAlphabetSize letterPosition = getCharacterPosition(prefixSumTable, letter);
    TSize sp = getPrefixSum(prefixSumTable, letterPosition);
    TSize ep = getPrefixSum(prefixSumTable, letterPosition + 1) - 1;

    // The search as proposed by Ferragina and Manzini.
    while ((sp <= ep) && (i > 0))
    {
        --i;
        letter = pattern[i];
        letterPosition = getCharacterPosition(prefixSumTable, letter);
        // TODO(esiragusa): Change this to getValue(lfTable) or something similar.
        TSize prefixSum = getPrefixSum(prefixSumTable, letterPosition);
        sp = prefixSum + _getRank(lfTable, sp - 1, letter);
        ep = prefixSum + _getRank(lfTable, ep, letter) - 1;
    }

    setPosition(range.i1, sp);
    setPosition(range.i2, ep + 1);
}

// ----------------------------------------------------------------------------
// Function _findFirstIndex()
// ----------------------------------------------------------------------------

// This function is used by the finder interface. It initializes the range of the finder.
template <typename TText, typename TPattern, typename TOccSpec, typename TSpec>
inline void
_findFirstIndex(Finder<Index<TText, FMIndex<TOccSpec, TSpec> >, FinderFMIndex> & finder,
		        TPattern const & pattern, FinderFMIndex const &)
{
    typedef Index<TText, FMIndex<TOccSpec, TSpec> >    TIndex;

    TIndex & index = haystack(finder);

    indexRequire(index, FibreSALF());
    setContainer(finder.range.i1, getFibre(container(finder), FibreSA()));
    setContainer(finder.range.i2, getFibre(container(finder), FibreSA()));

    _range(index, pattern, finder.range);
}

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

/**
.Function.Index#open
..class:Class.Index
..summary:This functions loads a dictionary from disk.
..signature:open(dictionary, fileName [, openMode])
..param.dictionary:The dictionary.
...type:Class.RankDictionary
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

// This function can be used to open a previously saved index.
template <typename TText, typename TOccSpec, typename TSpec>
inline bool open(Index<TText, FMIndex<TOccSpec, TSpec> > & index, const char * fileName, int openMode)
{
    String<char> name;

    String<FmIndexInfo_> infoString;

    name = fileName;    append(name, ".txt");
    if (!open(getFibre(index, FibreText()), toCString(name), openMode)) return false;

    name = fileName;    append(name, ".sa");
    if (!open(getFibre(index, FibreSA()), toCString(name), openMode)) return false;

    name = fileName;    append(name, ".lf");
    if (!open(getFibre(index, FibreLF()), toCString(name), openMode)) return false;

    name = fileName;    append(name, ".fma");
    if (!open(infoString, toCString(name), openMode)) return false;

    // Initialize index private members from info string.
    index.compressionFactor = infoString[0].compressionFactor;
    index.bwtLength = infoString[0].bwtLength;

    setLfTable(getFibre(index, FibreSA()), getFibre(index, FibreLF()));

    return true;
}

// This function can be used to open a previously saved index.
template <typename TText, typename TOccSpec, typename TSpec>
inline bool open(Index<TText, FMIndex<TOccSpec, TSpec> > & index, const char * fileName)
{
    return open(index, fileName, DefaultOpenMode<Index<TText, FMIndex<TOccSpec, TSpec> > >::VALUE);
}

// ----------------------------------------------------------------------------
// Function save()
// ----------------------------------------------------------------------------

/**
.Function.Index#save
..class:Class.Index
..summary:This functions saves an index to disk.
..signature:save(index, fileName [, openMode])
..param.index:The index.
...type:Class.RankDictionary
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

template <typename TText, typename TOccSpec, typename TSpec>
inline bool save(Index<TText, FMIndex<TOccSpec, TSpec> > const & index, const char * fileName, int openMode)
{
    String<char> name;

    typedef Index<TText, FMIndex<TOccSpec, TSpec> > TIndex;
    typedef typename Fibre<TIndex, FibreSA>::Type TSAFibre;
    typedef typename Value<TSAFibre>::Type TSAValue;

    String<FmIndexInfo_> infoString;
    FmIndexInfo_ info = { index.compressionFactor, sizeof(TSAValue), index.bwtLength };
    appendValue(infoString, info);

    name = fileName;    append(name, ".txt");
    if (!save(getFibre(index, FibreText()), toCString(name), openMode)) return false;

    name = fileName;    append(name, ".sa");
    if (!save(getFibre(index, FibreSA()), toCString(name), openMode)) return false;

    name = fileName;    append(name, ".lf");
    if (!save(getFibre(index, FibreLF()), toCString(name), openMode)) return false;

    name = fileName;    append(name, ".fma");
    if (!save(infoString, toCString(name), openMode)) return false;

    return true;
}

// This function can be used to save an index on disk.
template <typename TText, typename TOccSpec, typename TSpec>
inline bool save(Index<TText, FMIndex<TOccSpec, TSpec> > const & index, const char * fileName)
{
    return save(index, fileName, DefaultOpenMode<Index<TText, FMIndex<TOccSpec, TSpec> > >::VALUE);
}

}
#endif // INDEX_FM_H
