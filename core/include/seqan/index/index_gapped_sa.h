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
// Author: Sascha Meiers <meiers@inf.fu-berlin.de>
// ==========================================================================
// The gapped suffix array
// ==========================================================================

#ifndef SEQAN_HEADER_INDEX_GAPPEDSA_H
#define SEQAN_HEADER_INDEX_GAPPEDSA_H


namespace SEQAN_NAMESPACE_MAIN
{

    // ============================================================================
    // THings that should go somewhere else
    // ============================================================================
    
    
    // ----------------------------------------------------------------------------
    // Metafunction Suffix                                                for Index
    // ----------------------------------------------------------------------------
    
    // general Index
    template <typename TText, typename TSpec>
    struct Suffix<Index<TText, TSpec> >
    {
        typedef typename Suffix<TText>::Type Type;
    };
    
    // general Index; const variant
    template <typename TText, typename TSpec>
    struct Suffix<Index<TText, TSpec> const>
    {
        typedef typename Suffix<TText const>::Type Type;
    };

    // ----------------------------------------------------------------------------
    // Function suffix()                                                  for Index
    // ----------------------------------------------------------------------------
    
    // general Index
    template <typename TText, typename TSpec, typename TPosBegin>
    inline typename Suffix<Index<TText, TSpec> >::Type
    suffix(Index<TText, TSpec> & t, TPosBegin pos_begin)
    {
        return suffix(indexText(t), pos_begin);
    }
    
    // general Index; const variant
    template <typename TText, typename TSpec, typename TPosBegin>
    inline typename Suffix<Index<TText, TSpec> const>::Type
    suffix(Index<TText, TSpec> const &t, TPosBegin pos_begin)
    {
        return suffix(indexText(t), pos_begin);
    }

    

    
// ============================================================================
// Forwards
// ============================================================================

template <typename T = void> struct IndexSa {};

    
// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// TODO(sascha): Remove this definition as soon as Ticket #1096 has been solved: http://trac.seqan.de/ticket/1096
template < typename TText, typename TSpec >
struct DefaultFinder< Index<TText, IndexSa<TSpec> > > {
    typedef EsaFindMlr Type;	// standard suffix array finder is mlr-heuristic
};
    
template <typename TShape = GenericShape, typename TSpec = void>
struct WithGaps {};


template <typename TText, typename TShape, typename TSpec>
class Index<TText, IndexSa< WithGaps<CyclicShape<TShape>, TSpec > > > :
    public Index<TText, IndexSa<TSpec> >
{
public:
    typedef Index<TText, IndexSa<TSpec> > TSuper;
    
    // derive           text, cargo, sa
    CyclicShape<TShape> shape;
    
    Index() {}
    
    Index(Index & other) : TSuper(other), shape(other.shape)
    {}
    
    Index(Index const & other) : TSuper(other), shape(other.shape)
    {}
    
    template <typename TText_>
    Index(TText_ & _text) : TSuper(_text)
    {}
    
    template <typename TText_>
    Index(TText_ const & _text) : TSuper(_text)
    {}
    
};


// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// FibreShape
// ----------------------------------------------------------------------------

template < typename TText, typename TShape, typename TSpec >
struct Fibre< Index<TText, IndexSa<WithGaps<CyclicShape<TShape>, TSpec> > >, FibreShape>
{
    typedef CyclicShape<TShape>	Type;
}; 

// ----------------------------------------------------------------------------
// DefaultIndexCreator
// ----------------------------------------------------------------------------

template < typename TText, typename TShape, typename TSpec>
struct DefaultIndexCreator<Index<TText, IndexSa<WithGaps<CyclicShape<TShape>, TSpec> > >, FibreSA> {
    typedef SAQSort Type;
};

// ----------------------------------------------------------------------------
// Metafunction Suffix                                          gapped SA Index
// ----------------------------------------------------------------------------

// gapped Index SA
template <typename TText, typename TCycShape, typename TSpec>
struct Suffix<Index<TText, IndexSa<WithGaps<TCycShape, TSpec> > > >
{
    typedef ModifiedString<typename Suffix<TText>::Type, ModCyclicShape<TCycShape> > Type;
};

// gapped Index SA; const variant
template <typename TText, typename TCycShape, typename TSpec>
struct Suffix<Index<TText, IndexSa<WithGaps<TCycShape, TSpec> > > const>
{
    typedef ModifiedString<typename Suffix<TText const>::Type, ModCyclicShape<TCycShape> > Type;
};


// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function indexCreate (FibreSA)
// ----------------------------------------------------------------------------
template <typename TText, typename TShape, typename TSpec, typename TSpecAlg>
inline bool indexCreate(Index<TText, IndexSa<WithGaps<CyclicShape<TShape>, TSpec> > > &index, FibreSA, TSpecAlg const alg)
{
    resize(indexSA(index), length(indexRawText(index)), Exact());
    createGappedSuffixArray(indexSA(index), indexText(index), indexShape(index), alg);
    return true;
}

// ----------------------------------------------------------------------------
// Function suffix()                                            gapped SA Index
// ----------------------------------------------------------------------------

// gapped Index SA
template <typename TText, typename TShape, typename TSpec, typename TPosBegin>
inline typename Suffix<Index<TText, IndexSa<WithGaps<CyclicShape<TShape>, TSpec> > > >::Type
suffix(Index<TText, IndexSa<WithGaps<CyclicShape<TShape>, TSpec> > > & t, TPosBegin pos_begin)
{
    typedef typename Suffix<Index<TText, IndexSa<WithGaps<CyclicShape<TShape>, TSpec> > > >::Type TModifiedSuffix;
    return TModifiedSuffix(suffix(indexText(t), pos_begin), indexShape(t));
}

// gapped Index SA; const variant
template <typename TText, typename TShape, typename TSpec, typename TPosBegin>
inline typename Suffix<Index<TText, IndexSa<WithGaps<CyclicShape<TShape>, TSpec> > > const>::Type
suffix(Index<TText, IndexSa<WithGaps<CyclicShape<TShape>, TSpec> > > const &t, TPosBegin pos_begin)
{
    typedef typename Suffix<Index<TText, IndexSa<WithGaps<CyclicShape<TShape>, TSpec> > > const>::Type TModifiedSuffix;
    return TModifiedSuffix(suffix(indexText(t), pos_begin), indexShape(t));
}

// ----------------------------------------------------------------------------
// Function infix()                                             gapped SA Index
// ----------------------------------------------------------------------------

// gapped Index SA
template <typename TText, typename TShape, typename TSpec, typename TPosBegin, typename TPosEnd>
inline typename Prefix<typename Suffix<Index<TText, IndexSa<WithGaps<CyclicShape<TShape>, TSpec> > > >::Type>::Type
infix(Index<TText, IndexSa<WithGaps<CyclicShape<TShape>, TSpec> > > &index, TPosBegin pos_begin, TPosEnd pos_end)
{
    return prefix(suffix(index, pos_begin), pos_end - pos_begin);
}

// gapped Index SA; const variant
template <typename TText, typename TShape, typename TSpec, typename TPosBegin, typename TPosEnd>
inline typename Prefix<typename Suffix<Index<TText, IndexSa<WithGaps<CyclicShape<TShape>, TSpec> > > const>::Type>::Type
infix(Index<TText, IndexSa<WithGaps<CyclicShape<TShape>, TSpec> > > const &index, TPosBegin pos_begin, TPosEnd pos_end)
{
    return prefix(suffix(index, pos_begin), pos_end - pos_begin);
}

    

// ----------------------------------------------------------------------------
// Function _findFirstIndex()                                    used in find()
// ----------------------------------------------------------------------------

template < typename TText, typename TCycShape, typename TSpec, typename TSpecFinder, typename TPattern >
inline void
_findFirstIndex(Finder< Index<TText, IndexSa<WithGaps<TCycShape, TSpec> > >, TSpecFinder > &finder,
                TPattern const &pattern,
                EsaFindMlr const)
{
    typedef Index<TText, IndexSa<WithGaps<TCycShape, TSpec> > > TIndex;
    TIndex &index = haystack(finder);
    indexRequire(index, EsaSA());
    finder.range = equalRangeSAIterator(index, indexSA(index), pattern); // handle the index, not the text!
}

// TODO(sascha): _findFirstIndex for Lcp and STree variant, see index_find_esa.h







///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                             experimental stuff                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

    // TODO(sascha): Whole BottomUp iterator!
    // TODO(sascha): Missing functions for Iterator, e.g. repLength
    

// copy 'n paste from index_sa_stree.h
template <typename TText, typename TShapeSpec, typename TIndexSpec, typename TSpec, typename TString, typename TSize>
inline bool _goDownString(Iter<Index<TText, IndexSa<WithGaps<TShapeSpec, TIndexSpec> > >, VSTree<TopDown<TSpec> > > & it,
                          TString const & pattern, TSize & lcp)
{
    typedef Index<TText, IndexSa<WithGaps<TShapeSpec, TIndexSpec> > >   TIndex;
    typedef typename Fibre<TIndex, FibreSA>::Type                       TSA;
    typedef typename Size<TSA const>::Type                              TSASize;
    typedef typename Iterator<TSA const, Standard>::Type                TSAIterator;
    typedef SearchTreeIterator<TSA const, SortedList>                   TSearchTreeIterator;
    
    if (_isLeaf(it, HideEmptyEdges()))
        return false;
    
    TIndex const & index = container(it);
    TSA const & sa = indexSA(index);
    // TText const & text = indexText(index);                                        // CHANGE(sascha): handle index instead of text
        
    TSAIterator saBegin = begin(sa, Standard()) + value(it).range.i1;
    TSASize saLen = isRoot(it) ? length(sa) : value(it).range.i2 - value(it).range.i1;
    TSearchTreeIterator node(saBegin, saLen);
    Pair<TSAIterator> range = _equalRangeSA(index, node, pattern, value(it).repLen); // CHANGE(sascha): handle index instead of text
    
    if (range.i1 >= range.i2)
        return false;
    
    // Save vertex descriptor.
    _historyPush(it);
    
    // Update range, lastChar and repLen.
    value(it).range.i1 = range.i1 - begin(sa, Standard());
    value(it).range.i2 = range.i2 - begin(sa, Standard());
    value(it).lastChar = back(pattern);
    value(it).repLen += length(pattern);
        
    lcp = length(pattern);
    
    return true;
}


template < typename TText, typename TShapeSpec, typename TIndexSpec, typename TSpec >
inline typename Prefix<typename Suffix<Index<TText, IndexSa<WithGaps<TShapeSpec, TIndexSpec> > > >::Type>::Type
representative(Iter<Index<TText, IndexSa<WithGaps<TShapeSpec, TIndexSpec> > >, VSTree<TSpec> > &it)
{
    return prefix(suffix(container(it), getOccurrence(it)), repLength(it));
}

template < typename TText, typename TShapeSpec, typename TIndexSpec, typename TSpec >
inline typename Prefix<typename Suffix<Index<TText, IndexSa<WithGaps<TShapeSpec, TIndexSpec> > > const>::Type>::Type
representative(Iter<Index<TText, IndexSa<WithGaps<TShapeSpec, TIndexSpec> > > const, VSTree<TSpec> > const &it)
{
    return prefix(suffix(container(it), getOccurrence(it)), repLength(it));
}


// seperate specializations for String and StringSet needed to avoid ambiguities
template <typename TPos, typename TText, typename TShapeSpec, typename TIndexSpec>
inline typename Size<Index<TText, IndexSa<WithGaps<TShapeSpec, TIndexSpec> > > const>::Type
suffixLength(TPos pos, Index<TText, IndexSa<WithGaps<TShapeSpec, TIndexSpec> > > const &index)
{
    return length(suffix(index, pos));
}
    
template <typename TPos, typename TString, typename TSSetSpec, typename TShapeSpec, typename TIndexSpec>
inline typename Size<Index<StringSet<TString, TSSetSpec>, IndexSa<WithGaps<TShapeSpec, TIndexSpec> > > const>::Type
suffixLength(TPos pos, Index<StringSet<TString, TSSetSpec>, IndexSa<WithGaps<TShapeSpec, TIndexSpec> > > const &index)
{
    return length(suffix(index, pos));
}

template <typename TPos, typename TText, typename TShapeSpec, typename TIndexSpec>
inline typename Size<Index<TText, IndexSa<WithGaps<TShapeSpec, TIndexSpec> > > >::Type
suffixLength(TPos pos, Index<TText, IndexSa<WithGaps<TShapeSpec, TIndexSpec> > > &index)
{
    return length(suffix(index, pos));
}
template <typename TPos, typename TString, typename TSSetSpec, typename TShapeSpec, typename TIndexSpec>
inline typename Size<Index<StringSet<TString, TSSetSpec>, IndexSa<WithGaps<TShapeSpec, TIndexSpec> > > >::Type
suffixLength(TPos pos, Index<StringSet<TString, TSSetSpec>, IndexSa<WithGaps<TShapeSpec, TIndexSpec> > > &index)
{
    return length(suffix(index, pos));
}


template <typename TText, typename TShapeSpec, typename TIndexSpec, typename TSpec, typename TDfsOrder>
inline bool _goDown(Iter<Index<TText, IndexSa<WithGaps<TShapeSpec, TIndexSpec> > >, VSTree<TopDown<TSpec> > > & it,
                    VSTreeIteratorTraits<TDfsOrder, True> const)
{
    typedef Index<TText, IndexSa<WithGaps<TShapeSpec, TIndexSpec> > >              TIndex;
    typedef typename Fibre<TIndex, FibreSA>::Type           TSA;
    typedef typename Size<TIndex>::Type                     TSASize;
    typedef typename Iterator<TSA const, Standard>::Type    TSAIterator;
    typedef SearchTreeIterator<TSA const, SortedList>       TSearchTreeIterator;
    typedef typename Value<TIndex>::Type                    TAlphabet;
    
    // TODO(esiragusa): use HideEmptyEdges()
    if (_isLeaf(it, HideEmptyEdges()))
        return false;
    
    TIndex const & index = container(it);
    TSA const & sa = indexSA(index);
    
    // TODO(esiragusa): check nodeHullPredicate
    
    Pair<typename Size<TIndex>::Type> saRange = range(it);
    
    // TODO(esiragusa): remove this check.
    if (saRange.i1 >= saRange.i2) return false;
    
    // Skip $-edges.
    while (suffixLength(saAt(saRange.i1, index), index) <= value(it).repLen)
    {
        // TODO(esiragusa): remove this check and ++saRange.i1 in loop.
        // Interval contains only $-edges.
        if (++saRange.i1 >= saRange.i2)
            return false;
    }
    
    // Get first and last characters in interval.
    // CHANGED(sascha): get char via the suffix instead of textAt
    TAlphabet cLeft  = value(suffix(index, saAt(saRange.i1,     index)), value(it).repLen);
    TAlphabet cRight = value(suffix(index, saAt(saRange.i2 - 1, index)), value(it).repLen);
    
    // Save vertex descriptor.
    _historyPush(it);
    
    // Update left range.
    value(it).range.i1 = saRange.i1;
    
    // Update right range.
    // NOTE(esiragusa): I should use ordLess(cLeft, cRight) but masai redefines it for Ns.
    if (ordValue(cLeft) != ordValue(cRight))
    {
        TSAIterator saBegin = begin(sa, Standard()) + saRange.i1;
        TSASize saLen = saRange.i2 - saRange.i1;
        TSearchTreeIterator node(saBegin, saLen);
        
        // CHANGED(sascha): handle index instead of text
        TSAIterator upperBound = _upperBoundSA(index, node, cLeft, value(it).repLen);
        
        value(it).range.i2 = upperBound - begin(sa, Standard());
    }
    else
    {
        value(it).range.i2 = saRange.i2;
    }
    
    // Update child repLen, lastChar.
    value(it).repLen++;
    value(it).lastChar = cLeft;
        
    return true;
}
    
  
    
    
template <typename TText, typename TShapeSpec, typename TIndexSpec, typename TSpec, typename TDfsOrder>
inline bool _goRight(Iter<Index<TText, IndexSa<WithGaps<TShapeSpec, TIndexSpec> > >, VSTree<TopDown<TSpec> > > & it,
                    VSTreeIteratorTraits<TDfsOrder, True> const)
{
    typedef Index<TText, IndexSa<WithGaps<TShapeSpec, TIndexSpec> > >              TIndex;
    typedef typename Fibre<TIndex, FibreSA>::Type           TSA;
    typedef typename Size<TIndex>::Type                     TSASize;
    typedef typename Iterator<TSA const, Standard>::Type    TSAIterator;
    typedef SearchTreeIterator<TSA const, SortedList>       TSearchTreeIterator;
    typedef typename Value<TIndex>::Type                    TAlphabet;
            
    if (isRoot(it))
        return false;

    TIndex const & index = container(it);
    TSA const & sa = indexSA(index);
            
    Pair<typename Size<TIndex>::Type> saRange;
    saRange.i1 = value(it).range.i2;
    saRange.i2 = (_isSizeInval(value(it).parentRight)) ? length(sa) : value(it).parentRight;

    if (saRange.i1 >= saRange.i2) return false;

    // Change repLen to parent repLen.
    value(it).repLen--;

    // TODO(esiragusa): don't check for empty edges (do it in goDown)
    // Skip $-edges.
    while (suffixLength(saAt(saRange.i1, index), index) <= value(it).repLen)
    {
        // Interval contains only $-edges.
        if (++saRange.i1 >= saRange.i2)
            return false;
    }

    // Get first and last characters in interval.
    // CHANGED(sascha): get char via the suffix instead of textAt
    TAlphabet cLeft  = value(suffix(index, saAt(saRange.i1,     index)), value(it).repLen);
    TAlphabet cRight = value(suffix(index, saAt(saRange.i2 - 1, index)), value(it).repLen);

    SEQAN_ASSERT_NEQ(ordValue(cLeft), ordValue(value(it).lastChar));

    // Update left range.
    value(it).range.i1 = saRange.i1;

    // Update right range.
    // NOTE(esiragusa): I should use ordLess(cLeft, cRight) but masai redefines it for Ns.
    if (ordValue(cLeft) != ordValue(cRight))
    {
        TSAIterator saBegin = begin(sa, Standard()) + saRange.i1;
        TSASize saLen = saRange.i2 - saRange.i1;
        TSearchTreeIterator node(saBegin, saLen);
        
        // CHANGED(sascha): handle index instead of text
        TSAIterator upperBound = _upperBoundSA(index, node, cLeft, value(it).repLen);
        
        value(it).range.i2 = upperBound - begin(sa, Standard());
    }
    else
    {
        value(it).range.i2 = saRange.i2;
    }

    // Update repLen, lastChar.
    value(it).repLen++;
    value(it).lastChar = cLeft;

    return true;
}
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
#ifndef SEQAN_HEADER_INDEX_BASE_H
//////////////////////////////////////////////////////////////////////////////
// needful forward declarations

	// suffix array construction specs
	struct Skew3;
	struct Skew7;
	struct LarssonSadakane;
	struct ManberMyers;
	struct SAQSort;
	struct QGramAlg;
	struct GappedSAQSort;

	// lcp table construction algorithms
	struct Kasai;
	struct KasaiOriginal;	// original, but more space-consuming algorithm

	// enhanced suffix array construction algorithms
	struct Childtab;
	struct Bwt;

/*
.Tag.Index Find Algorithm
..summary:Tag to specify the index search algorithm.
..remarks:These tag can be used to specify the @Function.find@ algorithm 
for @Class.Index@ based substring searches.
..cat:Index

..tag.EsaFindMlr:Binary search with mlr-heuristic.
...remarks:Exact string matching using a suffix array binary search with the mlr-heuristic.

..tag.EsaFindLcpe:Binary search using lcp values.
...remarks:Exact string matching using a suffix array binary search and a lcp-interval tree.

..tag.FinderSTree:Suffix tree search.
...remarks:Exact string matching using a suffix tree.

..see:Class.Finder
..see:Spec.IndexEsa
..see:Spec.IndexQGram
..include:seqan/index.h
*/

/*
 * @defgroup IndexFindAlgorithm Index Find Algorithm
 * 
 * @brief Tag to specify the index search algorithm.
 * 
 * @section Remarks
 * 
 * These tags can be used to specify the @link find @endlink algorithm for @link
 * Index @endlink based substring searches.
 * 
 * @see Finder
 * 
 * @tag IndexFindAlgorithm#FinderSTree
 * 
 * @brief Suffix tree search.
 * 
 * @section Remarks
 * 
 * Exact string matching using a suffix tree.
 * 
 * @tag IndexFindAlgorithm#PizzaChiliFinder
 * 
 * @brief Finds an occurrence in a @link Pizza & Chili Index @endlink index.
 * 
 * @section Remarks
 * 
 * The actual algorithm used for searching depends on the @link Pizza & Chili
 * Index Tags @endlink used.
 * 
 * @tag IndexFindAlgorithm#QGramFindLookup
 * 
 * @brief q-gram search. Finds q-grams in a @link IndexQGram @endlink index
 *        using the hash table.
 * 
 * @tag IndexFindAlgorithm#EsaFindLcpe
 * 
 * @brief Binary search using lcp values.
 * 
 * @section Remarks
 * 
 * Exact string matching using a suffix array binary search and a lcp-interval
 * tree.
 * 
 * @tag IndexFindAlgorithm#EsaFindMlr
 * 
 * @brief Binary search with mlr-heuristic.
 * 
 * @section Remarks
 * 
 * Exact string matching using a suffix array binary search with the mlr-
 * heuristic.
 */
	// finder tags
    struct FinderMlr_;     // simple Suffix Array finder with mlr-heuristic
    struct FinderLcpe_;    // Suffix Array finder using an enhanced LCP-Table
    struct FinderSTree_;    // Suffix Array finder using an enhanced LCP-Table

    typedef Tag<FinderMlr_> const EsaFindMlr;
    typedef Tag<FinderLcpe_> const EsaFindLcpe;
    typedef Tag<FinderSTree_> const FinderSTree;

	template <typename TSpec = void>
	struct IndexEsa {};




//////////////////////////////////////////////////////////////////////////////
// default table type

	template < typename TObject, typename TSpec, typename TFibre >
	struct Fibre< Index<TObject, TSpec>, Tag<TFibre> const > {
		typedef String< 
			typename Size< Index<TObject, TSpec> >::Type,
			typename DefaultIndexStringSpec< Index<TObject, TSpec> >::Type 
		> Type;
	};

//////////////////////////////////////////////////////////////////////////////
// original text

	template < typename TText, typename TSpec >
	struct Fibre< Index<TText, TSpec>, FibreText> {
		typedef TText Type;
	};

//////////////////////////////////////////////////////////////////////////////
// concatenated text

	template < typename TText, typename TSpec >
	struct Fibre< Index<TText, TSpec>, FibreRawText> {
		typedef typename Concatenator<TText>::Type Type;
	};

//////////////////////////////////////////////////////////////////////////////
// suffix array type

	template < typename TText, typename TSpec >
	struct Fibre< Index<TText, TSpec>, FibreSA> {
		typedef String<
			typename SAValue< Index<TText, TSpec> >::Type,
			typename DefaultIndexStringSpec< Index<TText, TSpec> >::Type 
		> Type;
	};

//////////////////////////////////////////////////////////////////////////////
// globalize functor

	template <typename InType, typename TLimitsString, typename Result = typename Value<TLimitsString>::Type>
	struct FunctorGlobalize : public ::std::unary_function<InType,Result>
	{
		TLimitsString const *limits;
		FunctorGlobalize() {}
		FunctorGlobalize(TLimitsString const &_limits) : limits(&_limits) {}

		inline Result operator()(InType const &x) const
		{
			return posGlobalize(x, *limits);
		}
    };

	template <typename InType, typename Result>
	struct FunctorGlobalize<InType, Nothing, Result> : public ::std::unary_function<InType,InType>
	{
		FunctorGlobalize() {}
		FunctorGlobalize(Nothing const &) {}

        inline InType operator()(InType const &x) const
        {
			return x;
		}
    };

//////////////////////////////////////////////////////////////////////////////
// raw suffix array contains integer offsets relative to raw text
/*
	template < typename TString, typename TSpec >
	struct Fibre< Index<TString, TSpec>, FibreRawSA>:
		public Fibre< Index<TString, TSpec> const, FibreSA> {};

	template < typename TString, typename TSSetSpec, typename TSpec >
	struct Fibre< Index<StringSet<TString, TSSetSpec>, TSpec>, FibreRawSA> 
	{
		typedef Index< StringSet<TString, TSSetSpec>, TSpec> TIndex;
		typedef ModifiedString<
			typename Fibre<TIndex, FibreSA>::Type,
			ModView< FunctorGlobalize< 
				typename Value< typename Fibre<TIndex, FibreSA>::Type >::Type,
				typename StringSetLimits<TString>::Type >
			>
		> Type;
	};
*/
	template < typename TText, typename TSpec >
	struct Fibre< Index<TText, TSpec>, FibreRawSA> 
	{
		typedef Index<TText, TSpec> TIndex;
		typedef ModifiedString<
			typename Fibre<TIndex, FibreSA>::Type,
			ModView< FunctorGlobalize< 
				typename Value< typename Fibre<TIndex, FibreSA>::Type >::Type,
				typename StringSetLimits<TText>::Type >
			>
		> Type;
	};

//////////////////////////////////////////////////////////////////////////////
// default burrows-wheeler table

	template < typename TText, typename TSpec >
	struct Fibre< Index<TText, TSpec>, FibreBwt> {
		typedef String <
			typename Value< Index<TText, TSpec> >::Type,
			typename DefaultIndexStringSpec< Index<TText, TSpec> >::Type
		> Type;
	};


//////////////////////////////////////////////////////////////////////////////
// default fibre creators

	template < typename TText, typename TSpec >
	struct DefaultIndexCreator<Index<TText, TSpec>, FibreSA> {
        typedef Skew7 Type;							// standard suffix array creator is skew7
    };

	template < typename TText, typename TSpec >
	struct DefaultIndexCreator<Index<TText, IndexEsa<CyclicShape<TSpec> > >, FibreSA> {
        typedef GappedSAQSort Type;
    };

	template < typename TText, typename TSpec >
	struct DefaultIndexCreator<Index<TText, TSpec>, FibreLcp> {
        typedef Kasai Type;
    };

	template < typename TText, typename TSpec >
	struct DefaultIndexCreator<Index<TText, TSpec>, FibreBwt> {
        typedef Bwt Type;
    };

	template < typename TText, typename TSpec >
	struct DefaultIndexCreator<Index<TText, TSpec>, FibreChildtab> {
        typedef Childtab Type;
    };


	template <typename TText, typename TSpec>
	inline Holder<TText> & _dataHost(Index<TText, TSpec> &index) {
		return index.text;
	}
	template <typename TText, typename TSpec>
	inline Holder<TText> const & _dataHost(Index<TText, TSpec> const &index) {
		return index.text;
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.getFibre:
..summary:Returns a specific fibre of a container.
..signature:getFibre(index, fibreTag)
..class:Class.Index
..cat:Index
..param.index:The index holding the fibre.
...type:Class.Index
..param.fibreTag:A tag that identifies the @Metafunction.Fibre@.
...type:Tag.ESA Index Fibres
...type:Tag.QGram Index Fibres
...type:Tag.WOTD Index Fibres
...type:Tag.FM Index Fibres
..returns:A reference to the @Metafunction.Fibre@ object.
..include:seqan/index.h
..example
...text:The following code shows a simple example how the function @Function.getFibre@ is used.
...file:demos/index/index_begin_range_goDown_representative_repLength.cpp
...output:The string ISSI occurs 2 times in MISSISSIPPI and has 4 characters.
The string ISSI occurs 2 times in MISSISSIPPI and has 4 characters.
*/

/*
 * @fn Index#getFibre
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns a specific fibre of a container.
 * 
 * @signature getFibre(container, fibreTag)
 * 
 * @param fibreTag A tag that identifies the @link Index#Fibre @endlink. Types: @link IndexEsaFibres Index Esa Fibres @endlink, 
 * @link FMIndexFibres FM Index Fibres @endlink, @link IndexSaFibres Index SA Fibres @endlink, @link IndexWotdFibres Index Wotd Fibres @endlink, @link IndexDfiFibres Index Dfi Fibres @endlink and @link IndexQGramFibres Index QGram Fibres @endlink.
 *
 * @param container The container holding the fibre. Types: @link Index Index @endlink
 * 
 * @return TReturn A reference to the @link Fibre @endlink object.
 * 
 * @section Examples
 * 
 * @code{.cpp}
 * Index< String<char> > index_esa("tobeornottobe");
 *  
 * String<char> & text = getFibre(indexEsa, EsaText());
 * @endcode
 *
 * @see Index#Fibre
 */

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreText>::Type & 
	getFibre(Index<TText, TSpec> &index, FibreText) {
		return value(index.text);
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreText>::Type & 
	getFibre(Index<TText, TSpec> const &index, FibreText) {
		return value(index.text);
	}

//////////////////////////////////////////////////////////////////////////////

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreRawText>::Type & 
	getFibre(Index<TText, TSpec> &index, FibreRawText) {
		return concat(value(index.text));
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreRawText>::Type & 
	getFibre(Index<TText, TSpec> const &index, FibreRawText) {
		return concat(value(index.text));
	}

//////////////////////////////////////////////////////////////////////////////

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreSA>::Type & 
	getFibre(Index<TText, TSpec> &index, FibreSA) {
		return index.sa;
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreSA>::Type & 
	getFibre(Index<TText, TSpec> const &index, FibreSA) {
		return index.sa;
	}

//////////////////////////////////////////////////////////////////////////////
/*
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreRawSA>::Type const & 
	getFibre(Index<TText, TSpec> &index, FibreRawSA) {
		return indexSA(index);
	}

	template <typename TString, typename TSSetSpec, typename TSpec>
	inline typename Fibre<Index<StringSet<TString, TSSetSpec>, TSpec>, FibreRawSA>::Type
	getFibre(Index<StringSet<TString, TSSetSpec>, TSpec> &index, FibreRawSA) 
	{
		typedef Index< StringSet<TString, TSSetSpec>, TSpec> TIndex;
		
		typedef FunctorGlobalize<
			typename Value< typename Fibre<TIndex, FibreSA>::Type >::Type,
			typename StringSetLimits<StringSet<TString, TSSetSpec> >::Type
		> TFunctor;
		
		typedef ModifiedString<
			typename Fibre<Index<StringSet<TString, TSSetSpec>, TSpec>, FibreSA>::Type,
			ModView< TFunctor >
		> ModString;

		return ModString(indexSA(index), TFunctor(stringSetLimits(indexText(index))));
	}
*/

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreRawSA>::Type
	getFibre(Index<TText, TSpec> &index, FibreRawSA) 
	{
		typedef Index<TText, TSpec> TIndex;
		
		typedef FunctorGlobalize<
			typename Value< typename Fibre<TIndex, FibreSA>::Type >::Type,
			typename StringSetLimits<TText>::Type
		> TFunctor;
		
		typedef ModifiedString<
			typename Fibre<Index<TText, TSpec>, FibreSA>::Type,
			ModView< TFunctor >
		> ModString;

		return ModString(indexSA(index), TFunctor(stringSetLimits(indexText(index))));
	}
//////////////////////////////////////////////////////////////////////////////

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreLcp>::Type & 
	getFibre(Index<TText, TSpec> &index, FibreLcp) {
		return index.lcp;
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreLcp>::Type & 
	getFibre(Index<TText, TSpec> const &index, FibreLcp) {
		return index.lcp;
	}

//////////////////////////////////////////////////////////////////////////////

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreLcpe>::Type & 
	getFibre(Index<TText, TSpec> &index, FibreLcpe) {
		return index.lcpe;
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreLcpe>::Type & 
	getFibre(Index<TText, TSpec> const &index, FibreLcpe) {
		return index.lcpe;
	}

//////////////////////////////////////////////////////////////////////////////

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreChildtab>::Type & 
	getFibre(Index<TText, TSpec> &index, FibreChildtab) {
		return index.childtab;
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreChildtab>::Type & 
	getFibre(Index<TText, TSpec> const &index, FibreChildtab) {
		return index.childtab;
	}

//////////////////////////////////////////////////////////////////////////////

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreBwt>::Type & 
	getFibre(Index<TText, TSpec> &index, FibreBwt) {
		return index.bwt;
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreBwt>::Type & 
	getFibre(Index<TText, TSpec> const &index, FibreBwt) {
		return index.bwt;
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.Index#length
..cat:Index
..summary:The number of characters in the underlying text of the index is returned.
..signature:length(index)
..class:Class.Index
..param.index:The index to return the number of characters of.
...type:Class.Index
..returns:The number of characters in the raw underlying text of the index is returned.
...metafunction:Metafunction.Size
..remarks:If the underlying text is a @Class.StringSet@ then the sum of all characters of the sequneces in the string
set is returned.
..include:seqan/index.h
..example
...text:The following code shows how @Function.length@ can be used on an index in order to determine the number of characters in the underlying text.
...file:demos/index/index_length_countSequences.cpp
...output:Hit at position: < 1 , 2 >
Hit at position: < 0 , 0 >
 */

/*
 * @fn Index#length
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns the number of characters in the raw underlying text of the
 *        index.
 * 
 * @signature length(index)
 * 
 * @param index An index of a text. Types: @link Index @endlink
 * 
 * @return TSize Returns the number of characters in the raw underlying text of the
 *        index with TSize being the result of the @link Size @endlink metafunction
 *        of @link Index @endlink.
 */

	template <typename TText, typename TSpec>
	inline typename Size<Index<TText, TSpec> >::Type 
	length(Index<TText, TSpec> const &index) {
		return length(indexRawText(index));
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.countSequences
..cat:Index
..summary:Return the number of sequences in an index' underlying text.
..signature:countSequences(index)
..class:Class.Index
..param.index:The index to return the number of sequences of.
...type:Class.Index
..returns:The number of sequences in the index' underlying text.
...metafunction:Metafunction.Size
..include:seqan/index.h
..example
...text:The following code shows how @Function.countSequences@ can be used on an index in order to determine the number of sequences in the underlying text (which can be a @Class.StringSet@).
...file:demos/index/index_length_countSequences.cpp
...output:Hit at position: < 1 , 2 >
Hit at position: < 0 , 0 >

 */

/*
 * @fn Index#countSequences
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Return the number of sequences in an index' underlying text.
 * 
 * @signature countSequences(index)
 * 
 * @param index The index to return the number of sequences of. Types: @link Index @endlink
 * 
 * @return TReturn The number of sequences in the index' underlying text.
 *                 Metafunctions: Metafunction.Size
 */

	template <typename TText, typename TSpec>
	inline typename Size<TText>::Type 
	countSequences(Index<TText, TSpec> const &index) {
		return countSequences(indexText(index));
	}

//////////////////////////////////////////////////////////////////////////////
// TODO(singer): Since this is a public function it should be documented
	template <typename TText, typename TSpec>
	struct GetSequenceByNo< Index<TText, TSpec> >
	{
		typedef typename GetSequenceByNo<TText>::Type Type;
	};

	template <typename TText, typename TSpec>
	struct GetSequenceByNo< Index<TText, TSpec> const>
	{
		typedef typename GetSequenceByNo<TText const>::Type Type;
	};

//////////////////////////////////////////////////////////////////////////////
// TODO(singer): Since this is a public function it should be documented
	template <typename TSeqNo, typename TText, typename TSpec>
	inline typename GetSequenceByNo< Index<TText, TSpec> >::Type
	getSequenceByNo(TSeqNo seqNo, Index<TText, TSpec> &index)
	{
		return getSequenceByNo(seqNo, indexText(index));
	}

	template <typename TSeqNo, typename TText, typename TSpec>
	inline typename GetSequenceByNo< Index<TText, TSpec> const>::Type
	getSequenceByNo(TSeqNo seqNo, Index<TText, TSpec> const &index)
	{
		return getSequenceByNo(seqNo, indexText(index));
	}

//////////////////////////////////////////////////////////////////////////////
// TODO(singer): Since this is a public function it should be documented
	template <typename TSeqNo, typename TText, typename TSpec>
	inline typename Size<Index<TText, TSpec> >::Type 
	sequenceLength(TSeqNo seqNo, Index<TText, TSpec> const &index) {
		return sequenceLength(seqNo, indexText(index));
	}

//////////////////////////////////////////////////////////////////////////////
// TODO(singer): Since this is a public function it should be documented
	template <typename TPos, typename TText, typename TSpec>
	inline typename Size<Index<TText, TSpec> >::Type 
	suffixLength(TPos pos, Index<TText, TSpec> const &index)
    {
		return length(indexText(index)) - pos;
	}

	template <typename TPos, typename TString, typename TSSetSpec, typename TSpec>
	inline typename Size<Index<StringSet<TString, TSSetSpec>, TSpec> >::Type 
	suffixLength(TPos pos, Index<StringSet<TString, TSSetSpec>, TSpec> const &index)
    {
        typename StringSetLimits<StringSet<TString, TSSetSpec> >::Type const &limits = stringSetLimits(index);
		return sequenceLength(getSeqNo(pos, limits), index) - getSeqOffset(pos, limits);
	}


//////////////////////////////////////////////////////////////////////////////
/**
.Function.textAt:
..summary:Shortcut for $value(indexText(..), ..)$.
..cat:Index
..signature:textAt(position, index)
..class:Class.Index
..param.position:A position in the array on which the value should be accessed.
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference or proxy to the value.
..include:seqan/index.h
..example
...text:The following code shows how the BWT of a text can be computed.
...file:demos/index/index_textAt_indexText_saAt_indexRequire.cpp
...output:BWT	Suffices
P	PI
S	SIPPI
S	SISSIPPI
M	MISSISSIPPI
I	I
P	PPI
I	IPPI
S	SSIPPI
S	SSISSIPPI
I	ISSIPPI
I	ISSISSIPPI*/

/*
 * @fn Index#textAt
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Shortcut for <tt>value(indexText(..), ..)</tt>.
 * 
 * @signature textAt(position, index)
 * 
 * @param index The @link Index @endlink object. Types: @link Index @endlink
 *
 * @param position A position in the array on which the value should be
 *                 accessed.
 * 
 * @return TReturn A reference or proxy to the value.
 *
 * @section Note The result of this function when used on an Index<TText, FMIndex<TOccSpec, Compress> > is not defined.
 */

	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, FibreRawText>::Type>::Type 
	textAt(TPos i, TIndex &index) {
		return value(getFibre(index, FibreRawText()), i);
	}
	template <typename TPos, typename TString, typename TSSetSpec, typename TSpec>
	inline typename Reference<typename Fibre< Index< StringSet<TString, TSSetSpec>, TSpec>, FibreRawText>::Type>::Type 
	textAt(TPos i, Index< StringSet<TString, TSSetSpec>, TSpec> &index) {
		return value(getFibre(index, FibreRawText()), posGlobalize(i, stringSetLimits(index)));
	}
	template <typename TPos, typename TString, typename TSSetSpec, typename TSpec>
	inline typename Reference<typename Fibre< Index< StringSet<TString, TSSetSpec>, TSpec> const, FibreRawText>::Type>::Type 
	textAt(TPos i, Index< StringSet<TString, TSSetSpec>, TSpec> const &index) {
		return value(getFibre(index, FibreRawText()), posGlobalize(i, stringSetLimits(index)));
	}
	template <typename TPos, typename TString, typename TSpec>
	inline typename Reference<typename Fibre< Index< StringSet<TString, Owner<Default> >, TSpec>, FibreRawText>::Type>::Type 
	textAt(TPos i, Index< StringSet<TString, Owner<Default> >, TSpec> &index) {
		Pair <
			typename Size< StringSet<TString, Owner<Default> > >::Type,
			typename Size< TString >::Type > locPos;
		posLocalize(locPos, i, stringSetLimits(index));
		return value(value(getFibre(index, FibreText()), getValueI1(locPos)), getValueI2(locPos));
	}
	template <typename TPos, typename TString, typename TSpec>
	inline typename Reference<typename Fibre< Index< StringSet<TString, Owner<Default> >, TSpec> const, FibreRawText>::Type>::Type 
	textAt(TPos i, Index< StringSet<TString, Owner<Default> >, TSpec> const &index) {
		Pair <
        typename Size< StringSet<TString, Owner<Default> > >::Type,
        typename Size< TString >::Type > locPos;
		posLocalize(locPos, i, stringSetLimits(index));
		return value(value(getFibre(index, FibreText()), getValueI1(locPos)), getValueI2(locPos));
	}

//////////////////////////////////////////////////////////////////////////////
// infix

	template <typename TText, typename TSpec, typename TPosBegin, typename TPosEnd>
	inline typename Infix<TText>::Type
	infix(Index<TText, TSpec> &index, TPosBegin pos_begin, TPosEnd pos_end)
	{
		return infix(indexText(index), pos_begin, pos_end);
	}

	template <typename TText, typename TSpec, typename TPosBegin, typename TPosEnd>
	inline typename Infix<TText>::Type
	infix(Index<TText, TSpec> const &index, TPosBegin pos_begin, TPosEnd pos_end)
	{
		return infix(indexText(index), pos_begin, pos_end);
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.rawtextAt:
..summary:Shortcut for $value(indexRawText(..), ..)$.
..cat:Index
..signature:rawtextAt(position, index)
..class:Class.Index
..param.position:A position in the array on which the value should be accessed.
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference or proxy to the value.
..include:seqan/index.h
*/
/*
 * @fn Index#rawtextAt
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Shortcut for <tt>value(indexRawText(..), ..)</tt>.
 * 
 * @signature rawtextAt(position, index)
 * 
 * @param index The @link Index @endlink object. Types: @link Index @endlink
 *
 * @param position A position in the array on which the value should be
 *                 accessed.
 * 
 * @return TReturn A reference or proxy to the value.
 *
 * @section Note The result of this function when used on an Index<TText, FMIndex<TOccSpec, Compress> > is not defined.
 */

	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, FibreRawText>::Type>::Type rawtextAt(TPos i, TIndex &index) {
		return value(getFibre(index, FibreRawText()), i);
	}
	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex const, FibreRawText>::Type>::Type rawtextAt(TPos i, TIndex const &index) {
		return value(getFibre(index, FibreRawText()), i);
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.saAt:
..summary:Shortcut for $value(indexSA(..), ..)$.
..cat:Index
..signature:saAt(position, index)
..class:Class.Index
..param.position:A position in the array on which the value should be accessed.
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference or proxy to the value.
..include:seqan/index.h
..example
...text:The following code shows how the BWT of a text can be computed.
...file:demos/index/index_textAt_indexText_saAt_indexRequire.cpp
...output:BWT	Suffices
P	PI
S	SIPPI
S	SISSIPPI
M	MISSISSIPPI
I	I
P	PPI
I	IPPI
S	SSIPPI
S	SSISSIPPI
I	ISSIPPI
I	ISSISSIPPI*/
/*
 * @fn IndexEsa#saAt
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Shortcut for <tt>value(indexSA(..), ..)</tt>.
 *
 * @deprecated advanced
 * 
 * @signature saAt(position, index)
 * 
 * @param index The @link Index @endlink object holding the fibre. Types:
 *              @link IndexEsa @endlink
 * @param position A position in the array on which the value should be
 *                 accessed.
 * 
 * @return TReturn A reference or proxy to the value.
 */

	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, FibreSA>::Type>::Type saAt(TPos i, TIndex &index) {
		return value(getFibre(index, FibreSA()), i);
	}
	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex const, FibreSA>::Type>::Type saAt(TPos i, TIndex const &index) {
		return value(getFibre(index, FibreSA()), i);
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.rawsaAt:
..summary:Shortcut for $value(indexRawSA(..), ..)$.
..cat:Index
..signature:saAt(position, index)
..class:Class.Index
..param.position:A position in the array on which the value should be accessed.
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference or proxy to the value.
..include:seqan/index.h
*/

/*
 * @fn IndexEsa#rawsaAt
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Shortcut for <tt>value(indexRawSA(..), ..)</tt>.
 * 
 * @deprecated. advanced
 *
 * @signature rawsaAt(position, index)
 * 
 * @param index The @link Index @endlink object holding the fibre. Types:
 *              @link IndexEsa @endlink
 * @param position A position in the array on which the value should be
 *                 accessed.
 */

	template <typename TPos, typename TIndex>
	inline typename Value<typename Fibre<TIndex const, FibreRawSA>::Type>::Type rawsaAt(TPos i, TIndex const &index) {
		return posGlobalize(saAt(i, index), stringSetLimits(indexText(index)));
	}


//////////////////////////////////////////////////////////////////////////////
/**
.Function.IndexEsa#lcpAt:
..summary:Shortcut for $value(indexLcp(..), ..)$.
..cat:Index
..signature:lcpAt(position, index)
..class:Spec.IndexEsa
..param.position:A position in the array on which the value should be accessed.
..param.index:The @Spec.IndexEsa@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference or proxy to the value.
..include:seqan/index.h
*/
/*
 * @fn IndexEsa#lcpAt
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Shortcut for <tt>value(indexLcp(..), ..)</tt>.
 * 
 * @signature lcpAt(position, index)
 * 
 * @param index The @link Index @endlink object holding the fibre. Types:
 *              @link IndexEsa @endlink
 * @param position A position in the array on which the value should be
 *                 accessed.
 * 
 * @return TReturn A reference or proxy to the value.
 */

	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, FibreLcp>::Type>::Type lcpAt(TPos i, TIndex &index) {
		return value(getFibre(index, FibreLcp()), i);
	}
	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex const, FibreLcp>::Type>::Type lcpAt(TPos i, TIndex const &index) {
		return value(getFibre(index, FibreLcp()), i);
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.IndexEsa#lcpeAt:
..summary:Shortcut for $value(indexLcpe(..), ..)$.
..cat:Index
..signature:lcpeAt(position, index)
..class:Spec.IndexEsa
..param.position:A position in the array on which the value should be accessed.
..param.index:The @Spec.IndexEsa@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference or proxy to the value.
..include:seqan/index.h
*/
/*
 * @fn IndexEsa#lcpeAt
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Shortcut for <tt>value(indexLcpe(..), ..)</tt>.
 * 
 * @signature lcpeAt(position, index)
 * 
 * @param index The @link Index @endlink object holding the fibre. Types:
 *              @link IndexEsa @endlink
 * @param position A position in the array on which the value should be
 *                 accessed.
 * 
 * @return TReturn A reference or proxy to the value.
 */

	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, FibreLcpe>::Type>::Type lcpeAt(TPos i, TIndex &index) {
		return value(getFibre(index, FibreLcpe()), i);
	}
	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex const, FibreLcpe>::Type>::Type lcpeAt(TPos i, TIndex const &index) {
		return value(getFibre(index, FibreLcpe()), i);
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.IndexEsa#childAt:
..summary:Shortcut for $value(indexChildtab(..), ..)$.
..cat:Index
..signature:childAt(position, index)
..class:Spec.IndexEsa
..param.position:A position in the array on which the value should be accessed.
..param.index:The @Spec.IndexEsa@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference or proxy to the value.
..include:seqan/index.h
*/

	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, FibreChildtab>::Type>::Type childAt(TPos i, TIndex &index) {
		return value(getFibre(index, FibreChildtab()), i);
	}
	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex const, FibreChildtab>::Type>::Type childAt(TPos i, TIndex const &index) {
		return value(getFibre(index, FibreChildtab()), i);
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.IndexEsa#bwtAt:
..summary:Shortcut for $value(indexBwt(..), ..)$.
..cat:Index
..signature:bwtAt(position, index)
..class:Spec.IndexEsa
..param.position:A position in the array on which the value should be accessed.
..param.index:The @Spec.IndexEsa@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference or proxy to the value.
..include:seqan/index.h
*/
/*
 * @fn IndexEsa#childAt
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Shortcut for <tt>value(indexChildtab(..), ..)</tt>.
 * 
 * @signature childAt(position, index)
 * 
 * @param index The @link Index @endlink object holding the fibre. Types:
 *              @link IndexEsa @endlink
 * @param position A position in the array on which the value should be
 *                 accessed.
 * 
 * @return TReturn A reference or proxy to the value.
 */

	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, FibreBwt>::Type>::Type bwtAt(TPos i, TIndex &index) {
		return value(getFibre(index, FibreBwt()), i);
	}
	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex const, FibreBwt>::Type>::Type bwtAt(TPos i, TIndex const &index) {
		return value(getFibre(index, FibreBwt()), i);
	}

//////////////////////////////////////////////////////////////////////////////

    template <typename TIndex, typename TPos, typename TSize>
	inline typename SAValue<TIndex>::Type toSuffixPosition(TIndex &, TPos i, TSize) {
        return i;
	}
    template <typename TIndex, typename TPos, typename TSize>
	inline typename SAValue<TIndex const>::Type toSuffixPosition(TIndex const &, TPos i, TSize) {
        return i;
	}

//////////////////////////////////////////////////////////////////////////////
// interface for infinity/invalid values

	template <typename TValue>
	inline void _setSizeInval(TValue &v) {
		v = MaxValue<TValue>::VALUE;
	}

	template <typename TValue>
	inline bool _isSizeInval(TValue const &v) {
//IOREV _notio_
		return v == MaxValue<TValue>::VALUE;
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.indexText:
..summary:Shortcut for $getFibre(index, FibreText())$.
..cat:Index
..signature:indexText(index)
..class:Class.Index
..param.index:The @Class.Index@ object holding the fibre.
...type:Class.Index
..returns:A reference to the text fibre (original text).
..include:seqan/index.h
..example
...text:The following code shows how the BWT of a text can be computed.
...file:demos/index/index_textAt_indexText_saAt_indexRequire.cpp
...output:BWT	Suffices
P	PI
S	SIPPI
S	SISSIPPI
M	MISSISSIPPI
I	I
P	PPI
I	IPPI
S	SSIPPI
S	SSISSIPPI
I	ISSIPPI
I	ISSISSIPPI*/
/*
 * @fn Index#indexText
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Shortcut for <tt>getFibre(.., EsaText)</tt>.
 * 
 * @signature indexText(index)
 * 
 * @param index The @link Index @endlink object holding the fibre. Types:
 *              @link Index @endlink
 * 
 * @return TReturn A reference to the text of the index.
 *
 * @section Note The result of this function when used on an Index<TText, FMIndex<TOccSpec, Compress> > is not defined.
 *
 */

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreText>::Type & indexText(Index<TText, TSpec> &index) { return getFibre(index, FibreText()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreText>::Type & indexText(Index<TText, TSpec> const &index) { return getFibre(index, FibreText()); }

//////////////////////////////////////////////////////////////////////////////

	template <typename TText, typename TSpec>
	inline typename StringSetLimits<TText const>::Type
	stringSetLimits(Index<TText, TSpec> &) { 
		return Nothing(); 
	}

	template <typename TText, typename TSpec>
	inline typename StringSetLimits<TText const>::Type
	stringSetLimits(Index<TText, TSpec> const &) { 
		return Nothing(); 
	}

	template <typename TString, typename TSSetSpec, typename TSpec>
	inline typename StringSetLimits< StringSet<TString, TSSetSpec> const >::Type & 
	stringSetLimits(Index<StringSet<TString, TSSetSpec>, TSpec> &index) {
		return stringSetLimits(indexText(index)); 
	}

	template <typename TString, typename TSSetSpec, typename TSpec>
	inline typename StringSetLimits< StringSet<TString, TSSetSpec> const >::Type & 
	stringSetLimits(Index<StringSet<TString, TSSetSpec>, TSpec> const &index) {
		return stringSetLimits(indexText(index)); 
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.indexRawText:
..summary:Shortcut for $getFibre(.., EsaRawText)$.
..cat:Index
..signature:indexRawText(index)
..class:Class.Index
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference to the @Tag.ESA Index Fibres.EsaRawText@ fibre (concatenated input text).
..include:seqan/index.h
*/
//TODO(singer) The RawText Fibre exist for more then the Esa index
/*
 * @fn IndexEsa#indexRawText
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Shortcut for <tt>$getFibre(.., EsaRawText)</tt>.
 * 
 * @signature rawtextAt(position, index)
 * 
 * @param index The @link Index @endlink object. Types: @link Index @endlink
 *
 * @param position A position in the array on which the value should be
 *                 accessed.
 * 
 * @return TReturn A reference or proxy to the value.
 */

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreRawText>::Type & indexRawText(Index<TText, TSpec> &index) { return getFibre(index, FibreRawText()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreRawText>::Type & indexRawText(Index<TText, TSpec> const &index) { return getFibre(index, FibreRawText()); }

//////////////////////////////////////////////////////////////////////////////
/**
.Function.indexSA:
..summary:Shortcut for $getFibre(.., FibreSA)$.
..cat:Index
..signature:indexSA(index)
..class:Class.Index
..param.index:The @Class.Index@ object holding the fibre.
...type:Class.Index
..returns:A reference to the suffix array fibre.
..include:seqan/index.h
..example
...text:The following code shows how the BWT of a text can be computed.
...file:demos/index/index_textAt_indexText_saAt_indexRequire.cpp
...output:BWT	Suffices
P	PI
S	SIPPI
S	SISSIPPI
M	MISSISSIPPI
I	I
P	PPI
I	IPPI
S	SSIPPI
S	SSISSIPPI
I	ISSIPPI
I	ISSISSIPPI*/
//TODO(singer) The function in not only defined for the esa index
/*
 * @fn IndexEsa#indexSA
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Shortcut for <tt>getFibre(.., EsaSA)</tt>.
 * 
 * @signature indexSA(index)
 * 
 * @param index The @link Index @endlink object holding the fibre. Types:
 *              @link IndexEsa @endlink
 * 
 * @return TReturn A reference to the @link ESA Index Fibres.EsaSA @endlink
 *                 fibre (suffix array).
 */

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreSA>::Type & indexSA(Index<TText, TSpec> &index) { return getFibre(index, FibreSA()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreSA>::Type & indexSA(Index<TText, TSpec> const &index) { return getFibre(index, FibreSA()); }

//////////////////////////////////////////////////////////////////////////////
/**
.Function.indexRawSA:
..summary:Shortcut for $getFibre(.., EsaRawSA)$.
..cat:Index
..signature:indexRawSA(index)
..class:Class.Index
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference to the @Tag.ESA Index Fibres.EsaRawSA@ fibre (suffix array).
..include:seqan/index.h
*/
/*
 * @fn Index#indexRawSA
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Shortcut for <tt>getFibre(.., EsaRawSA)</tt>.
 * 
 * @signature indexRawSA(index)
 * 
 * @param index The @link Index @endlink object holding the fibre. Types:
 *              @link IndexEsa @endlink
 * 
 * @return TReturn A reference to the @link ESA Index Fibres.EsaRawSA @endlink
 *                 fibre (suffix array).
 */

/*
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreRawSA>::Type const & indexRawSA(Index<TText, TSpec> &index) { return getFibre(index, FibreRawSA()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreRawSA>::Type const & indexRawSA(Index<TText, TSpec> const &index) { return getFibre(index, FibreRawSA()); }
*/
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreRawSA>::Type indexRawSA(Index<TText, TSpec> &index) { return getFibre(index, FibreRawSA()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreRawSA>::Type indexRawSA(Index<TText, TSpec> const &index) { return getFibre(index, FibreRawSA()); }

//////////////////////////////////////////////////////////////////////////////
/**
.Function.IndexEsa#indexLcp:
..summary:Shortcut for $getFibre(.., EsaLcp)$.
..cat:Index
..signature:indexLcp(index)
..class:Spec.IndexEsa
..param.index:The @Spec.IndexEsa@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference to the @Tag.ESA Index Fibres.EsaLcp@ fibre (lcp table).
..include:seqan/index.h
*/
/*
 * @fn IndexEsa#indexLcp
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Shortcut for <tt>getFibre(.., EsaLcp)</tt>.
 *
 * deprecated. advanced
 * 
 * @signature indexLcp(index)
 * 
 * @param index The @link Index @endlink object holding the fibre. Types:
 *              @link IndexEsa @endlink
 * 
 * @return TReturn A reference to the @link ESA Index Fibres.EsaLcp @endlink
 *                 fibre (lcp table).
 */

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreLcp>::Type & indexLcp(Index<TText, TSpec> &index) { return getFibre(index, FibreLcp()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreLcp>::Type & indexLcp(Index<TText, TSpec> const &index) { return getFibre(index, FibreLcp()); }

//////////////////////////////////////////////////////////////////////////////
/**
.Function.IndexEsa#indexLcpe:
..summary:Shortcut for $getFibre(.., EsaLcpe)$.
..cat:Index
..signature:indexLcpe(index)
..class:Spec.IndexEsa
..param.index:The @Spec.IndexEsa@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference to the @Tag.ESA Index Fibres.EsaLcpe@ fibre (enhanced lcp table).
..include:seqan/index.h
*/
/*
 * @fn IndexEsa#indexLcpe
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Shortcut for <tt>getFibre(.., EsaLcpe)</tt>.
 * 
 * @signature indexLcpe(index)
 * 
 * @param index The @link Index @endlink object holding the fibre. Types:
 *              @link IndexEsa @endlink
 * 
 * @return TReturn A reference to the @link ESA Index Fibres.EsaLcpe @endlink
 *                 fibre (enhanced lcp table).
 */

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreLcpe>::Type & indexLcpe(Index<TText, TSpec> &index) { return getFibre(index, FibreLcpe()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreLcpe>::Type & indexLcpe(Index<TText, TSpec> const &index) { return getFibre(index, FibreLcpe()); }

//////////////////////////////////////////////////////////////////////////////
/**
.Function.IndexEsa#indexBwt:
..summary:Shortcut for $getFibre(.., EsaBwt)$.
..cat:Index
..signature:indexBwt(index)
..class:Spec.IndexEsa
..param.index:The @Spec.IndexEsa@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference to the @Tag.ESA Index Fibres.EsaBwt@ fibre (Burrows-Wheeler table).
..include:seqan/index.h
*/
/*
 * @fn IndexEsa#indexBwt
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Shortcut for <tt>getFibre(.., EsaBwt)</tt>.
 * 
 * @signature indexBwt(index)
 * 
 * @param index The @link IndexEsa @endlink object holding the fibre. Types:
 *              @link IndexEsa @endlink
 * 
 * @return TReturn A reference to the @link ESA Index Fibres.EsaBwt @endlink
 *                 fibre (Burrows-Wheeler table).
 */

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreBwt>::Type & indexBwt(Index<TText, TSpec> &index) { return getFibre(index, FibreBwt()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreBwt>::Type & indexBwt(Index<TText, TSpec> const &index) { return getFibre(index, FibreBwt()); }

//////////////////////////////////////////////////////////////////////////////
/**
.Function.IndexEsa#indexChildtab:
..summary:Shortcut for $getFibre(.., EsaChildtab)$.
..cat:Index
..signature:indexChildtab(index)
..class:Spec.IndexEsa
..param.index:The @Spec.IndexEsa@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference to the @Tag.ESA Index Fibres.EsaChildtab@ fibre (child table).
..include:seqan/index.h
*/
/*
 * @fn IndexEsa#indexChildtab
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Shortcut for <tt>getFibre(.., EsaChildtab)</tt>.
 * 
 * @signature indexChildtab(index)
 * 
 * @param index The @link Index @endlink object holding the fibre. Types:
 *              @link IndexEsa @endlink
 * 
 * @return TReturn A reference to the @link ESA Index Fibres.EsaChildtab
 *                 @endlink fibre (child table).
 */

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreChildtab>::Type & indexChildtab(Index<TText, TSpec> &index) { return getFibre(index, FibreChildtab()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreChildtab>::Type & indexChildtab(Index<TText, TSpec> const &index) { return getFibre(index, FibreChildtab()); }


// ----------------------------------------------------------------------------
// Function open
// ----------------------------------------------------------------------------
/**
.Function.Index#open
..class:Class.Index
..summary:This functions opens an index from disk.
..signature:open(index, fileName [, mode])
..param.index:The index to be opened.
...type:Class.Index
..param.fileName:C-style character string containing the file name.
..param.mode:The combination of flags defining how the file should be opened.
...remarks:To open a file read-only, write-only or to read and write use $OPEN_RDONLY$, $OPEN_WRONLY$, or $OPEN_RDWR$.
...remarks:To create or overwrite a file add $OPEN_CREATE$.
...remarks:To append a file if existing add $OPEN_APPEND$.
...remarks:To circumvent problems, files are always opened in binary mode.
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A $bool$ which is $true$ on success.
..include:seqan/index.h
..example
...text:The following code shows how the function @Function.open@ is used with indices.
...file:demos/index/index_open_save.cpp
...output:1
1
*/
/*
 * @fn Index#open
 * 
 * @headerfile seqan/index.h
 * 
 * @brief This functions opens an index from disk.
 * 
 * @signature open(index, fileName [, mode])
 * 
 * @param mode The combination of flags defining how the file should be
 *                 opened.To open a file read-only, write-only or to read and
 *                 write use <tt>OPEN_RDONLY</tt>, <tt>OPEN_WRONLY</tt>, or
 *                 <tt>OPEN_RDWR</tt>.To create or overwrite a file add
 *                 <tt>OPEN_CREATE</tt>.To append a file if existing add
 *                 <tt>OPEN_APPEND</tt>.To circumvent problems, files are always
 *                 opened in binary mode. Default: <tt>OPEN_RDWR | OPEN_CREATE |
 *                 OPEN_APPEND</tt>
 *
 * @param index The index to be opened. Types: @link Index @endlink
 *
 * @param fileName C-style character string containing the file name.
 * 
 * @return TReturn A <tt>bool</tt> which is <tt>true</tt> on success.
 */

// ----------------------------------------------------------------------------
// Function save
// ----------------------------------------------------------------------------

/* 
 * @fn Index#save
 * 
 * @headerfile seqan/index.h
 * 
 * @brief This functions saves an index to disk.
 * 
 * @signature save(index, fileName [, mode])
 * 
 * @param mode The combination of flags defining how the file should be
 *                 opened.To open a file read-only, write-only or to read and
 *                 write use <tt>OPEN_RDONLY</tt>, <tt>OPEN_WRONLY</tt>, or
 *                 <tt>OPEN_RDWR</tt>.To create or overwrite a file add
 *                 <tt>OPEN_CREATE</tt>.To append a file if existing add
 *                 <tt>OPEN_APPEND</tt>.To circumvent problems, files are always
 *                 opened in binary mode. Default: <tt>OPEN_RDWR | OPEN_CREATE |
 *                 OPEN_APPEND</tt>
 *
 * @param index The index to be saved to disk. Types: @link Index @endlink
 *
 * @param fileName C-style character string containing the file name.
 * 
 * @return TReturn A <tt>bool</tt> which is <tt>true</tt> on success.
 */
/**
.Function.Index#save
..class:Class.Index
..summary:This functions saves an index to disk.
..signature:save(index, fileName [, mode])
..param.index:The index to be saved to disk.
...type:Class.Index
..param.fileName:C-style character string containing the file name.
..param.mode:The combination of flags defining how the file should be opened.
...remarks:To open a file read-only, write-only or to read and write use $OPEN_RDONLY$, $OPEN_WRONLY$, or $OPEN_RDWR$.
...remarks:To create or overwrite a file add $OPEN_CREATE$.
...remarks:To append a file if existing add $OPEN_APPEND$.
...remarks:To circumvent problems, files are always opened in binary mode.
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A $bool$ which is $true$ on success.
..include:seqan/index.h
..example
...text:The following code shows how the function @Function.open@ is used with indices.
...file:demos/index/index_open_save.cpp
...output:1
1
*/
#endif // SEQAN_HEADER_INDEX_BASE_H
    
}

#endif

