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
// This file contains Index specializations for Views.
// ==========================================================================

#ifndef SEQAN_EXTRAS_INDEX_VIEW_H_
#define SEQAN_EXTRAS_INDEX_VIEW_H_

namespace seqan {

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction View                                                    [Index]
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec>
struct View<Index<TText, TSpec> >
{
    typedef Index<typename View<TText>::Type, TSpec>        Type;
};

template <typename TText, typename TSpec>
struct View<Index<TText, TSpec> const>
{
    typedef Index<typename View<TText>::Type, TSpec> const  Type;
};

// ----------------------------------------------------------------------------
// Metafunction View                                                  [LfTable]
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec>
struct View<LfTable<TText, TSpec> >
{
    typedef LfTable<typename View<TText>::Type, TSpec>          Type;
};

template <typename TText, typename TSpec>
struct View<LfTable<TText, TSpec> const>
{
    typedef LfTable<typename View<TText>::Type, TSpec> const    Type;
};

// ----------------------------------------------------------------------------
// Metafunction View                                             [CompressedSA]
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec>
struct View<CompressedSA<TText, TSpec> >
{
    typedef CompressedSA<typename View<TText>::Type, TSpec>         Type;
};

template <typename TText, typename TSpec>
struct View<CompressedSA<TText, TSpec> const>
{
    typedef CompressedSA<typename View<TText>::Type, TSpec> const   Type;
};

// ----------------------------------------------------------------------------
// Metafunction View                                           [PrefixSumTable]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct View<PrefixSumTable<TValue, TSpec> >
{
    typedef PrefixSumTable<TValue, View<TSpec> >        Type;
};

template <typename TValue, typename TSpec>
struct View<PrefixSumTable<TValue, TSpec> const>
{
    typedef PrefixSumTable<TValue, View<TSpec> > const  Type;
};

// ----------------------------------------------------------------------------
// Metafunction View                                           [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct View<RankDictionary<TwoLevels<TValue, TSpec> > >
{
    typedef RankDictionary<TwoLevels<TValue, View<TSpec> > >        Type;
};

template <typename TValue, typename TSpec>
struct View<RankDictionary<TwoLevels<TValue, TSpec> > const>
{
    typedef RankDictionary<TwoLevels<TValue, View<TSpec> > > const  Type;
};

// ----------------------------------------------------------------------------
// Metafunction View                                             [SparseString]
// ----------------------------------------------------------------------------

template <typename TString, typename TSpec>
struct View<SparseString<TString, TSpec> >
{
    typedef SparseString<ContainerView<TString>, TSpec>         Type;
};

template <typename TString, typename TSpec>
struct View<SparseString<TString, TSpec> const>
{
    typedef SparseString<ContainerView<TString>, TSpec> const   Type;
};

// ----------------------------------------------------------------------------
// Metafunction Fibre                                              [Index View]
// ----------------------------------------------------------------------------

// NOTE(esiragusa): This causes template resolution ambiguity.
//template <typename TText, typename TViewSpec, typename TSpec, typename TFibre>
//struct Fibre<Index<ContainerView<TText, TViewSpec>, TSpec>, Tag<TFibre> const>
//{
//    typedef typename View<typename Fibre<Index<TText, TSpec>, Tag<TFibre> const>::Type>::Type     Type;
//};

// ----------------------------------------------------------------------------
// Metafunction FibreText                                          [Index View]
// ----------------------------------------------------------------------------

template <typename TText, typename TViewSpec, typename TSpec>
struct Fibre<Index<ContainerView<TText, TViewSpec>, TSpec>, FibreText>
{
    typedef typename View<TText>::Type  Type;
};

template <typename TText, typename TViewSpec, typename TSpec>
struct Fibre<Index<ContainerView<TText, TViewSpec>, TSpec> const, FibreText> :
    public Fibre<Index<ContainerView<TText, TViewSpec>, TSpec>, FibreText> {};

// ----------------------------------------------------------------------------
// Metafunction FibreRawText                                       [Index View]
// ----------------------------------------------------------------------------

template <typename TText, typename TViewSpec, typename TSpec>
struct Fibre<Index<ContainerView<TText, TViewSpec>, TSpec>, FibreRawText>
{
    typedef typename View<typename Fibre<Index<TText, TSpec>, FibreRawText>::Type>::Type    Type;
};

template <typename TText, typename TViewSpec, typename TSpec>
struct Fibre<Index<ContainerView<TText, TViewSpec>, TSpec> const, FibreRawText> :
    public Fibre<Index<ContainerView<TText, TViewSpec>, TSpec>, FibreRawText> {};

// ----------------------------------------------------------------------------
// Metafunction FibreSA                                            [Index View]
// ----------------------------------------------------------------------------

template <typename TText, typename TViewSpec, typename TSpec>
struct Fibre<Index<ContainerView<TText, TViewSpec>, TSpec>, FibreSA>
{
    typedef typename View<typename Fibre<Index<TText, TSpec>, FibreSA>::Type>::Type         Type;
};

template <typename TText, typename TViewSpec, typename TSpec>
struct Fibre<Index<ContainerView<TText, TViewSpec>, TSpec> const, FibreSA>
{
    typedef typename View<typename Fibre<Index<TText, TSpec> const, FibreSA>::Type>::Type   Type;
};

// ----------------------------------------------------------------------------
// Metafunction FibreLcp                                           [Index View]
// ----------------------------------------------------------------------------

template <typename TText, typename TViewSpec, typename TSpec>
struct Fibre<Index<ContainerView<TText, TViewSpec>, TSpec>, FibreLcp>
{
    typedef typename View<typename Fibre<Index<TText, TSpec>, FibreLcp>::Type>::Type        Type;
};

template <typename TText, typename TViewSpec, typename TSpec>
struct Fibre<Index<ContainerView<TText, TViewSpec>, TSpec> const, FibreLcp>
{
    typedef typename View<typename Fibre<Index<TText, TSpec> const, FibreLcp>::Type>::Type  Type;
};

// ----------------------------------------------------------------------------
// Metafunction FibreChildtab                                      [Index View]
// ----------------------------------------------------------------------------

template <typename TText, typename TViewSpec, typename TSpec>
struct Fibre<Index<ContainerView<TText, TViewSpec>, TSpec>, FibreChildtab>
{
    typedef typename View<typename Fibre<Index<TText, TSpec>, FibreChildtab>::Type>::Type   Type;
};

template <typename TText, typename TViewSpec, typename TSpec>
struct Fibre<Index<ContainerView<TText, TViewSpec>, TSpec> const, FibreChildtab>
{
    typedef typename View<typename Fibre<Index<TText, TSpec> const, FibreChildtab>::Type>::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction FibreBwt                                           [Index View]
// ----------------------------------------------------------------------------

template <typename TText, typename TViewSpec, typename TSpec>
struct Fibre<Index<ContainerView<TText, TViewSpec>, TSpec>, FibreBwt>
{
    typedef typename View<typename Fibre<Index<TText, TSpec>, FibreBwt>::Type>::Type        Type;
};

template <typename TText, typename TViewSpec, typename TSpec>
struct Fibre<Index<ContainerView<TText, TViewSpec>, TSpec> const, FibreBwt>
{
    typedef typename View<typename Fibre<Index<TText, TSpec> const, FibreBwt>::Type>::Type  Type;
};

// ----------------------------------------------------------------------------
// Metafunction Member                                             [Index View]
// ----------------------------------------------------------------------------

template <typename TText, typename TViewSpec, typename TSpec>
struct Member<Index<ContainerView<TText, TViewSpec>, TSpec>, FibreText>
{
    typedef Index<ContainerView<TText, TViewSpec>, TSpec>       TIndex_;
    typedef typename Fibre<TIndex_, FibreText>::Type            Type;
};

template <typename TText, typename TViewSpec, typename TSpec>
struct Member<Index<ContainerView<TText, TViewSpec>, TSpec> const, FibreText>
{
    typedef Index<ContainerView<TText, TViewSpec>, TSpec> const TIndex_;
    typedef typename Fibre<TIndex_, FibreText>::Type            Type;
};

// ----------------------------------------------------------------------------
// Metafunction FibreSA                                          [FMIndex View]
// ----------------------------------------------------------------------------

template <typename TText, typename TViewSpec, typename TOccSpec, typename TSpec>
struct Fibre<Index<ContainerView<TText, TViewSpec>, FMIndex<TOccSpec, TSpec> >, FibreSA>
{
    typedef typename View<typename Fibre<Index<TText, FMIndex<TOccSpec, TSpec> >, FibreSA>::Type>::Type    Type;
};

template <typename TText, typename TViewSpec, typename TOccSpec, typename TSpec>
struct Fibre<Index<ContainerView<TText, TViewSpec>, FMIndex<TOccSpec, TSpec> > const, FibreSA>
{
    typedef typename View<typename Fibre<Index<TText, FMIndex<TOccSpec, TSpec> > const, FibreSA>::Type>::Type   Type;
};

// ----------------------------------------------------------------------------
// Metafunction FibrePrefixSum                                   [LfTable View]
// ----------------------------------------------------------------------------

template <typename TText, typename TViewSpec, typename TSpec>
struct Fibre<LfTable<ContainerView<TText, TViewSpec>, TSpec>, FibrePrefixSum>
{
    typedef typename View<typename Fibre<LfTable<TText, TSpec>, FibrePrefixSum>::Type>::Type   Type;
};

template <typename TText, typename TViewSpec, typename TSpec>
struct Fibre<LfTable<ContainerView<TText, TViewSpec>, TSpec> const, FibrePrefixSum>
{
    typedef typename View<typename Fibre<LfTable<TText, TSpec> const, FibrePrefixSum>::Type>::Type  Type;
};

// ----------------------------------------------------------------------------
// Metafunction FibreValues                                      [LfTable View]
// ----------------------------------------------------------------------------

template <typename TText, typename TViewSpec, typename TSpec>
struct Fibre<LfTable<ContainerView<TText, TViewSpec>, TSpec>, FibreValues>
{
    typedef typename View<typename Fibre<LfTable<TText, TSpec>, FibreValues>::Type>::Type   Type;
};

template <typename TText, typename TViewSpec, typename TSpec>
struct Fibre<LfTable<ContainerView<TText, TViewSpec>, TSpec> const, FibreValues>
{
    typedef typename View<typename Fibre<LfTable<TText, TSpec> const, FibreValues>::Type>::Type     Type;
};

// ----------------------------------------------------------------------------
// Metafunction FibreSentinels                                   [LfTable View]
// ----------------------------------------------------------------------------

//template <typename TText, typename TViewSpec, typename TSpec>
//struct Fibre<LfTable<ContainerView<TText, TViewSpec>, TSpec>, FibreSentinels>
//{
//    typedef typename Fibre<LfTable<TText, TSpec>, FibreSentinels>::Type   Type;
//};

// NOTE(esiragusa): StringSet specialization
//template <typename TText, typename TSSetSpec, typename TViewSpec, typename TSpec>
//struct Fibre<LfTable<ContainerView<StringSet<TText, TSSetSpec>, TViewSpec>, TSpec>, FibreSentinels>
//{
//    typedef typename View<typename Fibre<LfTable<StringSet<TText, TSSetSpec>, TSpec>, FibreSentinels>::Type>::Type    Type;
//};

// ----------------------------------------------------------------------------
// Metafunction FibreEntries                              [PrefixSumTable View]
// ----------------------------------------------------------------------------

template <typename TChar, typename TSpec>
struct Fibre<PrefixSumTable<TChar, View<TSpec> >, FibreEntries>
{
    typedef typename View<typename Fibre<PrefixSumTable<TChar, TSpec>, FibreEntries>::Type>::Type   Type;
};

template <typename TChar, typename TSpec>
struct Fibre<PrefixSumTable<TChar, View<TSpec> > const, FibreEntries>
{
    typedef typename View<typename Fibre<PrefixSumTable<TChar, TSpec> const, FibreEntries>::Type>::Type     Type;
};

// ----------------------------------------------------------------------------
// Metafunction Fibre                                     [RankDictionary View]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct Fibre<RankDictionary<TwoLevels<TValue, View<TSpec> > >, FibreRanks>
{
    typedef typename View<typename Fibre<RankDictionary<TwoLevels<TValue, TSpec> >, FibreRanks>::Type>::Type    Type;
};

template <typename TValue, typename TSpec>
struct Fibre<RankDictionary<TwoLevels<TValue, View<TSpec> > > const, FibreRanks>
{
    typedef typename View<typename Fibre<RankDictionary<TwoLevels<TValue, TSpec> > const, FibreRanks>::Type>::Type  Type;
};

// ----------------------------------------------------------------------------
// Metafunction Fibre                                       [CompressedSA View]
// ----------------------------------------------------------------------------

template <typename TText, typename TViewSpec, typename TSpec>
struct Fibre<CompressedSA<ContainerView<TText, TViewSpec>, TSpec>, FibreSparseString>
{
    typedef typename View<typename Fibre<CompressedSA<TText, TSpec>, FibreSparseString>::Type>::Type    Type;
};

template <typename TText, typename TViewSpec, typename TSpec>
struct Fibre<CompressedSA<ContainerView<TText, TViewSpec>, TSpec> const, FibreSparseString>
{
    typedef typename View<typename Fibre<CompressedSA<TText, TSpec> const, FibreSparseString>::Type>::Type  Type;
};

// ----------------------------------------------------------------------------
// Metafunction Member                                      [SparseString View]
// ----------------------------------------------------------------------------

template <typename TText, typename TViewSpec, typename TSpec>
struct Member<CompressedSA<ContainerView<TText, TViewSpec>, TSpec>, FibreLF>
{
    typedef typename Fibre<CompressedSA<ContainerView<TText, TViewSpec>, TSpec>, FibreLF>::Type     Type;
};

// ----------------------------------------------------------------------------
// Metafunction FibreIndicators                             [SparseString View]
// ----------------------------------------------------------------------------

template <typename TString, typename TSpec>
struct Fibre<SparseString<ContainerView<TString>, TSpec>, FibreIndicators>
{
    typedef typename View<typename Fibre<SparseString<TString, TSpec>, FibreIndicators>::Type>::Type        Type;
};

template <typename TString, typename TSpec>
struct Fibre<SparseString<ContainerView<TString>, TSpec> const, FibreIndicators>
{
    typedef typename View<typename Fibre<SparseString<TString, TSpec> const, FibreIndicators>::Type>::Type  Type;
};

// ----------------------------------------------------------------------------
// Metafunction FibreValues                                 [SparseString View]
// ----------------------------------------------------------------------------

template <typename TString, typename TSpec>
struct Fibre<SparseString<ContainerView<TString>, TSpec>, FibreValues>
{
    typedef typename View<typename Fibre<SparseString<TString, TSpec>, FibreValues>::Type>::Type        Type;
};

template <typename TString, typename TSpec>
struct Fibre<SparseString<ContainerView<TString>, TSpec> const, FibreValues>
{
    typedef typename View<typename Fibre<SparseString<TString, TSpec> const, FibreValues>::Type>::Type  Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function getFibre()                                             [Index View]
// ----------------------------------------------------------------------------

template <typename TText, typename TViewSpec, typename TSpec>
SEQAN_FUNC typename Fibre<Index<ContainerView<TText, TViewSpec>, TSpec>, FibreText>::Type &
getFibre(Index<ContainerView<TText, TViewSpec>, TSpec> & index, FibreText)
{
    return index.text;
}

template <typename TText, typename TViewSpec, typename TSpec>
SEQAN_FUNC typename Fibre<Index<ContainerView<TText, TViewSpec>, TSpec> const, FibreText>::Type &
getFibre(Index<ContainerView<TText, TViewSpec>, TSpec> const & index, FibreText)
{
    return index.text;
}

template <typename TText, typename TViewSpec, typename TSpec>
SEQAN_FUNC typename Fibre<Index<ContainerView<TText, TViewSpec>, TSpec>, FibreRawText>::Type &
getFibre(Index<ContainerView<TText, TViewSpec>, TSpec> & index, FibreRawText)
{
    return concat(getFibre(index, FibreText()));
}

template <typename TText, typename TViewSpec, typename TSpec>
SEQAN_FUNC typename Fibre<Index<ContainerView<TText, TViewSpec>, TSpec> const, FibreRawText>::Type &
getFibre(Index<ContainerView<TText, TViewSpec>, TSpec> const & index, FibreRawText)
{
    return concat(getFibre(index, FibreText()));
}

// ----------------------------------------------------------------------------
// Function indexText()                                            [Index View]
// ----------------------------------------------------------------------------

template <typename TText, typename TViewSpec, typename TSpec>
SEQAN_FUNC typename Fibre<Index<ContainerView<TText, TViewSpec>, TSpec>, FibreText>::Type &
indexText(Index<ContainerView<TText, TViewSpec>, TSpec> & index)
{
    return getFibre(index, FibreText());
}

template <typename TText, typename TViewSpec, typename TSpec>
SEQAN_FUNC typename Fibre<Index<ContainerView<TText, TViewSpec>, TSpec> const, FibreText>::Type &
indexText(Index<ContainerView<TText, TViewSpec>, TSpec> const & index)
{
    return getFibre(index, FibreText());
}

// ----------------------------------------------------------------------------
// Function indexRawText()                                         [Index View]
// ----------------------------------------------------------------------------

template <typename TText, typename TViewSpec, typename TSpec>
SEQAN_FUNC typename Fibre<Index<ContainerView<TText, TViewSpec>, TSpec>, FibreRawText>::Type &
indexRawText(Index<ContainerView<TText, TViewSpec>, TSpec> & index)
{
    return getFibre(index, FibreRawText());
}

template <typename TText, typename TViewSpec, typename TSpec>
SEQAN_FUNC typename Fibre<Index<ContainerView<TText, TViewSpec>, TSpec> const, FibreRawText>::Type &
indexRawText(Index<ContainerView<TText, TViewSpec>, TSpec> const & index)
{
    return getFibre(index, FibreRawText());
}

// ----------------------------------------------------------------------------
// Function indexRequire()                                         [Index View]
// ----------------------------------------------------------------------------

template <typename TText, typename TViewSpec, typename TSpec, typename TFibre>
SEQAN_FUNC bool indexRequire(Index<ContainerView<TText, TViewSpec>, TSpec> & index, Tag<TFibre> const fibre)
{
    bool supplied = indexSupplied(index, fibre);
    SEQAN_ASSERT_MSG(supplied, "Fibre must be supplied on a view.");
    return supplied;
}

// ----------------------------------------------------------------------------
// Function indexRequire()                                       [FMIndex View]
// ----------------------------------------------------------------------------

template <typename TText, typename TViewSpec, typename TOccSpec, typename TSpec, typename TFibre>
SEQAN_FUNC bool
indexRequire(Index<ContainerView<TText, TViewSpec>, FMIndex<TOccSpec, TSpec> > & index, Tag<TFibre> const fibre)
{
    bool supplied = indexSupplied(index, fibre);
    SEQAN_ASSERT_MSG(supplied, "Fibre must be supplied on a view.");
    return supplied;
}

// ----------------------------------------------------------------------------
// Function indexCreate()                                          [Index View]
// ----------------------------------------------------------------------------

template <typename TText, typename TViewSpec, typename TSpec, typename TFibre>
SEQAN_FUNC bool indexCreate(Index<ContainerView<TText, TViewSpec>, TSpec> & /* index */, Tag<TFibre> const /* fibre */)
{
    SEQAN_ASSERT_MSG(false, "Fibre cannot be created on a view.");
    return false;
}

// ----------------------------------------------------------------------------
// Function indexCreate()                                        [FMIndex View]
// ----------------------------------------------------------------------------

template <typename TText, typename TViewSpec, typename TOccSpec, typename TSpec, typename TFibre>
SEQAN_FUNC bool
indexCreate(Index<ContainerView<TText, TViewSpec>, FMIndex<TOccSpec, TSpec> > & /* index */, Tag<TFibre> const /* fibre */)
{
    SEQAN_ASSERT_MSG(false, "Fibre cannot be created on a view.");
    return false;
}

// ----------------------------------------------------------------------------
// Function view()                                                      [Index]
// ----------------------------------------------------------------------------
// NOTE(esiragusa): view() of a generic Index is not useful and potentially dangerous.

//template <typename TText, typename TSpec>
//Index<ContainerView<TText>, TSpec>
//view(Index<TText, TSpec> & index)
//{
//    Index<ContainerView<TText>, TSpec> indexView;
//
//    indexText(indexView) = view(indexText(index));
//
//    return indexView;
//}

// ----------------------------------------------------------------------------
// Function view()                                                    [IndexSa]
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec>
typename View<Index<TText, IndexSa<TSpec> > >::Type
view(Index<TText, IndexSa<TSpec> > & index)
{
    typename View<Index<TText, IndexSa<TSpec> > >::Type indexView;

    indexText(indexView) = view(indexText(index));
    indexSA(indexView) = view(indexSA(index));

    return indexView;
}

// ----------------------------------------------------------------------------
// Function view()                                                   [IndexEsa]
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec>
typename View<Index<TText, IndexEsa<TSpec> > >::Type
view(Index<TText, IndexEsa<TSpec> > & index)
{
    typename View<Index<TText, IndexEsa<TSpec> > >::Type indexView;

    indexText(indexView) = view(indexText(index));
    indexSA(indexView) = view(indexSA(index));
    indexLcp(indexView) = view(indexLcp(index));
    indexChildtab(indexView) = view(indexChildtab(index));
    // TODO(esiragusa): View of cargo?

    return indexView;
}

// ----------------------------------------------------------------------------
// Function view()                                                    [FMIndex]
// ----------------------------------------------------------------------------

template <typename TText, typename TOccSpec, typename TSpec>
typename View<Index<TText, FMIndex<TOccSpec, TSpec> > >::Type
view(Index<TText, FMIndex<TOccSpec, TSpec> > & index)
{
    typename View<Index<TText, FMIndex<TOccSpec, TSpec> > >::Type indexView;

    indexText(indexView) = view(indexText(index));
    indexLF(indexView) = view(indexLF(index));
    indexSA(indexView) = view(indexSA(index));

    return indexView;
}

// ----------------------------------------------------------------------------
// Function view()                                                    [LfTable]
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec>
typename View<LfTable<TText, TSpec> >::Type
view(LfTable<TText, TSpec> & lfTable)
{
    typename View<LfTable<TText, TSpec> >::Type lfTableView;

    getFibre(lfTableView, FibrePrefixSum()) = view(getFibre(lfTable, FibrePrefixSum()));
    getFibre(lfTableView, FibreValues()) = view(getFibre(lfTable, FibreValues()));
    getFibre(lfTableView, FibreSentinels()) = getFibre(lfTable, FibreSentinels());
    lfTableView.sentinelSubstitute = lfTable.sentinelSubstitute;

    return lfTableView;
}

// NOTE(esiragusa): StringSet specialization
//template <typename TText, typename TSSetSpec, typename TSpec>
//typename View<LfTable<StringSet<TText, TSSetSpec>, TSpec> >::Type
//view(LfTable<StringSet<TText, TSSetSpec>, TSpec> & lfTable)
//{
//    typename View<LfTable<StringSet<TText, TSSetSpec>, TSpec> >::Type lfTableView;
//
//    getFibre(lfTableView, FibrePrefixSum()) = view(getFibre(lfTable, FibrePrefixSum()));
//    getFibre(lfTableView, FibreValues()) = view(getFibre(lfTable, FibreValues()));
//    getFibre(lfTableView, FibreSentinels()) = view(getFibre(lfTable, FibreSentinels()));
//    lfTableView.sentinelSubstitute = lfTable.sentinelSubstitute;
//
//    return lfTableView;
//}

// ----------------------------------------------------------------------------
// Function view()                                             [PrefixSumTable]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
typename View<PrefixSumTable<TValue, TSpec> >::Type
view(PrefixSumTable<TValue, TSpec> & pst)
{
    typename View<PrefixSumTable<TValue, TSpec> >::Type pstView;

    getFibre(pstView, FibreEntries()) = view(getFibre(pst, FibreEntries()));

    return pstView;
}

// ----------------------------------------------------------------------------
// Function view()                                             [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
typename View<RankDictionary<TwoLevels<TValue, TSpec> > >::Type
view(RankDictionary<TwoLevels<TValue, TSpec> > & dict)
{
    typename View<RankDictionary<TwoLevels<TValue, TSpec> > >::Type dictView;

    getFibre(dictView, FibreRanks()) = view(getFibre(dict, FibreRanks()));
    dictView._length = dict._length;

    return dictView;
}

// ----------------------------------------------------------------------------
// Function view()                                               [CompressedSA]
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec>
typename View<CompressedSA<TText, TSpec> >::Type
view(CompressedSA<TText, TSpec> & sa)
{
    typename View<CompressedSA<TText, TSpec> >::Type saView;

    getFibre(saView, FibreLF()) = view(getFibre(sa, FibreLF()));
    getFibre(saView, FibreSparseString()) = view(getFibre(sa, FibreSparseString()));

    return saView;
}

// ----------------------------------------------------------------------------
// Function view()                                               [SparseString]
// ----------------------------------------------------------------------------

template <typename TString, typename TSpec>
typename View<SparseString<TString, TSpec> >::Type
view(SparseString<TString, TSpec> & sparseString)
{
    typename View<SparseString<TString, TSpec> >::Type sparseStringView;

    getFibre(sparseStringView, FibreValues()) = view(getFibre(sparseString, FibreValues()));
    getFibre(sparseStringView, FibreIndicators()) = view(getFibre(sparseString, FibreIndicators()));
    sparseStringView._length = sparseString._length;

    return sparseStringView;
}

}

#endif  // #ifndef SEQAN_EXTRAS_INDEX_VIEW_H_
