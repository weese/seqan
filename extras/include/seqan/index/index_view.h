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
    typedef Index<typename View<TText>::Type, TSpec>   Type;
};

// ----------------------------------------------------------------------------
// Metafunction View                                                  [LfTable]
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec>
struct View<LfTable<TText, TSpec> >
{
    typedef LfTable<typename View<TText>::Type, TSpec>   Type;
};

// ----------------------------------------------------------------------------
// Metafunction View                                             [CompressedSA]
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec>
struct View<CompressedSA<TText, TSpec> >
{
    typedef CompressedSA<typename View<TText>::Type, TSpec>   Type;
};

// ----------------------------------------------------------------------------
// Metafunction Device                                                  [Index]
// ----------------------------------------------------------------------------

#ifdef __CUDACC__
template <typename TText, typename TSpec>
struct Device<Index<TText, TSpec> >
{
    typedef Index<typename Device<TText>::Type, TSpec>  Type;
};
#endif

// ----------------------------------------------------------------------------
// Metafunction Fibre                                            [Device Index]
// ----------------------------------------------------------------------------

#ifdef __CUDACC__
// NOTE(esiragusa): This causes template resolution ambiguity.
//template <typename TValue, typename TAlloc, typename TSpec, typename TFibre>
//struct Fibre<Index<thrust::device_vector<TValue, TAlloc>, TSpec>, Tag<TFibre> const>
//{
//    typedef thrust::device_vector<typename Size<Index<thrust::device_vector<TValue, TAlloc>, TSpec> >::Type>  Type;
//};

template <typename TValue, typename TAlloc, typename TSpec>
struct Fibre<Index<thrust::device_vector<TValue, TAlloc>, TSpec>, FibreText>
{
    typedef thrust::device_vector<TValue, TAlloc>    Type;
};

template <typename TValue, typename TAlloc, typename TSpec>
struct Fibre<Index<thrust::device_vector<TValue, TAlloc>, TSpec>, FibreRawText>
{
    typedef typename Concatenator<thrust::device_vector<TValue, TAlloc> >::Type  Type;
};

template <typename TValue, typename TAlloc, typename TSpec>
struct Fibre<Index<thrust::device_vector<TValue, TAlloc>, TSpec>, FibreSA>
{
    typedef thrust::device_vector<typename SAValue<Index<thrust::device_vector<TValue, TAlloc>, TSpec> >::Type> Type;
};

template <typename TValue, typename TAlloc, typename TSpec>
struct Fibre<Index<thrust::device_vector<TValue, TAlloc>, TSpec>, FibreLcp>
{
    typedef thrust::device_vector<typename Size<Index<thrust::device_vector<TValue, TAlloc>, TSpec> >::Type>    Type;
};

template <typename TValue, typename TAlloc, typename TSpec>
struct Fibre<Index<thrust::device_vector<TValue, TAlloc>, TSpec>, FibreChildtab>
{
    typedef thrust::device_vector<typename Size<Index<thrust::device_vector<TValue, TAlloc>, TSpec> >::Type>    Type;
};

template <typename TValue, typename TAlloc, typename TSpec>
struct Fibre<Index<thrust::device_vector<TValue, TAlloc>, TSpec>, FibreBwt>
{
    typedef thrust::device_vector<typename Value<Index<thrust::device_vector<TValue, TAlloc>, TSpec> >::Type>   Type;
};
#endif

// ----------------------------------------------------------------------------
// Metafunction Fibre                                          [Device FMIndex]
// ----------------------------------------------------------------------------

#ifdef __CUDACC__
template <typename TValue, typename TAlloc, typename TOccSpec, typename TSpec>
struct Fibre<Index<thrust::device_vector<TValue, TAlloc>, FMIndex<TOccSpec, TSpec> >, FibreSA>
{
    typedef CompressedSA<thrust::device_vector<TValue, TAlloc>, TSpec>      Type;
};
#endif

// ----------------------------------------------------------------------------
// Metafunction Fibre                                          [Device LfTable]
// ----------------------------------------------------------------------------

#ifdef __CUDACC__
template <typename TValue, typename TAlloc, typename TSpec>
struct Fibre<LfTable<thrust::device_vector<TValue, TAlloc>, TSpec>, FibrePrefixSum>
{
    typedef thrust::device_vector<TValue, TAlloc>               TText_;
    typedef typename Value<LfTable<TText_, TSpec> >::Type       TValue_;
    typedef typename MakeUnsigned<TValue_>::Type                TUValue_;

    typedef PrefixSumTable<TUValue_, Device<TSpec> >            Type;
};

template <typename TValue, typename TAlloc, typename TSpec>
struct Fibre<LfTable<thrust::device_vector<TValue, TAlloc>, TSpec>, FibreValues>
{
    typedef thrust::device_vector<TValue, TAlloc>               TText_;
    typedef typename Value<LfTable<TText_, TSpec> >::Type       TValue_;

    typedef RankDictionary<TwoLevels<TValue_, Device<TSpec> > > Type;
};

//template <typename TText, typename TViewSpec, typename TSpec>
//struct Fibre<LfTable<thrust::device_vector<TValue, TAlloc>, TSpec>, FibreSentinels>
//{
//    typedef typename Fibre<LfTable<TText, TSpec>, FibreSentinels>::Type   Type;
//};
#endif

// ----------------------------------------------------------------------------
// Metafunction Fibre                                   [Device PrefixSumTable]
// ----------------------------------------------------------------------------

#ifdef __CUDACC__
template <typename TChar, typename TSpec>
struct Fibre<PrefixSumTable<TChar, Device<TSpec> >, FibreEntries>
{
    typedef thrust::device_vector<unsigned>         Type;
};
#endif

// ----------------------------------------------------------------------------
// Metafunction Fibre                                   [Device RankDictionary]
// ----------------------------------------------------------------------------

#ifdef __CUDACC__
template <typename TValue, typename TSpec>
struct Fibre<RankDictionary<TwoLevels<TValue, Device<TSpec> > >, FibreRanks>
{
    typedef TwoLevels<TValue, TSpec>        TRankDictionarySpec_;

    typedef thrust::device_vector<RankDictionaryEntry_<TRankDictionarySpec_> >  Type;
};
#endif

// ----------------------------------------------------------------------------
// Metafunction Fibre                                     [Device CompressedSA]
// ----------------------------------------------------------------------------

#ifdef __CUDACC__
template <typename TValue, typename TAlloc, typename TSpec>
struct Fibre<CompressedSA<thrust::device_vector<TValue, TAlloc>, TSpec>, FibreSparseString>
{
    typedef thrust::device_vector<TValue/*, TAlloc*/>       TText_;
    typedef typename SAValue<TText_>::Type                  TSAValue_;
    typedef thrust::device_vector<TSAValue_/*, TAlloc*/>    TString_;

    typedef SparseString<TString_, TSpec>                   Type;
};
#endif

// ----------------------------------------------------------------------------
// Metafunction Fibre                                     [Device SparseString]
// ----------------------------------------------------------------------------

//#ifdef __CUDACC__
//template <typename TValue, typename TAlloc, typename TSpec>
//struct Fibre<SparseString<thrust::device_vector<TValue, TAlloc>, TSpec>, FibreValues>
//{
//    typedef thrust::device_vector<TValue, TAlloc>       Type;
//};
//#endif

// ----------------------------------------------------------------------------
// Metafunction Fibre                                              [Index View]
// ----------------------------------------------------------------------------

// NOTE(esiragusa): This causes template resolution ambiguity.
//template <typename TText, typename TViewSpec, typename TSpec, typename TFibre>
//struct Fibre<Index<ContainerView<TText, TViewSpec>, TSpec>, Tag<TFibre> const>
//{
//    typedef typename View<typename Fibre<Index<TText, TSpec>, Tag<TFibre> const>::Type>::Type     Type;
//};

template <typename TText, typename TViewSpec, typename TSpec>
struct Fibre<Index<ContainerView<TText, TViewSpec>, TSpec>, FibreText>
{
    typedef typename View<TText>::Type  Type;
};

template <typename TText, typename TViewSpec, typename TSpec>
struct Fibre<Index<ContainerView<TText, TViewSpec>, TSpec>, FibreRawText>
{
    typedef typename View<typename Fibre<Index<TText, TSpec>, FibreRawText>::Type>::Type    Type;
};

template <typename TText, typename TViewSpec, typename TSpec>
struct Fibre<Index<ContainerView<TText, TViewSpec>, TSpec>, FibreSA>
{
    typedef typename View<typename Fibre<Index<TText, TSpec>, FibreSA>::Type>::Type         Type;
};

template <typename TText, typename TViewSpec, typename TSpec>
struct Fibre<Index<ContainerView<TText, TViewSpec>, TSpec>, FibreLcp>
{
    typedef typename View<typename Fibre<Index<TText, TSpec>, FibreLcp>::Type>::Type        Type;
};

template <typename TText, typename TViewSpec, typename TSpec>
struct Fibre<Index<ContainerView<TText, TViewSpec>, TSpec>, FibreChildtab>
{
    typedef typename View<typename Fibre<Index<TText, TSpec>, FibreChildtab>::Type>::Type   Type;
};

template <typename TText, typename TViewSpec, typename TSpec>
struct Fibre<Index<ContainerView<TText, TViewSpec>, TSpec>, FibreBwt>
{
    typedef typename View<typename Fibre<Index<TText, TSpec>, FibreBwt>::Type>::Type        Type;
};

// ----------------------------------------------------------------------------
// Metafunction FibreTextMember_                                   [Index View]
// ----------------------------------------------------------------------------

template <typename TText, typename TViewSpec, typename TSpec>
struct FibreTextMember_<Index<ContainerView<TText, TViewSpec>, TSpec> >
{
    typedef Index<ContainerView<TText, TViewSpec>, TSpec>       TIndex_;
    typedef typename Fibre<TIndex_, FibreText>::Type            Type;
};

// ----------------------------------------------------------------------------
// Metafunction Fibre                                            [FMIndex View]
// ----------------------------------------------------------------------------

template <typename TText, typename TViewSpec, typename TOccSpec, typename TSpec>
struct Fibre<Index<ContainerView<TText, TViewSpec>, FMIndex<TOccSpec, TSpec> >, FibreSA>
{
    typedef typename View<typename Fibre<Index<TText, FMIndex<TOccSpec, TSpec> >, FibreSA>::Type>::Type    Type;
};

// ----------------------------------------------------------------------------
// Metafunction Fibre                                            [LfTable View]
// ----------------------------------------------------------------------------

template <typename TText, typename TViewSpec, typename TSpec>
struct Fibre<LfTable<ContainerView<TText, TViewSpec>, TSpec>, FibrePrefixSum>
{
    typedef typename Fibre<LfTable<TText, ContainerView<TSpec, TViewSpec> >, FibrePrefixSum>::Type   Type;
};

template <typename TText, typename TViewSpec, typename TSpec>
struct Fibre<LfTable<ContainerView<TText, TViewSpec>, TSpec>, FibreValues>
{
    typedef typename Fibre<LfTable<TText, ContainerView<TSpec, TViewSpec> >, FibreValues>::Type      Type;
};

//template <typename TText, typename TViewSpec, typename TSpec>
//struct Fibre<LfTable<ContainerView<TText, TViewSpec>, TSpec>, FibreSentinels>
//{
//    typedef typename Fibre<LfTable<TText, TSpec>, FibreSentinels>::Type   Type;
//};

// NOTE(esiragusa): StringSet specialization
//template <typename TText, typename TSSetSpec, typename TViewSpec, typename TSpec>
//struct Fibre<LfTable<ContainerView<StringSet<TText, TSSetSpec>, TViewSpec>, TSpec>, FibreSentinels>
//{
//    typedef typename Fibre<LfTable<StringSet<TText, TSSetSpec>, ContainerView<TSpec, TViewSpec> >, FibreSentinels>::Type   Type;
//};

// ----------------------------------------------------------------------------
// Metafunction Fibre                                     [PrefixSumTable View]
// ----------------------------------------------------------------------------

template <typename TChar, typename TSpec, typename TViewSpec>
struct Fibre<PrefixSumTable<TChar, ContainerView<TSpec, TViewSpec> >, FibreEntries>
{
    typedef typename View<typename Fibre<PrefixSumTable<TChar, TSpec>, FibreEntries>::Type>::Type   Type;
};

// ----------------------------------------------------------------------------
// Metafunction Fibre                                     [RankDictionary View]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TViewSpec>
struct Fibre<RankDictionary<TwoLevels<TValue, ContainerView<TSpec, TViewSpec> > >, FibreRanks>
{
    typedef typename View<typename Fibre<RankDictionary<TwoLevels<TValue, TSpec> >, FibreRanks>::Type>::Type    Type;
};

// ----------------------------------------------------------------------------
// Metafunction Fibre                                       [CompressedSA View]
// ----------------------------------------------------------------------------

template <typename TText, typename TViewSpec, typename TSpec>
struct Fibre<CompressedSA<ContainerView<TText, TViewSpec>, TSpec>, FibreSparseString>
{
    typedef typename Fibre<CompressedSA<TText, ContainerView<TSpec, TViewSpec> >, FibreSparseString>::Type  Type;
};

// ----------------------------------------------------------------------------
// Metafunction Fibre                                       [SparseString View]
// ----------------------------------------------------------------------------

template <typename TFibreValues, typename TSpec>
struct Fibre<SparseString<TFibreValues, ContainerView<TSpec> >, FibreValues>
{
    typedef typename View<typename Fibre<SparseString<TFibreValues, TSpec>, FibreValues>::Type>::Type       Type;
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
PrefixSumTable<TValue, ContainerView<TSpec> >
view(PrefixSumTable<TValue, TSpec> & pst)
{
    PrefixSumTable<TValue, ContainerView<TSpec> > pstView;

    getFibre(pstView, FibreEntries()) = view(getFibre(pst, FibreEntries()));

    return pstView;
}

// ----------------------------------------------------------------------------
// Function view()                                             [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
RankDictionary<TwoLevels<TValue, ContainerView<TSpec> > >
view(RankDictionary<TwoLevels<TValue, TSpec> > & dict)
{
    RankDictionary<TwoLevels<TValue, ContainerView<TSpec> > > dictView;

    getFibre(dictView, FibreRanks()) = view(getFibre(dict, FibreRanks()));

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
SparseString<TString, ContainerView<TSpec> >
view(SparseString<TString, TSpec> & sparseString)
{
    SparseString<TString, ContainerView<TSpec> > sparseStringView;

    getFibre(sparseStringView, FibreValues()) = view(getFibre(sparseString, FibreValues()));
    getFibre(sparseStringView, FibreIndicators()) = view(getFibre(sparseString, FibreIndicators()));

    return sparseStringView;
}



// ----------------------------------------------------------------------------
// Function assign()                                                  [FMIndex]
// ----------------------------------------------------------------------------

template <typename TText, typename TOccSpec, typename TSpec, typename TText2, typename TOccSpec2, typename TSpec2>
inline void
assign(Index<TText, FMIndex<TOccSpec, TSpec> > & index, Index<TText2, FMIndex<TOccSpec2, TSpec2> > & source)
{
    assign(indexText(index), indexText(source));
    assign(indexSA(index), indexSA(source));
    assign(indexLF(index), indexLF(source));
}

template <typename TText, typename TSpec, typename TText2, typename TSpec2>
inline void
assign(LfTable<TText, TSpec> & lfTable, LfTable<TText2, TSpec2> & source)
{
    assign(getFibre(lfTable, FibrePrefixSum()), getFibre(source, FibrePrefixSum()));
    assign(getFibre(lfTable, FibreValues()), getFibre(source, FibreValues()));
    assign(getFibre(lfTable, FibreSentinels()), getFibre(source, FibreSentinels()));
    assign(lfTable.sentinelSubstitute, source.sentinelSubstitute);
}

template <typename TChar, typename TSpec, typename TChar2, typename TSpec2>
inline void
assign(PrefixSumTable<TChar, TSpec> & pst, PrefixSumTable<TChar2, TSpec2> & source)
{
    assign(getFibre(pst, FibreEntries()), getFibre(source, FibreEntries()));
}

template <typename TValue, typename TSpec, typename TValue2, typename TSpec2>
inline void
assign(RankDictionary<TwoLevels<TValue, TSpec> > & dict, RankDictionary<TwoLevels<TValue2, TSpec2> > & source)
{
    assign(getFibre(dict, FibreRanks()), getFibre(source, FibreRanks()));
}

template <typename TText, typename TSpec, typename TText2, typename TSpec2>
inline void
assign(CompressedSA<TText, TSpec> & sa, CompressedSA<TText2, TSpec2> & source)
{
    assign(getFibre(sa, FibreSparseString()), getFibre(source, FibreSparseString()));
    assign(getFibre(sa, FibreLF()), getFibre(source, FibreLF()));
}

template <typename TValueString, typename TSpec, typename TValueString2, typename TSpec2>
inline void
assign(SparseString<TValueString, TSpec> & sparseString, SparseString<TValueString2, TSpec2> & source)
{
    assign(getFibre(sparseString, FibreValues()), getFibre(source, FibreValues()));
    assign(getFibre(sparseString, FibreIndicators()), getFibre(source, FibreIndicators()));
}

}

#endif  // #ifndef SEQAN_EXTRAS_INDEX_VIEW_H_
