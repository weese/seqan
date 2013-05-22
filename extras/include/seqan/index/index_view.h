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
// This file contains Index class specializations for Views.
// ==========================================================================

#ifndef SEQAN_EXTRAS_INDEX_VIEW_H_
#define SEQAN_EXTRAS_INDEX_VIEW_H_

namespace seqan {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Classes
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Fibre                                            [Device Index]
// ----------------------------------------------------------------------------

#ifdef __CUDACC__
template <typename TValue, typename TAlloc, typename TSpec, typename TFibre>
struct Fibre<Index<thrust::device_vector<TValue, TAlloc>, TSpec>, Tag<TFibre> const>
{
    typedef thrust::device_vector<typename Size<Index<thrust::device_vector<TValue, TAlloc>, TSpec> >::Type>        Type;
};

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
    typedef thrust::device_vector<typename SAValue<Index<thrust::device_vector<TValue, TAlloc>, TSpec> >::Type>     Type;
};

template <typename TValue, typename TAlloc, typename TSpec>
struct Fibre<Index<thrust::device_vector<TValue, TAlloc>, TSpec>, FibreBwt>
{
    typedef thrust::device_vector<typename Value<Index<thrust::device_vector<TValue, TAlloc>, TSpec> >::Type>       Type;
};
#endif

// ----------------------------------------------------------------------------
// Metafunction Fibre                                          [Device FMIndex]
// ----------------------------------------------------------------------------

#ifdef __CUDACC__
template <typename TValue, typename TAlloc, typename TOccSpec, typename TSpec>
struct Fibre<Index<thrust::device_vector<TValue, TAlloc>, FMIndex<TOccSpec, TSpec> >, FibreSA>
{
    typedef Index<thrust::device_vector<TValue, TAlloc>, FMIndex<TOccSpec, TSpec> >     TIndex_;
    typedef typename SAValue<TIndex_>::Type                                             TSAValue_;
    typedef SparseString<thrust::device_vector<TSAValue_, TAlloc>, TSpec>               TSparseString_;
    typedef typename Fibre<TIndex_, FibreLF>::Type                                      TLF_;

    typedef CompressedSA<TSparseString_, TLF_, TSpec>                                   Type;
};
//
//template <typename TValue, typename TAlloc, typename TOccSpec, typename TSpec>
//struct Fibre<Index<thrust::device_vector<TValue, TAlloc>, FMIndex<TOccSpec, TSpec> >, FibreLF>
//{
//    typedef thrust::device_vector<TValue, TAlloc>                                       TText_;
//    typedef LfTable<TText_, TSpec>                                                      Type;
//};
//
//// TODO(esiragusa): Remove this TL specialization.
//template <typename TValue, typename TAlloc, typename TSpec>
//struct Fibre<Index<thrust::device_vector<TValue, TAlloc>, FMIndex<TL<TSpec>, TSpec> >, FibreLF>
//{
//    typedef thrust::device_vector<TValue, TAlloc>                                       TText_;
//    typedef LfTable<TText_, TSpec>                                                      Type;
//};
#endif

// ----------------------------------------------------------------------------
// Metafunction Fibre                                              [Index View]
// ----------------------------------------------------------------------------

//template <typename TText, typename TViewSpec, typename TSpec, typename TFibre>
//struct Fibre<Index<View<TText, TViewSpec>, TSpec>, Tag<TFibre> const>
//{
//    typedef View<typename Fibre<Index<TText, TSpec>, Tag<TFibre> const>::Type, TViewSpec> Type;
//};

template <typename TText, typename TViewSpec, typename TSpec>
struct Fibre<Index<View<TText, TViewSpec>, TSpec>, FibreText>
{
    typedef View<TText, TViewSpec>  Type;
};

template <typename TText, typename TViewSpec, typename TSpec>
struct Fibre<Index<View<TText, TViewSpec>, TSpec>, FibreRawText>
{
    typedef View<typename Fibre<Index<TText, TSpec>, FibreRawText>::Type, TViewSpec> Type;
};

template <typename TText, typename TViewSpec, typename TSpec>
struct Fibre<Index<View<TText, TViewSpec>, TSpec>, FibreSA>
{
    typedef View<typename Fibre<Index<TText, TSpec>, FibreSA>::Type, TViewSpec> Type;
};

template <typename TText, typename TViewSpec, typename TSpec>
struct Fibre<Index<View<TText, TViewSpec>, TSpec>, FibreBwt>
{
    typedef View<typename Fibre<Index<TText, TSpec>, FibreBwt>::Type, TViewSpec> Type;
};


// ----------------------------------------------------------------------------
// Metafunction FibreTextMember_                                   [Index View]
// ----------------------------------------------------------------------------

template <typename TText, typename TViewSpec, typename TSpec>
struct FibreTextMember_<Index<View<TText, TViewSpec>, TSpec> >
{
    typedef Index<View<TText, TViewSpec>, TSpec>        TIndex_;
    typedef typename Fibre<TIndex_, FibreText>::Type    Type;
};

// ----------------------------------------------------------------------------
// Metafunction Fibre                                            [FMIndex View]
// ----------------------------------------------------------------------------

template <typename TText, typename TViewSpec, typename TOccSpec, typename TSpec>
struct Fibre<Index<View<TText, TViewSpec>, FMIndex<TOccSpec, TSpec> >, FibreSA>
{
    typedef Index<View<TText, TViewSpec>, FMIndex<TOccSpec, TSpec> >    TIndex_;

    typedef typename SAValue<TIndex_>::Type                             TSAValue_;
    typedef SparseString<View<String<TSAValue_> >, TSpec>               TSparseString_;

    typedef typename Fibre<TIndex_, FibreLF>::Type                      TLF_;

    typedef CompressedSA<TSparseString_, TLF_, TSpec>                   Type;
};

//template <typename TText, typename TViewSpec, typename TOccSpec, typename TSpec>
//struct Fibre<Index<View<TText, TViewSpec>, FMIndex<TOccSpec, TSpec> >, FibreLF>
//{
//    typedef View<typename Fibre<Index<TText, FMIndex<TOccSpec, TSpec> >, FibreLF>::Type, TViewSpec> Type;
//};
//
//// TODO(esiragusa): Remove this TL specialization.
//template <typename TText, typename TViewSpec, typename TSpec>
//struct Fibre<Index<View<TText, TViewSpec>, FMIndex<TL<TSpec>, TSpec> >, FibreLF>
//{
//    typedef View<typename Fibre<Index<TText, FMIndex<TL<TSpec>, TSpec> >, FibreLF>::Type, TViewSpec> Type;
//};

// ----------------------------------------------------------------------------
// Metafunction Fibre                                            [LfTable View]
// ----------------------------------------------------------------------------

template <typename TText, typename TViewSpec, typename TSpec>
struct Fibre<LfTable<View<TText, TViewSpec>, TSpec>, FibrePrefixSum>
{
    typedef View<typename Fibre<LfTable<TText, TSpec>, FibrePrefixSum>::Type>    Type;
};

template <typename TText, typename TViewSpec, typename TSpec>
struct Fibre<LfTable<View<TText, TViewSpec>, TSpec>, FibreValues>
{
    typedef View<typename Fibre<LfTable<TText, TSpec>, FibreValues>::Type>    Type;
};

template <typename TText, typename TViewSpec, typename TSpec>
struct Fibre<LfTable<View<TText, TViewSpec>, TSpec>, FibreSentinels>
{
    typedef View<typename Fibre<LfTable<TText, TSpec>, FibreSentinels>::Type>    Type;
};

// ----------------------------------------------------------------------------
// Metafunction Fibre                                     [PrefixSumTable View]
// ----------------------------------------------------------------------------

template <typename TChar, typename TSpec, typename TViewSpec>
struct Fibre<View<PrefixSumTable<TChar, TSpec>, TViewSpec>, FibreEntries>
{
    typedef View<typename Fibre<PrefixSumTable<TChar, TSpec>, FibreEntries>::Type>   Type;
};

// ----------------------------------------------------------------------------
// Metafunction Fibre                                     [RankDictionary View]
// ----------------------------------------------------------------------------

template <typename TSpec>
struct Fibre<View<RankDictionary<TSpec> >, FibreRanks>
{
    typedef View<typename Fibre<RankDictionary<TSpec>, FibreRanks>::Type>   Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function getFibre()
// ----------------------------------------------------------------------------

template <typename TText, typename TViewSpec, typename TSpec>
SEQAN_FUNC typename Fibre<Index<View<TText, TViewSpec>, TSpec>, FibreText>::Type &
getFibre(Index<View<TText, TViewSpec>, TSpec> & index, FibreText)
{
    return index.text;
}

template <typename TText, typename TViewSpec, typename TSpec>
SEQAN_FUNC typename Fibre<Index<View<TText, TViewSpec>, TSpec> const, FibreText>::Type &
getFibre(Index<View<TText, TViewSpec>, TSpec> const & index, FibreText)
{
    return index.text;
}

template <typename TText, typename TViewSpec, typename TSpec>
SEQAN_FUNC typename Fibre<Index<View<TText, TViewSpec>, TSpec>, FibreRawText>::Type &
getFibre(Index<View<TText, TViewSpec>, TSpec> & index, FibreRawText)
{
    return concat(getFibre(index, FibreText()));
}

template <typename TText, typename TViewSpec, typename TSpec>
SEQAN_FUNC typename Fibre<Index<View<TText, TViewSpec>, TSpec> const, FibreRawText>::Type &
getFibre(Index<View<TText, TViewSpec>, TSpec> const & index, FibreRawText)
{
    return concat(getFibre(index, FibreText()));
}

// ----------------------------------------------------------------------------
// Function indexText()
// ----------------------------------------------------------------------------

template <typename TText, typename TViewSpec, typename TSpec>
SEQAN_FUNC typename Fibre<Index<View<TText, TViewSpec>, TSpec>, FibreText>::Type &
indexText(Index<View<TText, TViewSpec>, TSpec> & index)
{
    return getFibre(index, FibreText());
}

template <typename TText, typename TViewSpec, typename TSpec>
SEQAN_FUNC typename Fibre<Index<View<TText, TViewSpec>, TSpec> const, FibreText>::Type &
indexText(Index<View<TText, TViewSpec>, TSpec> const & index)
{
    return getFibre(index, FibreText());
}

// ----------------------------------------------------------------------------
// Function indexRawText()
// ----------------------------------------------------------------------------

template <typename TText, typename TViewSpec, typename TSpec>
SEQAN_FUNC typename Fibre<Index<View<TText, TViewSpec>, TSpec>, FibreRawText>::Type &
indexRawText(Index<View<TText, TViewSpec>, TSpec> & index)
{
    return getFibre(index, FibreRawText());
}

template <typename TText, typename TViewSpec, typename TSpec>
SEQAN_FUNC typename Fibre<Index<View<TText, TViewSpec>, TSpec> const, FibreRawText>::Type &
indexRawText(Index<View<TText, TViewSpec>, TSpec> const & index)
{
    return getFibre(index, FibreRawText());
}

// ----------------------------------------------------------------------------
// Function indexRequire()
// ----------------------------------------------------------------------------

template <typename TText, typename TViewSpec, typename TSpec, typename TFibre>
SEQAN_FUNC bool indexRequire(Index<View<TText, TViewSpec>, TSpec> & index, Tag<TFibre> const fibre)
{
    bool supplied = indexSupplied(index, fibre);
    SEQAN_ASSERT_MSG(supplied, "Fibre must be supplied on a view.");
    return supplied;
}

template <typename TText, typename TViewSpec, typename TOccSpec, typename TSpec, typename TFibre>
SEQAN_FUNC bool indexRequire(Index<View<TText, TViewSpec>, FMIndex<TOccSpec, TSpec> > & index, Tag<TFibre> const fibre)
{
    bool supplied = indexSupplied(index, fibre);
    SEQAN_ASSERT_MSG(supplied, "Fibre must be supplied on a view.");
    return supplied;
}

// ----------------------------------------------------------------------------
// Function indexCreate()
// ----------------------------------------------------------------------------

template <typename TText, typename TViewSpec, typename TSpec, typename TFibre>
SEQAN_FUNC bool indexCreate(Index<View<TText, TViewSpec>, TSpec> & /* index */, Tag<TFibre> const /* fibre */)
{
    SEQAN_ASSERT_MSG(false, "Fibre cannot be created on a view.");
    return false;
}

template <typename TText, typename TViewSpec, typename TOccSpec, typename TSpec, typename TFibre>
SEQAN_FUNC bool indexCreate(Index<View<TText, TViewSpec>, FMIndex<TOccSpec, TSpec> > & /* index */, Tag<TFibre> const /* fibre */)
{
    SEQAN_ASSERT_MSG(false, "Fibre cannot be created on a view.");
    return false;
}

// ----------------------------------------------------------------------------
// Function toView()
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec>
Index<View<TText>, TSpec>
toView(Index<TText, TSpec> & index)
{
    Index<View<TText>, TSpec> indexView;

    indexText(indexView) = toView(indexText(index));

    return indexView;
}

template <typename TText, typename TSpec>
Index<View<TText>, IndexSa<TSpec> >
toView(Index<TText, IndexSa<TSpec> > & index)
{
    Index<View<TText>, IndexSa<TSpec> > indexView;

    indexText(indexView) = toView(indexText(index));
    indexSA(indexView) = toView(indexSA(index));

    return indexView;
}

template <typename TText, typename TSpec>
Index<View<TText>, IndexEsa<TSpec> >
toView(Index<TText, IndexEsa<TSpec> > & index)
{
    Index<View<TText>, IndexEsa<TSpec> > indexView;

    indexText(indexView) = toView(indexText(index));
    indexSA(indexView) = toView(indexSA(index));
    indexLcp(indexView) = toView(indexLcp(index));
    indexChildtab(indexView) = toView(indexChildtab(index));
    // TODO(esiragusa): View of cargo?

    return indexView;
}

}

#endif  // #ifndef SEQAN_EXTRAS_INDEX_VIEW_H_
