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

// ==========================================================================
// Classes
// ==========================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Fibre
// ----------------------------------------------------------------------------

#ifdef __CUDACC__
template <typename TText, typename TAlloc, typename TSpec, typename TFibre>
struct Fibre<Index<thrust::device_vector<TText, TAlloc>, TSpec>, Tag<TFibre> const>
{
    typedef thrust::device_vector<typename Size<Index<thrust::device_vector<TText, TAlloc>, TSpec> >::Type>     Type;
};

template <typename TText, typename TAlloc, typename TSpec>
struct Fibre<Index<thrust::device_vector<TText, TAlloc>, TSpec>, FibreText>
{
    typedef thrust::device_vector<TText, TAlloc>    Type;
};

template <typename TText, typename TAlloc, typename TSpec>
struct Fibre<Index<thrust::device_vector<TText, TAlloc>, TSpec>, FibreRawText>
{
    typedef typename Concatenator<thrust::device_vector<TText, TAlloc> >::Type  Type;
};

template <typename TText, typename TAlloc, typename TSpec>
struct Fibre<Index<thrust::device_vector<TText, TAlloc>, TSpec>, FibreSA>
{
    typedef thrust::device_vector<typename SAValue<Index<thrust::device_vector<TText, TAlloc>, TSpec> >::Type>  Type;
};

template <typename TText, typename TAlloc, typename TSpec>
struct Fibre<Index<thrust::device_vector<TText, TAlloc>, TSpec>, FibreBwt>
{
    typedef thrust::device_vector<typename Value<Index<thrust::device_vector<TText, TAlloc>, TSpec> >::Type>    Type;
};
#endif

// ----------------------------------------------------------------------------
// Metafunction Fibre
// ----------------------------------------------------------------------------

template <typename TText, typename TViewSpec, typename TSpec, typename TFibre>
struct Fibre<Index<View<TText, TViewSpec>, TSpec>, Tag<TFibre> const>
{
    typedef View<typename Fibre<Index<TText, TSpec>, Tag<TFibre> const>::Type, TViewSpec> Type;
};

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
// Metafunction FibreTextMember_
// ----------------------------------------------------------------------------

template <typename TText, typename TViewSpec, typename TSpec>
struct FibreTextMember_<Index<View<TText, TViewSpec>, TSpec> >
{
    typedef Index<View<TText, TViewSpec>, TSpec>        TIndex_;
    typedef typename Fibre<TIndex_, FibreText>::Type    Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function getFibre()
// ----------------------------------------------------------------------------

template <typename TText, typename TViewSpec, typename TSpec>
SEQAN_FUNC typename Fibre<Index<View<TText, TViewSpec>, TSpec>, FibreText>::Type &
getFibre(Index<View<TText, TViewSpec>, TSpec> &index, FibreText) {
    return index.text;
}

template <typename TText, typename TViewSpec, typename TSpec>
SEQAN_FUNC typename Fibre<Index<View<TText, TViewSpec>, TSpec> const, FibreText>::Type &
getFibre(Index<View<TText, TViewSpec>, TSpec> const &index, FibreText) {
    return index.text;
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
// Function indexRequire()
// ----------------------------------------------------------------------------

template <typename TText, typename TViewSpec, typename TSpec, typename TFibre>
SEQAN_FUNC bool indexRequire(Index<View<TText, TViewSpec>, TSpec> & index, Tag<TFibre> const fibre)
{
    bool supplied = indexSupplied(index, fibre);
//    SEQAN_ASSERT_MSG(supplied, "Fibre must be supplied on a view.");
    return supplied;
}

// ----------------------------------------------------------------------------
// Function indexCreate()
// ----------------------------------------------------------------------------

template <typename TText, typename TViewSpec, typename TSpec, typename TFibre>
SEQAN_FUNC bool indexCreate(Index<View<TText, TViewSpec>, TSpec> & /* index */, Tag<TFibre> const /* fibre */)
{
//    SEQAN_ASSERT_MSG(false, "Fibre cannot be created on a view.");
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
    indexSA(indexView) = toView(indexSA(index));

    return indexView;
}

}

#endif  // #ifndef SEQAN_EXTRAS_INDEX_VIEW_H_
