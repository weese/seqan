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
template <typename TText, typename TAlloc, typename TSpec>
struct Fibre<Index<thrust::device_vector<TText, TAlloc>, TSpec>, FibreSA>
{
    typedef thrust::device_vector<typename SAValue<Index<thrust::device_vector<TText, TAlloc>, TSpec> >::Type>  Type;
};
#endif

template <typename TText, typename TViewSpec, typename TSpec>
struct Fibre<Index<View<TText, TViewSpec>, TSpec>, FibreSA>
{
    typedef View<typename Fibre<Index<TText, TSpec>, FibreSA>::Type, TViewSpec> Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function indexRequire()
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TFibre>
inline bool indexRequire(Index<View<TText, ContainerView>, TSpec> & index, Tag<TFibre> const fibre)
{
    bool supplied = indexSupplied(index, fibre);
    SEQAN_ASSERT_MSG(supplied, "Fibre must be supplied on a view.");
    return supplied;
}

// ----------------------------------------------------------------------------
// Function indexCreate()
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TFibre>
inline bool indexCreate(Index<View<TText, ContainerView>, TSpec> & /* index */, Tag<TFibre> const /* fibre */)
{
    SEQAN_ASSERT_MSG(false, "Fibre cannot be created on a view.");
    return false;
}

// ----------------------------------------------------------------------------
// Function toView()
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec>
Index<View<TText, ContainerView>, TSpec>
toView(Index<TText, TSpec> & index)
{
    Index<View<TText, ContainerView>, TSpec> indexView;

    indexText(indexView) = toView(indexText(index));
    indexSA(indexView) = toView(indexSA(index));

    return indexView;
}

// ----------------------------------------------------------------------------
// toRawView()
// ----------------------------------------------------------------------------

#ifdef __CUDACC__
template <typename TText, typename TSpec>
Index<View<typename thrust::detail::pointer_traits<typename TText::pointer>::raw_pointer, IteratorView>, TSpec>
toRawView(Index<TText, TSpec> & index)
{
	typedef typename TText::pointer TIterator;
	typedef typename thrust::detail::pointer_traits<TIterator>::raw_pointer TRawPointer;

    Index<View<TRawPointer, IteratorView>, TSpec> indexRawView;
    indexText(indexRawView) = toRawView(indexText(index));
    indexSA(indexRawView) = toRawView(indexSA(index));

    return indexRawView;
}
#endif

}

#endif  // #ifndef SEQAN_EXTRAS_INDEX_VIEW_H_
