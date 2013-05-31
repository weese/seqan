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

#ifndef SEQAN_HEADER_MISC_VIEW_CONCAT_DIRECT_H
#define SEQAN_HEADER_MISC_VIEW_CONCAT_DIRECT_H

namespace seqan {

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction View                                                [StringSet]
// ----------------------------------------------------------------------------
// TODO(esiragusa): Generic const version?

template <typename TString, typename TSpec>
struct View<StringSet<TString, TSpec> >
{
    typedef StringSet<typename View<TString>::Type, TSpec>      Type;
};

// ----------------------------------------------------------------------------
// Metafunction IsView                                         [StringSet View]
// ----------------------------------------------------------------------------

template <typename TString, typename TSpec>
struct IsView<StringSet<ContainerView<TString>, TSpec> > : public True {};

// ----------------------------------------------------------------------------
// Metafunction StringSetLimits                                [StringSet View]
// ----------------------------------------------------------------------------
// NOTE(esiragusa): Generic const version refers to this non-const one.

template <typename TString, typename TSpec>
struct StringSetLimits<StringSet<ContainerView<TString>, TSpec> >
{
    typedef typename View<typename StringSetLimits<StringSet<TString, TSpec> >::Type>::Type     Type;
};

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function _initStringSetLimits                [ConcatDirect StringSet View]
// --------------------------------------------------------------------------

template <typename TString, typename TSpec>
inline void _initStringSetLimits(StringSet<ContainerView<TString>, Owner<ConcatDirect<TSpec> > > & /* me */) {}

// ----------------------------------------------------------------------------
// Function view()                                     [ConcatDirect StringSet]
// ----------------------------------------------------------------------------

template <typename TString, typename TSpec>
typename View<StringSet<TString, Owner<ConcatDirect<TSpec> > > >::Type
view(StringSet<TString, Owner<ConcatDirect<TSpec> > > & stringSet)
{
    typename View<StringSet<TString, Owner<ConcatDirect<TSpec> > > >::Type stringSetView;

    concat(stringSetView) = view(concat(stringSet));
    stringSetLimits(stringSetView) = view(stringSetLimits(stringSet));

    return stringSetView;
}

// ----------------------------------------------------------------------------
// Function assign()                                   [ConcatDirect StringSet]
// ----------------------------------------------------------------------------

template <typename TString, typename TSpec, typename TString2, typename TSpec2>
void assign(StringSet<TString, Owner<ConcatDirect<TSpec> > > & stringSet,
            StringSet<TString2, Owner<ConcatDirect<TSpec2> > > & other)
{
    assign(concat(stringSet), concat(other));
    assign(stringSetLimits(stringSet), stringSetLimits(other));
}

}  // namespace seqan

#endif  // #ifndef SEQAN_HEADER_MISC_VIEW_CONCAT_DIRECT_H
