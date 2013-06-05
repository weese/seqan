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
// This file contains index type definitions.
// ==========================================================================

#ifndef SEQAN_EXTRAS_CUDAMAPPER_INDEX_H_
#define SEQAN_EXTRAS_CUDAMAPPER_INDEX_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index_extras.h>

using namespace seqan;

// ============================================================================
// Global Types
// ============================================================================

// ----------------------------------------------------------------------------
// Fragment Store Configuration
// ----------------------------------------------------------------------------

struct CUDAStoreConfig
{
    typedef String<Dna>             TReadSeq;
    typedef String<Dna>             TContigSeq;
    typedef Owner<ConcatDirect<> >  TContigSpec;

    typedef double                  TMean;
    typedef double                  TStd;
    typedef signed char             TMappingQuality;

    typedef void                    TReadStoreElementSpec;
    typedef Owner<ConcatDirect<> >  TReadSeqStoreSpec;
    typedef Alloc<>                 TReadNameSpec;
    typedef Owner<ConcatDirect<> >  TReadNameStoreSpec;
    typedef void                    TMatePairStoreElementSpec;
    typedef void                    TLibraryStoreElementSpec;
    typedef void                    TContigStoreElementSpec;
    typedef void                    TContigFileSpec;
    typedef void                    TAlignedReadStoreElementSpec;
    typedef Owner<ConcatDirect<> >  TAlignedReadTagStoreSpec;
    typedef void                    TAnnotationStoreElementSpec;
};

// ============================================================================
// Other Store Types
// ============================================================================

typedef StringSet<CUDAStoreConfig::TContigSeq, CUDAStoreConfig::TContigSpec>    TContigs;

// ============================================================================
// Index Types
// ============================================================================

// ----------------------------------------------------------------------------
// Suffix Array Value Type
// ----------------------------------------------------------------------------

namespace seqan {
template <>
struct SAValue<TContigs>
{
    typedef Pair<unsigned char, unsigned int, Pack> Type;
};
}

// ----------------------------------------------------------------------------
// FM Index Fibres
// ----------------------------------------------------------------------------

typedef FMIndex<>                           TGenomeIndexSpec;
typedef Index<TContigs, TGenomeIndexSpec>   TGenomeIndex;

// ----------------------------------------------------------------------------
// FM Index Size
// ----------------------------------------------------------------------------

namespace seqan {
template <typename TSpec>
struct Size<RankDictionary<TwoLevels<Dna, TSpec> > >
{
    typedef unsigned Type;
};

template <typename TSpec>
struct Size<RankDictionary<TwoLevels<bool, TSpec> > >
{
    typedef unsigned Type;
};
}

#endif  // #ifndef SEQAN_EXTRAS_CUDAMAPPER_INDEX_H_
