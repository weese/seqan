// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2013 NVIDIA Corporation
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
//     * Neither the name of NVIDIA Corporation nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL NVIDIA CORPORATION BE LIABLE
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
// This file contains the Genome class.
// ==========================================================================

#ifndef SEQAN_EXTRAS_MASAI_GENOME_H_
#define SEQAN_EXTRAS_MASAI_GENOME_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

template <typename TObject>
struct Contigs;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class GenomeConfig
// ----------------------------------------------------------------------------

template <typename TFragStoreConfig_ = FragmentStoreConfig<> >
struct GenomeConfig
{
    typedef TFragStoreConfig_   TFragStoreConfig;
};

// ----------------------------------------------------------------------------
// Class Genome
// ----------------------------------------------------------------------------

template <typename TSpec = void, typename TConfig = GenomeConfig<> >
struct Genome
{
    typedef typename TConfig::TFragStoreConfig          TFragStoreConfig_;
    typedef FragmentStore<void, TFragStoreConfig_>      TFragmentStore_;

    Holder<TFragmentStore_>                 _store;
    typename Contigs<Genome>::Type          contigs;
    String<typename Size<Genome>::Type>     contigsLength;

    Genome() :
        _store()
    {}

    template <typename TFragmentStore>
    Genome(TFragmentStore & store) :
        _store(store)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction GenomeHost<T>::Type                                   [TObject]
// ----------------------------------------------------------------------------

template <typename TObject>
struct GenomeHost {};

template <typename TObject>
struct GenomeHost<TObject const>
{
    typedef typename GenomeHost<TObject>::Type const    Type;
};

// ----------------------------------------------------------------------------
// Metafunction Contigs                                               [TObject]
// ----------------------------------------------------------------------------

template <typename TObject>
struct Contigs {};

template <typename TObject>
struct Contigs<TObject const>
{
    typedef typename Contigs<TObject>::Type const   Type;
};

// ----------------------------------------------------------------------------
// Metafunction Contigs                                                [Genome]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
struct Contigs<Genome<TSpec, TConfig> >
{
    typedef typename TConfig::TFragStoreConfig          TFragStoreConfig_;
    typedef FragmentStore<TSpec, TFragStoreConfig_>     TFragmentStore_;
    typedef typename TFragStoreConfig_::TContigSeq      TContigSeq_;

    typedef StringSet<TContigSeq_, Dependent<> >        Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size                                                   [Genome]
// ----------------------------------------------------------------------------

namespace seqan {
template <typename TSpec, typename TConfig>
struct Size<Genome<TSpec, TConfig> > :
    Size<typename Contigs<Genome<TSpec, TConfig> >::Type> {};
}

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function setGenome()                                               [TObject]
// ----------------------------------------------------------------------------

template <typename TObject, typename TSpec, typename TConfig>
inline void
setGenome(TObject & object, Genome<TSpec, TConfig> const & genome)
{
    setValue(object.genome, genome);
}

// ----------------------------------------------------------------------------
// Function getGenome()                                               [TObject]
// ----------------------------------------------------------------------------

template <typename TObject>
inline typename GenomeHost<TObject>::Type &
getGenome(TObject & object)
{
    return value(object.genome);
}

// ----------------------------------------------------------------------------
// Function load()                                                     [Genome]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TString>
bool load(Genome<TSpec, TConfig> & genome, TString const & genomeFile)
{
    // TODO(esiragusa): Use record reader instead of loadContigs() from FragmentStore.
    if (!loadContigs(value(genome._store), genomeFile))
        return false;

    // Shrink contigs.
    for (unsigned contigId = 0; contigId < length(value(genome._store).contigStore); ++contigId)
        shrinkToFit(value(genome._store).contigStore[contigId].seq);

    _updateContigs(genome);
    _updateContigsLength(genome);

    return true;
}

// ----------------------------------------------------------------------------
// Function _updateContigs()                                           [Genome]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
void _updateContigs(Genome<TSpec, TConfig> & genome)
{
    clear(genome.contigs);
    reserve(genome.contigs, length(value(genome._store).contigStore));

    for (unsigned contigId = 0; contigId < length(value(genome._store).contigStore); ++contigId)
        appendValue(genome.contigs, value(genome._store).contigStore[contigId].seq);
}

// ----------------------------------------------------------------------------
// Function _updateContigsLength()                                     [Genome]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
void _updateContigsLength(Genome<TSpec, TConfig> & genome)
{
    clear(genome.contigsLength);
    reserve(genome.contigsLength, length(genome.contigs), Exact());

    for (unsigned contigId = 0; contigId < length(genome.contigs); ++contigId)
        appendValue(genome.contigsLength, length(genome.contigs[contigId]));
}

// ----------------------------------------------------------------------------
// Function reverse()                                                  [Genome]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
void reverse(Genome<TSpec, TConfig> & genome)
{
    for (unsigned contigId = 0; contigId < length(value(genome._store).contigStore); ++contigId)
        reverse(value(genome._store).contigStore[contigId].seq);
}

// ----------------------------------------------------------------------------
// Function getContigs()                                               [Genome]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline typename Contigs<Genome<TSpec, TConfig> >::Type &
getContigs(Genome<TSpec, TConfig> & genome)
{
    return genome.contigs;
}

template <typename TSpec, typename TConfig>
inline typename Contigs<Genome<TSpec, TConfig> const>::Type &
getContigs(Genome<TSpec, TConfig> const & genome)
{
    return genome.contigs;
}

// ----------------------------------------------------------------------------
// Function contigLength()                                             [Genome]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TContigId>
inline typename Size<Genome<TSpec, TConfig> const>::Type
contigLength(Genome<TSpec, TConfig> const & genome, TContigId contigId)
{
    return genome.contigsLength[contigId];
}

// ----------------------------------------------------------------------------
// Function clear()                                                    [Genome]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
void clear(Genome<TSpec, TConfig> & genome)
{
    clear(genome.contigs);
    clearContigs(value(genome._store));
}

#endif  // #ifndef SEQAN_EXTRAS_MASAI_GENOME_H_
