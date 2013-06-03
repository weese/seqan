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
// This file contains the Genome class.
// ==========================================================================

#ifndef SEQAN_EXTRAS_MASAI_GENOME_H_
#define SEQAN_EXTRAS_MASAI_GENOME_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

using namespace seqan;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class GenomeConfig
// ----------------------------------------------------------------------------

template <typename TFragStoreConfig_    = FragmentStoreConfig<> >
struct GenomeConfig
{
    typedef TFragStoreConfig_   TFragStoreConfig;
};

// ----------------------------------------------------------------------------
// Class Genome
// ----------------------------------------------------------------------------

template <typename TSpec = void>
struct Genome
{
    Holder<TFragmentStore>  _store;
    TContigs                contigs;
    String<TContigSeqSize>  contigsLength;

    Genome() :
        _store()
    {}

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
struct GenomeHost
{};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function setGenome()                                               [TObject]
// ----------------------------------------------------------------------------

template <typename TObject, typename TSpec>
inline void
setGenome(TObject & object, Genome<TSpec> const & genome)
{
    setValue(object.genome, genome);
}

// ----------------------------------------------------------------------------
// Function getGenome()                                               [TObject]
// ----------------------------------------------------------------------------

template <typename TObject>
inline typename GenomeHost<TObject>::Type &
getGenome(TObject const & object)
{
    return value(object.genome);
}

// ----------------------------------------------------------------------------
// Function load()                                                     [Genome]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TString>
bool load(Genome<TSpec> & genome, TString const & genomeFile)
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

template <typename TSpec>
void _updateContigs(Genome<TSpec> & genome)
{
    clear(genome.contigs);
    reserve(genome.contigs, length(value(genome._store).contigStore));

    for (TContigStoreSize contigId = 0; contigId < length(value(genome._store).contigStore); ++contigId)
        appendValue(genome.contigs, value(genome._store).contigStore[contigId].seq);
}

// ----------------------------------------------------------------------------
// Function _updateContigsLength()                                     [Genome]
// ----------------------------------------------------------------------------

template <typename TSpec>
void _updateContigsLength(Genome<TSpec> & genome)
{
    clear(genome.contigsLength);
    reserve(genome.contigsLength, length(genome.contigs), Exact());

    for (TContigStoreSize contigId = 0; contigId < length(genome.contigs); ++contigId)
        appendValue(genome.contigsLength, length(genome.contigs[contigId]));
}

// ----------------------------------------------------------------------------
// Function reverse()                                                  [Genome]
// ----------------------------------------------------------------------------

template <typename TSpec>
void reverse(Genome<TSpec> & genome)
{
    for (TContigStoreSize contigId = 0; contigId < length(value(genome._store).contigStore); ++contigId)
        reverse(value(genome._store).contigStore[contigId].seq);
}

// ----------------------------------------------------------------------------
// Function getContigs()                                               [Genome]
// ----------------------------------------------------------------------------

template <typename TSpec>
inline TContigs &
getContigs(Genome<TSpec> & genome)
{
    return genome.contigs;
}

// ----------------------------------------------------------------------------
// Function contigLength()                                             [Genome]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TContigId>
inline TContigSeqSize
contigLength(Genome<TSpec> const & genome, TContigId contigId)
{
    return genome.contigsLength[contigId];
}

// ----------------------------------------------------------------------------
// Function clear()                                                    [Genome]
// ----------------------------------------------------------------------------

template <typename TSpec>
void clear(Genome<TSpec> & genome)
{
    clear(genome.contigs);
    clearContigs(value(genome._store));
}

#endif  // #ifndef SEQAN_EXTRAS_MASAI_GENOME_H_
