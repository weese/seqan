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
// This file contains the Seeder class.
// ==========================================================================

#ifndef SEQAN_EXTRAS_MASAI_SEEDER_H_
#define SEQAN_EXTRAS_MASAI_SEEDER_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/find.h>

#include "index.h"

using namespace seqan;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class SeederConfig
// ----------------------------------------------------------------------------

template <typename TDistance_       = HammingDistance,
          typename TAlgorithm_      = MultipleBacktracking,
          typename TReadsIndexSpec_ = TReadsWotdSpec>
struct SeederConfig
{
    typedef TDistance_          TDistance;
    typedef TAlgorithm_         TAlgorithm;
    typedef TReadsIndexSpec_    TReadsIndexSpec;
};

// ----------------------------------------------------------------------------
// Class Seeder
// ----------------------------------------------------------------------------

template <typename TReads, typename TManager, typename TDelegate, typename TSpec = void, typename TConfig = SeederConfig<> >
struct Seeder
{
    typedef ReadsIndex<TReads, typename TConfig::TReadsIndexSpec>   TReadsIndex;

    Holder<TReads>      reads;
    TReadsIndex         readsIndex;
    TManager            & manager;
    TDelegate           & delegate;
    unsigned long       hitsCount;

    Seeder(TManager & manager, TDelegate & delegate) :
        manager(manager),
        delegate(delegate),
        hitsCount(0)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction ReadsHost<T>::Type                                     [Seeder]
// ----------------------------------------------------------------------------

template <typename TReads, typename TManager, typename TDelegate, typename TSpec, typename TConfig>
struct ReadsHost<Seeder<TReads, TManager, TDelegate, TSpec, TConfig> >
{
    typedef TReads  Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function clear()                                                    [Seeder]
// ----------------------------------------------------------------------------

template <typename TReads, typename TManager, typename TDelegate, typename TSpec, typename TConfig>
void clear(Seeder<TReads, TManager, TDelegate, TSpec, TConfig> & seeder)
{
    clear(seeder.reads);
    clear(seeder.readsIndex);
    seeder.hitsCount = 0;
}

// ----------------------------------------------------------------------------
// Function setReads()                                                 [Seeder]
// ----------------------------------------------------------------------------

template <typename TReads, typename TManager, typename TDelegate, typename TSpec, typename TConfig>
void setReads(Seeder<TReads, TManager, TDelegate, TSpec, TConfig> & seeder, TReads & reads)
{
    clear(seeder);
    setValue(seeder.reads, reads);
    setReads(seeder.readsIndex, reads);
}

// ----------------------------------------------------------------------------
// Function find()                                                     [Seeder]
// ----------------------------------------------------------------------------

template <typename TReads, typename TManager, typename TDelegate, typename TSpec, typename TConfig, typename TGenomeIndex>
void find(Seeder<TReads, TManager, TDelegate, TSpec, TConfig> & seeder,
          TGenomeIndex & genomeIndex,
          TReadSeqSize seedsLength,
          TReadSeqSize errorsPerSeed,
          TReadSeqSize firstSeed,
          TReadSeqSize lastSeed)
{
    find(seeder, genomeIndex, seedsLength, firstSeed, lastSeed, errorsPerSeed, typename TConfig::TDistance(), typename TConfig::TAlgorithm());
}

// ----------------------------------------------------------------------------
// Function find()                               [Seeder<MultipleBacktracking>]
// ----------------------------------------------------------------------------

template <typename TReads, typename TManager, typename TDelegate, typename TSpec, typename TConfig, typename TGenomeIndex, typename TDistance>
void find(Seeder<TReads, TManager, TDelegate, TSpec, TConfig> & seeder,
          TGenomeIndex & genomeIndex,
          TReadSeqSize seedsLength,
          TReadSeqSize firstSeed,
          TReadSeqSize lastSeed,
          TReadSeqSize errorsPerSeed,
          TDistance const & /* tag */,
          MultipleBacktracking const & /* tag */)
{
    typedef Seeder<TReads, TManager, TDelegate, TSpec, TConfig>     TSeeder;
    typedef typename TSeeder::TReadsIndex::TIndex                   TReadsIndex;
    typedef Backtracking<TDistance>                                 TBacktracking;
    typedef Finder<TGenomeIndex, TBacktracking>                     TFinder;
    typedef Pattern<TReadsIndex, TBacktracking>                     TPattern;

    build(seeder.readsIndex, seeder.manager, seedsLength, firstSeed, lastSeed);

    TFinder finder(genomeIndex);
    TPattern pattern(seeder.readsIndex.index, seedsLength);

    while (find(finder, pattern, errorsPerSeed))
    {
        TReadSeqStoreSize seqId = getValueI1(position(pattern));

        // Skip disabled reads.
        if (isDisabled(seeder.manager, getReadId(getReads(seeder), seqId)))
            continue;

        ++seeder.hitsCount;

//        seeder.delegate.minErrorsPerRead = minErrors(seeder.manager, position(pattern).i1);
//        seeder.delegate.maxErrorsPerRead = maxErrors(seeder.manager, position(pattern).i1);

        onSeedHit(seeder.delegate,
                  getValueI1(position(finder)),
                  getValueI2(toSuffixPosition(host(finder), beginPosition(finder), length(finder))),
                  seqId,
                  getValueI2(beginPosition(pattern)),
                  pattern.prefix_aligner.errors);
    }
}

// ----------------------------------------------------------------------------
// Function find()                                 [Seeder<SingleBacktracking>]
// ----------------------------------------------------------------------------

template <typename TReads, typename TManager, typename TDelegate, typename TSpec, typename TConfig, typename TGenomeIndex, typename TDistance>
void find(Seeder<TReads, TManager, TDelegate, TSpec, TConfig> & seeder,
          TGenomeIndex & genomeIndex,
          TReadSeqSize seedsLength,
          TReadSeqSize firstSeed,
          TReadSeqSize lastSeed,
          TReadSeqSize errorsPerSeed,
          TDistance const & /* tag */,
          SingleBacktracking const & /* tag */)
{
    typedef Backtracking<TDistance>                 TBacktracking;
    typedef Finder<TGenomeIndex, TBacktracking>     TFinder;
    typedef Pattern<TReadSeq, TBacktracking>        TPattern;

    TFinder finder(genomeIndex);
    TPattern pattern;

    TReadSeqStoreSize seqsCount = length(getSeqs(getReads(seeder)));

    for (TReadSeqStoreSize seqId = 0; seqId < seqsCount; ++seqId)
    {
        // Skip disabled reads.
        if (isDisabled(seeder.manager, getReadId(getReads(seeder), seqId)))
            continue;

        TReadSeq & read = getSeqs(getReads(seeder))[seqId];

        for (TReadSeqSize seed = firstSeed; seed < lastSeed; ++seed)
        {
            setHost(pattern, infix(read, seedsLength * seed, seedsLength * (seed + 1)));

            clear(finder);

            while (find(finder, pattern, errorsPerSeed))
            {
                ++seeder.hitsCount;

//                seeder.delegate.minErrorsPerRead = minErrors(seeder.manager, seqId);
//                seeder.delegate.maxErrorsPerRead = maxErrors(seeder.manager, seqId);

                onSeedHit(seeder.delegate,
                          getValueI1(position(finder)),
                          getValueI2(toSuffixPosition(host(finder), beginPosition(finder), length(finder))),
                          seqId,
                          seedsLength * seed,
                          pattern.prefix_aligner.errors);
            }
        }
    }
}

// ----------------------------------------------------------------------------
// Function find()                         [Seeder<MultipleBacktracking,Exact>]
// ----------------------------------------------------------------------------

template <typename TReads, typename TManager, typename TDelegate, typename TSpec, typename TConfig, typename TGenomeIndex>
void find(Seeder<TReads, TManager, TDelegate, TSpec, TConfig> & seeder,
          TGenomeIndex & genomeIndex,
          TReadSeqSize seedsLength,
          TReadSeqSize firstSeed,
          TReadSeqSize lastSeed,
          TReadSeqSize /* errorsPerSeed */,
          Nothing const & /* tag */,
          MultipleBacktracking const & /* tag */)
{
    typedef Seeder<TReads, TManager, TDelegate, TSpec, TConfig>     TSeeder;
    typedef typename TSeeder::TReadsIndex::TIndex                   TReadsIndex;
    typedef Backtracking<HammingDistance, Stretched<> >             TBacktracking;
    typedef Finder<TGenomeIndex, TBacktracking>                     TFinder;
    typedef Pattern<TReadsIndex, TBacktracking>                     TPattern;

    build(seeder.readsIndex, seeder.manager, seedsLength, firstSeed, lastSeed);

    TFinder finder(genomeIndex);
    TPattern pattern(seeder.readsIndex.index, seedsLength);

    while (find(finder, pattern, 0u))
    {
        TReadSeqStoreSize seqId = getValueI1(position(pattern));

        // Skip disabled reads.
        if (isDisabled(seeder.manager, getReadId(getReads(seeder), seqId)))
            continue;

        ++seeder.hitsCount;

        onSeedHit(seeder.delegate,
                  getValueI1(position(finder)),
                  getValueI2(toSuffixPosition(host(finder), beginPosition(finder), length(finder))),
                  seqId,
                  getValueI2(beginPosition(pattern)),
                  0u);
    }
}

// ----------------------------------------------------------------------------
// Function find()                           [Seeder<SingleBacktracking,Exact>]
// ----------------------------------------------------------------------------

template <typename TReads, typename TManager, typename TDelegate, typename TSpec, typename TConfig, typename TGenomeIndex>
void find(Seeder<TReads, TManager, TDelegate, TSpec, TConfig> & seeder,
          TGenomeIndex & genomeIndex,
          TReadSeqSize seedsLength,
          TReadSeqSize firstSeed,
          TReadSeqSize lastSeed,
          TReadSeqSize /* errorsPerSeed */,
          Nothing const & /* tag */,
          SingleBacktracking const & /* tag */)
{
    typedef Finder<TGenomeIndex, FinderSTree>       TFinder;
    typedef Pattern<TReadSeq>                       TPattern;

    TFinder finder(genomeIndex);

    TReadSeqStoreSize seqsCount = length(getSeqs(getReads(seeder)));

    for (TReadSeqStoreSize seqId = 0; seqId < seqsCount; ++seqId)
    {
        // Skip disabled reads.
        if (isDisabled(seeder.manager, getReadId(getReads(seeder), seqId)))
            continue;

        TReadSeq & read = getSeqs(getReads(seeder))[seqId];

        for (TReadSeqSize seed = firstSeed; seed < lastSeed; ++seed)
        {
            TPattern pattern(infix(read, seedsLength * seed, seedsLength * (seed + 1)));

            clear(finder);

            while (find(finder, pattern))
            {
                ++seeder.hitsCount;

                onSeedHit(seeder.delegate,
                          getValueI1(position(finder)),
                          getValueI2(toSuffixPosition(host(finder), beginPosition(finder), length(finder))),
                          seqId,
                          seedsLength * seed,
                          0u);
            }
        }
    }
}

#endif  // #ifndef SEQAN_EXTRAS_MASAI_SEEDER_H_
