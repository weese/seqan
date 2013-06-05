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
// This file contains the Mapper class.
// ==========================================================================

#ifndef SEQAN_EXTRAS_MASAI_MAPPER_H_
#define SEQAN_EXTRAS_MASAI_MAPPER_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include "tags.h"
#include "store.h"
#include "seeder.h"
#include "extender.h"
#include "manager.h"

using namespace seqan;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class ReadMapperConfig
// ----------------------------------------------------------------------------

template <typename TDistance_       = EditDistance,
          typename TStrategy_       = AnyBest,
          typename TBacktracking_   = MultipleBacktracking,
          typename TGenomeConfig_   = GenomeConfig<> >
struct ReadMapperConfig
{
	typedef TDistance_          TDistance;
	typedef TStrategy_          TStrategy;
    typedef TBacktracking_      TBacktracking;
    typedef TGenomeConfig_      TGenomeConfig;
};

// ----------------------------------------------------------------------------
// Class Mapper
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec = void, typename TConfig = ReadMapperConfig<> >
struct Mapper
{
    typedef typename GenomeHost<Mapper>::Type                                                       TGenome;
    typedef Manager<TReads, TDelegate, typename TConfig::TStrategy>                                 TManager;
    typedef Extender<TGenome, TReads, TManager, typename TConfig::TDistance>                        TExtender;
    typedef SeederConfig<Nothing, typename TConfig::TBacktracking, TReadsQGramSpec>                 TSeederConfigExt;
    typedef SeederConfig<HammingDistance, typename TConfig::TBacktracking, TReadsWotdSpec>          TSeederConfigApx;
    typedef Seeder<TReads, TManager, TExtender, void, TSeederConfigExt>                             TSeederExt;
    typedef Seeder<TReads, TManager, TExtender, void, TSeederConfigApx>                             TSeederApx;

    Holder<TReads>      reads;
    TDelegate           & delegate;
    TManager            _manager;
    TExtender           _extender;
    TSeederExt          _seederExt;
    TSeederApx          _seederApx;
    TReadSeqSize        _seedLength;

    Mapper(TDelegate & delegate, bool disableExtender = false) :
        delegate(delegate),
        _manager(delegate),
        _extender(_manager, disableExtender),
        _seederExt(_manager, _extender),
        _seederApx(_manager, _extender),
        _seedLength(0)
    {}
};

// ----------------------------------------------------------------------------
// Class Seeding_
// ----------------------------------------------------------------------------

// TODO(esiragusa): Remove class Seeding_
template <typename TErrors = unsigned char, typename TSpec = void>
struct Seeding_
{
    typedef Pair<TReadSeqSize, TErrors> TSeed;
    typedef String<TSeed>               TSeeds;

    TSeeds          seeds;
    TReadSeqSize    seedLength;

    Seeding_(TReadSeqSize readLength, TErrors errors, TReadSeqSize seedLength) :
        seedLength(seedLength)
    {
        _computeSeeds(*this, readLength, errors);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction ReadsHost<T>::Type                                     [Mapper]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec, typename TConfig>
struct ReadsHost<Mapper<TReads, TDelegate, TSpec, TConfig> >
{
    typedef TReads  Type;
};

// ----------------------------------------------------------------------------
// Metafunction GenomeHost<T>::Type                                    [Mapper]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec, typename TConfig>
struct GenomeHost<Mapper<TReads, TDelegate, TSpec, TConfig> >
{
    typedef Genome<TSpec, typename TConfig::TGenomeConfig>  Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _computeSeeds()                                          [Seeding_]
// ----------------------------------------------------------------------------

template <typename TErrors, typename TSpec>
void _computeSeeds(Seeding_<TErrors, TSpec> & seeding, TReadSeqSize readLength, TErrors errors)
{
    typedef Seeding_<TErrors, TSpec>   TSeeding_;
    typedef typename TSeeding_::TSeed  TSeed;

    clear(seeding.seeds);

    TReadSeqSize seedsCount = readLength / seeding.seedLength;
    TErrors errorsPerSeed = errors / seedsCount;

    TReadSeqSize firstSeeds = (errors % seedsCount) + 1;

    std::cout << "Seeds:\t\t\t\t";

    if (errorsPerSeed > 0)
    {
        // Remaining seeds get errorsPerSeed - 1 errors.
        for (unsigned seed = 0; seed < seedsCount - firstSeeds; ++seed)
        {
            std::cout << "(" << seeding.seedLength << "," << (unsigned)(errorsPerSeed - 1) << ") ";
            appendValue(seeding.seeds, TSeed(seeding.seedLength, errorsPerSeed - 1));
        }
    }

    // First seeds get errorsPerSeed errors.
    for (unsigned seed = seedsCount - firstSeeds; seed < seedsCount; ++seed)
    {
        std::cout << "(" << seeding.seedLength << "," << (unsigned)errorsPerSeed << ") ";
        appendValue(seeding.seeds, TSeed(seeding.seedLength, errorsPerSeed));
    }

    std::cout << std::endl;
}

// ----------------------------------------------------------------------------
// Function setSeedLength()                                            [Mapper]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec, typename TConfig, typename TSize>
void setSeedLength(Mapper<TReads, TDelegate, TSpec, TConfig> & mapper, TSize seedLength)
{
    mapper._seedLength = seedLength;
    mapper._extender.seedLength = seedLength;
}

// ----------------------------------------------------------------------------
// Function setReads()                                                 [Mapper]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec, typename TConfig>
void setReads(Mapper<TReads, TDelegate, TSpec, TConfig> & mapper, TReads & reads)
{
    clear(mapper);

    setValue(mapper.reads, reads);
    setReads(mapper._seederExt, reads);
    setReads(mapper._seederApx, reads);
    setReads(mapper._extender, reads);
    setReads(mapper._manager, reads);
}

// ----------------------------------------------------------------------------
// Function clear()                                                    [Mapper]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec, typename TConfig>
void clear(Mapper<TReads, TDelegate, TSpec, TConfig> & mapper)
{
    clear(mapper.reads);
    clear(mapper._seederExt);
    clear(mapper._seederApx);
//    clear(mapper._extender);
    clear(mapper._manager);
}

// ----------------------------------------------------------------------------
// Function unmappedReads()                                            [Mapper]
// ----------------------------------------------------------------------------

//template <typename TReads, typename TDelegate, typename TSpec, typename TConfig>
//void unmappedReads(Mapper<TReads, TDelegate, TSpec, TConfig> & mapper)
//{
//}

// ----------------------------------------------------------------------------
// Function mapReads()                                                 [Mapper]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec, typename TConfig, typename TGenomeIndex, typename TErrors>
void mapReads(Mapper<TReads, TDelegate, TSpec, TConfig> & mapper, TGenomeIndex & genomeIndex, TErrors errors)
{
    setGenome(mapper._extender, getGenome(genomeIndex));
    _mapReads(mapper, genomeIndex, errors, typename TConfig::TStrategy());
}

// ----------------------------------------------------------------------------
// Function _mapReads()                                           [Mapper<All>]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec, typename TConfig, typename TGenomeIndex, typename TErrors>
void _mapReads(Mapper<TReads, TDelegate, TSpec, TConfig> & mapper, TGenomeIndex & genomeIndex, TErrors errors,
               All const & /*tag*/)
{
    typedef Seeding_<>                          TSeeding_;
    typedef TSeeding_::TSeeds                   TSeeds;
    typedef Iterator<TSeeds, Standard>::Type    TSeedsIterator;

    TSeeding_ seeding(avgSeqLength(getReads(mapper)), errors, mapper._seedLength);

    mapper._extender.minErrorsPerRead = 0;
    mapper._extender.maxErrorsPerRead = errors;

    unsigned position = 0;
    TSeedsIterator seedsEnd = end(seeding.seeds, Standard());

    std::cout << "Errors:\t\t\t\t" << errors << std::endl;

    for (TSeedsIterator seedsIt = begin(seeding.seeds, Standard()); seedsIt != seedsEnd; ++seedsIt)
    {
        _findSeeds(mapper, genomeIndex, *seedsIt, position, position + 1);

        // TODO(esiragusa):Compute minErrorsPerRead from seeds.
//        mapper._extender.minErrorsPerRead += getValueI2(*seedsIt);
        mapper._extender.minErrorsPerRead++;

        ++position;
    }

    std::cout << "Hits:\t\t\t\t" << hitsCount(mapper) << std::endl;
    std::cout << "Matches:\t\t\t" << matchesCount(mapper) << std::endl;
}

// ----------------------------------------------------------------------------
// Function _mapReads()                                       [Mapper<AllBest>]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec, typename TConfig, typename TGenomeIndex, typename TErrors>
void _mapReads(Mapper<TReads, TDelegate, TSpec, TConfig> & mapper, TGenomeIndex & genomeIndex, TErrors errors,
               AllBest const & /*tag*/)
{
    _mapReadsByStratum(mapper, genomeIndex, std::min(errors, (TErrors)1), AllBest());

    if (errors > 1)
    {
        mapper._manager.errors = 2;
        _mapReadsBySeed(mapper, genomeIndex, errors, AllBest());
    }
}

template <typename TReads, typename TDelegate, typename TSpec, typename TConfig, typename TGenomeIndex, typename TErrors>
void _mapReadsBySeed(Mapper<TReads, TDelegate, TSpec, TConfig> & mapper, TGenomeIndex & genomeIndex, TErrors errors,
                     AllBest const & /*tag*/)
{
    typedef Seeding_<>                          TSeeding_;
    typedef TSeeding_::TSeeds                   TSeeds;
    typedef Iterator<TSeeds, Standard>::Type    TSeedsIterator;

    TSeeding_ seeding(avgSeqLength(getReads(mapper)), errors, mapper._seedLength);

    mapper._extender.minErrorsPerRead = 0;
    mapper._extender.maxErrorsPerRead = errors;

    unsigned position = 0;
    TSeedsIterator seedsEnd = end(seeding.seeds, Standard());

    for (TSeedsIterator seedsIt = begin(seeding.seeds, Standard()); seedsIt != seedsEnd; ++seedsIt)
    {
        std::cout << "Errors:\t\t\t\t" << mapper._extender.minErrorsPerRead << std::endl;

        _findSeeds(mapper, genomeIndex, *seedsIt, position, position + 1);

        ++position;

        // TODO(esiragusa):Compute minErrorsPerRead from seeds.
//        mapper._extender.minErrorsPerRead += getValueI2(*seedsIt);
        mapper._extender.minErrorsPerRead++;

        std::cout << "Hits:\t\t\t\t" << hitsCount(mapper) << std::endl;
        std::cout << "Matches:\t\t\t" << matchesCount(mapper) << std::endl;

//        mapper._manager.errors += getValueI2(*seedsIt);
//        raiseErrorThreshold(mapper._manager);
        mapper._manager.errors = std::max((TErrors)(mapper._manager.errors), (TErrors)position);
    }
}

template <typename TReads, typename TDelegate, typename TSpec, typename TConfig, typename TGenomeIndex, typename TErrors>
void _mapReadsByStratum(Mapper<TReads, TDelegate, TSpec, TConfig> & mapper, TGenomeIndex & genomeIndex, TErrors errors,
                        AllBest const & /*tag*/)
{
    typedef Seeding_<>                          TSeeding_;
    typedef TSeeding_::TSeeds                   TSeeds;
    typedef Iterator<TSeeds, Standard>::Type    TSeedsIterator;

    mapper._extender.minErrorsPerRead = 0;
    mapper._extender.maxErrorsPerRead = 0;

    for (TErrors errors_ = 0; errors_ <= errors; ++errors_)
    {
        std::cout << "Errors:\t\t\t\t" << errors_ << std::endl;

        TReadSeqSize readsLength_ = avgSeqLength(getReads(mapper));
        TReadSeqSize seedLength_ = std::max(mapper._seedLength, readsLength_ / (errors_ + 1));

        TSeeding_ seeding(readsLength_, errors_, seedLength_);

        unsigned position = 0;
        TSeedsIterator seedsEnd = end(seeding.seeds, Standard());

        for (TSeedsIterator seedsIt = begin(seeding.seeds, Standard()); seedsIt != seedsEnd; ++seedsIt)
        {
            _findSeeds(mapper, genomeIndex, *seedsIt, position, position + 1);

            ++position;
        }

        mapper._extender.minErrorsPerRead++;
        mapper._extender.maxErrorsPerRead++;

        raiseErrorThreshold(mapper._manager);

        std::cout << "Hits:\t\t\t\t" << hitsCount(mapper) << std::endl;
        std::cout << "Matches:\t\t\t" << matchesCount(mapper) << std::endl;
    }
}

// ----------------------------------------------------------------------------
// Function _mapReads()                                         [Mapper<KBest>]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec, typename TConfig, typename TGenomeIndex, typename TErrors>
void _mapReads(Mapper<TReads, TDelegate, TSpec, TConfig> & mapper, TGenomeIndex & genomeIndex, TErrors errors,
               KBest const & /*tag*/)
{
    typedef Seeding_<>                         TSeeding_;
    typedef TSeeding_::TSeed                   TSeed;

    TReadSeqSize readsLength = avgSeqLength(getReads(mapper));
    TReadSeqSize seedLength  = mapper._seedLength;
    TReadSeqSize seedCount   = readsLength / seedLength;
    TErrors stratumDelta     = seedCount - 1;

    for (TErrors errors_ = 0; errors_ <= errors; ++errors_)
    {
        std::cout << "Errors:\t\t\t\t" << errors_ << std::endl;

        TErrors seedErrors_ = errors_ / seedCount;
        TReadSeqSize seed = errors_ % seedCount;

        std::cout << "Seed:\t\t\t\t(" << seedLength << "," << seedErrors_ << ")" << std::endl;

        mapper._extender.minErrorsPerRead = seed;
        mapper._extender.maxErrorsPerRead = std::min(errors, errors_ + stratumDelta);

        _findSeeds(mapper, genomeIndex, TSeed(seedLength, seedErrors_), seed, seed + 1);

        std::cout << "Hits:\t\t\t\t" << hitsCount(mapper) << std::endl;
        std::cout << "Matches:\t\t\t" << matchesCount(mapper) << std::endl;

        raiseErrorThreshold(mapper._manager);
    }
}

// ----------------------------------------------------------------------------
// Function _mapReads()                                       [Mapper<AnyBest>]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec, typename TConfig, typename TGenomeIndex, typename TErrors>
void _mapReads(Mapper<TReads, TDelegate, TSpec, TConfig> & mapper, TGenomeIndex & genomeIndex, TErrors errors,
               AnyBest const & /*tag*/)
{
    typedef Seeding_<>                          TSeeding_;
    typedef TSeeding_::TSeed                    TSeed;

    TReadSeqSize readsLength = avgSeqLength(getReads(mapper));
    TReadSeqSize seedLength  = mapper._seedLength;
    TReadSeqSize seedCount   = readsLength / seedLength;
    TErrors stratumDelta     = seedCount - 1;

    for (TErrors errors_ = 0; errors_ <= errors; ++errors_)
    {
        std::cout << "Errors:\t\t\t\t" << errors_ << std::endl;

        TErrors seedErrors_ = errors_ / seedCount;
        TReadSeqSize seed = errors_ % seedCount;

        std::cout << "Seed:\t\t\t\t(" << seedLength << "," << seedErrors_ << ")" << std::endl;

        mapper._extender.minErrorsPerRead = seed;
        mapper._extender.maxErrorsPerRead = std::min(errors, errors_ + stratumDelta);

//        if (mapper._extender.maxErrorsPerRead == errors)
//            mapper._extender.maxErrorsPerRead = errorsLossy;

        _findSeeds(mapper, genomeIndex, TSeed(seedLength, seedErrors_), seed, seed + 1);

        std::cout << "Hits:\t\t\t\t" << hitsCount(mapper) << std::endl;
        std::cout << "Matches:\t\t\t" << matchesCount(mapper) << std::endl;

        raiseErrorThreshold(mapper._manager);
    }

//    for (TErrors errors_ = errors + 1; errors_ <= errorsLossy; ++errors_)
//    {
//        std::cout << "Errors:\t\t\t\t" << errors_ << std::endl;
//        std::cout << "Matches:\t\t\t" << matchesCount(mapper) << std::endl;
//
//        raiseErrorThreshold(mapper._manager);
//    }
}

// ----------------------------------------------------------------------------
// Function _findSeeds()                                               [Mapper]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec, typename TConfig, typename TGenomeIndex, typename TSeed, typename TPosition>
void _findSeeds(Mapper<TReads, TDelegate, TSpec, TConfig> & mapper, TGenomeIndex & genomeIndex, TSeed seed, TPosition firstSeed, TPosition lastSeed)
{
    if (getValueI2(seed) > 0u)
        find(mapper._seederApx, genomeIndex.index, getValueI1(seed), getValueI2(seed), firstSeed, lastSeed);
    else
        find(mapper._seederExt, genomeIndex.index, getValueI1(seed), getValueI2(seed), firstSeed, lastSeed);
}

// ----------------------------------------------------------------------------
// Function hitsCount()                                                [Mapper]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec, typename TConfig>
unsigned long hitsCount(Mapper<TReads, TDelegate, TSpec, TConfig> const & mapper)
{
    return mapper._seederExt.hitsCount + mapper._seederApx.hitsCount;
}

// ----------------------------------------------------------------------------
// Function matchesCount()                                             [Mapper]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec, typename TConfig>
unsigned long matchesCount(Mapper<TReads, TDelegate, TSpec, TConfig> const & mapper)
{
    return mapper._manager.matchesCount;
}

#endif  // #ifndef SEQAN_EXTRAS_MASAI_MAPPER_H_
