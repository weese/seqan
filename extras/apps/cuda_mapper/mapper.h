// ==========================================================================
//                                cuda_mapper
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

// ============================================================================
// Prerequisites
// ============================================================================

#include <seqan/basic_extras.h>
#include <seqan/sequence_extras.h>
#include <seqan/store.h>

#include "../masai/tags.h"
#include "../masai/options.h"
#include "../masai/store/reads.h"
#include "../masai/store/genome.h"
#include "../masai/index/genome_index.h"

#include "index.h"
#include "kernels.h"

using namespace seqan;

// ============================================================================
// Tags
// ============================================================================

struct CPU_;
typedef Tag<CPU_>     CPU;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class Options
// ----------------------------------------------------------------------------

struct Options
{
    CharString  genomeFile;
    CharString  genomeIndexFile;
    CharString  readsFile;

    int         mappingBlock;
    unsigned    seedLength;

    Options() :
        mappingBlock(MaxValue<int>::VALUE),
        seedLength(33)
    {}
};

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function mapReads()                                                  [CPU]
// --------------------------------------------------------------------------

template <typename TIndex, typename TReadSeqs>
inline void
mapReads(TIndex & index, TReadSeqs & readSeqs, CPU const & /* tag */)
{
    for (unsigned i = 0; i < length(readSeqs); i++)
    {
        typename Iterator<TIndex, TopDown<> >::Type it(index);

        unsigned occurrences = goDown(it, readSeqs[i]) ? countOccurrences(it) : 0;

        printf("occurrences=%d\n", occurrences);
    }
}

// ----------------------------------------------------------------------------
// Function runMapper()
// ----------------------------------------------------------------------------

template <typename TCPUGPU>
int runMapper(Options & options)
{
    typedef Genome<void, CUDAStoreConfig>                           TGenome;
    typedef GenomeLoader<void, CUDAStoreConfig>                     TGenomeLoader;
    typedef GenomeIndex<TGenome, TGenomeIndexSpec, void>            TGenomeIndex;

    typedef FragmentStore<void, CUDAStoreConfig>                    TStore;
    typedef ReadsConfig<False, False, True, True, CUDAStoreConfig>  TReadsConfig;
    typedef Reads<void, TReadsConfig>                               TReads;
    typedef ReadsLoader<void, TReadsConfig>                         TReadsLoader;

    TGenome             genome;
    TGenomeLoader       genomeLoader(genome);
    TGenomeIndex        genomeIndex(genome);

    TStore              store;
    TReads              reads(store);
    TReadsLoader        readsLoader(reads);

    double start, finish;

    // Load genome.
    if (!open(genomeLoader, options.genomeFile))
    {
        std::cerr << "Error while loading genome" << std::endl;
        return 1;
    }

    std::cout << "Loading genome:\t\t\t" << std::flush;
    start = sysTime();
    if (!load(genomeLoader))
    {
        std::cerr << "Error while loading genome" << std::endl;
        return 1;
    }
    finish = sysTime();
    std::cout << finish - start << " sec" << std::endl;

    // Load genome index.
    std::cout << "Loading genome index:\t\t" << std::flush;
    start = sysTime();
    if (!load(genomeIndex, options.genomeIndexFile))
    {
        std::cout << "Error while loading genome index" << std::endl;
//        return 1;
        std::cout << "Building genome index:\t\t" << std::flush;
        build(genomeIndex);
    }
    finish = sysTime();
    std::cout << finish - start << " sec" << std::endl;

    // Open reads file.
    start = sysTime();
    if (!open(readsLoader, options.readsFile))
    {
        std::cerr << "Error while opening reads file" << std::endl;
        return 1;
    }

    // Reserve space for reads.
    if (options.mappingBlock < MaxValue<int>::VALUE)
        reserve(reads, options.mappingBlock);
    else
        reserve(reads);

    // Process reads in blocks.
    while (!atEnd(readsLoader))
    {
        // Load reads.
        std::cout << "Loading reads:\t\t\t" << std::flush;
        if (!load(readsLoader, options.mappingBlock))
        {
            std::cerr << "Error while loading reads" << std::endl;
            return 1;
        }
        finish = sysTime();
        std::cout << finish - start << " sec" << std::endl;
        std::cout << "Reads count:\t\t\t" << reads.readsCount << std::endl;

        // Map reads.
        start = sysTime();
        mapReads(genomeIndex.index, getSeqs(reads), TCPUGPU());
        finish = sysTime();
        std::cout << "Mapping time:\t\t\t" << std::flush;
        std::cout << finish - start << " sec" << std::endl;

        // Clear mapped reads.
        clear(reads);
    }

    // Close reads file.
    close(readsLoader);

    return 0;
}
