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

#include <seqan/basic_extras.h>
#include <seqan/sequence_extras.h>
#include <seqan/store.h>

//#include "../masai/options.h"
#include "../masai/tags.h"
#include "../masai/store/reads.h"
#include "../masai/store/genome.h"
#include "../masai/index/genome_index.h"

#include "index.h"

using namespace seqan;

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
        genomeFile("genome.fasta"),
        genomeIndexFile("genome"),
        readsFile("reads.fastq"),
        mappingBlock(MaxValue<int>::VALUE),
        seedLength(33)
    {}
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function setupArgumentParser()                              [ArgumentParser]
// ----------------------------------------------------------------------------

//void setupArgumentParser(ArgumentParser & parser, Options const & options)
//{
//    setAppName(parser, "cuda_mapper");
//    setShortDescription(parser, "CUDA Mapper");
//    setCategory(parser, "Read Mapping");
//
//    setDateAndVersion(parser);
//    setDescription(parser);
//
//    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIGENOME FILE\\fP> <\\fIREADS FILE\\fP>");
//
//    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE));
//    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE));
//    setValidValues(parser, 0, "fasta fa");
//    setValidValues(parser, 1, "fastq fasta fa");
//
//    addSection(parser, "Mapping Options");
//
//    addOption(parser, ArgParseOption("mb", "mapping-block", "Maximum number of reads to be mapped at once.", ArgParseOption::INTEGER));
//    setMinValue(parser, "mapping-block", "10000");
//    setDefaultValue(parser, "mapping-block", options.mappingBlock);
//
//    addOption(parser, ArgParseOption("sl", "seed-length", "Minimum seed length.", ArgParseOption::INTEGER));
//    setMinValue(parser, "seed-length", "10");
//    setMaxValue(parser, "seed-length", "100");
//    setDefaultValue(parser, "seed-length", options.seedLength);
//
//
//    addSection(parser, "Genome Index Options");
//
//    setIndexPrefix(parser);
//}

// ----------------------------------------------------------------------------
// Function parseCommandLine()                                        [Options]
// ----------------------------------------------------------------------------

//ArgumentParser::ParseResult
//parseCommandLine(Options & options, ArgumentParser & parser, int argc, char const ** argv)
//{
//    ArgumentParser::ParseResult res = parse(parser, argc, argv);
//
//    if (res != seqan::ArgumentParser::PARSE_OK)
//        return res;
//
//    // Parse genome input file.
//    getArgumentValue(options.genomeFile, parser, 0);
//
//    // Parse reads input file.
//    getArgumentValue(options.readsFile, parser, 1);
//
//    // Parse mapping block.
//    getOptionValue(options.mappingBlock, parser, "mapping-block");
//
//    // Parse mapping options.
//    getOptionValue(options.seedLength, parser, "seed-length");
//
//    // Parse genome index prefix.
//    getIndexPrefix(options, parser);
//
//    return seqan::ArgumentParser::PARSE_OK;
//}

// --------------------------------------------------------------------------
// Function mapReadsKernel()
// --------------------------------------------------------------------------

#ifdef __CUDACC__
template <typename TIndex, typename TReads>
__global__ void
mapReadsKernel(TIndex index, TReads reads)
{
    typedef typename Iterator<TIndex, TopDown<> >::Type TIterator;

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    printf("index=%i\n", idx);

//    TIterator it(index);
//    goDown(it, reads[idx]);
}
#endif

// ----------------------------------------------------------------------------
// Function runMapper()
// ----------------------------------------------------------------------------

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
//    build(genomeIndex);
    if (!load(genomeIndex, options.genomeIndexFile))
    {
        std::cout << "Error while loading genome index" << std::endl;
        return 1;
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
        mapReadsKernel<<<32,512>>>(view(genomeIndex.index), view(getSeqs(reads)));
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

// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
//    ArgumentParser parser;
    Options options;
//    setupArgumentParser(parser, options);
//
//    ArgumentParser::ParseResult res = parseCommandLine(options, parser, argc, argv);
//
//    if (res != seqan::ArgumentParser::PARSE_OK)
//        return res == seqan::ArgumentParser::PARSE_ERROR;

    return runMapper(options);
}
