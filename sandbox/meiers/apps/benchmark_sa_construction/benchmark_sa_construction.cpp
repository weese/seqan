// ==========================================================================
//                         benchmark_sa_construction
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
// Author: Sascha Meiers <meiers@inf.fu-berlin.de>
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/arg_parse.h>
#include <seqan/index.h>

#include <time.h>

using namespace seqan;


// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class AppOptions
// --------------------------------------------------------------------------

struct AppOptions
{
    seqan::String<char> infile;
    seqan::String<char> outfile;
    bool qsort;
};


// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(AppOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("seqanBuildSuffixArray");
    // Set short description, version, and date.
    setShortDescription(parser, "Builds an ESA-index.");
    setVersion(parser, "Sascha.0.1");
    setDate(parser, "July 2013");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fITEXT\\fP\"");
    addDescription(parser, "Builds an index (ESA) of the specified genome and writes it to disk.");

    // We require one argument.
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUTFILE, "IN"));

    addOption(parser, seqan::ArgParseOption("o",
                                            "output",
                                            "Output file (<IN>.index if not specified)",
                                            seqan::ArgParseArgument::OUTPUTFILE, "OUT"));

    addOption(parser, seqan::ArgParseOption("n", "noqsort", "Disable QuickSort (for large input files)"));

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
    return res;

    // disable qicksort
    options.qsort = ! isSet(parser, "noqsort");

    // set outfile to infile.db if not outfile was specified
    seqan::String<char> outfile;
    seqan::getArgumentValue(options.infile, parser, 0);

    outfile = options.	infile;
    outfile += ".index";
    getOptionValue(outfile, parser, "output");
    options.outfile = outfile;

    return seqan::ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function benchmark()
// --------------------------------------------------------------------------

template <typename TIndex, typename TAlgo, typename TSuffAr, typename TLabel>
void benchmarkSACA(
                   TIndex & index,
                   TAlgo const & algo,
                   TSuffAr const & correctSA,
                   TLabel const & labelPattern,
                   TLabel const & labelAlgorithm,
                   TLabel const & labelText)
{
    double start = sysTime();
    indexCreate(index, EsaSA(), algo);

    std::cout << labelAlgorithm << "\t" << labelPattern << "\t" << labelText << "\t" << sysTime() - start << "\t";

    if(length(correctSA) == length(indexSA(index)))
    {
        typename Size<TIndex>::Type errs =0;
        typename Iterator<typename Fibre<TIndex, EsaSA>::Type>::Type it2 = begin(indexSA(index));
        for(typename Iterator<TSuffAr const>::Type it=begin(correctSA); it < end(correctSA); ++it, ++it2)
        if(*it != *it2) std::cout << std::endl << "   |-> err " << ++errs << " at pos " << it - begin(correctSA) << ": correct " << *it << " , seen " << *it2;
        //if(*it != *it2) ++errs;
        std::cout << "errors:" << errs << std::endl;
    } else {
        std::cout << "/" << std::endl;
    }
}

/*
template <typename TPair>
struct _V1comparator
{
    bool operator()(TPair const & a, TPair const & b)
    {
        return a.i1 < b.i1;
    }
};

template <typename TString, typename TMod>
String<Pair<Dna5String, unsigned> >
___buildListOfSuffixes(TString & str, TMod const &)
{
    String<unsigned> sa;
    resize(sa, length(str));
    _initializeSA(sa, str);

    typedef ModifiedString<typename Suffix<TString>::Type, TMod> TModStr;
    typedef Pair<Dna5String, unsigned> TPair;

    String<TPair> list;

    for (typename Iterator<String<unsigned> >::Type saIt=begin(sa); saIt != end(sa); ++saIt)
    {
        Dna5String s = TModStr(suffix(str, *saIt));
        TPair pa(s, *saIt);
        appendValue(list, pa);
        std::cout << *saIt << "\t" << s << std::endl;
    }

    return list;
}

template <typename TString, typename TMod, typename TSpec>
String<Pair<Dna5String, Pair<unsigned> > >
___buildListOfSuffixes(StringSet<TString, TSpec> & str, TMod const &)
{
    String<Pair<unsigned> > sa;
    resize(sa, lengthSum(str));
    _initializeSA(sa, str);

    typedef ModifiedString<typename Suffix<TString>::Type, TMod> TModStr;
    typedef Pair<Dna5String, Pair<unsigned> > TPair;

    String<TPair> list;

    for (typename Iterator<String<Pair<unsigned> > >::Type saIt=begin(sa); saIt != end(sa); ++saIt)
    {
        Dna5String s = TModStr(suffix(str, *saIt));
        TPair pa(s, *saIt);
        appendValue(list, pa);
        std::cout << *saIt << "\t" << s << std::endl;
    }

    return list;
}


template <typename TList, typename TSA>
void ___generateCorrectSA(TSA & sa, TList & list)
{
    _V1comparator<typename Value<TList>::Type > comp;
    std::sort(begin(list), end(list), comp);

    typedef typename Value<TSA>::Type SAValue;
    resize(sa, length(list));

    typename Iterator<TList>::Type listIt = begin(list);
    for (typename Iterator<TSA>::Type saIt=begin(sa); saIt != end(sa); ++saIt, ++listIt)
        *saIt = listIt->i2;
}
*/


// --------------------------------------------------------------------------
// Function callBenchmark()
// --------------------------------------------------------------------------

template <typename TString, typename TSpec>
void callBenchmarks(StringSet<TString, TSpec> const & set, bool qsort) {

    // header
    std::cout << "# alphabet size: " << (int)ValueSize<typename Value<TString>::Type>::VALUE << std::endl;
    std::cout << "# dataset: " << length(set) << " strings, total length " << length(concat(set)) << std::endl;
    std::cout <<  "algor.\tpattern_______\ttext\truntime\tcheck" << std::endl;

    // labels for output
    CharString ungapped  = "ungapped______";
    CharString shape101  = "101___________";
    CharString shape2    = "11000100101110";
    CharString shape3    = "0001010_______";
    CharString skew7     = "Skew7";
    CharString skew3     = "Skew3";
    CharString saqsort   = "QSort";
    CharString dislex    = "Dislex";
    CharString dislexExt = "Extern";
    CharString radix     = "Radix";
    CharString string    = "String";
    CharString stringset = "StrSet";


    {   //- Ungapped Indices ----------------------------------------------------------------

        String<Pair<typename Size<TString>::Type> > correctSA1;
        String<typename Size<TString>::Type> correctSA2;

        {
            typedef Index<StringSet<TString, TSpec> const, IndexSa<> > TIndex;
            TIndex index(set);
            benchmarkSACA(index, Skew7(), correctSA1, ungapped, skew7, stringset);

            correctSA1 = indexSA(index);
        }
        {
            typedef Index<TString, IndexSa<> > TIndex;
            TIndex index (concat(set));
            benchmarkSACA(index, Skew7(), correctSA2, ungapped, skew7, string);

            correctSA2 = indexSA(index);
        }
        {
            typedef Index<StringSet<TString, TSpec> const, IndexSa<> > TIndex;
            TIndex index(set);
            benchmarkSACA(index, Skew3(), correctSA1, ungapped, skew3, stringset);
        }
        {
            typedef Index<TString, IndexSa<> > TIndex;
            TIndex index (concat(set));
            benchmarkSACA(index, Skew3(), correctSA2, ungapped, skew3, string);
        }
        if (qsort)
        {
            typedef Index<StringSet<TString, TSpec> const, IndexSa<> > TIndex;
            TIndex index(set);
            benchmarkSACA(index, SAQSort(), correctSA1, ungapped, saqsort, stringset);
        }
        if (qsort)
        {
            typedef Index<TString, IndexSa<> > TIndex;
            TIndex index (concat(set));
            benchmarkSACA(index, SAQSort(), correctSA2, ungapped, saqsort, string);
        }
        {
            typedef Index<StringSet<TString, TSpec> const, IndexSa<> > TIndex;
            TIndex index(set);
            benchmarkSACA(index, InplaceRadixSort(), correctSA1, ungapped, radix, stringset);
        }
        {
            typedef Index<TString, IndexSa<> > TIndex;
            TIndex index(concat(set));
            benchmarkSACA(index, InplaceRadixSort(), correctSA2, ungapped, radix, string);
        }
    }





    {   //- Gapped Indices: 101 ----------------------------------------------------------------

        typedef CyclicShape<FixedShape<0,GappedShape<HardwiredShape<2> >, 0> > TShape;
        String<Pair<typename Size<TString>::Type> > correctSA1;
        String<typename Size<TString>::Type> correctSA2;

        {
            typedef Index<StringSet<TString, TSpec>, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set);
            benchmarkSACA(index, InplaceRadixSort(), correctSA1, shape101, radix, stringset);

            correctSA1 = indexSA(index);
        }
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(concat(set));
            benchmarkSACA(index, InplaceRadixSort(), correctSA2, shape101, radix, string);

            correctSA2 = indexSA(index);
        }
        if (qsort)
        {
            typedef Index<StringSet<TString, TSpec> const, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set);
            benchmarkSACA(index, SAQSort(), correctSA1, shape101, saqsort, stringset);
        }
        if(qsort)
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(concat(set));
            benchmarkSACA(index, SAQSort(), correctSA2, shape101, saqsort, string);
        }
        {
            typedef Index<StringSet<TString,TSpec>, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set);
            benchmarkSACA(index, Dislex<Skew7>(), correctSA1, shape101, dislex, stringset);
        }
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(concat(set));
            benchmarkSACA(index, Dislex<Skew7>(), correctSA2, shape101, dislex, string);
        }
        {
            typedef Index<StringSet<TString,TSpec>, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set);
            benchmarkSACA(index, DislexExternal<TShape>(), correctSA1, shape101, dislexExt, stringset);
        }
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(concat(set));
            benchmarkSACA(index, DislexExternal<TShape>(), correctSA2, shape101, dislexExt, string);
        }
    }



    { //----- Gapped Indices: 11000100101110 ---------------------------------------------------------

        typedef CyclicShape<FixedShape<0,GappedShape<HardwiredShape<1,4,3,2,1,1> >, 1> > TShape;
        String<Pair<typename Size<TString>::Type> > correctSA1;
        String<typename Size<TString>::Type> correctSA2;

        {
            typedef Index<StringSet<TString, TSpec>, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set);
            benchmarkSACA(index, InplaceRadixSort(), correctSA1, shape2, radix, stringset);

            correctSA1 = indexSA(index);
        }
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(concat(set));
            benchmarkSACA(index, InplaceRadixSort(), correctSA2, shape2, radix, string);

            correctSA2 = indexSA(index);
        }
        if (qsort)
        {
            typedef Index<StringSet<TString, TSpec> const, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set);
            benchmarkSACA(index, SAQSort(), correctSA1, shape2, saqsort, stringset);
        }
        if (qsort)
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(concat(set));
            benchmarkSACA(index, SAQSort(), correctSA2, shape2, saqsort, string);
        }
        {
            typedef Index<StringSet<TString,TSpec>, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set);
            benchmarkSACA(index, Dislex<Skew7>(), correctSA1, shape2, dislex, stringset);
        }
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(concat(set));
            benchmarkSACA(index, Dislex<Skew7>(), correctSA2, shape2, dislex, string);
        }
        {
            typedef Index<StringSet<TString,TSpec>, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set);
            benchmarkSACA(index, DislexExternal<TShape>(), correctSA1, shape2, dislexExt, stringset);
        }
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(concat(set));
            benchmarkSACA(index, DislexExternal<TShape>(), correctSA2, shape2, dislexExt, string);
        }
    }

    { //----- Gapped Indices: 0001010 ---------------------------------------------------------

        typedef CyclicShape<FixedShape<3,GappedShape<HardwiredShape<2> >, 1> > TShape;
        String<Pair<typename Size<TString>::Type> > correctSA1;
        String<typename Size<TString>::Type> correctSA2;

        {
            typedef Index<StringSet<TString, TSpec>, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set);
            benchmarkSACA(index, InplaceRadixSort(), correctSA1, shape3, radix, stringset);

            correctSA1 = indexSA(index);
        }
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(concat(set));
            benchmarkSACA(index, InplaceRadixSort(), correctSA2, shape3, radix, string);

            correctSA2 = indexSA(index);
        }
        if (qsort)
        {
            typedef Index<StringSet<TString, TSpec> const, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set);
            benchmarkSACA(index, SAQSort(), correctSA1, shape3, saqsort, stringset);
        }
        if (qsort)
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(concat(set));
            benchmarkSACA(index, SAQSort(), correctSA2, shape3, saqsort, string);
        }
        {
            typedef Index<StringSet<TString,TSpec>, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set);
            benchmarkSACA(index, Dislex<Skew7>(), correctSA1, shape3, dislex, stringset);
        }
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(concat(set));
            benchmarkSACA(index, Dislex<Skew7>(), correctSA2, shape3, dislex, string);
        }
        {
            typedef Index<StringSet<TString,TSpec>, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set);
            benchmarkSACA(index, DislexExternal<TShape>(), correctSA1, shape3, dislexExt, stringset);
        }
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(concat(set));
            benchmarkSACA(index, DislexExternal<TShape>(), correctSA2, shape3, dislexExt, string);
        }
    }



}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{
    // Parse the command line.
    seqan::ArgumentParser parser;
    AppOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != seqan::ArgumentParser::PARSE_OK)
    return res == seqan::ArgumentParser::PARSE_ERROR;


    // start main program ///////////////////////////////////////////////////

    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::Dna5String> seqs;



    // New reading:
    // use the sequenceStream
    seqan::SequenceStream seqStream(seqan::toCString(options.infile));
    if (!isGood(seqStream))
    {
        std::cerr << "ERROR: Could not open the file." << std::endl;
        return 1;
    }
    if (readAll(ids, seqs, seqStream) != 0)
    {
        std::cerr << "ERROR: Could not read from " << options.infile << std::endl;
        return 1;
    }
    
    
    // Old reading:
    // This fails if the fasta file is too small !
    //    std::fstream in;
    //    in.open( seqan::toCString(options.infile));
    //    seqan::RecordReader<std::fstream, seqan::SinglePass<> > reader(in);
    //    if (!in.good())
    //    {
    //        std::cerr << "Couldn't open " << options.infile << std::endl;
    //        return 1;
    //    }
    //    if (read2(ids, seqs, reader, seqan::Fasta()) != 0)
    //        return 2;  // Could not record from file.
    
    
    
    time_t t; time(&t);
    std::cout << "# input file : " << options.infile << std::endl;
    std::cout << "# timestamp  : " << ctime(&t);
    
    callBenchmarks(seqs, options.qsort);
    
    return 0;
}