// ==========================================================================
//                                 uselessApp
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


#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>

using namespace seqan;

// --------------------------------------------------------------------------
// Class AppOptions
// --------------------------------------------------------------------------

struct AppOptions
{
    seqan::String<char> infile;
    seqan::String<char> outfile;
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

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

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
// Function callBenchmark()
// --------------------------------------------------------------------------

template <typename TString, typename TSpec>
void callBenchmarks(StringSet<TString, TSpec> const & set) {

    CharString dislexExt = "Extern";
    CharString shapeStr  = "10011001100";
    CharString string    = "String";
    CharString stringset = "StrSet";

    //- Gapped Indices: 10011001100 ----------------------------------------------------------------

    typedef CyclicShape<FixedShape<0,GappedShape<HardwiredShape<2,1,3,1> >, 2> > TShape;

    {
        typedef typename Concatenator<StringSet<TString, TSpec> >::Type			TConcat;
        typedef Multi<GappedTupler<TShape, false, BitPacked<> >,
                Pair<typename Size<TString>::Type>,
                typename StringSetLimits<StringSet<TString, TSpec> >::Type>     TMulti;

        // specialization
		typedef Pipe< TConcat, Source<> >                                       TSource;
        typedef Pipe< TSource, TMulti>                                          TTupler;

        TSource source(concat(set));
        TTupler tupler(source, stringSetLimits(set));
        beginRead(tupler);
        while (!eof(tupler))
            pop(tupler);
        endRead(tupler);
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


    time_t t; time(&t);
    std::cout << "# input file : " << options.infile << std::endl;
    std::cout << "# timestamp  : " << ctime(&t);

    callBenchmarks(seqs);

    return 0;
}

