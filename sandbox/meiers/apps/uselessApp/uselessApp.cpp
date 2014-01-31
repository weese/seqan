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
#include <seqan/seq_io.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{

    StringSet<CharString> ids;
    StringSet<Dna5String> seqs;

    /*
    // use the sequenceStream
    SequenceStream seqStream(argv[1]);
    if (!isGood(seqStream))
    {
        std::cerr << "ERROR: Could not open the file." << std::endl;
        return 1;
    }
    if (readAll(ids, seqs, seqStream) != 0)
    {
        std::cerr << "ERROR: Could not read from file" << std::endl;
        return 1;
    }
     */


    // variable typedefs:
    typedef CyclicShape<FixedShape<0, GappedShape<HardwiredShape<2> >, 0> > TShape;
    typedef GappedTupler<TShape, false> GappedTupler_;

    CharString shapeString;
    cyclicShapeToString(shapeString, TShape());
    std::cout << "Shape: " << shapeString << std::endl << std::endl;

    {
        // StringSet
        std::cout << "StringSet" << std::endl;
        StringSet<CharString> set;
        appendValue(set, "0123456789");
        appendValue(set, "012");
        appendValue(set, "0");
        appendValue(set, "0123456789");
        appendValue(set, "012");



        typedef Concatenator<StringSet<CharString> >::Type         TConcat;
        typedef StringSetLimits<StringSet<CharString> >::Type      TLimits;
;
    }

    {
        // String
        std::cout << "String" << std::endl;
        CharString s = "banaaanaaa";

        typedef Pipe< CharString, Source<> >                    TPipeSource;
        typedef Pipe< TPipeSource, DislexExternal<TShape> >     TPipeDislex;
        TPipeSource source(s);
        TPipeDislex dislex(source);

        String<unsigned> sa;
        sa << dislex;
        std::cout << sa;
    }

    
}

