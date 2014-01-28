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

int main(int argc, char * argv[])
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
    typedef CyclicShape<FixedShape<1, GappedShape<HardwiredShape<1,1,2> >, 1> > TShape;
    typedef GappedTupler<TShape, false> GappedTupler_;

    CharString shape;
    cyclicShapeToString(shape, TShape());
    std::cout << "Shape: " << shape << std::endl << std::endl;



    // Strings
    CharString str = "0123456789012345678901234567890123456789";
    std::cout << "String: " << str << std::endl;

    typedef Pipe< CharString, Source<> >            src_t;
    typedef Pipe< src_t, GappedTupler_ >	        tupler_t;
    src_t		src(str);
    tupler_t	tupler(src);
    std::cout << tupler;





    // StringSet
    std::cout << "StringSet" << std::endl;
    StringSet<CharString> set;
    appendValue(set, "0123456789");
    appendValue(set, "012");
    appendValue(set, "0");
    appendValue(set, "0123456789");

    typedef typename Concatenator<StringSet<CharString> >::Type         TConcat;
    typedef typename StringSetLimits<StringSet<CharString> >::Type      TLimits;
    typedef Multi< GappedTupler_, Pair<typename Size<CharString>::Type>, TLimits>             TMultiGappedTupler;
    typedef Pipe< TConcat, Source<> >                                   src_t2;
    typedef Pipe< src_t2, TMultiGappedTupler >                          tupler_t2;
    src_t2		src2(concat(set));
    tupler_t2	tupler2(src2,stringSetLimits(set));
    std::cout << tupler2;




    /*
     // test for aldabi praktikum
     // write SA as binary file

     using namespace std;
     ifstream f(argv[1], ios::in);
     CharString text;
     string line;
     while(getline(f, line))
     append(text,line);

     cout << "length = " << length(text) << endl;
     toUpper(text);

     typedef Index<CharString, IndexEsa<> > TIndex;
     TIndex index(text);

     indexRequire(index, FibreSA());
     cout << "sizeof(SAValue) = " << sizeof(SAValue<TIndex>::Type) << endl;

     save(index, "/Users/mauli87/Downloads/chr1_noN.sa");
     f.close();
     
     */
    
}

