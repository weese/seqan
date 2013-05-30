// ==========================================================================
//                                cuda_finder
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
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
#include <seqan/index_extras.h>

using namespace seqan;

// ==========================================================================
// Metafunctions
// ==========================================================================

// --------------------------------------------------------------------------
// Metafunction Size
// --------------------------------------------------------------------------

namespace seqan {
template <typename TSpec>
struct Size<RankDictionary<TwoLevels<Dna, TSpec> > >
{
    typedef unsigned Type;
};

template <typename TSpec>
struct Size<RankDictionary<TwoLevels<bool, TSpec> > >
{
    typedef unsigned Type;
};
}

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function findCUDA()
// --------------------------------------------------------------------------

#ifdef __CUDACC__
template <typename TIndex, typename TPattern>
__global__ void
findCUDA(TIndex index, TPattern pattern)
{
    typedef typename Iterator<TIndex, TopDown<> >::Type TIterator;
    typedef typename EdgeLabel<TIterator>::Type         TEdgeLabel;
    typedef typename Fibre<TIndex, FibreText>::Type     TTextFibre;
    typedef typename Infix<TTextFibre const>::Type      TRepresentative;

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    printf("index=%i\n", idx);

    // Print the Compressed SA values.
    printf("lengthSA=%ld\n", length(indexSA(index)));
    for (unsigned i = 0; i < length(indexSA(index)); ++i)
        printf("%ld\n", indexSA(index)[i]);

    // Instantiate a virtual suffix tree iterator.
    TIterator it(index);

    // At root.
    printf("isRoot()=%d\n", isRoot(it));
    printf("repLength()=%ld\n", repLength(it));
    printf("countOccurrences()=%ld\n", countOccurrences(it));

    // Visit the leftmost children of the root.
    if (goDown(it))
    {
        // Visit all the siblings at depth one.
        do
        {
            printf("repLength()=%ld\n", repLength(it));
            printf("parentEdgeLabel()=%c\n", static_cast<char>(parentEdgeLabel(it)));
            printf("countOccurrences()=%ld\n", countOccurrences(it));
            printf("isLeaf()=%d\n", isLeaf(it));
        }
        while (goRight(it));
    }

    // Restart from root.
    goRoot(it);
    printf("goRoot()\n");

    // Search the pattern.
    printf("goDown(pattern)=%d\n", goDown(it, pattern));
    printf("repLength()=%ld\n", repLength(it));
    printf("countOccurrences()=%ld\n", countOccurrences(it));
    printf("isLeaf()=%d\n", isLeaf(it));
}
#endif

// --------------------------------------------------------------------------
// Function testInfix()
// --------------------------------------------------------------------------

template <typename TAlphabet>
void testInfix()
{
    typedef String<TAlphabet>                           TString;
    typedef typename View<TString>::Type                TStringView;
    typedef typename Infix<TString>::Type               TStringInfix;
    typedef typename Infix<TStringView>::Type           TStringViewInfix;

    TString s = "AAACCCGGGTTT";
    TStringInfix sInfix = infix(s, 3, 6);

    TStringView sView = view(s);
    TStringViewInfix sViewInfix = infix(sView, 3, 6);

    SEQAN_ASSERT(isEqual(sInfix, sViewInfix));
}

// --------------------------------------------------------------------------
// Function testStringSet()
// --------------------------------------------------------------------------

template <typename TAlphabet>
void testStringSet()
{
    typedef String<TAlphabet>                           TString;
    typedef StringSet<TString, Owner<ConcatDirect<> > > TStringSet;
    typedef typename View<TString>::Type                TStringView;
    typedef typename View<TStringSet>::Type             TStringSetView;
    typedef typename Device<TString>::Type              TDeviceString;
    typedef typename Device<TStringSet>::Type           TDeviceStringSet;

    TStringSet ss;
    appendValue(ss, "AAAAAAAA");
    appendValue(ss, "CCCCCCC");
    appendValue(ss, "GGGGGGGGGGGGGG");
    appendValue(ss, "T");

    TStringSetView ssView = view(ss);

    SEQAN_ASSERT_EQ(length(ss), length(ssView));
    for (unsigned i = 0; i < length(ss); ++i)
        SEQAN_ASSERT(isEqual(ss[i], ssView[i]));
}

// --------------------------------------------------------------------------
// Function testIndex()
// --------------------------------------------------------------------------

template <typename TAlphabet, typename TIndexSpec>
void testIndex()
{
    typedef String<TAlphabet>                           TString;
    typedef StringSet<TString, Owner<ConcatDirect<> > > TStringSet;
    typedef Index<TString, TIndexSpec>                  TIndex;
    typedef typename Device<TString>::Type              TDeviceString;
    typedef typename Device<TStringSet>::Type           TDeviceStringSet;
    typedef typename Device<TIndex>::Type               TDeviceIndex;

    // Instantiate an index over a text.
    TString text("ACGTACGTACGT");
    TIndex index(text);

    // Create the index on the reversed text.
    reverse(text);
    indexCreate(index);
    reverse(text);

    // Copy index to device.
    TDeviceIndex deviceIndex;
    assign(deviceIndex, index);

    // Create a pattern.
    TString pattern("TA");
    TDeviceString devicePattern;
    assign(devicePattern, pattern);

    // Find on GPU.
    findCUDA<<< 1,1 >>>(view(deviceIndex), view(devicePattern));
    cudaDeviceSynchronize();
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    testInfix<Dna>();
    testStringSet<Dna>();
    testIndex<Dna, FMIndex<> >();
};
