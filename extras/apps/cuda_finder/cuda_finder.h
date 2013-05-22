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
// Classes
// ==========================================================================

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function findCUDA()
// --------------------------------------------------------------------------

#ifdef __CUDACC__
template <typename TText, typename TViewSpec, typename TSpec, typename TPattern>
__global__ void
findCUDA(Index<View<TText, TViewSpec>, TSpec> index, View<TPattern, TViewSpec> pattern)
{
    typedef Index<View<TText, TViewSpec>, TSpec>        TIndex;
    typedef typename Iterator<TIndex, TopDown<> >::Type TIterator;
    typedef typename EdgeLabel<TIterator>::Type         TEdgeLabel;
    typedef typename Fibre<TIndex, FibreText>::Type     TTextFibre;
    typedef typename Infix<TTextFibre const>::Type      TRepresentative;

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    printf("index=%i\n", idx);

    printf("lengthSA=%ld\n", length(indexSA(index)));
    for (unsigned i = 0; i < length(indexSA(index)); ++i)
        printf("%ld\n", indexSA(index)[i]);

    printf("lengthLcp=%ld\n", length(indexLcp(index)));
    for (unsigned i = 0; i < length(indexLcp(index)); ++i)
        printf("%ld\n", indexLcp(index)[i]);

    printf("lengthChildtab=%ld\n", length(indexChildtab(index)));
    for (unsigned i = 0; i < length(indexChildtab(index)); ++i)
        printf("%ld\n", indexChildtab(index)[i]);

    TIterator it(index);

    printf("isRoot()=%d\n", isRoot(it));
    printf("countOccurrences()=%ld\n", countOccurrences(it));
    printf("isLeaf()=%d\n", isLeaf(it));
    printf("repLength()=%ld\n", repLength(it));
    printf("goDown()=%d\n", goDown(it));
    printf("repLength()=%ld\n", repLength(it));
    printf("goRight()=%d\n", goRight(it));

//    TEdgeLabel edgeLabel = parentEdgeLabel(it);
//    TRepresentative repr = representative(it);

//    goRoot(it);
//    printf("goDown(pattern)=%d\n", goDown(it, pattern));

//    begin(repr, Standard());

//    typedef typename Iterator<TRepresentative, Standard>::Type     TReprIt;
//    TReprIt reprIt = begin(repr, Standard());
}
#endif

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    typedef char                                        TAlphabet;
    typedef String<TAlphabet>                           TText;
    typedef thrust::device_vector<TAlphabet>            TDeviceText;
    typedef Index<TText, IndexEsa<> >                   TIndex;
    typedef Index<TDeviceText, IndexEsa<> >             TDeviceIndex;

    // Create text.
    TText text("text");

    // Create index.
    TIndex index(text);
    indexCreate(index, FibreSA());
    indexCreate(index, FibreLcp());
    indexCreate(index, FibreChildtab());

    // Copy index to device.
    TDeviceIndex deviceIndex;
    assign(deviceIndex, index);

    // Create a pattern.
    TText pattern("ex");
    TDeviceText devicePattern;
    assign(devicePattern, pattern);

    // Find on GPU.
    findCUDA<<< 1,1 >>>(toView(deviceIndex), toView(devicePattern));
    cudaDeviceSynchronize();

    return 0;
}
