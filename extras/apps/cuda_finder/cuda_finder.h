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

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/misc/cuda.h>
#include <seqan/index_extras.h>
#include <seqan/misc/misc_view.h>
#include <seqan/index/index_view.h>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>

using namespace seqan;

// ==========================================================================
// Classes
// ==========================================================================

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function findOnGPU()
// --------------------------------------------------------------------------

#ifdef __CUDACC__
template <typename TText, typename TViewSpec, typename TSpec>
__global__ void
findOnGPU(Index<View<TText, TViewSpec>, TSpec> index)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    cuPrintf("index=%i", idx);

    cuPrintf("lengthSA=%i", length(indexSA(index)));
}
#endif

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    typedef CharString                                  TText;
    typedef thrust::host_vector<char>                   THostText;
    typedef thrust::device_vector<char>                 TDeviceText;
    typedef Index<TDeviceText, IndexSa<> >              TDeviceIndex;

    TText text("text");
    THostText hostText(length(text));
    thrust::copy(begin(text, Standard()), end(text, Standard()), hostText.begin());
    TDeviceText deviceText(hostText);
    TDeviceIndex deviceIndex(deviceText);

    int block_size = 1;
    int n_blocks = 1;

    cudaPrintfInit();
    findOnGPU<<< n_blocks, block_size >>>(toView(deviceIndex));
    cudaDeviceSynchronize();
    cudaPrintfDisplay(stdout, true);
    cudaPrintfEnd();

    return 0;
}
