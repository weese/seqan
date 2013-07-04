// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2013 NVIDIA Corporation
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
//     * Neither the name of NVIDIA Corporation nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL NVIDIA CORPORATION BE LIABLE
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

#include "cuda_demo.h"

using namespace seqan;

int main()
{
    // ----------------------------------------------------------------------
    // Create input data on the CPU.
    // ----------------------------------------------------------------------

    // Select the input types.
    typedef DnaString                                       THaystack;
    typedef StringSet<DnaString, Owner<ConcatDirect<> > >   TNeedles;

    // Create a haystack.
    THaystack haystack = "ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA";

    // Create a set of needles.
    TNeedles needles;
    appendValue(needles, "GTTG");
    appendValue(needles, "TAC");
    appendValue(needles, "CAAC");

    // ----------------------------------------------------------------------
    // Build the FM-index on the CPU.
    // ----------------------------------------------------------------------

    // Select the index type.
    typedef Index<THaystack, FMIndex<> >    TIndex;

    // Build the index over the reversed haystack.
    TIndex index(haystack);
    reverse(haystack);
    indexCreate(index);
    reverse(haystack);

    // ----------------------------------------------------------------------
    // Count on the CPU.
    // ----------------------------------------------------------------------

    omp_set_num_threads(8);
    std::cout << "CPU Occurrences: " << count(index, needles, ExecSpace<TIndex>::Type()) << std::endl;

#ifndef CUDA_DISABLED
    // ----------------------------------------------------------------------
    // Copy data to the GPU.
    // ----------------------------------------------------------------------

    // Select the GPU types.
    typedef Device<TNeedles>::Type     TDeviceNeedles;
    typedef Device<TIndex>::Type       TDeviceIndex;

    // Copy the needles to the GPU.
    TDeviceNeedles deviceNeedles;
    assign(deviceNeedles, needles);

    // Copy the index to the GPU.
    TDeviceIndex deviceIndex;
    assign(deviceIndex, index);

    // ----------------------------------------------------------------------
    // Count on the GPU.
    // ----------------------------------------------------------------------

    std::cout << "GPU Occurrences: " << count(deviceIndex, deviceNeedles, ExecSpace<TDeviceIndex>::Type()) << std::endl;
#endif

    return 0;
}
