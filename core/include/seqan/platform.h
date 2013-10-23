// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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

#ifndef SEQAN_PLATFORM_H
#define SEQAN_PLATFORM_H

#ifdef __MINGW32__
	#include "platform/platform_mingw.h"
#elif _MSC_VER
	#include "platform/platform_windows.h"
#elif __SUNPRO_C
	#include "platform/platform_solaris.h"
#elif __ICC
	#include "platform/platform_icc.h"
#elif __PGI
	#include "platform/platform_pgi.h"
#else
	#include "platform/platform_gcc.h"
#endif

// NOTE(esiragusa): nvcc header must be included even if __CUDACC__ is not defined.
#include "platform/platform_nvcc.h"

// This macro defines an alias for the auto_ptr class template deprecated in C++11.
#ifdef SEQAN_CXX11_STANDARD
    #define SEQAN_AUTO_PTR_NAME unique_ptr
#else  // #ifdef SEQAN_CXX11_STANDARD
    #define SEQAN_AUTO_PTR_NAME auto_ptr
#endif  // #ifdef SEQAN_CXX11_STANDARD

// branch prediction macros, see platform_gcc.h
#ifndef SEQAN_LIKELY
#define SEQAN_LIKELY(expr) (expr)
#endif

#ifndef SEQAN_UNLIKELY
#define SEQAN_UNLIKELY(expr) (expr)
#endif

#endif
