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

#ifndef SEQAN_EXTRAS_CUDAMAPPER_MAPPER_H_
#define SEQAN_EXTRAS_CUDAMAPPER_MAPPER_H_

// ============================================================================
// Prerequisites
// ============================================================================

#include <seqan/basic_extras.h>
#include <seqan/sequence_extras.h>

#include "index.h"


namespace seqan {

// ============================================================================
// ============================================================================
// ============================================================================

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class Delegator
// ----------------------------------------------------------------------------

template <typename TObject, typename TDelegate, typename TSpec = void>
struct Delegator
{
    TObject & object;
    TDelegate & delegate;

    SEQAN_FUNC
    Delegator(TObject & object, TDelegate & delegate) :
        object(object),
        delegate(delegate)
    {}

    template <typename TOther>
    SEQAN_FUNC void
    operator()(TOther & /* other */)
    {
        delegate(object);
    }
};

// ============================================================================
// ============================================================================
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction View
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TSpec>
struct View<Iter<Index<TText, TIndexSpec>, VSTree<TSpec> > >
{
    typedef Iter<typename View<Index<TText, TIndexSpec> >::Type, VSTree<TSpec> >    Type;
};

// ----------------------------------------------------------------------------
// Function view()
// ----------------------------------------------------------------------------

//#ifdef __CUDACC__
//template <typename TText, typename TIndexSpec, typename TSpec>
//inline typename View<Iter<Index<TText, TIndexSpec>, VSTree<TSpec> > >::Type
//view(Iter<Index<TText, TIndexSpec>, VSTree<TSpec> > & it)
//{
//    typename View<Iter<Index<TText, TIndexSpec>, VSTree<TSpec> > >::Type itView;
//
//    // TODO(esiragusa): Take device pointer only for TText being a device_vector
////    itView.index = thrust::raw_pointer_cast(it.index);
//
//    return itView;
//}
//#endif

// ============================================================================
// ============================================================================
// ============================================================================

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class Finder
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec = void>
struct Finder2
{};

// ----------------------------------------------------------------------------
// Class Finder
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec>
struct Finder2<Index<TText, TIndexSpec>, TPattern, TSpec>
{
    typedef Index<TText, TIndexSpec>                        TIndex;
    typedef typename Iterator<TIndex, TopDown<> >::Type     TIterator;

    TIterator   _textIt;

    SEQAN_FUNC
    Finder2() {}

    SEQAN_FUNC
    Finder2(TIndex & index) :
        _textIt(index)
    {}
};

// ============================================================================
// Metafunction
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction View
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec>
struct View<Finder2<TText, TPattern, TSpec> >
{
    typedef Finder2<typename View<TText>::Type, typename View<TPattern>::Type, TSpec>   Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function view()
// ----------------------------------------------------------------------------

//template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec>
//inline typename View<Finder2<Index<TText, TIndexSpec>, TPattern, TSpec> >::Type
//view(Finder2<Index<TText, TIndexSpec>, TPattern, TSpec> & finder)
//{
//    typename View<Finder2<Index<TText, TIndexSpec>, TPattern, TSpec> >::Type finderView;
//
//    finderView._textIt = view(finder._textIt);
//
//    return finderView;
//}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec>
SEQAN_FUNC void
clear(Finder2<Index<TText, TIndexSpec>, TPattern, TSpec> & finder)
{
    goRoot(finder._textIt);
}

// ----------------------------------------------------------------------------
// Function find()
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPattern, typename TDelegate>
SEQAN_FUNC void
find(Finder2<Index<TText, TIndexSpec>, TPattern, FinderSTree> & finder,
     TPattern const & pattern,
     TDelegate & delegate)
{
    if (goDown(finder._textIt, pattern)) delegate(finder);
}

// ============================================================================
// ============================================================================
// ============================================================================

// ============================================================================
// Tags
// ============================================================================

// ----------------------------------------------------------------------------
// Tag Multiple
// ----------------------------------------------------------------------------

template <typename TSpec>
struct Multiple;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Finder_
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec>
struct Finder_;

template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec>
struct Finder_<Index<TText, TIndexSpec>, TPattern, Multiple<TSpec> >
{
    typedef Index<TText, TIndexSpec>            TIndex_;
    typedef typename Value<TPattern>::Type      TPattern_;
    typedef Finder2<TIndex_, TPattern_, TSpec>  Type;
};

#ifdef __CUDACC__
template <typename TValue, typename TAlloc, typename TIndexSpec, typename TPattern, typename TSpec>
struct Finder_<Index<thrust::device_vector<TValue, TAlloc>, TIndexSpec>, TPattern, Multiple<TSpec> >
{
    typedef thrust::device_vector<TValue, TAlloc>   TText_;
    typedef Index<TText_, TIndexSpec>               TIndex_;
    typedef typename View<TPattern>::Type           TPatternView_;
    typedef typename Value<TPatternView_>::Type     TPattern_;
    typedef Finder2<TIndex_, TPattern_, TSpec>      Type;
};

template <typename TValue, typename TAlloc, typename TSSetSpec, typename TIndexSpec, typename TPattern, typename TSpec>
struct Finder_<Index<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, TIndexSpec>, TPattern, Multiple<TSpec> >
{
    typedef StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec> TText_;
    typedef Index<TText_, TIndexSpec>                                   TIndex_;
    typedef typename View<TPattern>::Type                               TPatternView_;
    typedef typename Value<TPatternView_>::Type                         TPattern_;
    typedef Finder2<TIndex_, TPattern_, TSpec>                          Type;
};
#endif

// ----------------------------------------------------------------------------
// Metafunction FinderPool_
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec>
struct FinderPool_;

template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec>
struct FinderPool_<Index<TText, TIndexSpec>, TPattern, Multiple<TSpec> >
{
    typedef String<typename Finder_<Index<TText, TIndexSpec>, TPattern, Multiple<TSpec> >::Type>    Type;
};

template <typename TText, typename TViewSpec, typename TIndexSpec, typename TPattern, typename TSpec>
struct FinderPool_<Index<ContainerView<TText, TViewSpec>, TIndexSpec>, TPattern, Multiple<TSpec> >
{
    typedef typename View<typename FinderPool_<Index<TText, TIndexSpec>, TPattern, Multiple<TSpec> >::Type>::Type Type;
};

template <typename TText, typename TViewSpec, typename TSSetSpec, typename TIndexSpec, typename TPattern, typename TSpec>
struct FinderPool_<Index<StringSet<ContainerView<TText, TViewSpec>, TSSetSpec>, TIndexSpec>, TPattern, Multiple<TSpec> >
{
    typedef typename View<typename FinderPool_<Index<StringSet<TText, TSSetSpec>, TIndexSpec>, TPattern, Multiple<TSpec> >::Type>::Type Type;
};

#ifdef __CUDACC__
template <typename TValue, typename TAlloc, typename TIndexSpec, typename TPattern, typename TSpec>
struct FinderPool_<Index<thrust::device_vector<TValue, TAlloc>, TIndexSpec>, TPattern, Multiple<TSpec> >
{
    typedef thrust::device_vector<TValue, TAlloc>                       TText_;
    typedef Index<TText_, TIndexSpec>                                   TIndex_;
    typedef typename Finder_<TIndex_, TPattern, Multiple<TSpec> >::Type TFinder_;
    typedef typename View<TFinder_>::Type                               TFinderView_;
    typedef thrust::device_vector<TFinderView_>                         Type;
};

template <typename TValue, typename TAlloc, typename TSSetSpec, typename TIndexSpec, typename TPattern, typename TSpec>
struct FinderPool_<Index<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, TIndexSpec>, TPattern, Multiple<TSpec> >
{
    typedef StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec> TText_;
    typedef Index<TText_, TIndexSpec>                                   TIndex_;
    typedef typename Finder_<TIndex_, TPattern, Multiple<TSpec> >::Type TFinder_;
    typedef typename View<TFinder_>::Type                               TFinderView_;
    typedef thrust::device_vector<TFinderView_>                         Type;
};
#endif

// ----------------------------------------------------------------------------
// Metafunction FinderBlockSize_
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec>
struct FinderBlockSize_ {};

template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec>
struct FinderBlockSize_<Index<TText, TIndexSpec>, TPattern, Multiple<TSpec> >
{
    static const unsigned VALUE = 256;
};

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class Finder
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec>
struct Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<TSpec> >
{
    typedef Index<TText, TIndexSpec>                                        TIndex;
    typedef typename FinderPool_<TIndex, TPattern, Multiple<TSpec> >::Type  TPool;

    TIndex      * index;
    TPool       _pool;
    unsigned    _poolSize;

//    typename Member<Finder2, FinderPattern_>::Type    _patternsIt;
//    typename Member<Finder2, FinderIndices_>::Type    _idxs;
//    typename Member<Finder2, FinderHashes_>::Type     _hashes;

    SEQAN_FUNC
    Finder2() :
        index(),
        _poolSize(0)
    {}

    Finder2(TIndex & index) :
        index(&index),
        _poolSize(0)
    {}
};

// ============================================================================
// Kernels
// ============================================================================

// ----------------------------------------------------------------------------
// Kernel _computeHashesKernel()
// ----------------------------------------------------------------------------

#ifdef __CUDACC__
template <typename TFinderView, typename TPatternsView>
__global__ void
_computeHashesKernel(TFinderView finder, TPatternsView patterns)
{
    typedef typename Value<TPatternsView>::Type         TPatternView;
    typedef typename Value<TPatternView>::Type          TAlphabet;
    typedef Shape<TAlphabet, UngappedShape<10> >        TShape;

    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    // Return silently if there is no job left.
    if (idx >= length(patterns)) return;

    // Compute the hash of a read.
    TShape shape;
    finder._hashes[idx] = hash(shape, begin(patterns[idx], Standard()));

    // Fill idxs with the identity permutation.
    finder._idxs[idx] = idx;
}
#endif

// ----------------------------------------------------------------------------
// Kernel _findKernel()
// ----------------------------------------------------------------------------

#ifdef __CUDACC__
template <typename TFinderView, typename TPatternsView, typename TDelegateView>
__global__ void
_findKernel(TFinderView finder, TPatternsView patterns, TDelegateView delegate)
{
    typedef Delegator<TFinderView, TDelegateView>       TDelegator;

    unsigned threadId = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned gridThreads = gridDim.x * blockDim.x;

    // Use a delegator object to delegate this finder instead of the pooled finders.
    TDelegator delegator(finder, delegate);

    unsigned patternsCount = length(patterns);
	for (unsigned idx = threadId; idx < patternsCount; idx += gridThreads)
    {
        // Find a single pattern.
    //    find(finder._pool[idx], patterns[finder._idxs[idx]], delegator);
        find(finder._pool[idx], patterns[idx], delegator);
    }
}
#endif

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _updatePoolSize()
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec, typename TDelegate>
inline void
_updatePoolSize(Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<TSpec> > & finder, TDelegate & /* delegate */)
{
    finder._poolSize = omp_get_max_threads();
}

#ifdef __CUDACC__
template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec, typename TDelegate>
inline void
_updatePoolSizeCuda(Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<TSpec> > & finder, TDelegate & delegate)
{
    typedef Index<TText, TIndexSpec>                    TIndex;
    typedef Finder2<TIndex, TPattern, Multiple<TSpec> > TFinder;
    typedef typename View<TFinder>::Type                TFinderView;
    typedef typename View<TPattern>::Type               TPatternView;

    unsigned blockSize = FinderBlockSize_<TIndex, TPattern, Multiple<TSpec> >::VALUE;

    finder._poolSize = cudaMaxActiveBlocks(_findKernel<TFinderView, TPatternView, TDelegate>, blockSize, 0);
}

template <typename TValue, typename TAlloc, typename TIndexSpec, typename TPattern, typename TSpec, typename TDelegate>
inline void
_updatePoolSize(Finder2<Index<thrust::device_vector<TValue, TAlloc>, TIndexSpec>, TPattern, Multiple<TSpec> > & finder,
                TDelegate & delegate)
{
    _updatePoolSizeCuda(finder, delegate);
}

template <typename TValue, typename TAlloc, typename TSSetSpec, typename TIndexSpec, typename TPattern, typename TSpec, typename TDelegate>
inline void
_updatePoolSize(Finder2<Index<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, TIndexSpec>, TPattern, Multiple<TSpec> > & finder,
                TDelegate & delegate)
{
    _updatePoolSizeCuda(finder, delegate);
}
#endif

// ----------------------------------------------------------------------------
// Function _preparePool()
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec, typename TDelegate>
inline void
_preparePool(Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<TSpec> > & finder, TDelegate & delegate)
{
    typedef Index<TText, TIndexSpec>                                    TIndex;
    typedef typename Finder_<TIndex, TPattern, Multiple<TSpec> >::Type  TFinder_;

    _updatePoolSize(finder, delegate);
    resize(finder._pool, finder._poolSize);
//    resize(finder._pool, finder._poolSize, TFinder_(value(finder.index)));

    std::cout << "Pool size:\t\t\t" << finder._poolSize << std::endl;
}

// ----------------------------------------------------------------------------
// Function _preprocess()
// ----------------------------------------------------------------------------

//template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec, typename TPatterns>
//SEQAN_FUNC void
//_preprocess(Finder2<Index<TText, TIndexSpec>, TPattern, TSpec> & finder, TPatterns const & patterns)
//{
//    _resize(finder, length(patterns));
//    _computeHashes(finder, patterns);
//    _fillIdxsWithIdentity(finder);
//    _sortIdxsByHashes(finder);
//    cudaDeviceSynchronize();
//}

// ----------------------------------------------------------------------------
// Function find()
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec, typename TDelegate>
inline void
find(Finder2<Index<TText, TIndexSpec>, TPattern,  Multiple<TSpec> > & finder, TPattern & pattern, TDelegate & delegate)
{
    typedef Index<TText, TIndexSpec>                                        TIndex;
    typedef Finder2<TIndex, TPattern,  Multiple<TSpec> >                    TFinder;
    typedef Delegator<TFinder, TDelegate>                                   TDelegator;
    typedef typename Iterator<TPattern, Standard>::Type                     TPatternIter;

    // Use a delegator object to delegate this finder instead of the pooled finders.
    TDelegator delegator(finder, delegate);

    // Resize the pool.
    _preparePool(finder, delegator);

    // Find all patterns in parallel.
    TPatternIter patternBegin = begin(pattern, Standard());
    TPatternIter patternEnd = end(pattern, Standard());
    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic))
    for (TPatternIter patternIt = patternBegin; patternIt != patternEnd; ++patternIt)
    {
//        finder._patternIt[omp_get_thread_num()] = patternIt;
        clear(finder._pool[omp_get_thread_num()]);
        find(finder._pool[omp_get_thread_num()], value(patternIt), delegator);
    }
}

// ----------------------------------------------------------------------------
// Function _find()
// ----------------------------------------------------------------------------

#ifdef __CUDACC__
template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec, typename TDelegate>
inline void
_find(Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<TSpec> > & finder,
      TPattern /* const */ & pattern,
      TDelegate & delegate)
{
    // Preprocess patterns.
//    _preprocess(finder, patterns);

    // Resize the pool.
    _preparePool(finder, delegate);

    unsigned blockSize = FinderBlockSize_<Index<TText, TIndexSpec>, TPattern, Multiple<TSpec> >::VALUE;

    std::cout << "Block size:\t\t\t" << blockSize << std::endl;

    _findKernel<<<finder._poolSize, blockSize>>>(view(finder), view(pattern), delegate);
}
#endif

// ----------------------------------------------------------------------------
// Function find()
// ----------------------------------------------------------------------------

#ifdef __CUDACC__
template <typename TValue, typename TAlloc, typename TIndexSpec, typename TPattern, typename TSpec, typename TDelegate>
inline void
find(Finder2<Index<thrust::device_vector<TValue, TAlloc>, TIndexSpec>, TPattern, Multiple<TSpec> > & finder,
     TPattern /* const */ & pattern,
     TDelegate & delegate)
{
    _find(finder, pattern, delegate);
}

template <typename TValue, typename TAlloc, typename TSSetSpec, typename TIndexSpec, typename TPattern, typename TSpec, typename TDelegate>
inline void
find(Finder2<Index<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, TIndexSpec>, TPattern, Multiple<TSpec> > & finder,
     TPattern /* const */ & pattern,
     TDelegate & delegate)
{
    _find(finder, pattern, delegate);
}
#endif

// ----------------------------------------------------------------------------
// Function view()
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec>
inline typename View<Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<TSpec> > >::Type
view(Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<TSpec> > & finder)
{
    typename View<Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<TSpec> > >::Type finderView;

    finderView._pool = view(finder._pool);

//    finderView._idxs = view(finder._idxs);
//    finderView._hashes = view(finder._hashes);

    return finderView;
}

}


using namespace seqan;

// ============================================================================
// Tags
// ============================================================================

struct CPU_;
typedef Tag<CPU_>     CPU;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class HitsCounter
// ----------------------------------------------------------------------------

template <typename TSpec = void>
struct HitsCounter
{
    unsigned long hits;

    HitsCounter() :
        hits(0)
    {}

    template <typename TFinder>
    SEQAN_FUNC void
    operator() (TFinder const & finder)
    {
        SEQAN_OMP_PRAGMA(atomic)
        hits += countOccurrences(finder._pool[omp_get_thread_num()]._textIt);
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function mapReads()                                                    [CPU]
// ----------------------------------------------------------------------------

template <typename TIndex, typename TReadSeqs>
inline void
mapReads(TIndex & index, TReadSeqs & readSeqs, CPU const & /* tag */)
{
    typedef Finder2<TIndex, TReadSeqs, Multiple<FinderSTree> >  TFinder;

    // Instantiate a multiple finder.
    TFinder finder(index);

    // Count hits.
    HitsCounter<> counter;
    find(finder, readSeqs, counter);
    std::cout << "Hits count:\t\t\t" << counter.hits << std::endl;
}


#endif  // #ifndef SEQAN_EXTRAS_CUDAMAPPER_MAPPER_H_
