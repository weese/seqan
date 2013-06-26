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
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function getThreadId()
// ----------------------------------------------------------------------------

//SEQAN_FUNC unsigned getThreadId()
//{
//#ifdef _OPENMP
//    return omp_get_thread_num();
//#elif __CUDA_ARCH__
//    return blockIdx.x * blockDim.x + threadIdx.x;
//#else
//    return 0;
//#endif
//}

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
// Metafunction HistoryStack_
// ----------------------------------------------------------------------------

template <typename TText, typename TViewSpec, typename TIndexSpec, typename TSpec>
struct HistoryStack_<Iter<Index<ContainerView<TText, TViewSpec>, TIndexSpec>,
                     VSTree<TopDown<ParentLinks<TSpec> > > > >
{
    typedef Index<TText, TIndexSpec>                                    TIndex_;
    typedef Iter<TIndex_, VSTree<TopDown<ParentLinks<TSpec> > > >       TIter_;
    typedef typename View<typename HistoryStack_<TIter_>::Type>::Type   Type;
};

template <typename TText, typename TViewSpec, typename TSSetSpec, typename TIndexSpec, typename TSpec>
struct HistoryStack_<Iter<Index<StringSet<ContainerView<TText, TViewSpec>, TSSetSpec>, TIndexSpec>,
                          VSTree<TopDown<ParentLinks<TSpec> > > > >
{
    typedef Index<StringSet<TText, TSSetSpec>, TIndexSpec>              TIndex_;
    typedef Iter<TIndex_, VSTree<TopDown<ParentLinks<TSpec> > > >       TIter_;
    typedef typename View<typename HistoryStack_<TIter_>::Type>::Type   Type;
};

#ifdef __CUDACC__
template <typename TValue, typename TAlloc, typename TIndexSpec, typename TSpec>
struct HistoryStack_<Iter<Index<thrust::device_vector<TValue, TAlloc>, TIndexSpec>,
                     VSTree<TopDown<ParentLinks<TSpec> > > > >
{
    typedef Index<thrust::device_vector<TValue, TAlloc>, TIndexSpec>            TIndex_;
    typedef Iter<TIndex_, VSTree<TopDown<ParentLinks<TSpec> > > >               TIter_;
    typedef thrust::device_vector<typename HistoryStackEntry_<TIter_>::Type>    Type;
};

template <typename TValue, typename TAlloc, typename TSSetSpec, typename TIndexSpec, typename TSpec>
struct HistoryStack_<Iter<Index<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, TIndexSpec>,
                          VSTree<TopDown<ParentLinks<TSpec> > > > >
{
    typedef Index<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, TIndexSpec>  TIndex_;
    typedef Iter<TIndex_, VSTree<TopDown<ParentLinks<TSpec> > > >                           TIter_;
    typedef thrust::device_vector<typename HistoryStackEntry_<TIter_>::Type>                Type;
};
#endif

// ============================================================================
// ============================================================================
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction TextIterator_
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec>
struct TextIterator_
{
    typedef typename Iterator<TText>::Type  Type;
};

template <typename TText, typename TIndexSpec, typename TSpec>
struct TextIterator_<Index<TText, TIndexSpec>, TSpec>
{
    typedef typename Iterator<Index<TText, TIndexSpec>, TopDown<> >::Type  Type;
};

template <typename TText, typename TIndexSpec, typename TDistance, typename TSpec>
struct TextIterator_<Index<TText, TIndexSpec>, Backtracking<TDistance, TSpec> >
{
    typedef typename Iterator<Index<TText, TIndexSpec>, TopDown<ParentLinks<> > >::Type  Type;
};

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class Finder
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec = void>
struct Finder2;

// ----------------------------------------------------------------------------
// Class Finder
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec>
struct Finder2<Index<TText, TIndexSpec>, TPattern, TSpec>
{
    typedef Index<TText, TIndexSpec>            TIndex;

    typename TextIterator_<TIndex, TSpec>::Type _textIt;

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

template <typename TText, typename TIndexSpec, typename TPattern, typename TDistance, typename TSpec, typename TDelegate>
SEQAN_FUNC void
find(Finder2<Index<TText, TIndexSpec>, TPattern, Backtracking<TDistance, TSpec> > & finder,
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
// Metafunction FinderFlyweight_
// ----------------------------------------------------------------------------

template <typename TFinder>
struct FinderFlyweight_;

template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec>
struct FinderFlyweight_<Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<TSpec> > >
{
    typedef Finder2<Index<TText, TIndexSpec>, typename Value<TPattern>::Type, TSpec>    Type;
};

// ----------------------------------------------------------------------------
// Metafunction FinderCTASize_
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec>
struct FinderCTASize_;

template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec>
struct FinderCTASize_<Index<TText, TIndexSpec>, TPattern, Multiple<TSpec> >
{
    static const unsigned VALUE = 256;
};

// ----------------------------------------------------------------------------
// Member FinderText_
// ----------------------------------------------------------------------------

struct FinderText_;

template <typename TText, typename TPattern, typename TSpec>
struct Member<Finder2<TText, TPattern, TSpec>, FinderText_>
{
    typedef Holder<TText>   Type;
};

template <typename TText, typename TViewSpec, typename TIndexSpec, typename TPattern, typename TSpec>
struct Member<Finder2<Index<ContainerView<TText, TViewSpec>, TIndexSpec>, TPattern, TSpec>, FinderText_>
{
    typedef typename View<Index<TText, TIndexSpec> >::Type  Type;
};

template <typename TText, typename TViewSpec, typename TSSetSpec, typename TIndexSpec, typename TPattern, typename TSpec>
struct Member<Finder2<Index<StringSet<ContainerView<TText, TViewSpec>, TSSetSpec>, TIndexSpec>, TPattern, TSpec>, FinderText_>
{
    typedef typename View<Index<StringSet<ContainerView<TText, TViewSpec>, TSSetSpec>, TIndexSpec> >::Type  Type;
};

// ----------------------------------------------------------------------------
// Member FinderHistory_
// ----------------------------------------------------------------------------

struct FinderHistory_;

template <typename TText, typename TPattern, typename TSpec>
struct Member<Finder2<TText, TPattern, Multiple<TSpec> >, FinderHistory_>
{
    typedef typename TextIterator_<TText, TSpec>::Type      TTextIterator_;
    typedef typename HistoryStack_<TTextIterator_>::Type    Type;
};

template <typename TText, typename TViewSpec, typename TIndexSpec, typename TPattern, typename TSpec>
struct Member<Finder2<Index<ContainerView<TText, TViewSpec>, TIndexSpec>, TPattern, Multiple<TSpec> >, FinderHistory_>
{
    typedef Index<TText, TIndexSpec>                        TIndex_;
    typedef Finder2<TIndex_, TPattern, Multiple<TSpec> >    TFinder_;
    typedef typename Member<TFinder_, FinderHistory_>::Type TFinderHistory_;
    typedef typename View<TFinderHistory_>::Type            Type;
};

template <typename TText, typename TViewSpec, typename TSSetSpec, typename TIndexSpec, typename TPattern, typename TSpec>
struct Member<Finder2<Index<StringSet<ContainerView<TText, TViewSpec>, TSSetSpec>, TIndexSpec>, TPattern, Multiple<TSpec> >, FinderHistory_>
{
    typedef Index<StringSet<TText, TSSetSpec>, TIndexSpec>  TIndex_;
    typedef Finder2<TIndex_, TPattern, Multiple<TSpec> >    TFinder_;
    typedef typename Member<TFinder_, FinderHistory_>::Type TFinderHistory_;
    typedef typename View<TFinderHistory_>::Type            Type;
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
    typedef Index<TText, TIndexSpec>                TIndex;

    typename Member<Finder2, FinderText_>::Type     _index;

    Finder2() :
        _index()
    {}

    Finder2(TIndex & index) :
        _index(index)
    {}
};

template <typename TText, typename TIndexSpec, typename TPattern, typename TDistance, typename TSpec>
struct Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<Backtracking<TDistance, TSpec> > >
{
    typedef Index<TText, TIndexSpec>                TIndex;

    typename Member<Finder2, FinderText_>::Type     _index;
    typename Member<Finder2, FinderHistory_>::Type  _history;
    typename Size<TPattern>::Type                   _historyLength;
    // TODO(esiragusa): change type to typename Size<typename Value<TPattern>::Type>::Type.

    Finder2() :
        _index()
    {}

    Finder2(TIndex & index) :
        _index(index)
    {}
};

// ----------------------------------------------------------------------------
// Class Proxy                                                         [Finder]
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec>
class Proxy<Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<TSpec> > >
{
public:
    typedef Index<TText, TIndexSpec>                    TIndex;
    typedef typename TextIterator_<TIndex, TSpec>::Type TTextIterator;
    typedef typename Size<TPattern>::Type               TSize;

    TTextIterator const & _textIt;
    TSize patternId;

    template <typename TFinder>
    SEQAN_FUNC
    Proxy(TFinder const & finder) :
        _textIt(finder._textIt)
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
    typedef Proxy<TFinderView>                              TFinderProxy;
    typedef Delegator<TFinderProxy, TDelegateView>          TDelegator;
    typedef typename FinderFlyweight_<TFinderView>::Type    TFlyweightFinder;

    unsigned threadId = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned gridThreads = gridDim.x * blockDim.x;

    // Instantiate a flyweight finder.
    TFlyweightFinder finderFlyweight = _getFlyweightFinder(finder, threadId);

    // Instantiate a finder proxy to be delegated.
    TFinderProxy finderProxy(finderFlyweight);

    // Instantiate a delegator object to delegate the finder proxy instead of the flyweight finder.
    TDelegator delegator(finderProxy, delegate);

    unsigned patternsCount = length(patterns);

	for (unsigned patternId = threadId; patternId < patternsCount; patternId += gridThreads)
    {
        finderProxy.patternId = patternId;

        // Find a single pattern.
//        find(finderFlyweight, patterns[patternId], delegator);
//        find(finderFlyweight, patterns[finder._idxs[patternId]], delegator);
    }
}
#endif

// ============================================================================
// Functions
// ============================================================================

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
// Function _computeHistoryLength()
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPattern, typename TDistance, typename TSpec>
inline void
_computeHistoryLength(Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<Backtracking<TDistance, TSpec> > > & finder,
                      TPattern const & pattern)
{
    typedef typename Iterator<TPattern const, Standard>::Type   TIterator;
    typedef typename Size<TPattern>::Type                       TPatternSize;

    TPatternSize maxPatternLength = MinValue<TPatternSize>::VALUE;

    // TODO(esiragusa): Add a kernel to get the max string length.
    TIterator patternEnd = end(pattern, Standard());
    for (TIterator patternIt = begin(pattern, Standard()); patternIt != patternEnd; ++patternIt)
        std::max(maxPatternLength, length(value(patternIt)));

    finder._historyLength = maxPatternLength;
}

// ----------------------------------------------------------------------------
// Function _resizeHistory()
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec, typename TSize>
inline void
_resizeHistory(Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<TSpec> > & /* finder */, TSize /* threadsCount */)
{}

template <typename TText, typename TIndexSpec, typename TPattern, typename TDistance, typename TSpec, typename TSize>
inline void
_resizeHistory(Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<Backtracking<TDistance, TSpec> > > & finder,
               TSize threadsCount)
{
    resize(finder._history, threadsCount * finder._historyLength, Exact());
}

// ----------------------------------------------------------------------------
// Function _getFlyweightFinder()
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec, typename TId>
SEQAN_FUNC typename FinderFlyweight_<Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<TSpec> > >::Type
_getFlyweightFinder(Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<TSpec> > & finder, TId /* threadId */)
{
    typedef Index<TText, TIndexSpec>                    TIndex;
    typedef Multiple<TSpec>                             TFinderSpec;
    typedef Finder2<TIndex, TPattern, TFinderSpec>      TFinder;

    return TFlyweightFinder(finder._index);
}

template <typename TText, typename TIndexSpec, typename TPattern, typename TDistance, typename TSpec, typename TId>
SEQAN_FUNC typename FinderFlyweight_<Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<Backtracking<TDistance, TSpec> > > >::Type
_getFlyweightFinder(Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<Backtracking<TDistance, TSpec> > > & finder,
                    TId threadId)
{
    typedef Index<TText, TIndexSpec>                    TIndex;
    typedef Multiple<Backtracking<TDistance, TSpec> >   TFinderSpec;
    typedef Finder2<TIndex, TPattern, TFinderSpec>      TFinder;
    typedef typename FinderFlyweight_<TFinder>::Type    TFlyweightFinder;

    TFlyweightFinder finderFlyweight(finder._index);

    finderFlyweight._textIt.history._begin = begin(finder._history, Standard()) + finder._historyLength * threadId;
    finderFlyweight._textIt.history._end = end(finder._history, Standard()) + finder._historyLength * (threadId + 1);

    return finderFlyweight;
}

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

template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec, typename TDelegate>
inline void
_find(Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<TSpec> > & finder,
      TPattern /* const */ & pattern,
      TDelegate & delegate)
{
    typedef Index<TText, TIndexSpec>                    TIndex;
    typedef Finder2<TIndex, TPattern, Multiple<TSpec> > TFinder;
    typedef typename View<TFinder>::Type                TFinderView;
    typedef typename View<TPattern>::Type               TPatternView;

    // Preprocess patterns.
//    _preprocess(finder, patterns);

    // Compute grid size.
    unsigned ctaSize = FinderCTASize_<TIndex, TPattern, Multiple<TSpec> >::VALUE;
    unsigned activeBlocks = cudaMaxActiveBlocks(_findKernel<TFinderView, TPatternView, TDelegate>, ctaSize, 0);

    std::cout << "CTA Size:\t\t\t" << ctaSize << std::endl;
    std::cout << "Active Blocks:\t\t\t" << activeBlocks << std::endl;

    // Initialize the history.
    _computeHistoryLength(finder, pattern);
    _resizeHistory(finder, activeBlocks * ctaSize);

    // Launch the find kernel.
    _findKernel<<<activeBlocks, ctaSize>>>(view(finder), view(pattern), delegate);
}
#endif

// ----------------------------------------------------------------------------
// Function find()
// ----------------------------------------------------------------------------

//template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec, typename TDelegate>
//inline void
//find(Finder2<Index<TText, TIndexSpec>, TPattern,  Multiple<TSpec> > & finder, TPattern & pattern, TDelegate & delegate)
//{
//    typedef Index<TText, TIndexSpec>                                        TIndex;
//    typedef Finder2<TIndex, TPattern,  Multiple<TSpec> >                    TFinder;
//    typedef Delegator<TFinder, TDelegate>                                   TDelegator;
//    typedef typename Iterator<TPattern, Standard>::Type                     TPatternIter;
//
//    // Use a delegator object to delegate this finder instead of the pooled finders.
//    TDelegator delegator(finder, delegate);
//
//    // Initialize the pool.
//    _initPool(finder, omp_get_max_threads());
//
//    // Find all patterns in parallel.
//    TPatternIter patternBegin = begin(pattern, Standard());
//    TPatternIter patternEnd = end(pattern, Standard());
//    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic))
//    for (TPatternIter patternIt = patternBegin; patternIt != patternEnd; ++patternIt)
//    {
////        finder._patternIt[omp_get_thread_num()] = patternIt;
//        clear(finder._pool[omp_get_thread_num()]);
//        find(finder._pool[omp_get_thread_num()], value(patternIt), delegator);
//    }
//}

// ----------------------------------------------------------------------------
// Function view()
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec>
inline typename View<Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<TSpec> > >::Type
view(Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<TSpec> > & finder)
{
    typename View<Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<TSpec> > >::Type finderView;

    finderView._index = view(value(finder._index));

    return finderView;
}

template <typename TText, typename TIndexSpec, typename TPattern, typename TDistance, typename TSpec>
inline typename View<Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<Backtracking<TDistance, TSpec> > > >::Type
view(Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<Backtracking<TDistance, TSpec> > > & finder)
{
    typename View<Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<Backtracking<TDistance, TSpec> > > >::Type finderView;

    finderView._index = view(value(finder._index));
    finderView._history = view(finder._history);
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
    operator() (TFinder const & /* finder */)
    {
//        SEQAN_OMP_PRAGMA(atomic)
//        hits += countOccurrences(finder._pool[omp_get_thread_num()]._textIt);
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
//    typedef Finder2<TIndex, TReadSeqs, Multiple<FinderSTree> >  TFinder;
//
//    // Instantiate a multiple finder.
//    TFinder finder(index);
//
//    // Count hits.
//    HitsCounter<> counter;
//    find(finder, readSeqs, counter);
//    std::cout << "Hits count:\t\t\t" << counter.hits << std::endl;
}


#endif  // #ifndef SEQAN_EXTRAS_CUDAMAPPER_MAPPER_H_
