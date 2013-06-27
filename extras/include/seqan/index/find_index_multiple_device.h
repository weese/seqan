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

#ifndef SEQAN_EXTRAS_INDEX_FIND_INDEX_MULTIPLE_DEVICE_H
#define SEQAN_EXTRAS_INDEX_FIND_INDEX_MULTIPLE_DEVICE_H

namespace seqan {

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

// TODO(esiragusa): FinderFlyweight_ for CPU ?

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


// ============================================================================
// Kernels
// ============================================================================

// ----------------------------------------------------------------------------
// Kernel _computeHashesKernel()
// ----------------------------------------------------------------------------

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

// ----------------------------------------------------------------------------
// Kernel _findKernel()
// ----------------------------------------------------------------------------

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
//        finderProxy._patternIt = begin(patterns, Standard()) + patternId;
        finderProxy._patternIt = patternId;

        // Find a single pattern.
        find(finderFlyweight, patterns[patternId], delegator);
//        find(finderFlyweight, patterns[finder._idxs[patternId]], delegator);
    }
}

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

template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec>
inline void
_computeHistoryLength(Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<TSpec> > & /* finder */,
                      TPattern const & /* pattern */)
{}

template <typename TText, typename TIndexSpec, typename TPattern, typename TDistance, typename TSpec>
inline void
_computeHistoryLength(Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<Backtracking<TDistance, TSpec> > > & finder,
                      TPattern const & /* pattern */)
{
    typedef typename Iterator<TPattern const, Standard>::Type   TIterator;
    typedef typename Size<TPattern>::Type                       TPatternSize;

    // TODO(esiragusa): Add a kernel to get the max string length.
    finder._historyLength = 100;

//    TPatternSize maxPatternLength = MinValue<TPatternSize>::VALUE;
//
//    TIterator patternEnd = end(pattern, Standard());
//    for (TIterator patternIt = begin(pattern, Standard()); patternIt != patternEnd; ++patternIt)
//        std::max(maxPatternLength, length(value(patternIt)));
//
//    finder._historyLength = maxPatternLength;
}

// ----------------------------------------------------------------------------
// Function _resizeHistory()
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec, typename TSize>
inline void
_resizeHistory(Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<TSpec> > & /* finder */,
               TSize /* threadsCount */)
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
    typedef typename FinderFlyweight_<TFinder>::Type    TFlyweightFinder;

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
    finderFlyweight._textIt.history._end = finderFlyweight._textIt.history._begin;
    finderFlyweight._textIt.history._capacity = finder._historyLength;

    // TODO(esiragusa): Maybe use a StringSet of Histories and take a view of single Strings?
//    finderFlyweight._textIt.history = finder._history[threadId];
//    clear(finderFlyweight._textIt.history);

    return finderFlyweight;
}

// ----------------------------------------------------------------------------
// Function find()
// ----------------------------------------------------------------------------

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
    typedef typename View<TDelegate>::Type              TDelegateView;

    // Preprocess patterns.
//    _preprocess(finder, patterns);

    // Compute grid size.
    unsigned ctaSize = FinderCTASize_<TIndex, TPattern, Multiple<TSpec> >::VALUE;
    unsigned activeBlocks = cudaMaxActiveBlocks(_findKernel<TFinderView, TPatternView, TDelegateView>, ctaSize, 0);

    std::cout << "CTA Size:\t\t\t" << ctaSize << std::endl;
    std::cout << "Active Blocks:\t\t\t" << activeBlocks << std::endl;

    // Initialize the history.
    _computeHistoryLength(finder, pattern);
    _resizeHistory(finder, activeBlocks * ctaSize);

    // Launch the find kernel.
    _findKernel<<<activeBlocks, ctaSize>>>(view(finder), view(pattern), view(delegate));
}

}

#endif  // #ifndef SEQAN_EXTRAS_INDEX_FIND_INDEX_MULTIPLE_DEVICE_H
