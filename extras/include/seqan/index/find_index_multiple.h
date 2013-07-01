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

#ifndef SEQAN_EXTRAS_INDEX_FIND_INDEX_MULTIPLE_H
#define SEQAN_EXTRAS_INDEX_FIND_INDEX_MULTIPLE_H

#ifdef PLATFORM_CUDA
#include <thrust/sort.h>
#include <thrust/reduce.h>
#endif

namespace seqan {

// ============================================================================
// Tags
// ============================================================================

// ----------------------------------------------------------------------------
// Tag Multiple
// ----------------------------------------------------------------------------

template <typename TSpec>
struct Multiple;

// ----------------------------------------------------------------------------
// Finder Member Tags
// ----------------------------------------------------------------------------

struct Factory_;

// ----------------------------------------------------------------------------
// Pattern Member Tags
// ----------------------------------------------------------------------------

struct Needles_;
struct Hashes_;
struct Permutation_;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class Finder; Multiple
// ----------------------------------------------------------------------------

template <typename TText, typename TPatterns, typename TSpec>
struct Finder2<TText, TPatterns, Multiple<TSpec> >
{
    typename Member<Finder2, Factory_>::Type    _factory;

    SEQAN_HOST_DEVICE
    Finder2() {}

    Finder2(TText & text) :
        _factory(text)
    {}
};

// ----------------------------------------------------------------------------
// Class Finder; Index, Multiple
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPatterns, typename TSpec>
struct Finder2<Index<TText, TIndexSpec>, TPatterns, Multiple<TSpec> >
{
    typedef Index<TText, TIndexSpec>            TIndex;

    typename Member<Finder2, Factory_>::Type    _factory;

    SEQAN_HOST_DEVICE
    Finder2() {}

    Finder2(TIndex & index) :
        _factory(index)
    {}
};

// ----------------------------------------------------------------------------
// Class Pattern; Multiple
// ----------------------------------------------------------------------------

template <typename TNeedles, typename TSpec>
class Pattern<TNeedles, Multiple<TSpec> >
{
public:
    typename Member<Pattern, Needles_>::Type        data_host;
    typename Member<Pattern, Hashes_>::Type         _hashes;
    typename Member<Pattern, Permutation_>::Type    _permutation;

    SEQAN_HOST_DEVICE
    Pattern() {}

    Pattern(TNeedles const & needles) :
        data_host(needles)
    {
        _preprocess(*this);
    }

//    Pattern(TNeedles const & needles)
//    {
//        setHost(*this, needles);
//    }
};

// ----------------------------------------------------------------------------
// Class Proxy                                                         [Finder]
// ----------------------------------------------------------------------------

template <typename TText, typename TPatterns, typename TSpec>
class Proxy<Finder2<TText, TPatterns, Multiple<TSpec> > >
{
public:
    typedef typename TextIterator_<TText, TSpec>::Type      TTextIterator;
    typedef typename Iterator<TPatterns, Standard>::Type    TPatternIterator;

    TTextIterator const &   _textIt;
    unsigned                _patternIt;

    template <typename TFinder>
    SEQAN_HOST_DEVICE
    Proxy(TFinder const & finder) :
        _textIt(textIterator(finder))
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Host                                                  [Pattern]
// ----------------------------------------------------------------------------

template <typename TNeedles, typename TSpec>
struct Host<Pattern<TNeedles, Multiple<TSpec> > >
{
    typedef TNeedles Type;
};

// ----------------------------------------------------------------------------
// Metafunction PatternShape_                                         [Pattern]
// ----------------------------------------------------------------------------
// TODO(esiragusa): Automatically select a good shape.

template <typename TPattern>
struct PatternShape_
{
    typedef typename Host<TPattern>::Type               TNeedles_;
    typedef typename Value<TNeedles_>::Type             TNeedle_;
    typedef typename Value<TNeedle_>::Type              TAlphabet_;

    typedef Shape<TAlphabet_, UngappedShape<10> >       Type;
};

// ----------------------------------------------------------------------------
// Member Needles_                                                    [Pattern]
// ----------------------------------------------------------------------------

template <typename TNeedles, typename TSpec>
struct Member<Pattern<TNeedles, Multiple<TSpec> >, Needles_>
{
    typedef typename IfView<TNeedles, TNeedles, Holder<TNeedles> >::Type    Type;
};

// ----------------------------------------------------------------------------
// Member Hashes_                                                     [Pattern]
// ----------------------------------------------------------------------------

template <typename TNeedles, typename TSpec>
struct Member<Pattern<TNeedles, Multiple<TSpec> >, Hashes_>
{
//    typedef Pattern<TNeedles, Multiple<TSpec> >     TPattern_;
//    typedef typename PatternShape_<TPattern_>::Type TShape_;
//    typedef typename Value<TShape_>::Type           THash_;
//
//    typedef String<THash_>                          Type;

    typedef Nothing                                 Type;
};

#ifdef PLATFORM_CUDA
template <typename TValue, typename TAlloc, typename TSSetSpec, typename TSpec>
struct Member<Pattern<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, Multiple<TSpec> >, Hashes_>
{
    typedef thrust::device_vector<TValue, TAlloc>   TNeedle_;
    typedef StringSet<TNeedle_, TSSetSpec>          TNeedles_;
    typedef Pattern<TNeedles_, Multiple<TSpec> >    TPattern_;
    typedef typename PatternShape_<TPattern_>::Type TShape_;
    typedef typename Value<TShape_>::Type           THash_;

    typedef thrust::device_vector<THash_>           Type;
};
#endif

template <typename TNeedle, typename TViewSpec, typename TSSetSpec, typename TSpec>
struct Member<Pattern<StringSet<ContainerView<TNeedle, TViewSpec>, TSSetSpec>, Multiple<TSpec> >, Hashes_>
{
    typedef Pattern<StringSet<TNeedle, TSSetSpec>, Multiple<TSpec> >        TPattern_;
    typedef typename View<typename Member<TPattern_, Hashes_>::Type>::Type  Type;
};

// ----------------------------------------------------------------------------
// Member Permutation_                                                [Pattern]
// ----------------------------------------------------------------------------

template <typename TNeedles, typename TSpec>
struct Member<Pattern<TNeedles, Multiple<TSpec> >, Permutation_>
{
    typedef Nothing Type;
};

#ifdef PLATFORM_CUDA
template <typename TValue, typename TAlloc, typename TSSetSpec, typename TSpec>
struct Member<Pattern<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, Multiple<TSpec> >, Permutation_>
{
    typedef thrust::device_vector<TValue, TAlloc>   TNeedle_;
    typedef StringSet<TNeedle_, TSSetSpec>          TNeedles_;
    typedef typename Size<TNeedles_>::Type          TSize_;

    typedef thrust::device_vector<TSize_>           Type;
};
#endif

template <typename TNeedle, typename TViewSpec, typename TSSetSpec, typename TSpec>
struct Member<Pattern<StringSet<ContainerView<TNeedle, TViewSpec>, TSSetSpec>, Multiple<TSpec> >, Permutation_>
{
    typedef Pattern<StringSet<TNeedle, TSSetSpec>, Multiple<TSpec> >            TPattern_;
    typedef typename View<typename Member<TPattern_, Permutation_>::Type>::Type Type;
};

// ----------------------------------------------------------------------------
// Member Factory_                                                     [Finder]
// ----------------------------------------------------------------------------

template <typename TText, typename TPatterns, typename TSpec>
struct Member<Finder2<TText, TPatterns, Multiple<TSpec> >, Factory_>
{
    typedef Factory<typename TextIterator_<TText, TSpec>::Type>  Type;
};

// ----------------------------------------------------------------------------
// Metafunction Delegated                                              [Finder]
// ----------------------------------------------------------------------------

template <typename TText, typename TPatterns, typename TSpec>
struct Delegated<Finder2<TText, TPatterns, Multiple<TSpec> > >
{
    typedef Proxy<Finder2<TText, TPatterns, Multiple<TSpec> > > Type;
};

// ----------------------------------------------------------------------------
// Metafunction FinderCTASize_                                         [Finder]
// ----------------------------------------------------------------------------

template <typename TFinder>
struct FinderCTASize_
{
    static const unsigned VALUE = 256;
};

// ============================================================================
// Kernels
// ============================================================================

// ----------------------------------------------------------------------------
// Kernel _preprocessKernel()                                         [Pattern]
// ----------------------------------------------------------------------------

#ifdef PLATFORM_CUDA
template <typename TNeedles, typename TSpec>
SEQAN_GLOBAL void
_preprocessKernel(Pattern<TNeedles, Multiple<TSpec> > pattern)
{
    typedef Pattern<TNeedles, Multiple<TSpec> >     TPattern;
    typedef typename PatternShape_<TPattern>::Type  TShape;

    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    // Return silently if there is no job left.
    if (idx >= length(pattern.data_host)) return;

    // Compute the hash of a needle.
    TShape shape;
    pattern._hashes[idx] = hash(shape, begin(pattern.data_host[idx], Standard()));

    // Fill with the identity permutation.
    pattern._permutation[idx] = idx;
}
#endif

// ----------------------------------------------------------------------------
// Kernel _findKernel()                                                [Finder]
// ----------------------------------------------------------------------------

#ifdef PLATFORM_CUDA
template <typename TText, typename TPatterns, typename TSpec, typename TDelegate>
SEQAN_GLOBAL void
_findKernel(Finder2<TText, TPatterns, Multiple<TSpec> > finder, TPatterns patterns, TDelegate delegate)
{
    typedef typename Value<TPatterns>::Type                 TPattern;
    typedef Finder2<TText, TPattern, TSpec>                 TFinderSimple;
    typedef Finder2<TText, TPatterns, Multiple<TSpec> >     TFinderView;
    typedef Proxy<TFinderView>                              TFinderProxy;
    typedef Delegator<TFinderProxy, TDelegate>              TDelegator;

    unsigned threadId = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned gridThreads = gridDim.x * blockDim.x;

    // Instantiate a simple finder.
    TFinderSimple simpleFinder = getObject(finder._factory, threadId);

    // Instantiate a finder proxy to be delegated.
    TFinderProxy finderProxy(simpleFinder);

    // Instantiate a delegator object to delegate the finder proxy instead of the serial finder.
    TDelegator delegator(finderProxy, delegate);

    unsigned patternsCount = length(patterns);

	for (unsigned patternId = threadId; patternId < patternsCount; patternId += gridThreads)
    {
        finderProxy._patternIt = patternId;
//        finderProxy._patternIt = threadId;

        // Find a single pattern.
        find(simpleFinder, patterns[finderProxy._patternIt], delegator);
    }
}
#endif

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _preprocess()                                   [Pattern; ExecHost]
// ----------------------------------------------------------------------------

template <typename TNeedles, typename TSpec>
inline void
_preprocess(Pattern<TNeedles, Multiple<TSpec> > & /* pattern */, ExecHost const & /* tag */) {}

// ----------------------------------------------------------------------------
// Function _preprocess()                                 [Pattern; ExecDevice]
// ----------------------------------------------------------------------------

#ifdef PLATFORM_CUDA
template <typename TNeedles, typename TSpec>
inline void
_preprocess(Pattern<TNeedles, Multiple<TSpec> > & pattern, ExecDevice const & /* tag */)
{
    typedef typename Size<TNeedles>::Type   TSize;

    TSize needlesCount = length(needle(pattern));

    resize(pattern._hashes, needlesCount, Exact());
    resize(pattern._permutation, needlesCount, Exact());

    // Compute grid size.
    unsigned ctaSize = 256;
    unsigned activeBlocks = (needlesCount + ctaSize - 1) / ctaSize;

    // Launch the preprocessing kernel.
    _preprocessKernel<<<activeBlocks, ctaSize>>>(view(pattern));

    // Sort the patterns.
    thrust::sort_by_key(pattern._hashes.begin(), pattern._hashes.end(), pattern._permutation.begin());

    cudaDeviceSynchronize();
}
#endif

// ----------------------------------------------------------------------------
// Function _preprocess()                                             [Pattern]
// ----------------------------------------------------------------------------

template <typename TNeedles, typename TSpec>
inline void
_preprocess(Pattern<TNeedles, Multiple<TSpec> > & pattern)
{
    typedef Pattern<TNeedles, Multiple<TSpec> > TPattern;

    _preprocess(pattern, typename ExecSpace<TPattern>::Type());
}

// ----------------------------------------------------------------------------
// Function setHost()                                                 [Pattern]
// ----------------------------------------------------------------------------

//template <typename TNeedles, typename TSpec, typename TOtherNeedles>
//SEQAN_FUNC void
//setHost(Pattern<TNeedles, Multiple<TSpec> > & pattern, TOtherNeedles const & needles)
//{
//}

// ----------------------------------------------------------------------------
// Function _find()                                          [Finder; ExecHost]
// ----------------------------------------------------------------------------

template <typename TText, typename TPatterns, typename TSpec, typename TDelegate>
inline void
_find(Finder2<TText, TPatterns, Multiple<TSpec> > & finder,
      TPatterns & patterns,
      TDelegate & delegate,
      ExecHost const & /* tag */)
{
    typedef Finder2<TText, TPatterns,  Multiple<TSpec> >    TFinder;
    typedef typename Value<TPatterns>::Type                 TPattern;
    typedef Finder2<TText, TPattern, TSpec>                 TFinderSimple;
    typedef typename View<TFinder>::Type                    TFinderView;
    typedef typename View<TFinderSimple>::Type              TFinderSimpleView;
    typedef Proxy<TFinderView>                              TFinderProxy;
    typedef Delegator<TFinderProxy, TDelegate>              TDelegator;
    typedef typename Iterator<TPatterns, Standard>::Type    TPatternsIter;

    // Initialize the iterator factory.
    setMaxHistoryLength(finder._factory, length(back(patterns)));
    setMaxObjects(finder._factory, omp_get_max_threads());
    build(finder._factory);

    // Use a finder view.
    TFinderView finderView = view(finder);

    // Instantiate a finder.
    TFinderSimpleView simpleFinder = getObject(finderView._factory, omp_get_thread_num());

    // Instantiate a finder proxy to be delegated.
    TFinderProxy finderProxy(simpleFinder);

    // Instantiate a delegator object to delegate the finder proxy instead of the simple finder.
    TDelegator delegator(finderProxy, delegate);

    // Find all patterns in parallel.
    TPatternsIter patternsBegin = begin(patterns, Standard());
    TPatternsIter patternsEnd = end(patterns, Standard());
//    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic))
    for (TPatternsIter patternsIt = patternsBegin; patternsIt != patternsEnd; ++patternsIt)
    {
        clear(simpleFinder);
        find(simpleFinder, value(patternsIt), delegator);
    }
}

// ----------------------------------------------------------------------------
// Function _find();                                       [Finder; ExecDevice]
// ----------------------------------------------------------------------------

#ifdef PLATFORM_CUDA
template <typename TText, typename TPatterns, typename TSpec, typename TDelegate>
inline void
_find(Finder2<TText, TPatterns, Multiple<TSpec> > & finder,
      TPatterns /* const */ & patterns,
      TDelegate & delegate,
      ExecDevice const & /* tag */)
{
    typedef Finder2<TText, TPatterns, Multiple<TSpec> > TFinder;
    typedef typename View<TFinder>::Type                TFinderView;
    typedef typename View<TText>::Type                  TTextView;
    typedef typename View<TPatterns>::Type              TPatternsView;
    typedef typename View<TDelegate>::Type              TDelegateView;

    // Preprocess patterns.
//    _preprocess(finder, patterns);

    // Compute grid size.
    unsigned ctaSize = FinderCTASize_<TFinderView>::VALUE;
    unsigned activeBlocks = cudaMaxActiveBlocks(_findKernel<TTextView, TPatternsView, TSpec, TDelegateView>, ctaSize, 0);
//    unsigned activeBlocks = (length(patterns) + ctaSize - 1) / ctaSize;

    std::cout << "CTA Size:\t\t\t" << ctaSize << std::endl;
    std::cout << "Active Blocks:\t\t\t" << activeBlocks << std::endl;

    // Initialize the iterator factory.
    setMaxHistoryLength(finder._factory, length(back(patterns)));
//    setMaxObjects(finder._factory, length(patterns));
    setMaxObjects(finder._factory, activeBlocks * ctaSize);
    build(finder._factory);

    // Launch the find kernel.
    _findKernel<<<activeBlocks, ctaSize>>>(view(finder), view(patterns), view(delegate));
}
#endif

// ----------------------------------------------------------------------------
// Function find()                                                     [Finder]
// ----------------------------------------------------------------------------

template <typename TText, typename TPatterns, typename TSpec, typename TDelegate>
inline void
find(Finder2<TText, TPatterns, Multiple<TSpec> > & finder, TPatterns & patterns, TDelegate & delegate)
{
    typedef Finder2<TText, TPatterns,  Multiple<TSpec> > TFinder;

    _find(finder, patterns, delegate, typename ExecSpace<TFinder>::Type());
}

// ----------------------------------------------------------------------------
// Function view()                                                    [Pattern]
// ----------------------------------------------------------------------------

template <typename TNeedles, typename TSpec>
inline typename View<Pattern<TNeedles, Multiple<TSpec> > >::Type
view(Pattern<TNeedles, Multiple<TSpec> > & pattern)
{
    typename View<Pattern<TNeedles, Multiple<TSpec> > >::Type   patternView;

    patternView.data_host = view(needle(pattern));
    patternView._hashes = view(pattern._hashes);
    patternView._permutation = view(pattern._permutation);

    return patternView;
}

// ----------------------------------------------------------------------------
// Function view()                                                     [Finder]
// ----------------------------------------------------------------------------

template <typename TText, typename TPatterns, typename TSpec>
inline typename View<Finder2<TText, TPatterns, Multiple<TSpec> > >::Type
view(Finder2<TText, TPatterns, Multiple<TSpec> > & finder)
{
    typename View<Finder2<TText, TPatterns, Multiple<TSpec> > >::Type finderView;

    finderView._factory = view(finder._factory);
//    finderView._idxs = view(finder._idxs);
//    finderView._hashes = view(finder._hashes);

    return finderView;
}

// ----------------------------------------------------------------------------
// Function textIterator()                                       [Finder Proxy]
// ----------------------------------------------------------------------------

template <typename TText, typename TPatterns, typename TSpec>
inline SEQAN_HOST_DEVICE typename TextIterator_<TText, Multiple<TSpec> >::Type &
textIterator(Proxy<Finder2<TText, TPatterns, Multiple<TSpec> > > & finder)
{
    return finder._textIt;
}

template <typename TText, typename TPatterns, typename TSpec>
inline SEQAN_HOST_DEVICE typename TextIterator_<TText, Multiple<TSpec> >::Type const &
textIterator(Proxy<Finder2<TText, TPatterns, Multiple<TSpec> > > const & finder)
{
    return finder._textIt;
}

}

#endif  // #ifndef SEQAN_EXTRAS_INDEX_FIND_INDEX_MULTIPLE_H
