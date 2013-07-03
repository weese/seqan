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

template <typename TText, typename TPattern, typename TSpec>
struct Finder2<TText, TPattern, Multiple<TSpec> >
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

template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec>
struct Finder2<Index<TText, TIndexSpec>, TPattern, Multiple<TSpec> >
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

template <typename TText, typename TPattern, typename TSpec>
class Proxy<Finder2<TText, TPattern, Multiple<TSpec> > >
{
public:
    typedef typename TextIterator_<TText, TSpec>::Type     TTextIterator;
    typedef typename Iterator<TPattern, Standard>::Type    TPatternIterator;

    TTextIterator const &   _textIt;
    unsigned                _patternIt;

    template <typename TFinder>
    SEQAN_HOST_DEVICE
    Proxy(TFinder const & finder) :
        _textIt(textIterator(finder))
    {}
};

// ----------------------------------------------------------------------------
// Class FinderContext_
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec, typename TDelegate>
struct FinderContext_;

template <typename TText, typename TPattern, typename TSpec, typename TDelegate>
struct FinderContext_<TText, TPattern, Multiple<TSpec>, TDelegate>
{
    typedef typename Needle<TPattern>::Type                 TNeedles;
    typedef typename Value<TPattern>::Type                  TNeedle;
    typedef Finder2<TText, TNeedle, TSpec>                  TSingleFinder;
    typedef Finder2<TText, TPattern, Multiple<TSpec> >      TFinder;
    typedef Proxy<TFinder>                                  TFinderProxy;
    typedef Delegator<TFinderProxy, TDelegate>              TDelegator;
    typedef typename Member<TFinder, Factory_>::Type        TFactory;

    TFactory &      factory;
    TSingleFinder   finder;
    TFinderProxy    proxy;
    TDelegator      delegator;

    explicit
    FinderContext_(TFactory & factory, TDelegate & delegate) :
        factory(factory),
        finder(getObject(factory, getThreadId())),
        proxy(finder),
        delegator(proxy, delegate)
    {}

    // NOTE(esiragusa): This is called on firstprivate.
    FinderContext_(FinderContext_ & ctx) :
        factory(ctx.factory),
        finder(getObject(factory, getThreadId())),
        proxy(finder),
        delegator(proxy, ctx.delegator.delegate)
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
    typedef Nothing Type;
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

#ifdef PLATFORM_CUDA
template <typename TValue, typename TAlloc, typename TViewSpec, typename TAllocOuter, typename TSpec>
struct Member<Pattern<thrust::device_vector<ContainerView<thrust::device_vector<TValue, TAlloc>, TViewSpec>, TAllocOuter>, Multiple<TSpec> >, Hashes_>
{
    typedef thrust::device_vector<TValue, TAlloc>               TNeedle_;
    typedef ContainerView<TNeedle_, TViewSpec>                  TNeedleView_;
    typedef thrust::device_vector<TNeedleView_, TAllocOuter>    TNeedles_;
    typedef Pattern<TNeedles_, Multiple<TSpec> >                TPattern_;
    typedef typename PatternShape_<TPattern_>::Type             TShape_;
    typedef typename Value<TShape_>::Type                       THash_;

    typedef thrust::device_vector<THash_>                       Type;
};
#endif

template <typename TNeedle, typename TViewSpec, typename TSSetSpec, typename TSpec>
struct Member<Pattern<StringSet<ContainerView<TNeedle, TViewSpec>, TSSetSpec>, Multiple<TSpec> >, Hashes_>
{
    typedef Pattern<StringSet<TNeedle, TSSetSpec>, Multiple<TSpec> >        TPattern_;
    typedef typename View<typename Member<TPattern_, Hashes_>::Type>::Type  Type;
};

template <typename TNeedles, typename TViewSpec, typename TSpec>
struct Member<Pattern<ContainerView<TNeedles, TViewSpec>, Multiple<TSpec> >, Hashes_>
{
    typedef Pattern<TNeedles, Multiple<TSpec> >                             TPattern_;
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

#ifdef PLATFORM_CUDA
template <typename TValue, typename TAlloc, typename TViewSpec, typename TAllocOuter, typename TSpec>
struct Member<Pattern<thrust::device_vector<ContainerView<thrust::device_vector<TValue, TAlloc>, TViewSpec>, TAllocOuter>, Multiple<TSpec> >, Permutation_>
{
    typedef thrust::device_vector<TValue, TAlloc>               TNeedle_;
    typedef ContainerView<TNeedle_, TViewSpec>                  TNeedleView_;
    typedef thrust::device_vector<TNeedleView_, TAllocOuter>    TNeedles_;
    typedef typename Size<TNeedles_>::Type                      TSize_;

    typedef thrust::device_vector<TSize_>                       Type;
};
#endif

template <typename TNeedle, typename TViewSpec, typename TSSetSpec, typename TSpec>
struct Member<Pattern<StringSet<ContainerView<TNeedle, TViewSpec>, TSSetSpec>, Multiple<TSpec> >, Permutation_>
{
    typedef Pattern<StringSet<TNeedle, TSSetSpec>, Multiple<TSpec> >            TPattern_;
    typedef typename View<typename Member<TPattern_, Permutation_>::Type>::Type Type;
};

template <typename TNeedles, typename TViewSpec, typename TSpec>
struct Member<Pattern<ContainerView<TNeedles, TViewSpec>, Multiple<TSpec> >, Permutation_>
{
    typedef Pattern<TNeedles, Multiple<TSpec> >                                 TPattern_;
    typedef typename View<typename Member<TPattern_, Permutation_>::Type>::Type Type;
};

// ----------------------------------------------------------------------------
// Member Factory_                                                     [Finder]
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec>
struct Member<Finder2<TText, TPattern, Multiple<TSpec> >, Factory_>
{
    typedef Factory<typename TextIterator_<TText, TSpec>::Type>  Type;
};

// ----------------------------------------------------------------------------
// Metafunction Delegated                                              [Finder]
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec>
struct Delegated<Finder2<TText, TPattern, Multiple<TSpec> > >
{
    typedef Proxy<Finder2<TText, TPattern, Multiple<TSpec> > > Type;
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

    unsigned threadId = getThreadId();

    // Return silently if there is no job left.
    if (threadId >= length(pattern.data_host)) return;

    // Compute the hash of a needle.
    TShape shape;
    pattern._hashes[threadId] = hash(shape, begin(pattern.data_host[threadId], Standard()));

    // Fill with the identity permutation.
    pattern._permutation[threadId] = threadId;
}
#endif

// ----------------------------------------------------------------------------
// Kernel _findKernel()                                                [Finder]
// ----------------------------------------------------------------------------

#ifdef PLATFORM_CUDA
template <typename TText, typename TPattern, typename TSpec, typename TDelegate>
SEQAN_GLOBAL void
_findKernel(Finder2<TText, TPattern, Multiple<TSpec> > finder, TPattern pattern, TDelegate delegate)
{
    typedef FinderContext_<TText, TPattern, Multiple<TSpec>, TDelegate> TFinderContext;

    unsigned threadId = getThreadId();

    // Return if there is no job left.
    if (threadId >= length(pattern.data_host)) return;

    // Instantiate a thread context.
    TFinderContext ctx(finder._factory, delegate);

    // Get the sorted needle id.
    ctx.proxy._patternIt = pattern._permutation[threadId];

    // Find a single needle.
    find(ctx.finder, pattern.data_host[ctx.proxy._patternIt], ctx.delegator);
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

    // Sort the pattern.
    thrust::sort_by_key(pattern._hashes.begin(), pattern._hashes.end(), pattern._permutation.begin());
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

template <typename TText, typename TPattern, typename TSpec, typename TDelegate>
inline void
_find(Finder2<TText, TPattern, Multiple<TSpec> > & finder,
      TPattern & pattern,
      TDelegate & delegate,
      ExecHost const & /* tag */)
{
    typedef Finder2<TText, TPattern,  Multiple<TSpec> >                         TFinder;
    typedef typename View<TText>::Type                                          TTextView;
    typedef typename View<TPattern>::Type                                       TPatternView;
    typedef typename View<TFinder>::Type                                        TFinderView;
    typedef FinderContext_<TTextView, TPatternView, Multiple<TSpec>, TDelegate> TFinderContext;
    typedef typename Needle<TPattern>::Type                                     TNeedles;
    typedef typename Size<TNeedles>::Type                                       TSize;

    // Initialize the iterator factory.
    setMaxHistoryLength(finder._factory, length(back(needle(pattern))));
    setMaxObjects(finder._factory, omp_get_max_threads());
    build(finder._factory);

    // Use a finder view.
    TFinderView finderView = view(finder);

    // Instantiate a thread context.
    // NOTE(esiragusa): Each thread initializes its private context on firstprivate.
    TFinderContext ctx(finderView._factory, delegate);

    TNeedles & needles = needle(pattern);
    TSize needlesCount = length(needles);

    // Find all needles in parallel.
    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic) firstprivate(ctx))
    for (TSize needleId = 0; needleId < needlesCount; ++needleId)
    {
        clear(ctx.finder);
        find(ctx.finder, needles[needleId], ctx.delegator);
    }
}

// ----------------------------------------------------------------------------
// Function _find();                                       [Finder; ExecDevice]
// ----------------------------------------------------------------------------

#ifdef PLATFORM_CUDA
template <typename TText, typename TPattern, typename TSpec, typename TDelegate>
inline void
_find(Finder2<TText, TPattern, Multiple<TSpec> > & finder,
      TPattern /* const */ & pattern,
      TDelegate & delegate,
      ExecDevice const & /* tag */)
{
    typedef Finder2<TText, TPattern, Multiple<TSpec> >  TFinder;
    typedef typename View<TFinder>::Type                TFinderView;
    typedef typename View<TText>::Type                  TTextView;
    typedef typename View<TPattern>::Type               TPatternsView;
    typedef typename View<TDelegate>::Type              TDelegateView;

    // NOTE(esiragusa): Use this to switch to persistent threads.
//    unsigned activeBlocks = cudaMaxActiveBlocks(_findKernel<TTextView, TPatternsView, TSpec, TDelegateView>, ctaSize, 0);
//    setMaxObjects(finder._factory, activeBlocks * ctaSize);

    // Compute grid size.
    unsigned ctaSize = FinderCTASize_<TFinderView>::VALUE;
    unsigned activeBlocks = (length(needle(pattern)) + ctaSize - 1) / ctaSize;

    // Initialize the iterator factory.
    setMaxHistoryLength(finder._factory, length(back(needle(pattern))));
    setMaxObjects(finder._factory, length(needle(pattern)));
    build(finder._factory);

    // Launch the find kernel.
    _findKernel<<<activeBlocks, ctaSize>>>(view(finder), view(pattern), view(delegate));
}
#endif

// ----------------------------------------------------------------------------
// Function find()                                                     [Finder]
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec, typename TDelegate>
inline void
find(Finder2<TText, TPattern, Multiple<TSpec> > & finder, TPattern & pattern, TDelegate & delegate)
{
    typedef Finder2<TText, TPattern,  Multiple<TSpec> > TFinder;

    _find(finder, pattern, delegate, typename ExecSpace<TFinder>::Type());
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

template <typename TText, typename TPattern, typename TSpec>
inline typename View<Finder2<TText, TPattern, Multiple<TSpec> > >::Type
view(Finder2<TText, TPattern, Multiple<TSpec> > & finder)
{
    typename View<Finder2<TText, TPattern, Multiple<TSpec> > >::Type finderView;

    finderView._factory = view(finder._factory);

    return finderView;
}

// ----------------------------------------------------------------------------
// Function textIterator()                                       [Finder Proxy]
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec>
inline SEQAN_HOST_DEVICE typename TextIterator_<TText, Multiple<TSpec> >::Type &
textIterator(Proxy<Finder2<TText, TPattern, Multiple<TSpec> > > & finder)
{
    return finder._textIt;
}

template <typename TText, typename TPattern, typename TSpec>
inline SEQAN_HOST_DEVICE typename TextIterator_<TText, Multiple<TSpec> >::Type const &
textIterator(Proxy<Finder2<TText, TPattern, Multiple<TSpec> > > const & finder)
{
    return finder._textIt;
}

}

#endif  // #ifndef SEQAN_EXTRAS_INDEX_FIND_INDEX_MULTIPLE_H
