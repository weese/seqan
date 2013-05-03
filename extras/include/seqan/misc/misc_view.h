// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_MISC_VIEW_H
#define SEQAN_HEADER_MISC_VIEW_H


// disable CUDA functionality for now...
//#undef SEQAN_FUNC
//#define SEQAN_FUNC inline

#ifdef __CUDACC__
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#endif

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// class View
// ----------------------------------------------------------------------------

template <typename TObject, typename TSpec = void>
class View
{
public:
    typedef typename Iterator<View<TObject, TSpec>, Standard>::Type TIterator;
    typedef typename Reference<TIterator>::Type                     TReference;
    typedef typename GetValue<TIterator>::Type                      TGetValue;

    TIterator _begin;
    TIterator _end;

    SEQAN_FUNC
    View() {}

    template <typename TContainer>
    SEQAN_FUNC
    View(TContainer &cont):
        _begin(begin(cont, Standard())),
        _end(end(cont, Standard())) {}

    template <typename TContainer>
    SEQAN_FUNC
    View(TContainer const &cont):
        _begin(begin(cont, Standard())),
        _end(end(cont, Standard())) {}

    SEQAN_FUNC
    View(TIterator const &begin, TIterator const &end):
        _begin(begin),
        _end(end) {}

    template <typename TContainer>
    SEQAN_FUNC
    View &
    operator= (TContainer &other)
    {
        assign(*this, other);
        return *this;
    }

    template <typename TPos>
    SEQAN_FUNC
    TReference
    operator[] (TPos pos)
    {
        return *(_begin + pos);
    }

    template <typename TPos>
    TGetValue
    operator[] (TPos pos) const
    {
        return getValue(_begin + pos);
    }

};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TObject, typename TSpec>
struct Value<View<TObject, TSpec> >
{
    typedef typename Value<TObject>::Type Type;
};

template <typename TObject, typename TSpec>
struct Value<View<TObject, TSpec> const>
{
    typedef typename Value<TObject>::Type Type;
};

template <typename TObject, typename TSpec>
struct Iterator<View<TObject, TSpec>, Standard>:
    public Iterator<TObject, Standard> {};

template <typename TObject, typename TSpec>
struct Iterator<View<TObject, TSpec> const, Standard>:
    public Iterator<TObject const, Standard> {};

#ifdef __CUDACC__
template <typename TObject, typename TAlloc, typename TSpec>
struct Iterator<View<thrust::device_vector<TObject, TAlloc>, TSpec>, Standard>
{
    typedef typename thrust::device_vector<TObject, TAlloc>::pointer            TIterator_;
    typedef typename thrust::detail::pointer_traits<TIterator_>::raw_pointer    Type;
};
#endif

template <typename TObject, typename TSpec>
struct Difference<View<TObject, TSpec> >
{
    typedef typename Difference<TObject>::Type Type;
};

template <typename TObject, typename TSpec>
struct Size<View<TObject, TSpec> >
{
    typedef typename Difference<TObject>::Type          TDifference;
    typedef typename MakeUnsigned<TDifference>::Type    Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// begin()
// ----------------------------------------------------------------------------

template <typename TObject, typename TSpec>
SEQAN_FUNC
typename Iterator<View<TObject, TSpec>, Standard>::Type
begin(View<TObject, TSpec> & view, Standard)
{
    return view._begin;
}

template <typename TObject, typename TSpec>
SEQAN_FUNC
typename Iterator<View<TObject, TSpec>, Standard>::Type
begin(View<TObject, TSpec> const & view, Standard)
{
    return view._begin;
}

// ----------------------------------------------------------------------------
// end()
// ----------------------------------------------------------------------------

template <typename TObject, typename TSpec>
SEQAN_FUNC
typename Iterator<View<TObject, TSpec>, Standard>::Type
end(View<TObject, TSpec> & view, Standard)
{
    return view._end;
}

template <typename TObject, typename TSpec>
SEQAN_FUNC
typename Iterator<View<TObject, TSpec>, Standard>::Type
end(View<TObject, TSpec> const & view, Standard)
{
    return view._end;
}

// ----------------------------------------------------------------------------
// value()
// ----------------------------------------------------------------------------

template <typename TObject, typename TSpec, typename TPos>
SEQAN_FUNC
typename Reference<View<TObject, TSpec> >::Type
value(View<TObject, TSpec> & view, TPos pos)
{
    return *(view._begin + pos);
}

template <typename TObject, typename TSpec, typename TPos>
SEQAN_FUNC
typename Reference<View<TObject, TSpec> const>::Type
value(View<TObject, TSpec> const & view, TPos pos)
{
    return *(view._begin + pos);
}

// ----------------------------------------------------------------------------
// length()
// ----------------------------------------------------------------------------

template <typename TObject, typename TSpec>
SEQAN_FUNC
typename Difference<View<TObject, TSpec> >::Type
length(View<TObject, TSpec> const & view)
{
    return view._end - view._begin;
}

// ----------------------------------------------------------------------------
// resize()
// ----------------------------------------------------------------------------

// this function doesn't do anything as we are not allowed to change the host (only its elements)
// it is, however, implemented for algorithms that get a sequence to work on
// and need to make sure that it has a certain length

template <typename TObject, typename TSpec, typename TSize, typename TExpand>
inline typename Size< View<TObject, TSpec> >::Type
resize(
    View<TObject, TSpec> & me,
    TSize new_length,
    Tag<TExpand>)
{
    ignoreUnusedVariableWarning(new_length);

    SEQAN_ASSERT_EQ(new_length, (TSize)length(me));
    return length(me);
}

// ----------------------------------------------------------------------------
// assign()
// ----------------------------------------------------------------------------

template <typename TObject, typename TSpec, typename TContainer>
void
assign(View<TObject, TSpec> &view, TContainer const & cont)
{
    view._begin = begin(cont, Standard());
    view._end = end(cont, Standard());
}

// ----------------------------------------------------------------------------
// operator<<
// ----------------------------------------------------------------------------

template <typename TStream, typename TObject, typename TSpec>
inline TStream &
operator<<(TStream & target,
           View<TObject, TSpec> const & source)
{
SEQAN_CHECKPOINT
    write(target, source);
    return target;
}

// ----------------------------------------------------------------------------
// toView()
// ----------------------------------------------------------------------------

template <typename TContainer>
View<TContainer>
toView(TContainer & container)
{
    return View<TContainer>(container);
}

#ifdef __CUDACC__
template <typename TObject, typename TAlloc>
View<thrust::device_vector<TObject, TAlloc> >
toView(thrust::device_vector<TObject, TAlloc> & container)
{
    typedef thrust::device_vector<TObject, TAlloc>    TContainer;
    return View<TContainer>(
            thrust::raw_pointer_cast(container.data()),
            thrust::raw_pointer_cast(&container[container.size()])
    );
}
#endif

// ----------------------------------------------------------------------------
// pipe interface
// ----------------------------------------------------------------------------

//template < typename TObject,
//           typename TSpec,
//           typename TInput,
//           typename TPipeSpec >
//inline void assign(View<TObject, TSpec> &dest, Pipe<TInput, TPipeSpec> &src)
//{
//    typedef typename Iterator<View<TObject, TSpec>, Standard>::Type TDestIter;
//    resize(dest, length(src));
//    beginRead(src);
//    for (TDestIter _cur = begin(dest, Standard()), _end = end(dest, Standard()); _cur != _end; ++_cur, ++src)
//        *_cur = *src;
//    endRead(src);
//}
//
//template < typename TObject,
//           typename TSpec,
//           typename TInput,
//           typename TPipeSpec >
//inline void operator << (View<TObject, TSpec> &dest, Pipe<TInput, TPipeSpec> &src)
//{
//    assign(dest, src);
//}



}  // namespace seqan

#endif  // #ifndef SEQAN_HEADER_MISC_VIEW_H
