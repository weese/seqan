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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_MISC_VIEW_H
#define SEQAN_HEADER_MISC_VIEW_H

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class ContainerView
// ----------------------------------------------------------------------------

template <typename TContainer, typename TSpec = void>
class ContainerView
{
public:
    typedef typename Iterator<ContainerView, Standard>::Type        TIterator;
    typedef typename Reference<TIterator>::Type                     TReference;
    typedef typename GetValue<TIterator>::Type                      TGetValue;

    TIterator _begin;
    TIterator _end;

    // ------------------------------------------------------------------------
    // ContainerView Constructors
    // ------------------------------------------------------------------------

    SEQAN_FUNC
    ContainerView() {}

    template <typename TOtherContainer>
    SEQAN_FUNC
    ContainerView(TOtherContainer & cont):
        _begin(begin(cont, Standard())),
        _end(end(cont, Standard())) {}

    template <typename TOtherContainer>
    SEQAN_FUNC
    ContainerView(TOtherContainer const & cont):
        _begin(begin(cont, Standard())),
        _end(end(cont, Standard())) {}

    SEQAN_FUNC
    ContainerView(TIterator const & begin, TIterator const & end):
        _begin(begin),
        _end(end) {}

    // ------------------------------------------------------------------------
    // Operator =
    // ------------------------------------------------------------------------

    template <typename TOtherContainer>
    SEQAN_FUNC
    ContainerView &
    operator= (TOtherContainer & other)
    {
        assign(*this, other);
        return *this;
    }

    // ------------------------------------------------------------------------
    // Operator []
    // ------------------------------------------------------------------------

    template <typename TPos>
    SEQAN_FUNC
    TReference
    operator[] (TPos pos)
    {
        return value(*this, pos);
    }

    template <typename TPos>
    TGetValue
    operator[] (TPos pos) const
    {
        // TODO(esiragusa): There should be a getValue(view, pos)
        return getValue(_begin + pos);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction View
// ----------------------------------------------------------------------------
// TODO(esiragusa): Move generic View metafunction somewhere else.

template <typename TObject>
struct View
{
    typedef typename If<typename IsSimple<TObject>::Type, TObject, ContainerView<TObject> >::Type    Type;
};

// ----------------------------------------------------------------------------
// Metafunction Device
// ----------------------------------------------------------------------------
// TODO(esiragusa): Move generic Device metafunction somewhere else.

#ifdef __CUDACC__
template <typename TObject>
struct Device
{
    typedef typename Value<TObject>::Type               TValue_;
//    typedef typename DefaultDeviceAlloc<TObject>::Type  TAlloc_;
    typedef thrust::device_vector<TValue_/*, TAlloc_*/>     TDevice_;

    typedef typename If<typename IsSimple<TObject>::Type, TObject, TDevice_>::Type      Type;
};
#endif

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TContainer, typename TSpec>
struct Value<ContainerView<TContainer, TSpec> >
{
    typedef typename Value<TContainer>::Type Type;
};

template <typename TContainer, typename TSpec>
struct Value<ContainerView<TContainer, TSpec> const> :
    public Value<ContainerView<TContainer const, TSpec> > {};

// ----------------------------------------------------------------------------
// Metafunction Iterator
// ----------------------------------------------------------------------------

template <typename TContainer, typename TSpec>
struct Iterator<ContainerView<TContainer, TSpec>, Standard>:
    public Iterator<TContainer, Standard> {};

template <typename TContainer, typename TSpec>
struct Iterator<ContainerView<TContainer, TSpec> const, Standard>:
    public Iterator<TContainer const, Standard> {};

#ifdef __CUDACC__
template <typename TContainer, typename TAlloc, typename TSpec>
struct Iterator<ContainerView<thrust::device_vector<TContainer, TAlloc>, TSpec>, Standard>
{
    typedef typename thrust::device_vector<TContainer, TAlloc>::pointer            TIterator_;
    typedef typename thrust::detail::pointer_traits<TIterator_>::raw_pointer    Type;
};

template <typename TContainer, typename TAlloc, typename TSpec>
struct Iterator<ContainerView<thrust::device_vector<TContainer, TAlloc>, TSpec> const, Standard>
{
    typedef typename thrust::device_vector<TContainer const, TAlloc>::pointer      TIterator_;
    typedef typename thrust::detail::pointer_traits<TIterator_>::raw_pointer    Type;
};
#endif

// ----------------------------------------------------------------------------
// Metafunction Difference
// ----------------------------------------------------------------------------

template <typename TContainer, typename TSpec>
struct Difference<ContainerView<TContainer, TSpec> >
{
    typedef typename Difference<TContainer>::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <typename TContainer, typename TSpec>
struct Size<ContainerView<TContainer, TSpec> >
{
    typedef typename Difference<TContainer>::Type          TDifference;
    typedef typename MakeUnsigned<TDifference>::Type    Type;
};

// ----------------------------------------------------------------------------
// Metafunction IsSequence
// ----------------------------------------------------------------------------

template <typename TContainer, typename TSpec>
struct IsSequence<ContainerView<TContainer, TSpec> >:
    public IsSequence<TContainer> {};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function begin()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TSpec>
SEQAN_FUNC
typename Iterator<ContainerView<TContainer, TSpec>, Standard>::Type
begin(ContainerView<TContainer, TSpec> & view, Standard)
{
    return view._begin;
}

template <typename TContainer, typename TSpec>
SEQAN_FUNC
typename Iterator<ContainerView<TContainer, TSpec> const, Standard>::Type
begin(ContainerView<TContainer, TSpec> const & view, Standard)
{
    return view._begin;
}

// ----------------------------------------------------------------------------
// Function end()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TSpec>
SEQAN_FUNC
typename Iterator<ContainerView<TContainer, TSpec>, Standard>::Type
end(ContainerView<TContainer, TSpec> & view, Standard)
{
    return view._end;
}

template <typename TContainer, typename TSpec>
SEQAN_FUNC
typename Iterator<ContainerView<TContainer, TSpec> const, Standard>::Type
end(ContainerView<TContainer, TSpec> const & view, Standard)
{
    return view._end;
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TSpec, typename TPos>
SEQAN_FUNC
typename Reference<ContainerView<TContainer, TSpec> >::Type
value(ContainerView<TContainer, TSpec> & view, TPos pos)
{
    return *(view._begin + pos);
}

template <typename TContainer, typename TSpec, typename TPos>
SEQAN_FUNC
typename Reference<ContainerView<TContainer, TSpec> const>::Type
value(ContainerView<TContainer, TSpec> const & view, TPos pos)
{
    return *(view._begin + pos);
}

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TSpec>
SEQAN_FUNC
typename Difference<ContainerView<TContainer, TSpec> >::Type
length(ContainerView<TContainer, TSpec> const & view)
{
    return view._end - view._begin;
}

// ----------------------------------------------------------------------------
// Function resize()
// ----------------------------------------------------------------------------

// this function doesn't do anything as we are not allowed to change the host (only its elements)
// it is, however, implemented for algorithms that get a sequence to work on
// and need to make sure that it has a certain length

template <typename TContainer, typename TSpec, typename TSize, typename TValue, typename TExpand>
inline typename Size< ContainerView<TContainer, TSpec> >::Type
resize(ContainerView<TContainer, TSpec> & me, TSize new_length, TValue /* val */, Tag<TExpand>)
{
    ignoreUnusedVariableWarning(new_length);

    SEQAN_ASSERT_EQ(new_length, (TSize)length(me));
    return length(me);
}

template <typename TContainer, typename TSpec, typename TSize, typename TExpand>
inline typename Size< ContainerView<TContainer, TSpec> >::Type
resize(ContainerView<TContainer, TSpec> & me, TSize new_length, Tag<TExpand> tag)
{
    return resize(me, new_length, Nothing(), tag);
}

// ----------------------------------------------------------------------------
// Function assign()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TSpec, typename TOtherContainer>
void
assign(ContainerView<TContainer, TSpec> & view, TOtherContainer const & cont)
{
    view._begin = begin(cont, Standard());
    view._end = end(cont, Standard());
}

// ----------------------------------------------------------------------------
// Function view()
// ----------------------------------------------------------------------------

template <typename TContainer>
ContainerView<TContainer>
view(TContainer & container)
{
    return ContainerView<TContainer>(container);
}

#ifdef __CUDACC__
template <typename TContainer, typename TAlloc>
ContainerView<thrust::device_vector<TContainer, TAlloc> >
view(thrust::device_vector<TContainer, TAlloc> & container)
{
    typedef thrust::device_vector<TContainer, TAlloc>    TDeviceContainer;
    return ContainerView<TDeviceContainer>(
            thrust::raw_pointer_cast(container.data()),
            thrust::raw_pointer_cast(&container[container.size()])
    );
}
#endif

// ----------------------------------------------------------------------------
// Operator<<
// ----------------------------------------------------------------------------

template <typename TStream, typename TContainer, typename TSpec>
inline TStream &
operator<<(TStream & target, ContainerView<TContainer, TSpec> const & source)
{
    write(target, source);
    return target;
}

// ----------------------------------------------------------------------------
// Pipe interface
// ----------------------------------------------------------------------------

//template < typename TContainer,
//           typename TSpec,
//           typename TInput,
//           typename TPipeSpec >
//inline void assign(ContainerView<TContainer, TSpec> &dest, Pipe<TInput, TPipeSpec> &src)
//{
//    typedef typename Iterator<ContainerView<TContainer, TSpec>, Standard>::Type TDestIter;
//    resize(dest, length(src));
//    beginRead(src);
//    for (TDestIter _cur = begin(dest, Standard()), _end = end(dest, Standard()); _cur != _end; ++_cur, ++src)
//        *_cur = *src;
//    endRead(src);
//}
//
//template < typename TContainer,
//           typename TSpec,
//           typename TInput,
//           typename TPipeSpec >
//inline void operator << (ContainerView<TContainer, TSpec> &dest, Pipe<TInput, TPipeSpec> &src)
//{
//    assign(dest, src);
//}

}  // namespace seqan

#endif  // #ifndef SEQAN_HEADER_MISC_VIEW_H
