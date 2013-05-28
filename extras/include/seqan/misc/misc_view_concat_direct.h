// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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

#ifndef SEQAN_HEADER_MISC_VIEW_CONCAT_DIRECT_H
#define SEQAN_HEADER_MISC_VIEW_CONCAT_DIRECT_H

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class ConcatDirect ContainerView
// ----------------------------------------------------------------------------

template <typename TContainer, typename TSpec>
class ContainerView<TContainer, ConcatDirect<TSpec> >
{
public:
    typedef typename Iterator<ContainerView, Standard>::Type        TIterator;

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
    typename Reference<ContainerView>::Type
    operator[] (TPos pos)
    {
        return value(*this, pos);
    }

    template <typename TPos>
    SEQAN_FUNC
    typename GetValue<ContainerView>::Type
    operator[] (TPos pos) const
    {
        return getValue(*this, pos);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Metafunction GetValue
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Metafunction Iterator
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Metafunction Difference
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Metafunction IsSequence
// ----------------------------------------------------------------------------

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function begin()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TSpec>
SEQAN_FUNC
typename Iterator<ContainerView<TContainer, ConcatDirect<TSpec> >, Standard>::Type
begin(ContainerView<TContainer, ConcatDirect<TSpec> > & view, Standard)
{
    return view._begin;
}

template <typename TContainer, typename TSpec>
SEQAN_FUNC
typename Iterator<ContainerView<TContainer, ConcatDirect<TSpec> > const, Standard>::Type
begin(ContainerView<TContainer, ConcatDirect<TSpec> > const & view, Standard)
{
    return view._begin;
}

// ----------------------------------------------------------------------------
// Function end()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TSpec>
SEQAN_FUNC
typename Iterator<ContainerView<TContainer, ConcatDirect<TSpec> >, Standard>::Type
end(ContainerView<TContainer, ConcatDirect<TSpec> > & view, Standard)
{
    return view._end;
}

template <typename TContainer, typename TSpec>
SEQAN_FUNC
typename Iterator<ContainerView<TContainer, ConcatDirect<TSpec> > const, Standard>::Type
end(ContainerView<TContainer, ConcatDirect<TSpec> > const & view, Standard)
{
    return view._end;
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TSpec, typename TPos>
SEQAN_FUNC
typename Reference<ContainerView<TContainer, ConcatDirect<TSpec> > >::Type
value(ContainerView<TContainer, ConcatDirect<TSpec> > & view, TPos pos)
{
    return *(view._begin + pos);
}

template <typename TContainer, typename TSpec, typename TPos>
SEQAN_FUNC
typename Reference<ContainerView<TContainer, ConcatDirect<TSpec> > const>::Type
value(ContainerView<TContainer, ConcatDirect<TSpec> > const & view, TPos pos)
{
    return *(view._begin + pos);
}

// ----------------------------------------------------------------------------
// Function getValue()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TSpec, typename TPos>
SEQAN_FUNC
typename GetValue<ContainerView<TContainer, ConcatDirect<TSpec> > >::Type
getValue(ContainerView<TContainer, ConcatDirect<TSpec> > & view, TPos pos)
{
    return getValue(view._begin + pos);
}

template <typename TContainer, typename TSpec, typename TPos>
SEQAN_FUNC
typename GetValue<ContainerView<TContainer, ConcatDirect<TSpec> > const>::Type
getValue(ContainerView<TContainer, ConcatDirect<TSpec> > const & view, TPos pos)
{
    return getValue(view._begin + pos);
}

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TSpec>
SEQAN_FUNC
typename Difference<ContainerView<TContainer, ConcatDirect<TSpec> > >::Type
length(ContainerView<TContainer, ConcatDirect<TSpec> > const & view)
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
inline typename Size< ContainerView<TContainer, ConcatDirect<TSpec> > >::Type
resize(ContainerView<TContainer, ConcatDirect<TSpec> > & me, TSize new_length, TValue /* val */, Tag<TExpand>)
{
    ignoreUnusedVariableWarning(new_length);

    SEQAN_ASSERT_EQ(new_length, (TSize)length(me));
    return length(me);
}

template <typename TContainer, typename TSpec, typename TSize, typename TExpand>
inline typename Size< ContainerView<TContainer, ConcatDirect<TSpec> > >::Type
resize(ContainerView<TContainer, ConcatDirect<TSpec> > & me, TSize new_length, Tag<TExpand> tag)
{
    return resize(me, new_length, Nothing(), tag);
}

// ----------------------------------------------------------------------------
// Function assign()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TSpec, typename TOtherContainer>
void
assign(ContainerView<TContainer, ConcatDirect<TSpec> > & view, TOtherContainer const & cont)
{
    view._begin = begin(cont, Standard());
    view._end = end(cont, Standard());
}

// ----------------------------------------------------------------------------
// Function view()
// ----------------------------------------------------------------------------

//template <typename TContainer>
//ContainerView<TContainer>
//view(TContainer & container)
//{
//    return ContainerView<TContainer>(container);
//}
//
//#ifdef __CUDACC__
//template <typename TContainer, typename TAlloc>
//ContainerView<thrust::device_vector<TContainer, TAlloc> >
//view(thrust::device_vector<TContainer, TAlloc> & container)
//{
//    typedef thrust::device_vector<TContainer, TAlloc>    TDeviceContainer;
//    return ContainerView<TDeviceContainer>(
//            thrust::raw_pointer_cast(container.data()),
//            thrust::raw_pointer_cast(&container[container.size()])
//    );
//}
//#endif

}  // namespace seqan

#endif  // #ifndef SEQAN_HEADER_MISC_VIEW_CONCAT_DIRECT_H
