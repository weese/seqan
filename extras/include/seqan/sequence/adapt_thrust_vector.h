// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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
// Author: Knut Reinert <knut.reinert@fu-berlin.de>
// ==========================================================================
// Adaptions for STL vectors to SeqAn strings.
// ==========================================================================

// SEQAN_NO_GENERATED_FORWARDS: No forwards are generated for this file.

#ifndef SEQAN_SEQUENCE_ADAPT_THRUST_VECTOR_H_
#define SEQAN_SEQUENCE_ADAPT_THRUST_VECTOR_H_

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

/**
.Adaption."thrust::device_vector"
..summary:Adaption for STL vector class.
 */

// ===========================================================================
// Metafunctions
// ===========================================================================

///.Metafunction.IsContiguous.param.T.type:Adaption.thrust::device_vector
///.Metafunction.IsContiguous.class:Adaption.thrust::device_vector

template <typename TChar, typename TAlloc>
struct IsContiguous< thrust::device_vector<TChar, TAlloc> >
{
    enum { VALUE = true };
};

template <typename  TChar, typename TAlloc>
struct IsContiguous< thrust::device_vector<TChar, TAlloc> const>
        : IsContiguous< thrust::device_vector<TChar, TAlloc> > {};

///.Metafunction.Value.param.T.type:Adaption.thrust::device_vector
///.Metafunction.Value.class:Adaption.thrust::device_vector

template <typename TChar, typename TAlloc>
struct Value< thrust::device_vector<TChar, TAlloc> >
{
    typedef typename thrust::device_vector<TChar, TAlloc>::value_type Type;
};

template <typename TChar, typename TAlloc>
struct Value< thrust::device_vector<TChar, TAlloc> const>
        : Value< thrust::device_vector<TChar, TAlloc> > {};

///.Metafunction.GetValue.param.T.type:Adaption.thrust::device_vector
// TODO(holtgrew): GetValue is a reference?! I thought the reverse was true in respect to Value<>.
///.Metafunction.GetValue.class:Adaption.thrust::device_vector

template <typename TChar, typename TAlloc>
struct GetValue< thrust::device_vector<TChar, TAlloc> >
{
    typedef typename thrust::device_vector<TChar, TAlloc>::reference Type;
};

template <typename TChar, typename TAlloc>
struct GetValue< thrust::device_vector<TChar,  TAlloc> const>
{
    typedef typename thrust::device_vector<TChar, TAlloc>::const_reference Type;
};

///.Metafunction.Reference.param.T.type:Adaption.thrust::device_vector
///.Metafunction.Reference.class:Adaption.thrust::device_vector

template <typename TChar, typename TAlloc>
struct Reference< thrust::device_vector<TChar, TAlloc> >
{
    typedef typename thrust::device_vector<TChar, TAlloc>::reference Type;
};

template <typename TChar,  typename TAlloc>
struct Reference< thrust::device_vector<TChar, TAlloc> const>
{
    typedef typename thrust::device_vector<TChar,  TAlloc>::const_reference Type;
};

///.Metafunction.Iterator.param.T.type:Adaption.thrust::device_vector
///.Metafunction.Iterator.class:Adaption.thrust::device_vector

template <typename TChar, typename TAlloc>
struct Iterator< thrust::device_vector<TChar, TAlloc>, Rooted>
{
    typedef thrust::device_vector<TChar, TAlloc> TVector_;
    typedef Iter<TVector_, StdIteratorAdaptor> TIterator_;
    typedef Iter<TVector_, AdaptorIterator<TIterator_> > Type;
};

template <typename TChar, typename TAlloc>
struct Iterator< thrust::device_vector<TChar, TAlloc> const, Rooted>
{
    typedef thrust::device_vector<TChar, TAlloc> const TVector_;
    typedef Iter<TVector_, StdIteratorAdaptor> TIterator_;
    typedef Iter<TVector_, AdaptorIterator<TIterator_> > Type;
};

template <typename TChar,  typename TAlloc>
struct Iterator< thrust::device_vector<TChar, TAlloc>, Standard >
{
	typedef typename thrust::device_vector<TChar,  TAlloc>::iterator Type;
    // typedef Iter< thrust::device_vector<TChar,  TAlloc>, StdIteratorAdaptor > Type;
};

template <typename TChar,  typename TAlloc>
struct Iterator< thrust::device_vector<TChar,  TAlloc> const, Standard>
{
	typedef typename thrust::device_vector<TChar,  TAlloc>::const_iterator Type;
    // typedef Iter< thrust::device_vector<TChar, TAlloc> const, StdIteratorAdaptor > Type;
};

template <typename TValue>
struct Reference< thrust::detail::normal_iterator<thrust::device_ptr<TValue> > > 
{
	typedef thrust::device_reference<TValue> Type;
};

template <typename TValue>
struct Reference< thrust::detail::normal_iterator<thrust::device_ptr<TValue> > const> 
{
	typedef thrust::device_reference<TValue> const Type;
};

///.Metafunction.Position.param.T.type:Adaption.thrust::device_vector
///.Metafunction.Position.class:Adaption.thrust::device_vector

template <typename TChar,  typename TAlloc>
struct Position< thrust::device_vector<TChar, TAlloc> >
{
    typedef typename thrust::device_vector<TChar,  TAlloc>::size_type Type;
};

template <typename TChar,  typename TAlloc>
struct Position< thrust::device_vector<TChar,  TAlloc> const>
        : Position< thrust::device_vector<TChar,  TAlloc> > {};

///.Metafunction.Position.param.T.type:Adaption.thrust::device_vector
///.Metafunction.Position.class:Adaption.thrust::device_vector

template <typename TChar,  typename TAlloc>
struct Size< thrust::device_vector<TChar, TAlloc> >
{
    typedef typename thrust::device_vector<TChar, TAlloc>::size_type Type;
};

template <typename TChar, typename TAlloc>
struct Size< thrust::device_vector<TChar, TAlloc> const>
        : Size< thrust::device_vector<TChar, TAlloc> > {};

///.Metafunction.Size.param.T.type:Adaption.thrust::device_vector
///.Metafunction.Size.class:Adaption.thrust::device_vector

template <typename TChar, typename TAlloc>
struct DefaultOverflowImplicit< thrust::device_vector<TChar, TAlloc> >
{
    typedef Generous Type;
};

///.Metafunction.StdContainerIterator.param.T.type:Adaption.thrust::device_vector
///.Metafunction.StdContainerIterator.class:Adaption.thrust::device_vector

template <typename TChar, typename TAlloc>
struct StdContainerIterator< thrust::device_vector<TChar, TAlloc> >
{
    typedef thrust::device_vector<TChar, TAlloc> TContainer_;
    typedef typename TContainer_::iterator Type;
};

template <typename TChar, typename TAlloc>
struct StdContainerIterator< thrust::device_vector<TChar, TAlloc> const>
{
    typedef thrust::device_vector<TChar, TAlloc> TContainer_;
    typedef typename TContainer_::const_iterator Type;
};

// ===========================================================================
// Functions
// ===========================================================================

///.Function.getObjectId.param.object.type:Adaption.thrust::device_vector
///.Function.getObjectId.class:Adaption.thrust::device_vector

template <typename TChar,  typename TAlloc>
SEQAN_FUNC void const *
getObjectId(thrust::device_vector<TChar, TAlloc> const & me)
{
    SEQAN_CHECKPOINT;
    if (me.empty())
        return NULL;
    else
        return (& *(me.end() - 1)) + 1;
}

///.Function.begin.param.object.type:Adaption.thrust::device_vector
///.Function.begin.class:Adaption.thrust::device_vector

template <typename TChar,  typename TAlloc>
SEQAN_FUNC typename Iterator< thrust::device_vector<TChar,  TAlloc>, Standard>::Type
begin(thrust::device_vector<TChar,  TAlloc> & me,
      Standard)
{
    SEQAN_CHECKPOINT;
    return typename Iterator< thrust::device_vector<TChar,  TAlloc>, Standard>::Type(me.begin());
}
template <typename TChar,  typename TAlloc>
SEQAN_FUNC typename Iterator< thrust::device_vector<TChar,  TAlloc> const, Standard>::Type
begin(thrust::device_vector<TChar, TAlloc> const & me,
      Standard)
{
    SEQAN_CHECKPOINT;
    return typename Iterator< thrust::device_vector<TChar,  TAlloc> const, Standard>::Type(me.begin());
}

///.Function.end.param.object.type:Adaption.thrust::device_vector
///.Function.end.class:Adaption.thrust::device_vector

template <typename TChar, typename TAlloc>
SEQAN_FUNC typename Iterator< thrust::device_vector<TChar, TAlloc>, Standard>::Type
end(thrust::device_vector<TChar,  TAlloc> & me,
    Standard)
{
    SEQAN_CHECKPOINT;
    return typename Iterator< thrust::device_vector<TChar, TAlloc>, Standard>::Type(me.end());
}
template <typename TChar,  typename TAlloc>
SEQAN_FUNC typename Iterator< thrust::device_vector<TChar,  TAlloc> const, Standard>::Type
end(thrust::device_vector<TChar,  TAlloc> const & me,
    Standard)
{
    SEQAN_CHECKPOINT;
    return typename Iterator< thrust::device_vector<TChar,  TAlloc> const, Standard>::Type(me.end());
}

///.Function.value.param.container.type:Adaption.thrust::device_vector
///.Function.value.class:Adaption.thrust::device_vector

template <typename TChar,  typename TAlloc, typename TPos>
SEQAN_FUNC typename GetValue< thrust::device_vector<TChar, TAlloc> >::Type
value(thrust::device_vector<TChar,  TAlloc> & me,
      TPos pos)
{
    SEQAN_CHECKPOINT;
    return me[pos];
}
template <typename TChar,  typename TAlloc, typename TPos>
SEQAN_FUNC typename GetValue< thrust::device_vector<TChar,  TAlloc> const>::Type
value(thrust::device_vector<TChar, TAlloc> const & me,
      TPos pos)
{
    SEQAN_CHECKPOINT;
    return me[pos];
}

///.Function.value.param.container.type:Adaption.thrust::device_vector
///.Function.value.class:Adaption.thrust::device_vector

template <typename TChar, typename TAlloc>
SEQAN_FUNC typename Size< thrust::device_vector<TChar, TAlloc> >::Type
length(thrust::device_vector<TChar, TAlloc> const & me)
{
    SEQAN_CHECKPOINT;
    return me.size();
}

///.Function.capacity.param.object.type:Adaption.thrust::device_vector
///.Function.capacity.class:Adaption.thrust::device_vector

template <typename TChar, typename TAlloc>
SEQAN_FUNC typename Size< thrust::device_vector<TChar, TAlloc> >::Type
capacity(thrust::device_vector<TChar, TAlloc> const & me)
{
    SEQAN_CHECKPOINT;
    return me.capacity();
}

///.Function.empty.param.object.type:Adaption.thrust::device_vector
///.Function.empty.class:Adaption.thrust::device_vector

template <typename TChar, typename TAlloc>
SEQAN_FUNC bool
empty(thrust::device_vector<TChar, TAlloc> const & me)
{
    SEQAN_CHECKPOINT;
    return me.empty();
}

///.Function.clear.param.object.type:Adaption.thrust::device_vector
///.Function.clear.class:Adaption.thrust::device_vector

template <typename TChar,  typename TAlloc>
SEQAN_FUNC void
clear(thrust::device_vector<TChar, TAlloc> & me)
{
    SEQAN_CHECKPOINT;
    me.clear();
}

///.Function.front.param.container.type:Adaption.thrust::device_vector
///.Function.front.class:Adaption.thrust::device_vector

template <typename TChar>
SEQAN_FUNC typename Reference<thrust::device_vector<TChar> >::Type
front(thrust::device_vector<TChar> & list)
{
    SEQAN_CHECKPOINT;
    return list.front();
}

template <typename TChar>
SEQAN_FUNC typename Reference<thrust::device_vector<TChar> const>::Type
front(thrust::device_vector<TChar> const & list)
{
    SEQAN_CHECKPOINT;
    return list.front();
}

///.Function.back.param.container.type:Adaption.thrust::device_vector
///.Function.back.class:Adaption.thrust::device_vector

template <typename TChar>
SEQAN_FUNC typename Reference<thrust::device_vector<TChar> >::Type
back(thrust::device_vector<TChar> & list)
{
    SEQAN_CHECKPOINT;
    return list.back();
}

template <typename TChar>
SEQAN_FUNC typename Reference<thrust::device_vector<TChar> const>::Type
back(thrust::device_vector<TChar> const & list)
{
    SEQAN_CHECKPOINT;
    return list.back();
}

//////////////////////////////////////////////////////////////////////////////
//assign to thrust::device_vector

///.Function.assign.param.target.type:Adaption.thrust::device_vector
///.Function.assign.class:Adaption.thrust::device_vector

template <typename TChar,  typename TAlloc, typename TSource>
SEQAN_FUNC void
assign(thrust::device_vector<TChar,  TAlloc> & target,
       TSource & source)
{
    SEQAN_CHECKPOINT;
    assign(target, source, Generous());
}
template <typename TChar,  typename TAlloc, typename TSource>
SEQAN_FUNC void
assign(thrust::device_vector<TChar,  TAlloc> & target,
       TSource const & source)
{
    SEQAN_CHECKPOINT;
    assign(target, source, Generous());
}

template <typename TChar,  typename TAlloc, typename TSource, typename TSize>
SEQAN_FUNC void
assign(thrust::device_vector<TChar,  TAlloc> & target,
       TSource & source,
       TSize limit)
{
    SEQAN_CHECKPOINT;
    assign(target, source, limit, Generous());
}
template <typename TChar,  typename TAlloc, typename TSource, typename TSize>
SEQAN_FUNC void
assign(thrust::device_vector<TChar,  TAlloc> & target,
       TSource const & source,
       TSize limit)
{
    SEQAN_CHECKPOINT;
    assign(target, source, limit, Generous());
}

//____________________________________________________________________________

template <typename TChar,  typename TAlloc, typename TSource>
SEQAN_FUNC void
assign(thrust::device_vector<TChar, TAlloc> & target,
       TSource & source,
       Generous)
{
    SEQAN_CHECKPOINT;
    target.assign(begin(source, Standard()), end(source, Standard()));
}
template <typename TChar, typename TAlloc, typename TSource>
SEQAN_FUNC void
assign(thrust::device_vector<TChar, TAlloc> & target,
       TSource const & source,
       Generous)
{
    SEQAN_CHECKPOINT;
    target.assign(begin(source, Standard()), end(source, Standard()));
}


template <typename TChar,  typename TAlloc, typename TSource>
SEQAN_FUNC void
assign_std_vector_Generous_impl(thrust::device_vector<TChar,  TAlloc> & target,
                                TSource & source,
                                typename Size< thrust::device_vector<TChar,  TAlloc> >::Type limit)
{
    SEQAN_CHECKPOINT;
    typename Iterator<TSource const, Standard>::Type source_begin = begin(source, Standard());
    typename Size<TSource const>::Type source_length = length(source);
    if (source_length > limit)
    {
        source_length = limit;
    }
    target.assign(source_begin, source_begin + source_length);
}
template <typename TChar,  typename TAlloc, typename TSource>
SEQAN_FUNC void
assign(thrust::device_vector<TChar,  TAlloc> & target,
       TSource & source,
       typename Size< thrust::device_vector<TChar, TAlloc> >::Type limit,
       Generous)
{
    SEQAN_CHECKPOINT;
    assign_std_vector_Generous_impl(target, source, limit);
}
template <typename TChar,  typename TAlloc, typename TSource>
SEQAN_FUNC void
assign(thrust::device_vector<TChar, TAlloc> & target,
       TSource const & source,
       typename Size< thrust::device_vector<TChar, TAlloc> >::Type limit,
       Generous)
{
    SEQAN_CHECKPOINT;
    assign_std_vector_Generous_impl(target, source, limit);
}

//____________________________________________________________________________

template <typename TChar, typename TAlloc, typename TSource>
SEQAN_FUNC void
assign(thrust::device_vector<TChar,  TAlloc> & target,
       TSource & source,
       Limit)
{
    SEQAN_CHECKPOINT;
    assign(target, source, target.capacity(), Generous());
}
template <typename TChar, typename TAlloc, typename TSource>
SEQAN_FUNC void
assign(thrust::device_vector<TChar,  TAlloc> & target,
       TSource const & source,
       Limit)
{
    SEQAN_CHECKPOINT;
    assign(target, source, target.capacity(), Generous());
}

template <typename TChar, typename TAlloc, typename TSource>
SEQAN_FUNC void
assign(thrust::device_vector<TChar, TAlloc> & target,
       TSource & source,
       typename Size< thrust::device_vector<TChar, TAlloc> >::Type limit,
       Limit)
{
    SEQAN_CHECKPOINT;
    if (limit > target.capacity())
    {
        limit = target.capacity();
    }

    assign(target, source, limit, Generous());
}
template <typename TChar, typename TAlloc, typename TSource>
SEQAN_FUNC void
assign(thrust::device_vector<TChar,  TAlloc> & target,
       TSource const & source,
       typename Size< thrust::device_vector<TChar,  TAlloc> >::Type limit,
       Limit)
{
    SEQAN_CHECKPOINT;
    if (limit > target.capacity())
    {
        limit = target.capacity();
    }

    assign(target, source, limit, Generous());
}

//////////////////////////////////////////////////////////////////////////////
//append to thrust::device_vector

///.Function.append.param.target.type:Adaption.thrust::device_vector
///.Function.append.class:Adaption.thrust::device_vector

template <typename TChar, typename TAlloc, typename TSource>
SEQAN_FUNC void
append(thrust::device_vector<TChar,  TAlloc> & target,
       TSource const & source,
       Generous)
{
    SEQAN_CHECKPOINT;
    target.insert(target.end(), begin(source, Standard()), end(source, Standard()));
}

template <typename TChar,  typename TAlloc, typename TSource>
SEQAN_FUNC void
append(thrust::device_vector<TChar, TAlloc> & target,
       TSource const & source,
       typename Size< thrust::device_vector<TChar, TAlloc> >::Type limit,
       Generous)
{
    SEQAN_CHECKPOINT;
    typename Size< thrust::device_vector<TChar, TAlloc> >::Type target_length = target.length();
    if (target_length > limit)
    {
        target.resize(limit);
    }
    else
    {
        limit -= target_length;
        typename Iterator<TSource const, Standard>::Type source_begin = begin(source, Standard());
        typename Size<TSource const>::Type source_length = length(source);
        if (source_length > limit)
        {
            source_length = limit;
        }

        target.insert(target.end(), source_begin, source_begin + source_length);
    }
}

//____________________________________________________________________________

template <typename TChar, typename TAlloc, typename TSource>
SEQAN_FUNC void
append(thrust::device_vector<TChar,  TAlloc> & target,
       TSource const & source,
       Limit)
{
    SEQAN_CHECKPOINT;
    append(target, source, target.capacity(), Generous());
}

template <typename TChar, typename TAlloc, typename TSource>
SEQAN_FUNC void
append(thrust::device_vector<TChar, TAlloc> & target,
       TSource const & source,
       typename Size< thrust::device_vector<TChar, TAlloc> >::Type limit,
       Limit)
{
    SEQAN_CHECKPOINT;
    if (limit > target.capacity())
    {
        limit = target.capacity();
    }

    append(target, source, limit, Generous());
}

//////////////////////////////////////////////////////////////////////////////
///.Function.appendValue.param.target.type:Adaption.thrust::device_vector
///.Function.appendValue.class:Adaption.thrust::device_vector

template <typename TChar, typename TAlloc, typename TValue, typename TTag>
SEQAN_FUNC void
appendValue(thrust::device_vector<TChar, TAlloc> & me,
            TValue const & _value,
            TTag)
{
    SEQAN_CHECKPOINT;
    me.push_back(_value);
}

template <typename TChar, typename TAlloc, typename TValue>
SEQAN_FUNC void
appendValue(thrust::device_vector<TChar,  TAlloc> & me,
            TValue const & _value,
            Limit)
{
    SEQAN_CHECKPOINT;
    if (capacity(me) > length(me)) me.push_back(_value);
}

//////////////////////////////////////////////////////////////////////////////
//replace to thrust::device_vector

///.Function.replace.param.target.type:Adaption.thrust::device_vector
///.Function.replace.param.source.type:Adaption.thrust::device_vector
///.Function.replace.class:Adaption.thrust::device_vector

template <typename TChar,  typename TAlloc, typename TSource>
SEQAN_FUNC void
replace(thrust::device_vector<TChar, TAlloc> & target,
        typename Position< thrust::device_vector<TChar, TAlloc> >::Type pos_begin,
        typename Position< thrust::device_vector<TChar, TAlloc> >::Type pos_end,
        TSource const & source,
        Generous)
{
    SEQAN_CHECKPOINT;
    typename Size< thrust::device_vector<TChar, TAlloc> >::Type target_size = pos_end-pos_begin;
    typename Size< thrust::device_vector<TChar, TAlloc> >::Type source_size =length(source);

    if(target_size >= source_size)
        {
            copy(source.begin(),source.end(), target.begin()+pos_begin);
            if( target_size > source_size )
                target.erase(target.begin()+pos_begin+source_size,target.begin()+pos_end);
        }
    else
        {
            copy(source.begin(),source.begin()+target_size,target.begin()+pos_begin);
            target.insert(target.begin()+pos_end,source.begin()+target_size,source.end());
        }
}

template <typename TChar, typename TAlloc, typename TSource>
SEQAN_FUNC void
replace(thrust::device_vector<TChar, TAlloc> & target,
        typename Position< thrust::device_vector<TChar, TAlloc> >::Type pos_begin,
        typename Position< thrust::device_vector<TChar, TAlloc> >::Type pos_end,
        TSource const & source,
        typename Size< thrust::device_vector<TChar, TAlloc> >::Type limit,
        Generous)
{
    SEQAN_CHECKPOINT;
    if (pos_begin >= limit)
    {
        target.resize(limit);
    }
    else
    {
        typename Size<TSource const>::Type source_length = length(source);
        typename Size< thrust::device_vector<TChar, TAlloc> >::Type pos_mid = pos_begin + source_length;
        typename Size< thrust::device_vector<TChar, TAlloc> >::Type pos_limit(limit);
        if (pos_mid > limit)
        {
            target.resize(limit);
            replace(target,pos_begin,pos_limit,source);
            target.resize(limit);
        }
        else
        {
            replace(target,pos_begin,pos_end,source);
            if (target.size() > limit)
            {
                target.resize(limit);
            }
        }
    }

}

template <typename TChar,  typename TAlloc, typename TSource>
SEQAN_FUNC void
replace(thrust::device_vector<TChar,  TAlloc> & target,
        typename Position< thrust::device_vector<TChar, TAlloc> >::Type pos_begin,
        typename Position< thrust::device_vector<TChar, TAlloc> >::Type pos_end,
        TSource const & source,
        Limit)
{
    SEQAN_CHECKPOINT;
    replace(target, pos_begin, pos_end, source, target.capacity(), Generous());
}

template <typename TChar, typename TAlloc, typename TSource>
SEQAN_FUNC void
replace(thrust::device_vector<TChar, TAlloc> & target,
        typename Position< thrust::device_vector<TChar,  TAlloc> >::Type pos_begin,
        typename Position< thrust::device_vector<TChar,  TAlloc> >::Type pos_end,
        TSource const & source,
        typename Size< thrust::device_vector<TChar, TAlloc> >::Type limit,
        Limit)
{
    SEQAN_CHECKPOINT;
    if (limit > target.capacity())
    {
        limit = target.capacity();
    }
    replace(target, pos_begin, pos_end, source, limit, Generous());
}


//////////////////////////////////////////////////////////////////////////////
// handling of iterators as begin and end

template<typename TChar, typename TCharTraits, typename TAlloc, typename TSource, typename TExpand>
SEQAN_FUNC void
replace(thrust::device_vector<TChar, TAlloc> & target,
        typename Iterator< thrust::device_vector<TChar, TAlloc>, Rooted>::Type pos_begin,
        typename Iterator< thrust::device_vector<TChar, TAlloc>, Rooted>::Type pos_end,
        TSource & source,
        Tag<TExpand> const tag)
{
    replace(target, position(pos_begin), position(pos_end), source, tag);
}

/*
template<typename TChar, typename TAlloc, typename TSource, typename TExpand>
SEQAN_FUNC void
replace(thrust::device_vector<TChar, TAlloc> & target,
        typename Iterator< thrust::device_vector<TChar, TAlloc>, Rooted>::Type pos_begin,
        typename Iterator< thrust::device_vector<TChar, TAlloc>, Rooted>::Type pos_end,
        TSource & source,
        typename Size< thrust::device_vector<TChar, TAlloc> >::Type limit,
        Tag<TExpand> const tag)
{
    replace(target,  position(pos_begin),  position(pos_end), source, tag);
}
*/


///.Function.reserve.param.object.type:Adaption.thrust::device_vector
///.Function.reserve.remarks:For @Adaption.thrust::device_vector|STL Adaptions@, $reserve$ is only guaranteed to have the specified behaviour with $Insist$ and $Generous$.
///.Function.reserve.class:Adaption.thrust::device_vector

template <typename TChar,  typename TAlloc, typename TSize, typename TExpand>
SEQAN_FUNC typename Size< thrust::device_vector<TChar, TAlloc> >::Type
reserve(
    thrust::device_vector<TChar, TAlloc> & seq,
    TSize new_capacity,
    Tag<TExpand> const & tag)
{
    SEQAN_CHECKPOINT;
    seq.reserve(new_capacity);
    return _capacityReturned(seq, new_capacity, tag);
}

template <typename TChar, typename TAlloc, typename TSize>
SEQAN_FUNC typename Size< thrust::device_vector<TChar, TAlloc> >::Type
reserve(
    thrust::device_vector<TChar, TAlloc> & seq,
    TSize new_capacity,
    Insist const &)
{
    SEQAN_CHECKPOINT;
    // do nothing
    return _capacityReturned(seq, new_capacity, Insist());
}

template <typename TChar,  typename TAlloc, typename TSize>
SEQAN_FUNC typename Size< thrust::device_vector<TChar, TAlloc> >::Type
reserve(
    thrust::device_vector<TChar,  TAlloc> & seq,
    TSize new_capacity,
    Limit const &)
{
    SEQAN_CHECKPOINT;
    // do nothing
    return _capacityReturned(seq, new_capacity, Limit());
}

///.Function.resize.param.object.type:Adaption.thrust::device_vector
template <typename TChar,  typename TAlloc, typename TSize, typename TExpand>
SEQAN_FUNC typename Size< thrust::device_vector<TChar,  TAlloc> >::Type
resize(
    thrust::device_vector<TChar, TAlloc> & me,
    TSize new_length,
    Tag<TExpand> const &)
{
    SEQAN_CHECKPOINT;
    me.resize(new_length);
    return me.size();
}

template <typename TChar, typename TAlloc, typename TSize, typename TExpand>
SEQAN_FUNC typename Size< thrust::device_vector<TChar,  TAlloc> >::Type
fill(
    thrust::device_vector<TChar, TAlloc> & me,
    TSize new_length,
    TChar const & val,
    Tag<TExpand> const &)
{
    SEQAN_CHECKPOINT;
    me.resize(new_length, val);
    return me.length();
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEQUENCE_ADAPT_THRUST_VECTOR_H_
