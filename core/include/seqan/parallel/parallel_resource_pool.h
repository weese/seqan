// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2014, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this pool of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this pool of conditions and the following disclaimer in the
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
// Thread-safe pool
// ==========================================================================

#ifndef SEQAN_PARALLEL_PARALLEL_LIST_H_
#define SEQAN_PARALLEL_PARALLEL_LIST_H_

namespace seqan {

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class ConcurrentResourcePool
// ----------------------------------------------------------------------------
/*!
 * @class ConcurrentResourcePool Concurrent Resource Pool
 * @headerfile <seqan/parallel.h>
 * @brief Thread-safe pool for multiple producers and multiple consumers.
 *
 * @signature template <typename TValue, typename TSpec>
 *            class ConcurrentResourcePool;
 *
 * @tparam TValue Element type of the pool.
 * @tparam TSpec  Tag for further specializing the Concurrent Pool. Default is <tt>void</tt>.
 *
 * The Concurrent Pool is a thread-safe FIFO pool that supports multiple producers and multiple consumers (MPMC).
 * Elements are inserted via @link ConcurrentResourcePool#appendValue @endlink and depoold with @link
 * ConcurrentResourcePool#tryPopFront @endlink or @link ConcurrentResourcePool#popFront @endlink.
 * Depending on the expansion tag of appendValue it can grow dynamically or have a fixed size.
 *
 * The implementation is lock-free and uses a @Class.AllocString@ as ring buffer.
 *
 * @section Examples
 *
 * Simple example for a single producer single consumer (SPSC) dynamic pool.
 *
 * @include demos/parallel/pool_example.cpp
 *
 * The output is:
 *
 * @include demos/parallel/pool_example.cpp.stdout
 */


template <typename TValue>
struct PoolElement
{
    ListElement *next;
    TValue      val;

    PoolElement() :
        next(NULL),
        val(TValue())
    {}
};

template <typename TValue, typename TSpec = void>
class ConcurrentResourcePool
{
public:
    typedef PoolElement<TValue> TPoolElement;

    TPoolElement    *head;          char pad1[SEQAN_CACHE_LINE_SIZE - sizeof(PoolElement)];
    TPoolElement    *headEmpty;     char pad2[SEQAN_CACHE_LINE_SIZE - sizeof(PoolElement)];
    TPoolElement    *unused;        char pad3[SEQAN_CACHE_LINE_SIZE - sizeof(PoolElement)];

    ~ConcurrentResourcePool()
    {
        clear(*this);
    }

private:
    ConcurrentResourcePool(ConcurrentResourcePool const &);
    void operator=(ConcurrentResourcePool const &);
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct Value<ConcurrentResourcePool<TValue, TSpec> >
{
    typedef TValue Type;
};

template <typename TValue, typename TSpec>
struct Host<ConcurrentResourcePool<TValue, TSpec> >
{
    typedef String<TValue> Type;
};

template <typename TValue, typename TSpec>
struct Size<ConcurrentResourcePool<TValue, TSpec> >:
    Size<Host<ConcurrentResourcePool<TValue, TSpec> > >
{};

// ============================================================================
// Functions
// ============================================================================

template <typename TValue, typename TSpec>
inline void
tryAquireSwap(TValue &val, ConcurrentResourcePool<TValue, TSpec> & me)
{
    PoolElement<TValue> *elem;
    while (true)
    {
        elem = me.head;
        if (elem == NULL)
            return;

        // try to take one element from stack of used elements
        if (atomicCasBool(&me.head, elem, elem->next))
            break;
    }

    std::swap(val, elem->val);

    while (true)
    {
        elem->next = me.headEmpty;
        // try to put that element on stack of empty elements
        if (atomicCasBool(&me.headEmpty, elem->next, elem))
            break;
    }
}

template <typename TValue, typename TSpec>
inline void
releaseSwap(ConcurrentResourcePool<TValue, TSpec> & me, TValue SEQAN_FORWARD_ARG val)
{
    PoolElement<TValue> *elem;
    while (true)
    {
        elem = me.headEmpty;
        if (elem == NULL)
        {
            elem = new PoolElement();
            break;
        }

        // try to take one element from stack of empty elements
        if (atomicCasBool(&me.headEmpty, elem, elem->next))
            break;
    }

    std::swap(elem->val, val);

    while (true)
    {
        elem->next = me.head;
        // try to put that element on stack of used elements
        if (atomicCasBool(&me.head, elem->next, elem))
            break;
    }
}

}  // namespace seqan

#endif  // #ifndef SEQAN_PARALLEL_PARALLEL_LIST_H_
