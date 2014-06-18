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

#ifndef SEQAN_STREAM_STREAM_CONCATTER_H_
#define SEQAN_STREAM_STREAM_CONCATTER_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Classes
// ============================================================================

// writes a sequence of pages linked by keys ... -> prevKey -> key -> ...
template <typename TAlgTag, typename TKey>
Pager<TTarget, Concatter<TKey> >
{
    typedef std::map<TKey, std::pair<TPage*, TKey> >    TPageMap;
    rtpedef ConcurrentResourcePool<TValue, TSpec>       TPagePool;

    TTarget         &target;
    TKey            waitForKey;
    TPageMap        pageMap;
    TPagePool       pagePool;
    ReadWriteLock   lock;

    Pager(TTarget &target):
        target(target),
        waitForKey(TKey())
    {}

    Pager(TTarget &target, TKey firstKey):
        target(target),
        waitForKey(firstKey)
    {}

    TPage getPage ()
    {
        TPage tmp;
        tryAquireSwap(tmp, pagePool);
        return tmp;
    }

    void insertPage (TPage &page, TKey key, TKey nextKey)
    {
        if (key == waitForKey)
        {
            writeN(target, page.buffer.begin, length(page.buffer));
            waitForKey = nextKey;
            ReleaseSwap(pagePool, page);
        }
        else
        {
            ScopedWriteLock writeLock(lock);
            pageMap.insert(std::make_pair(key, std::make_pair(&page, nextKey)));

            typename TPageMap::iterator it;
            while ((it = pageMap.find(waitForKey)) != pageMap.end())
            {
                TPage &writePage = it->second.first;
                // currently this write is synchronous
                writeN(target, writePage.buffer.begin, length(writePage.buffer));
                ReleaseSwap(pagePool, writePage);

                waitForKey = it->second.second;
                pageMap.erase(it);
            }
        }
    }
};

// ============================================================================
// Functions
// ============================================================================

#endif  // SEQAN_STREAM_STREAM_CONCATTER_H_

