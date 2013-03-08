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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// This file contains the Manager class.
// ==========================================================================

#ifndef SEQAN_EXTRAS_MASAI_MANAGER_H_
#define SEQAN_EXTRAS_MASAI_MANAGER_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/find.h>

#include "tags.h"
#include "store.h"
#include "matches.h"

using namespace seqan;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Manager
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec = void, typename TMatch = Match<TSpec> >
struct Manager
{
    Holder<TReads>      reads;
    TDelegate           & delegate;
    unsigned long       matchesCount;

    Manager(TDelegate & delegate) :
        delegate(delegate),
        matchesCount(0)
    {}
};

// ----------------------------------------------------------------------------
// Class Manager<AllBest>
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TMatch>
struct Manager<TReads, TDelegate, AllBest, TMatch> :
    public Manager<TReads, TDelegate, void, TMatch>
{
    typedef Manager<TReads, TDelegate, void, TMatch>    TBase;

    unsigned                    errors;
    String<unsigned char>       _minErrors;

    Manager(TDelegate & delegate) :
        TBase(delegate),
        errors(0)
    {}
};

// ----------------------------------------------------------------------------
// Class Manager<KBest>
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TMatch>
struct Manager<TReads, TDelegate, KBest, TMatch>:
    public Manager<TReads, TDelegate, void, TMatch>
{
    typedef Manager<TReads, TDelegate, void, TMatch>    TBase;
    typedef String<TMatch>                              TMatches;

    unsigned                    errors;
    String<unsigned char>       _minErrors;
    String<TMatches>            _matches;

    Manager(TDelegate & delegate, TReadSeqStoreSize readsCount) :
        TBase(delegate, readsCount),
        errors(0)
    {}
};

// ----------------------------------------------------------------------------
// Class Manager<AnyBest>
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TMatch>
struct Manager<TReads, TDelegate, AnyBest, TMatch>:
    public Manager<TReads, TDelegate, void, TMatch>
{
    typedef Manager<TReads, TDelegate, void, TMatch>    TBase;
    typedef String<TMatch>                              TMatches;

    unsigned                    errors;
    TMatches                    _matches;

    Manager(TDelegate & delegate) :
        TBase(delegate),
        errors(0)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction ReadsHost<T>::Type                                    [Manager]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec, typename TMatch>
struct ReadsHost<Manager<TReads, TDelegate, TSpec, TMatch> >
{
    typedef TReads  Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function clear()                                                   [Manager]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec, typename TMatch>
void clear(Manager<TReads, TDelegate, TSpec, TMatch> & manager)
{
    manager.matchesCount = 0;
}

// ----------------------------------------------------------------------------
// Function clear()                                          [Manager<AllBest>]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TMatch>
void clear(Manager<TReads, TDelegate, AllBest, TMatch> & manager)
{
    manager.matchesCount = 0;
    clear(manager._minErrors);
}

// ----------------------------------------------------------------------------
// Function clear()                                            [Manager<KBest>]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TMatch>
void clear(Manager<TReads, TDelegate, KBest, TMatch> & manager)
{
    manager.matchesCount = 0;
    clear(manager._minErrors);
    clear(manager._matches);
}

// ----------------------------------------------------------------------------
// Function clear()                                          [Manager<AnyBest>]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TMatch>
void clear(Manager<TReads, TDelegate, AnyBest, TMatch> & manager)
{
    manager.matchesCount = 0;
    clear(manager._matches);
}

// ----------------------------------------------------------------------------
// Function resize()                                                  [Manager]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec, typename TMatch, typename TSize>
void resize(Manager<TReads, TDelegate, TSpec, TMatch> & /* manager */, TSize /* count */)
{}

// ----------------------------------------------------------------------------
// Function resize()                                         [Manager<AllBest>]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TMatch, typename TSize>
void resize(Manager<TReads, TDelegate, AllBest, TMatch> & manager, TSize count)
{
    resize(manager._minErrors, count, MaxValue<unsigned char>::VALUE, Exact());
}

// ----------------------------------------------------------------------------
// Function resize()                                           [Manager<KBest>]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TMatch, typename TSize>
void resize(Manager<TReads, TDelegate, KBest, TMatch> & manager, TSize count)
{
    resize(manager._minErrors, count, MaxValue<unsigned char>::VALUE, Exact());
    // TODO(esiragusa): Change hardcoded size.
    resize(manager._matches, 32, Exact());
}

// ----------------------------------------------------------------------------
// Function resize()                                         [Manager<AnyBest>]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TMatch, typename TSize>
void resize(Manager<TReads, TDelegate, AnyBest, TMatch> & manager, TSize count)
{
    TMatch match;
    fill(match, 0, 0, 0, 0, MaxValue<unsigned char>::VALUE, false);
    resize(manager._matches, count, match, Exact());
}

// ----------------------------------------------------------------------------
// Function setReads()                                                [Manager]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec, typename TMatch>
void setReads(Manager<TReads, TDelegate, TSpec, TMatch> & manager, TReads & reads)
{
    setValue(manager.reads, reads);
    clear(manager);
    resize(manager, reads.readsCount);
}

// ----------------------------------------------------------------------------
// Function onMatch()                                                 [Manager]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec, typename TMatch,
          typename TContigId, typename TContigPos, typename TReadId, typename TErrors>
inline void onMatch(Manager<TReads, TDelegate, TSpec, TMatch> & manager,
                    TContigId contigId,
                    TContigPos beginPos,
                    TContigPos endPos,
                    TReadId readId,
                    TErrors errors,
                    bool reverseComplemented)
{
    // Call matches delegate.
    onMatch(manager.delegate, contigId, beginPos, endPos, readId, errors, reverseComplemented);

    // Increment matches counters.
    manager.matchesCount++;
}

// ----------------------------------------------------------------------------
// Function onMatch()                                        [Manager<AllBest>]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TMatch,
          typename TContigId, typename TContigPos, typename TReadId, typename TErrors>
inline void onMatch(Manager<TReads, TDelegate, AllBest, TMatch> & manager,
                    TContigId contigId,
                    TContigPos beginPos,
                    TContigPos endPos,
                    TReadId readId,
                    TErrors errors,
                    bool reverseComplemented)
{
    // This match is useful.
    if (errors <= manager._minErrors[readId])
    {
        // Call matches delegate.
        onMatch(manager.delegate, contigId, beginPos, endPos, readId, errors, reverseComplemented);

        // Increment matches counters.
        manager.matchesCount++;

        // Disable the read after current stratum.
        manager._minErrors[readId] = errors;
    }

    // Otherwise this match is superfluous.
}

// ----------------------------------------------------------------------------
// Function onMatch()                                          [Manager<KBest>]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TMatch,
          typename TContigId, typename TContigPos, typename TReadId, typename TErrors>
inline void onMatch(Manager<TReads, TDelegate, KBest, TMatch> & manager,
                    TContigId contigId,
                    TContigPos beginPos,
                    TContigPos endPos,
                    TReadId readId,
                    TErrors errors,
                    bool reverseComplemented)
{
    // This match is useful.
    if (errors < manager._minErrors[readId])
    {
        if (errors <= manager.errors)
        {
            // Call matches delegate.
            onMatch(manager.delegate, contigId, beginPos, endPos, readId, errors, reverseComplemented);

            // Increment matches counters.
            manager.matchesCount++;
        }
        else
        {
            // Store match.
            TMatch match;
            fill(match, contigId, beginPos, endPos, readId, errors, reverseComplemented);
            appendValue(manager._matches[errors], match, Generous());
        }
        // Disable the read.
        manager._minErrors[readId] = errors;
    }

    // Otherwise this match is superfluous.
}

// ----------------------------------------------------------------------------
// Function onMatch()                                        [Manager<AnyBest>]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TMatch,
typename TContigId, typename TContigPos, typename TReadId, typename TErrors>
inline void onMatch(Manager<TReads, TDelegate, AnyBest, TMatch> & manager,
                    TContigId contigId,
                    TContigPos beginPos,
                    TContigPos endPos,
                    TReadId readId,
                    TErrors errors,
                    bool reverseComplemented)
{
    // NOTE(esiragusa): This match should always be useful...
    if (errors < manager._matches[readId].errors)
    {
        // Store match.
        fill(manager._matches[readId], contigId, beginPos, endPos, readId, errors, reverseComplemented);

        if (errors == manager.errors)
        {
            // Call matches delegate.
            onMatch(manager.delegate, manager._matches[readId]);

            // Increment matches counters.
            manager.matchesCount++;
        }
    }
}

// ----------------------------------------------------------------------------
// Function isDisabled()                                              [Manager]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec, typename TMatch, typename TReadId>
inline bool isDisabled(Manager<TReads, TDelegate, TSpec, TMatch> const &, TReadId)
{
    return false;
}

// ----------------------------------------------------------------------------
// Function isDisabled()                                     [Manager<AllBest>]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TMatch, typename TReadId>
inline bool isDisabled(Manager<TReads, TDelegate, AllBest, TMatch> const & manager, TReadId readId)
{
    // Reads with at least one match in lower strata are disabled.
    return manager._minErrors[readId] < manager.errors;
}

// ----------------------------------------------------------------------------
// Function isDisabled()                                       [Manager<KBest>]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TMatch, typename TReadId>
inline bool isDisabled(Manager<TReads, TDelegate, KBest, TMatch> const & manager, TReadId readId)
{
    // Reads with at least one best match are disabled.
    return manager._minErrors[readId] <= manager.errors;
}

// ----------------------------------------------------------------------------
// Function isDisabled()                                     [Manager<AnyBest>]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TMatch, typename TReadId>
inline bool isDisabled(Manager<TReads, TDelegate, AnyBest, TMatch> const & manager, TReadId readId)
{
    // Reads with at least one best match are disabled.
    return manager._matches[readId].errors <= manager.errors;
}

// ----------------------------------------------------------------------------
// Function minErrors()                                               [Manager]
// ----------------------------------------------------------------------------

//template <typename TReads, typename TDelegate, typename TSpec, typename TMatch, typename TReadId>
//inline unsigned char minErrors(Manager<TReads, TDelegate, TSpec, TMatch> const &, TReadId)
//{
//    // TODO
//    return 0;
//}
//
//template <typename TReads, typename TDelegate, typename TMatch, typename TReadId>
//inline unsigned char minErrors(Manager<TReads, TDelegate, AnyBest, TMatch> const & manager, TReadId readId)
//{
//    return manager._matches[getReadId(manager, readId)].errors + 1;
//}

// ----------------------------------------------------------------------------
// Function maxErrors()                                               [Manager]
// ----------------------------------------------------------------------------

//template <typename TReads, typename TDelegate, typename TSpec, typename TMatch, typename TReadId>
//inline unsigned char maxErrors(Manager<TReads, TDelegate, TSpec, TMatch> const &, TReadId)
//{
//    // TODO
//    return 0;
//}
//
//template <typename TReads, typename TDelegate, typename TMatch, typename TReadId>
//inline unsigned char maxErrors(Manager<TReads, TDelegate, AnyBest, TMatch> const & manager, TReadId readId)
//{
//    // TODO
//    return 0;
//}

// ----------------------------------------------------------------------------
// Function raiseErrorThreshold()                                     [Manager]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec, typename TMatch>
inline void raiseErrorThreshold(Manager<TReads, TDelegate, TSpec, TMatch> &)
{}

// ----------------------------------------------------------------------------
// Function raiseErrorThreshold()                            [Manager<AllBest>]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TMatch>
inline void raiseErrorThreshold(Manager<TReads, TDelegate, AllBest, TMatch> & manager)
{
    // Increment error threshold.
    manager.errors++;
}

// ----------------------------------------------------------------------------
// Function raiseErrorThreshold()                              [Manager<KBest>]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TMatch>
inline void raiseErrorThreshold(Manager<TReads, TDelegate, KBest, TMatch> & manager)
{
    typedef String<TMatch>                                  TMatches;
    typedef typename Iterator<TMatches, Standard>::Type     TMatchesIterator;

    TMatchesIterator matchesEnd = end(manager._matches[manager.errors + 1], Standard());

    for (TMatchesIterator matchesIt = begin(manager._matches[manager.errors + 1], Standard()); matchesIt != matchesEnd; ++matchesIt)
    {
        // If a read is disabled we already wrote its best match.
        if (isDisabled(manager, (*matchesIt).readId))
            continue;

        // Call matches delegate.
        onMatch(manager.delegate, *matchesIt);

        // Increment matches counter.
        manager.matchesCount++;

        // Disable the read.
        manager._minErrors[(*matchesIt).readId] = manager.errors + 1;
    }

    // Increment error threshold.
    manager.errors++;

    // Forget matches below error threshold.
    clear(manager._matches[manager.errors]);
    shrinkToFit(manager._matches[manager.errors]);
}

// ----------------------------------------------------------------------------
// Function raiseErrorThreshold()                            [Manager<AnyBest>]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TMatch>
inline void raiseErrorThreshold(Manager<TReads, TDelegate, AnyBest, TMatch> & manager)
{
    typedef String<TMatch>                                  TMatches;
    typedef typename Iterator<TMatches, Standard>::Type     TMatchesIterator;

    // Increment error threshold.
    manager.errors++;

    TMatchesIterator matchesEnd = end(manager._matches, Standard());

    for (TMatchesIterator matchesIt = begin(manager._matches, Standard()); matchesIt != matchesEnd; ++matchesIt)
    {
        // Write current best matches.
        if ((*matchesIt).errors == manager.errors)
        {
            // Call matches delegate.
            onMatch(manager.delegate, *matchesIt);

            // Increment matches counter.
            manager.matchesCount++;

        }
    }
}

#endif  // #ifndef SEQAN_EXTRAS_MASAI_MATCHES_H_
