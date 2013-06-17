// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// This file contains the Extender class.
// ==========================================================================

#ifndef SEQAN_EXTRAS_MASAI_EXTENDER_H_
#define SEQAN_EXTRAS_MASAI_EXTENDER_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/find.h>

#include "store.h"

using namespace seqan;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Extender
// ----------------------------------------------------------------------------

template <typename TGenome, typename TReads, typename TDelegate, typename TDistance = HammingDistance, typename TSpec = void>
struct Extender
{
    Holder<TGenome>         genome;
    Holder<TReads>          reads;
    TDelegate               & delegate;

    TReadSeqSize            minErrorsPerRead;
    TReadSeqSize            maxErrorsPerRead;
    TReadSeqSize            seedLength;

    bool                    disabled;

    Extender(TDelegate & delegate, bool disabled = false) :
        delegate(delegate),
        minErrorsPerRead(0),
        maxErrorsPerRead(0),
        seedLength(0),
        disabled(disabled)
    {}
};

// ----------------------------------------------------------------------------
// Class Extender<EditDistance>
// ----------------------------------------------------------------------------

template <typename TGenome, typename TReads, typename TDelegate, typename TSpec>
struct Extender<TGenome, TReads, TDelegate, EditDistance, TSpec>:
    public Extender<TGenome, TReads, TDelegate>
{
    typedef Extender<TGenome, TReads, TDelegate>            TBase;
    typedef Myers<AlignTextBanded<FindPrefix, NMatchesNone_, NMatchesNone_>, True, void> TAlgorithmSpec;

    typedef Segment<TReadSeq, InfixSegment>                 TReadInfix;
    typedef ModifiedString<TReadInfix, ModReverse>          TReadInfixRev;

    typedef PatternState_<TReadInfix, TAlgorithmSpec>       TPatternState;
    typedef PatternState_<TReadInfixRev, TAlgorithmSpec>    TPatternStateRev;

    TPatternState       _patternState;
    TPatternStateRev    _patternStateRev;

    Extender(TDelegate & delegate, bool disabled = false) :
        TBase(delegate, disabled)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction GenomeHost<T>::Type                                  [Extender]
// ----------------------------------------------------------------------------

template <typename TGenome, typename TReads, typename TDelegate, typename TDistance, typename TSpec>
struct GenomeHost<Extender<TGenome, TReads, TDelegate, TDistance, TSpec> >
{
    typedef TGenome Type;
};

// ----------------------------------------------------------------------------
// Metafunction ReadsHost<T>::Type                                   [Extender]
// ----------------------------------------------------------------------------

template <typename TGenome, typename TReads, typename TDelegate, typename TDistance, typename TSpec>
struct ReadsHost<Extender<TGenome, TReads, TDelegate, TDistance, TSpec> >
{
    typedef TReads  Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function onSeedHit()                                              [Extender]
// ----------------------------------------------------------------------------

template <typename TGenome, typename TReads, typename TDelegate, typename TSpec>
inline bool onSeedHit(Extender<TGenome, TReads, TDelegate, HammingDistance, TSpec> & extender,
                      TContigStoreSize contigId,
                      TContigSeqSize contigBegin,
                      TReadSeqStoreSize seqId,
                      TContigSeqSize seedBegin,
                      TReadSeqSize seedErrors)
{
    typedef Segment<TReadSeq, InfixSegment>                 TReadInfix;

    if (extender.disabled)
        return false;

    TReadSeqSize errors = seedErrors;

    TContigSeq & contig = getContigs(getGenome(extender))[contigId];
    TReadSeq & read = getSeqs(getReads(extender))[seqId];
    TReadSeqSize readLength = length(read);

    // Extend left.
    TContigSeqSize matchBegin = contigBegin;

    if (seedBegin > 0)
    {
        TContigSeqSize contigLeftBegin = 0;
        if (contigBegin > seedBegin)
            contigLeftBegin = contigBegin - seedBegin;

        TContigInfix contigLeft(contig, contigLeftBegin, contigBegin);
        TReadInfix readLeft(read, 0, seedBegin);

        if (!_extend(extender, contigLeft, readLeft, errors))
            return false;

        matchBegin = contigLeftBegin;
    }

    // This removes some duplicates.
    if (errors - seedErrors < extender.minErrorsPerRead)
        return false;

    // Extend right.
    TContigSeqSize matchEnd = contigBegin + extender.seedLength;

    if (seedBegin + extender.seedLength < readLength)
    {
        TContigSeqSize contigRightEnd = contigLength(getGenome(extender), contigId);
        if (contigRightEnd > contigBegin + readLength - seedBegin)
            contigRightEnd = contigBegin + readLength - seedBegin;

        TContigInfix contigRight(contig, contigBegin + extender.seedLength, contigRightEnd);
        TReadInfix readRight(read, seedBegin + extender.seedLength, readLength);

        if (!_extend(extender, contigRight, readRight, errors))
            return false;

        matchEnd = contigRightEnd;
    }

    // This removes some duplicates.
    if (errors < extender.minErrorsPerRead)
        return false;

    bool isReverseComplemented = isReverse(getReads(extender), seqId);
    TReadSeqStoreSize readId = getReadId(getReads(extender), seqId);

    onMatch(extender.delegate, contigId, matchBegin, matchEnd, readId, errors, isReverseComplemented);

    return true;
}

template <typename TGenome, typename TReads, typename TDelegate, typename TSpec, typename TContigInfix, typename TReadInfix>
inline bool _extend(Extender<TGenome, TReads, TDelegate, HammingDistance, TSpec> & extender,
                    TContigInfix & contigInfix,
                    TReadInfix & readInfix,
                    TReadSeqSize & errors)
{
    typedef typename Iterator<TContigInfix, Standard>::Type          TContigInfixIterator;
    typedef typename Iterator<TReadInfix, Standard>::Type            TReadInfixIterator;

    if (length(contigInfix) != length(readInfix))
        return false;

    TContigInfixIterator contigIt = begin(contigInfix, Standard());
    TReadInfixIterator readBegin = begin(readInfix, Standard());
    TReadInfixIterator readEnd = end(readInfix, Standard());

    for (TReadInfixIterator readIt = readBegin; readIt != readEnd; ++readIt, ++contigIt)
        if (*readIt != *contigIt)
            if (++errors > extender.maxErrorsPerRead)
                return false;

    return true;
}

// ----------------------------------------------------------------------------
// Function onSeedHit()                                [Extender<EditDistance>]
// ----------------------------------------------------------------------------

template <typename TGenome, typename TReads, typename TDelegate, typename TSpec>
inline bool onSeedHit(Extender<TGenome, TReads, TDelegate, EditDistance, TSpec> & extender,
                      TContigStoreSize contigId,
                      TContigSeqSize contigBegin,
                      TReadSeqStoreSize seqId,
                      TContigSeqSize seedBegin,
                      TReadSeqSize seedErrors)
{
    typedef Segment<TReadSeq, InfixSegment> TReadInfix;

    if (extender.disabled)
        return false;

    TReadSeqSize errors = seedErrors;

    TContigSeq & contig = getContigs(getGenome(extender))[contigId];
    TReadSeq & read = getSeqs(getReads(extender))[seqId];
    TReadSeqSize readLength = length(read);

    // Extend left.
    TContigSeqSize matchBegin = contigBegin;

    if (seedBegin > 0)
    {
        TContigSeqSize contigLeftBegin = 0;
        if (contigBegin > seedBegin + extender.maxErrorsPerRead - errors)
            contigLeftBegin = contigBegin - (seedBegin + extender.maxErrorsPerRead - errors);

        TContigInfix contigLeft(contig, contigLeftBegin, contigBegin);
        TReadInfix readLeft(read, 0, seedBegin);

        if (!_extendLeft(extender, extender._patternStateRev, contigLeft, readLeft, errors, matchBegin))
            return false;
    }

    // This removes some duplicates.
    if (errors - seedErrors < extender.minErrorsPerRead)
        return false;

    // Extend right.
    TContigSeqSize matchEnd = contigBegin + extender.seedLength;

    if (seedBegin + extender.seedLength < readLength)
    {
        TContigSeqSize contigRightEnd = contigLength(getGenome(extender), contigId);
        if (contigRightEnd > contigBegin + readLength - seedBegin + extender.maxErrorsPerRead - errors)
            contigRightEnd = contigBegin + readLength - seedBegin + extender.maxErrorsPerRead - errors;

        if (contigBegin + extender.seedLength >= contigRightEnd)
            return false;

        TContigInfix contigRight(contig, contigBegin + extender.seedLength, contigRightEnd);
        TReadInfix readRight(read, seedBegin + extender.seedLength, readLength);

        if (!_extendRight(extender, extender._patternState, contigRight, readRight, errors, matchEnd))
            return false;
    }

    // This removes some duplicates.
    if (errors < extender.minErrorsPerRead)
        return false;

    bool isReverseComplemented = isReverse(getReads(extender), seqId);
    TReadSeqStoreSize readId = getReadId(getReads(extender), seqId);

    onMatch(extender.delegate, contigId, matchBegin, matchEnd, readId, errors, isReverseComplemented);

    return true;
}

template <typename TGenome, typename TReads, typename TDelegate, typename TSpec, typename TPatternState, typename TContigInfix, typename TReadInfix>
inline bool _extendLeft(Extender<TGenome, TReads, TDelegate, EditDistance, TSpec> & extender,
                        TPatternState & _patternState,
                        TContigInfix & contigInfix,
                        TReadInfix & readInfix,
                        TReadSeqSize & errors,
                        TContigSeqSize & matchBegin)
{
    typedef ModifiedString<TReadInfix, ModReverse>          TReadInfixRev;
    typedef ModifiedString<TContigInfix, ModReverse>        TContigInfixRev;
    typedef Finder<TContigInfixRev>                         TFinder;

    // Lcp trick.
    TContigSeqSize lcp = 0;
    {  // TODO(holtgrew): Workaround to storing and returning copies in host() for nested infixes/modified strings. This is ugly and should be fixed later.
        TReadInfixRev readInfixRev(readInfix);
        TContigInfixRev contigInfixRev(contigInfix);
        lcp = lcpLength(contigInfixRev, readInfixRev);
    }
    if (lcp == length(readInfix))
    {
        matchBegin -= lcp;
        return true;
    }
    setEndPosition(contigInfix, endPosition(contigInfix) - lcp);
    setEndPosition(readInfix, endPosition(readInfix) - lcp);

    TReadSeqSize remainingErrors = extender.maxErrorsPerRead - errors;
    TReadSeqSize minErrors = remainingErrors + 1;
    TContigSeqSize endPos = 0;

    // Stop seed extension.
    if (!remainingErrors)
        return false;

    // Align.
    TReadInfixRev readInfixRev(readInfix);
    TContigInfixRev contigInfixRev(contigInfix);
    TFinder finder(contigInfixRev);
    _patternState.leftClip = remainingErrors;

    // TODO(esiragusa): Use a generic type for errors.
    while (find(finder, readInfixRev, _patternState, -static_cast<int>(remainingErrors)))
    {
        TReadSeqSize currentErrors = -getScore(_patternState);

        if (currentErrors <= minErrors)
        {
            minErrors = currentErrors;
            endPos = position(finder) + 1;
        }
    }

    errors += minErrors;
    matchBegin -= endPos + lcp;

    return errors <= extender.maxErrorsPerRead;
}

template <typename TGenome, typename TReads, typename TDelegate, typename TSpec, typename TPatternState, typename TContigInfix, typename TReadInfix>
inline bool _extendRight(Extender<TGenome, TReads, TDelegate, EditDistance, TSpec> & extender,
                         TPatternState & _patternState,
                         TContigInfix & contigInfix,
                         TReadInfix & readInfix,
                         TReadSeqSize & errors,
                         TContigSeqSize & matchEnd)
{
    typedef Finder<TContigInfix>    TFinder;

    // Lcp trick.
    TContigSeqSize lcp = lcpLength(contigInfix, readInfix);
    if (lcp == length(readInfix))
    {
        matchEnd += lcp;
        return true;
    }
    else if (lcp == length(contigInfix))
    {
        errors += length(readInfix) - length(contigInfix);
        matchEnd += length(readInfix);
        return errors <= extender.maxErrorsPerRead;
    }
    setBeginPosition(contigInfix, beginPosition(contigInfix) + lcp);
    setBeginPosition(readInfix, beginPosition(readInfix) + lcp);

    // NOTE Uncomment this to disable lcp trick.
//    TContigSeqSize lcp = 0;

    TReadSeqSize remainingErrors = extender.maxErrorsPerRead - errors;
    TReadSeqSize minErrors = remainingErrors + 1;
    TContigSeqSize endPos = 0;

    // NOTE Comment this to disable lcp trick.
    // Stop seed extension.
    if (!remainingErrors)
        return false;

    // Remove last base.
    TContigInfix contigPrefix(contigInfix);
    TReadInfix readPrefix(readInfix);
    setEndPosition(contigPrefix, endPosition(contigPrefix) - 1);
    setEndPosition(readPrefix, endPosition(readPrefix) - 1);

    // Align.
    TFinder finder(contigPrefix);
    _patternState.leftClip = remainingErrors;

    // TODO(esiragusa): Use a generic type for errors.
    while (find(finder, readPrefix, _patternState, -static_cast<int>(remainingErrors)))
    {
        TContigSeqSize currentEnd = position(finder) + 1;
        TReadSeqSize currentErrors = -getScore(_patternState);

        // Compare last base.
        if (contigInfix[currentEnd] != back(readInfix))
            if (++currentErrors > remainingErrors)
                continue;

        if (currentErrors <= minErrors)
        {
            minErrors = currentErrors;
            endPos = currentEnd;
        }
    }

    errors += minErrors;
    matchEnd += endPos + lcp + 1;

    return errors <= extender.maxErrorsPerRead;
}

#endif  // #ifndef SEQAN_EXTRAS_MASAI_EXTENDER_H_
