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
// This file contains the ReadsIndex class.
// ==========================================================================

#ifndef SEQAN_EXTRAS_MASAI_READS_INDEX_H_
#define SEQAN_EXTRAS_MASAI_READS_INDEX_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index_extras.h>

using namespace seqan;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class ReadsIndex
// ----------------------------------------------------------------------------

template <typename TReads, typename TIndex, typename TSpec = void>
struct ReadsIndex
{
    Holder<TReads>      reads;
    TIndex              index;

    ReadsIndex() {}

    ReadsIndex(TReads const & reads) :
        reads(reads)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction ReadsHost<T>::Type                                 [ReadsIndex]
// ----------------------------------------------------------------------------

template <typename TReads, typename TIndex, typename TSpec>
struct ReadsHost<ReadsIndex<TReads, TIndex, TSpec> >
{
    typedef TReads  Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function build()                                                [ReadsIndex]
// ----------------------------------------------------------------------------

template <typename TReads, typename TSpec, typename TManager>
void build(ReadsIndex<TReads, TReadsQGram, TSpec> & readsIndex,
           TManager & manager,
           TReadSeqSize seedsLength,
           TReadSeqSize firstSeed,
           TReadSeqSize lastSeed)
{
    typedef typename Fibre<TReadsQGram, QGramSA>::Type          TReadsIndexSAFibre;
    typedef typename Fibre<TReadsQGram, QGramDir>::Type         TReadsIndexDirFibre;
    typedef typename Fibre<TReadsQGram, QGramShape>::Type       TReadsIndexShape;
    typedef typename Fibre<TReadsQGram, QGramBucketMap>::Type   TReadsIndexBucketMap;

    typedef typename Value<TReadsIndexDirFibre>::Type           TSize;
    typedef Iterator<TReadSeq, Standard>::Type                  TReadSeqIterator;

    TReadSeqStoreSize seqsCount = length(getSeqs(getReads(readsIndex)));

    readsIndex.index = TReadsQGram(getSeqs(getReads(readsIndex)));

    setStepSize(readsIndex.index, seedsLength);

    TReadsIndexSAFibre & sa = indexSA(readsIndex.index);
    TReadsIndexDirFibre & dir = indexDir(readsIndex.index);
    TReadsIndexShape & shape = indexShape(readsIndex.index);
    TReadsIndexBucketMap & bucketMap = indexBucketMap(readsIndex.index);

    // Resize suffix array and directory.
    resize(sa, (lastSeed - firstSeed) * seqsCount, Exact());
    resize(dir, _fullDirLength(readsIndex.index), Exact());

    // Clear directory.
    _qgramClearDir(dir, bucketMap);

    // Count qgrams.
    for (TReadSeqStoreSize seqId = 0; seqId < seqsCount; ++seqId)
    {
        // Skip disabled reads.
        if (isDisabled(manager, getReadId(getReads(readsIndex), seqId)))
            continue;

        TReadSeq & read = getSeqs(getReads(readsIndex))[seqId];
        TReadSeqIterator itText = begin(read, Standard());

        itText += seedsLength * firstSeed;
        for (TSize i = firstSeed; i < lastSeed; ++i)
        {
            ++dir[requestBucket(bucketMap, hash(shape, itText))];
            itText += seedsLength;
        }
    }

    // Compute cumulative sum.
    _qgramCummulativeSum(dir, False());

    // Fill suffix array.
    for (TReadSeqStoreSize seqId = 0; seqId < seqsCount; ++seqId)
    {
        // Skip disabled reads.
        if (isDisabled(manager, getReadId(getReads(readsIndex), seqId)))
            continue;

        TReadSeq & read = getSeqs(getReads(readsIndex))[seqId];
        TReadSeqIterator itText = begin(read, Standard());

        typename Value<TReadsIndexSAFibre>::Type localPos;
        assignValueI1(localPos, seqId);
        assignValueI2(localPos, 0);

        itText += seedsLength * firstSeed;
        for (TSize i = firstSeed; i < lastSeed; ++i)
        {
            assignValueI2(localPos, seedsLength * i);

            sa[dir[getBucket(bucketMap, hash(shape, itText)) + 1]++] = localPos;

            itText += seedsLength;
        }
    }

    // Refine buckets.
//    _refineQGramIndex(sa, dir, indexText(readsIndex.index), weight(shape), seedsLength);
//    _setHost(readsIndex.index);
}

// ----------------------------------------------------------------------------
// Function build()                                                [ReadsIndex]
// ----------------------------------------------------------------------------

template <typename TReads, typename TSpec, typename TManager>
void build(ReadsIndex<TReads, TReadsWotd, TSpec> & readsIndex,
           TManager & manager,
           TReadSeqSize seedsLength,
           TReadSeqSize firstSeed,
           TReadSeqSize lastSeed)
{
    typedef typename Fibre<TReadsWotd, FibreSA>::Type           TReadsIndexSAFibre;
    typedef typename Value<TReadsIndexSAFibre>::Type            TReadsIndexSAPos;

//    clear(readsIndex.index);
    readsIndex.index = TReadsWotd(getSeqs(getReads(readsIndex)));
    TReadSeqStoreSize seqsCount = length(getSeqs(getReads(readsIndex)));

    TReadsIndexSAFibre & sa = indexSA(readsIndex.index);

    reserve(sa, (lastSeed - firstSeed) * seqsCount, Exact());

    for (TReadSeqStoreSize seqId = 0; seqId < seqsCount; ++seqId)
    {
        // Skip disabled reads.
        if (isDisabled(manager, getReadId(getReads(readsIndex), seqId)))
            continue;

        for (TReadSeqSize seed = firstSeed; seed < lastSeed; ++seed)
        {
            TReadsIndexSAPos localPos;
            assignValueI1(localPos, seqId);
            assignValueI2(localPos, seed * seedsLength);
            appendValue(sa, localPos, Exact());
        }
    }
}

// ----------------------------------------------------------------------------
// Function visit()                                                [ReadsIndex]
// ----------------------------------------------------------------------------

template <typename TReads, typename TIndex, typename TSpec, typename TDepth>
void visit(ReadsIndex<TReads, TIndex, TSpec> & readsIndex, TDepth depth)
{
    typedef typename Iterator<TIndex, TopDown<ParentLinks<> > >::Type  TReadsIndexIterator;
    TReadsIndexIterator readsIt(readsIndex.index);

    do
    {
        std::cout << representative(readsIt) << std::endl;
        if (repLength(readsIt) >= depth || !goDown(readsIt))
            if (!goRight(readsIt))
                while (goUp(readsIt) && !goRight(readsIt)) ;
    }
    while (!isRoot(readsIt));
}

// ----------------------------------------------------------------------------
// Function clear()                                                [ReadsIndex]
// ----------------------------------------------------------------------------

template <typename TReads, typename TIndex, typename TSpec>
void clear(ReadsIndex<TReads, TIndex, TSpec> & readsIndex)
{
    clear(readsIndex.index);
}

#endif  // #ifndef SEQAN_EXTRAS_MASAI_READS_INDEX_H_
