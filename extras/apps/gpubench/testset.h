// File:   testset.cpp
// Author: David Weese <david.weese@fu-berlin.de>

#ifndef SEQAN_TESTSET_H
#define SEQAN_TESTSET_H

#include <iostream>
#include <seqan/index.h>
#include <seqan/random.h>
#include <seqan/misc/misc_cmdparser.h>
//#include "misc_resources.h"

namespace seqan {

struct BenchParams 
{
    CharString              refFileName;
    CharString              refName;
    String<__uint64>        textLengths;
    String<unsigned>        patternLengths;

//  READ GENERATION
    unsigned                patternAmount;
    double                  findPositiveFrac;
    Rng<MersenneTwister>    rng;

//  STRING SEARCH
    unsigned                maxErrors;

    BenchParams()
    {
        appendValue(patternLengths, 50);
        patternAmount = 100000;
        findPositiveFrac = 0.9;

        maxErrors = 0;
    }

};

template <typename TPattern>
inline void
generateRandomString(BenchParams & params, TPattern & pattern)
{
    typedef typename Value<TPattern>::Type TValue;
    typedef typename Size<TPattern>::Type TSize;

    Pdf<Uniform<int> > pdf(0, ValueSize<TValue>::VALUE - 1);
    for (TSize i = 0; i < length(pattern); ++i)
//        pattern[i] = (TValue)pickRandomNumber(params.rng, pdf);
        pattern[i] = (TValue)0;
}

// FIND BENCH
template <typename TPattern, typename TSpec, typename TPosString, typename TText, typename TPatternSize>
inline void
generatePatterns(BenchParams & params, StringSet<TPattern, TSpec> & patterns, TPosString & pos, TText & text, TPatternSize patternLength)
{
    typedef typename Size<TText>::Type TSize;
    TSize textLength = length(text);
    TPattern randomPattern;
    Pdf<Uniform<TSize> > pdfPos(0, textLength - patternLength);
    Pdf<Uniform<double> > pdfFrac(0.0, 1.0);

    resize(randomPattern, patternLength);
    reserve(patterns.concat, patternLength * params.patternAmount, Exact());
    resize(pos, params.patternAmount);
    for (unsigned i = 0; i < params.patternAmount; ++i)
    {
        assert(textLength >= patternLength);
        TSize beginPos = pickRandomNumber(params.rng, pdfPos);
        appendValue(pos, beginPos);

        if (pickRandomNumber(params.rng, pdfFrac) <= params.findPositiveFrac)
            appendValue(patterns, infix(text, beginPos, beginPos + patternLength));
        else
        {
            generateRandomString(params, randomPattern);
            appendValue(patterns, randomPattern);
        }
    }
}

// FIND BENCH
template <typename TPosString, typename TPattern, typename TText, typename TMatchCount>
inline void
generateMatches(BenchParams & params, TPosString & pos, TPattern & pattern, TText & text, TMatchCount matchCount)
{
    typedef typename Size<TText>::Type TSize;
    TSize textLength = length(text);
    TSize patternLength = length(pattern);

    Pdf<Uniform<TSize> > pdfPos(patternLength, textLength - 2 * patternLength);
    Pdf<Uniform<double> > pdfFrac(0.0, 1.0);

    for (unsigned i = 0; i < matchCount; ++i)
    {
        assert(textLength >= patternLength);
        TSize beginPos = pickRandomNumber(params.rng, pdfPos);
        appendValue(pos, beginPos);

        if (pickRandomNumber(params.rng, pdfFrac) <= params.findPositiveFrac)
            infix(text, beginPos, beginPos + patternLength) = pattern;
    }
}

} // namespace seqan

#endif