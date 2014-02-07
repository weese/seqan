// ==========================================================================
//                              radix_inplace.h
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
// Author: Sascha Meiers <meiers@inf.fu-berlin.de>
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_INDEX_RADIX_INPLACE_H_
#define CORE_INCLUDE_SEQAN_INDEX_RADIX_INPLACE_H_

namespace SEQAN_NAMESPACE_MAIN
{


// turn on debug output for radix sort (will spam the screen)
//#define CORE_INCLUDE_SEQAN_INDEX_RADIX_INPLACE_H_DEBUG_RADIX_SORT

// ============================================================================
// struct RadixTextAccessor
// ============================================================================

template <typename TSAValue,
typename TSize,
typename TString,
typename TSpec=void >
struct RadixTextAccessor;
/*
 * NOTE:
 * These accessors cannot resolve the correct order of out-of-bound-positions,
 * i.e. when suffixes are equal up to their last character.
 * All these cases get collected in a 0 bucket.
 * The InplaceRadixSorter takes care of that by calling a special
 * sort function on the 0 buckets.
 */

// ----------------------------------------------------------------------------
// struct RadixTextAccessor                                            [String]
// ----------------------------------------------------------------------------

template <typename TSAValue,
typename TSize,
typename TString>
struct RadixTextAccessor<TSAValue, TSize, TString, void> :
public std::unary_function<TSAValue, TSize>
{
    TString const & text;
    typename Size<TString>::Type const L;

    RadixTextAccessor(TString const &str) : text(str), L(length(str))
    {}

    inline TSize operator()(TSAValue x, TSize depth) const
    {
        typename Size<TString>::Type pos = x + depth;
        if (pos >= L)   return 0;
        TSize ret = ordValue(text[pos]);
        return ret+1;
    }
};

// ----------------------------------------------------------------------------
// struct RadixTextAccessor                                         [StringSet]
// ----------------------------------------------------------------------------

template <typename TSAValue,
typename TSize,
typename TString,
typename TSetSpec>
struct RadixTextAccessor<TSAValue, TSize, StringSet<TString, TSetSpec>, void> :
public std::unary_function<TSAValue, TSize>
{
    StringSet<TString, TSetSpec> const & text;
    String<typename Size<TString>::Type> L;

    RadixTextAccessor(StringSet<TString, TSetSpec> const &str) : text(str)
    {
        resize(L, length(text), Exact());
        for(typename Size<TString>::Type i = 0; i < length(text); ++i)
        L[i] = length(text[i]);
    }

    inline TSize operator()(TSAValue x, TSize depth) const
    {
        typename Size<TString>::Type pos = getSeqOffset(x) + depth;
        typename Size<TString>::Type seq = getSeqNo(x);
        if (pos >= L[seq])   return 0;
        TSize ret = ordValue(text[seq][pos]);
        return ret+1;
    }
};

// ----------------------------------------------------------------------------
// struct RadixTextAccessor                               [String, CyclicShape]
// ----------------------------------------------------------------------------

template <typename TSAValue,
typename TSize,
typename TString,
typename TShape>
struct RadixTextAccessor<TSAValue, TSize,TString,ModCyclicShape<CyclicShape<TShape> > > :
public std::unary_function<TSAValue, TSize> // in, out
{
    TString const & text;
    TSize const L,w,s;
    String<unsigned> positions;

    // Cargo < ModifiedString<ModCyclicShape> > = CyclicShape
    RadixTextAccessor(TString const &str, CyclicShape<TShape> const & shape) :
        text(str), L(length(str)), w(weight(shape)), s(shape.span)
    {
        carePositions(positions, shape);
    }

    inline TSize operator()(TSAValue x, TSize depth) const
    {
        typename Size<TString>::Type pos = x + depth/w * s + positions[ depth % w ];
        if (pos >= L) return 0;
        TSize ret = ordValue(text[pos]);
        return ret+1;
    }
};

// ----------------------------------------------------------------------------
// struct RadixTextAccessor                            [StringSet, CyclicShape]
// ----------------------------------------------------------------------------

// TODO: Maybe specialised version of hardwired shape that accesses the text even faster??

template <typename TSAValue,
typename TSize,
typename TString,
typename TSetSpec,
typename TShape>
struct RadixTextAccessor<TSAValue, TSize, StringSet<TString, TSetSpec>,
ModCyclicShape<CyclicShape<TShape> > >: public std::unary_function<TSAValue, TSize>
{
    StringSet<TString, TSetSpec> const & text;
    String<typename Size<TString>::Type> L;
    String<unsigned> positions;
    const TSize w, s;

    RadixTextAccessor(StringSet<TString, TSetSpec> const &str, CyclicShape<TShape> const & shape) :
    text(str),
    w(weight(shape)),
    s(shape.span)
    {
        carePositions(positions, shape);
        resize(L, length(text), Exact());
        for(typename Size<TString>::Type i = 0; i < length(text); ++i)
        L[i] = length(text[i]);
    }

    inline TSize operator()(TSAValue x, TSize depth) const
    {
        typename Size<TString>::Type pos = getSeqOffset(x) + depth/w * s
        + positions[ depth % w ];
        typename Size<TString>::Type seq = getSeqNo(x);
        if (pos >= L[seq])   return 0;
        TSize ret = ordValue(text[seq][pos]);
        return ret+1;
    }
};


// ============================================================================
// functor InplaceRadixSorter and accessories
//
// The following Radix Sort functions are adapted from Martin Frith's "last"
// tool (last.cbrc.jp), but he himself adapted the code from McIlroy, Bostic:
// "Engineering radix sort" as well as Kärkkäinen, Rantala: "Engineering radix
// sort for strings". Thanks to Martin for showing this to me.
// ============================================================================

// ----------------------------------------------------------------------------
// RecursionStack.
// ----------------------------------------------------------------------------
// Self written in the hope of being efficient. Note the hardcoded stack
// size. Maybe switch to some std::deque?

template <typename TSAValue, typename TSmallSize>
struct _RadixRecursionStackEntry
{
    TSAValue * from;
    TSAValue * to;
    TSmallSize depth;
    _RadixRecursionStackEntry()
    {}
    _RadixRecursionStackEntry(TSAValue *a, TSAValue *b, TSmallSize d) :
    from(a), to(b), depth(d)
    {}
};
/*
 template <typename TSAValue, typename TSmallSize=unsigned>
 struct RadixRecursionStack
 {
 typedef _RadixRecursionStackEntry<TSAValue, TSmallSize> TEntry;
 String<TEntry> stack;

 RadixRecursionStack()
 {}

 inline bool empty() { return length(stack) <=0; }

 inline void push(TSAValue *beg, TSAValue *end, TSmallSize depth)
 {
 appendValue(stack, TEntry(beg, end, depth), Generous());
 }
 inline void pop(TSAValue *& beg, TSAValue *& end, TSmallSize &depth)
 {
 TEntry & top = back(stack);
 beg = top.from;
 end = top.to;
 depth = top.depth;
 eraseBack(stack);
 }
 }; */

// TODO: this stack is not very generic, but so much faster than the one above :/

template <typename TSAValue, typename TSmallSize=unsigned>
struct RadixRecursionStack
{
    typedef _RadixRecursionStackEntry<TSAValue, TSmallSize> TEntry;
    TEntry stack[256*200 + 1]; // enough for pattern length 200 on char
    TEntry *top;

    RadixRecursionStack() : top(stack) {}

    inline bool empty() { return top <= stack; }

    inline void push(TSAValue *beg, TSAValue *end, TSmallSize depth)
    {
        top->from = beg;
        top->to = end;
        top->depth = depth;
        ++top;
    }
    inline void pop(TSAValue *& beg, TSAValue *& end, TSmallSize &depth)
    {
        --top;
        beg = top->from;
        end = top->to;
        depth = top->depth;
    }
};


// ----------------------------------------------------------------------------
// InplaceRadixSorter                                   general alphabet <= 256
// ----------------------------------------------------------------------------

template <typename TValue,
unsigned Q,                             // alph size = ValueSize + 1
typename TAccessFunctor,                // text accessor
typename TOrderFunctor,                 // For seperate sort of the 0 bucket.
typename TSize = unsigned,              // type of depth and bucketCount a.s.o
typename TBucketValue = unsigned>       // type the alphabet gets translated to
struct InplaceRadixSorter {

    // TODO: define type 'uchar' according to alphabet
    typedef unsigned char uchar;
    static const unsigned ORACLESIZE = 256;
    TAccessFunctor const & textAccess;
    TOrderFunctor const & comp;

    InplaceRadixSorter(TAccessFunctor const & f, TOrderFunctor const & c) : textAccess(f), comp(c)
    {}

    inline void operator()(TValue * beg,
                           TValue * end,
                           TSize depth,
                           RadixRecursionStack<TValue, TSize> & stack)
    {
        static TSize bucketSize[Q];  // initialized to zero at startup
        TValue* bucketEnd[Q];  // "static" makes little difference to speed

        // get bucket sizes (i.e. letter counts):
        // The intermediate oracle array makes it faster (see "Engineering
        // Radix Sort for Strings" by J Karkkainen & T Rantala)
        for( TValue* i = beg; i < end; /* noop */ )
        {
            uchar oracle [ORACLESIZE]; // buffer for the next chars
            uchar* oracleEnd = oracle + std::min( sizeof(oracle), std::size_t(end - i) );

            for( uchar* j = oracle; j < oracleEnd; ++j )
                *j = textAccess(*i++, depth);

            for( uchar* j = oracle; j < oracleEnd; ++j )
                ++bucketSize[ *j ];
        }

        // get bucket ends, and put buckets on the stack to sort within them later:
        // EDIT: 0 bucket is not sorted here !
        TSize zeroBucketSize = bucketSize[0];
        TValue* pos     = beg + bucketSize[0];
        bucketEnd[0] = pos;

        for( unsigned i = 1; i < Q; ++i )
        {
            TValue* nextPos = pos + bucketSize[i];
            if (nextPos - pos > 1)
            stack.push(pos, nextPos, depth+1);
            pos = nextPos;
            bucketEnd[i] = pos;
        }

        // permute items into the correct buckets:
        for( TValue* i = beg; i < end; ) {
            unsigned subset;  // unsigned is faster than uchar!
            TValue holdOut = *i;
            while( --bucketEnd[ subset = textAccess(holdOut, depth) ] > i )
            std::swap( *bucketEnd[subset], holdOut );
            *i = holdOut;
            i += bucketSize[subset];
            bucketSize[subset] = 0;  // reset it so we can reuse it
        }

        // sort the 0 bucket using std::sort
        if(zeroBucketSize > 1) {

#ifdef CORE_INCLUDE_SEQAN_INDEX_RADIX_INPLACE_H_DEBUG_RADIX_SORT
            std::cout << "Sorted zero bucket of size " << zeroBucketSize << std::endl;
#endif
            std::sort(beg, beg+zeroBucketSize, comp);
            //std::cout << "sort 0 bucket: " << zeroBucketSize << std::endl;
        }

    }
};

// ----------------------------------------------------------------------------
// InplaceRadixSorter                                                 Dna (4+1)
// ----------------------------------------------------------------------------
/*
 template <typename TValue, typename TFunctor, typename TSize>
 struct Radix<TValue, 5, TFunctor, TSize> {

 TFunctor const & textAccess;

 Radix(TFunctor const & f) : textAccess(f)
 {}

 inline void operator()(
 TValue *beg,
 TValue *end,
 TSize depth,
 RadixRecursionStack<TValue, TSize> & stack)
 {
 TValue *end0 = beg,
 *end1 = beg,
 *end2 = beg,
 *end3 = beg,
 *beg4 = end;
 while(end3 < beg4)
 {
 TValue x = *end3;
 switch(textAccess(x, depth) ) {
 case 0:
 *end3++ = *end2;
 *end2++ = *end1;
 *end1++ = *end0;
 *end0++ = x;
 break;
 case 1: // A
 *end3++ = *end2;
 *end2++ = *end1;
 *end1++ = x;
 break;
 case 2: // C
 *end3++ = *end2;
 *end2++ = x;
 break;
 case 3: // G
 ++end3;
 break;
 default: // T
 *end3 = *--beg4;
 *beg4 = x;
 break;
 }
 }
 stack.push(beg, end0, depth+1); // which order?
 stack.push(end3, end, depth+1);
 stack.push(end2, end3, depth+1);
 stack.push(end1, end2, depth+1);
 stack.push(end0, end1, depth+1);
 }
 }; */
template <typename TStr, typename TSA, typename TPos>
void __outputSA(TStr const & str, TSA const & sa, TPos from, TPos to)
{
    for(TPos x = from; x < to; ++x)
        if (length(suffix(str, sa[x])) > 20)
            std::cout << x << ": " << sa[x] << "  \t" << prefix(suffix(str, sa[x]),20) << "..." << std::endl;
        else
            std::cout << x << ": " << sa[x] << "  \t" << suffix(str, sa[x]) << std::endl;
}

// ----------------------------------------------------------------------------
// Functors to compare suffixes from 0 bucket (no proper suffixes)
// ----------------------------------------------------------------------------

template <typename TSAValue, typename TSuffixModifier=void>
struct _ZeroBucketComparator
{
    _ZeroBucketComparator(){}

    inline bool operator()(TSAValue const & a, TSAValue const & b)
    {
        return a > b;
    }
};

// StringSet
template <typename TV1, typename TV2, typename TSpec, typename TSuffixModifier>
struct _ZeroBucketComparator<Pair<TV1, TV2, TSpec>, TSuffixModifier >
{
    typedef Pair<TV1, TV2, TSpec> TSAValue;
    _ZeroBucketComparator(){}

    inline bool operator()(TSAValue const & a, TSAValue const & b)
    {
        if(getSeqNo(a) != getSeqNo(b))
        return getSeqNo(a) > getSeqNo(b);
        return getSeqOffset(a) > getSeqOffset(b);
    }
};

// ----------------------------------------------------------------------------
// Function inplaceRadixSort()                                        [default]
// ----------------------------------------------------------------------------

template <typename TSA, typename TString>
void inplaceRadixSort(
                      TSA & sa,
                      TString const & str,
                      typename Size<TString>::Type maxDepth)
{
    typedef typename Value<typename Concatenator<TString>::Type>::Type TAlphabet;
    typedef typename Value<TSA>::Type                               TSAValue;
    typedef typename Size<TAlphabet>::Type                          TOrdValue;
    typedef typename Size<TString>::Type                            TSize;
    typedef RadixTextAccessor<TSAValue, TOrdValue, TString>         TAccessFunctor;
    typedef _ZeroBucketComparator<TSAValue>                         TCompareFunctor;
    static const unsigned SIGMA = ValueSize<TAlphabet>::VALUE + 1;
    typedef InplaceRadixSorter<TSAValue, SIGMA, TAccessFunctor, TCompareFunctor, TSize>    TSorter;

    if (length(sa) < 1) return; // otherwise access sa[0] fails

    TAccessFunctor textAccess(str);
    TSorter radixSort(textAccess, TCompareFunctor());

    RadixRecursionStack<TSAValue, TSize> stack;
    stack.push(&sa[0], &sa[0]+length(sa), 0);

    while(!stack.empty())
    {
        TSAValue *from;
        TSAValue *to;
        TSize currDepth;
        stack.pop(from, to, currDepth);

        if(currDepth >= maxDepth)
        continue;

        radixSort(from, to, currDepth, stack);

#ifdef CORE_INCLUDE_SEQAN_INDEX_RADIX_INPLACE_H_DEBUG_RADIX_SORT
        unsigned dbg_from = from - &sa[0];
        unsigned dbg_to = to - &sa[0];
        std::cout << "Sorted from " << dbg_from << " to " << dbg_to << " in depth " << currDepth <<  std::endl;
        __outputSA(str, sa, dbg_from, dbg_to);
#endif

    }
}


// ----------------------------------------------------------------------------
// Function inplaceRadixSort()                              [modified Suffixes]
// ----------------------------------------------------------------------------
// NOTE: General for all cyclic suffix modifiers, as long as a corresponding
//      text accessor exists.

template <typename TSA, typename TString, typename TMod>
void inplaceRadixSort(
                      TSA & sa,
                      TString const & str,
                      typename Size<TString>::Type maxDepth,
                      typename Cargo<ModifiedString<TString, TMod> >::Type const & modiferCargo,
                      TMod const &)
{
    typedef typename Value<typename Concatenator<TString>::Type>::Type TAlphabet;
    typedef typename Value<TSA>::Type                               TSAValue;
    typedef typename Size<TAlphabet>::Type                          TOrdValue;
    typedef typename Size<TString>::Type                            TSize;
    typedef RadixTextAccessor<TSAValue, TOrdValue, TString, TMod>   TAccessFunctor;
    typedef _ZeroBucketComparator<TSAValue>                         TCompareFunctor;
    static const unsigned SIGMA = ValueSize<TAlphabet>::VALUE + 1;
    typedef InplaceRadixSorter<TSAValue, SIGMA, TAccessFunctor, TCompareFunctor, TSize>    TSorter;

    if (length(sa) < 1) return; // otherwise access sa[0] fails

    TAccessFunctor textAccess(str, modiferCargo);
    TSorter radixSort(textAccess, TCompareFunctor());

    RadixRecursionStack<TSAValue, TSize> stack;
    stack.push(&sa[0], &sa[0]+length(sa), 0);

    while(!stack.empty())
    {
        TSAValue *from;
        TSAValue *to;
        TSize currDepth;
        stack.pop(from, to, currDepth);

        if(currDepth >= maxDepth)
            continue;

        radixSort(from, to, currDepth, stack);

#ifdef CORE_INCLUDE_SEQAN_INDEX_RADIX_INPLACE_H_DEBUG_RADIX_SORT
        unsigned dbg_from = from - &sa[0];
        unsigned dbg_to = to - &sa[0];
        std::cout << "Sorted from " << dbg_from << " to " << dbg_to << " in depth " << currDepth <<  std::endl;
        __outputSA(str, sa, dbg_from, dbg_to);
#endif

    }
}


// ----------------------------------------------------------------------------
// Function inplaceFullRadixSort()                                    [default]
// ----------------------------------------------------------------------------

template <typename TSA, typename TString>
void inplaceFullRadixSort( TSA & sa, TString const & str)
{
    typedef typename Value<typename Concatenator<TString>::Type>::Type TAlphabet;
    typedef typename Value<TSA>::Type                               TSAValue;
    typedef typename Size<TAlphabet>::Type                          TOrdValue;
    typedef typename Size<TString>::Type                            TSize;
    typedef RadixTextAccessor<TSAValue, TOrdValue, TString>         TAccessFunctor;
    typedef _ZeroBucketComparator<TSAValue>                         TCompareFunctor;
    static const unsigned SIGMA = ValueSize<TAlphabet>::VALUE + 1;
    typedef InplaceRadixSorter<TSAValue, SIGMA, TAccessFunctor, TCompareFunctor, TSize>    TSorter;

    if (length(sa) < 1) return; // otherwise access sa[0] fails

    TAccessFunctor textAccess(str);
    TSorter radixSort(textAccess, TCompareFunctor());

    RadixRecursionStack<TSAValue, TSize> stack;
    stack.push(&sa[0], &sa[0]+length(sa), 0);

    while(!stack.empty())
    {
        // DEBUG
        //        for(int i=0; i< (20u < length(sa) ? 20u : length(sa)); ++i)
        //            std::cout << getSeqOffset(sa[i]) << ":" << suffix(str,sa[i]) << ", ";
        //        std::cout << std::endl;
        // END DEBUG


        TSAValue *from;
        TSAValue *to;
        TSize currDepth;
        stack.pop(from, to, currDepth);

        if(to - from < 2)
        continue;
        /*
        // other sort algorithm for small buckets:
        if(to - from < 9)
        {
            ::std::sort( from, to,
                        SuffixLess_<TSAValue, TString const, void>(str, currDepth));
#ifdef CORE_INCLUDE_SEQAN_INDEX_RADIX_INPLACE_H_DEBUG_RADIX_SORT
            unsigned dbg_from = from - &sa[0];
            unsigned dbg_to = to - &sa[0];
            std::cout << "Sorted from " << dbg_from << " to " << dbg_to << " in depth " << currDepth <<  " with Q Sort " << std::endl;
            __outputSA(str, sa, dbg_from, dbg_to);
#endif
            continue;
        }
        */
        radixSort(from, to, currDepth, stack);

#ifdef CORE_INCLUDE_SEQAN_INDEX_RADIX_INPLACE_H_DEBUG_RADIX_SORT
        unsigned dbg_from = from - &sa[0];
        unsigned dbg_to = to - &sa[0];
        std::cout << "Sorted from " << dbg_from << " to " << dbg_to << " in depth " << currDepth <<  std::endl;
        __outputSA(str, sa, dbg_from, dbg_to);
#endif

    }
}


// ----------------------------------------------------------------------------
// Function inplaceFullRadixSort()                          [modified Suffixes]
// ----------------------------------------------------------------------------
// NOTE: General for all cyclic suffix modifiers, as long as a corresponding
//      text accessor in radixSort exists.
template <typename TSA, typename TString, typename TMod>
void inplaceFullRadixSort(TSA & sa,
                          TString const & str,
                          typename Cargo<ModifiedString<TString, TMod> >::Type const & modiferCargo,
                          TMod const &)
{
    typedef typename Value<typename Concatenator<TString>::Type>::Type TAlphabet;
    typedef typename Value<TSA>::Type                               TSAValue;
    typedef typename Size<TAlphabet>::Type                          TOrdValue;
    typedef typename Size<TString>::Type                            TSize;
    typedef RadixTextAccessor<TSAValue, TOrdValue, TString, TMod>   TAccessFunctor;
    typedef _ZeroBucketComparator<TSAValue>                         TCompareFunctor;
    static const unsigned SIGMA = ValueSize<TAlphabet>::VALUE + 1;
    typedef InplaceRadixSorter<TSAValue, SIGMA, TAccessFunctor, TCompareFunctor, TSize>    TSorter;
    
    if (length(sa) < 1) return; // otherwise access sa[0] fails

    TAccessFunctor textAccess(str, modiferCargo);
    TSorter radixSort(textAccess, TCompareFunctor());
    
    RadixRecursionStack<TSAValue, TSize> stack;
    stack.push(&sa[0], &sa[0]+length(sa), 0);
    
    while(!stack.empty())
    {
        TSAValue *from;
        TSAValue *to;
        TSize currDepth;
        stack.pop(from, to, currDepth);
        
        if(to - from < 2)
        continue;
        /*
        // other sort algorithm for small buckets:
        if(to - from < 9)
        {
            ::std::sort( from, to,
                        SuffixLess_<TSAValue, TString const, TMod>(str, modiferCargo, currDepth));
#ifdef CORE_INCLUDE_SEQAN_INDEX_RADIX_INPLACE_H_DEBUG_RADIX_SORT
            unsigned dbg_from = from - &sa[0];
            unsigned dbg_to = to - &sa[0];
            std::cout << "Sorted from " << dbg_from << " to " << dbg_to << " in depth " << currDepth << " with Q Sort " <<  std::endl;
            __outputSA(str, sa, dbg_from, dbg_to);
#endif
            continue;
        }
        */
        radixSort(from, to, currDepth, stack);

#ifdef CORE_INCLUDE_SEQAN_INDEX_RADIX_INPLACE_H_DEBUG_RADIX_SORT
        unsigned dbg_from = from - &sa[0];
        unsigned dbg_to = to - &sa[0];
        std::cout << "Sorted from " << dbg_from << " to " << dbg_to << " in depth " << currDepth <<  std::endl;
        __outputSA(str, sa, dbg_from, dbg_to);
#endif
    }
}
    
    
    
}


#endif  // #ifndef CORE_INCLUDE_SEQAN_INDEX_RADIX_INPLACE_H_
