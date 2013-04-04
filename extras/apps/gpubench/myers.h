#ifndef SEQAN_MYERS_H
#define SEQAN_MYERS_H

namespace seqan {

template <typename TSimdVector_>
struct MyersStateSimd
{
    // SIMD version
    // 1 pattern x m texts
    //
    // patternBitmasks is a SIMD vector and
    // encodes all bitmasks
    
    typedef TSimdVector_ TSimdVector;
    typedef TSimdVector_ TPatternBitmasks;

    TSimdVector vp;
    TSimdVector vn;
    TSimdVector errors;
    TPatternBitmasks patternBitmasks;
    TPatternBitmasks patternOr[5];
    unsigned maxErrors;
    unsigned leftClip;

    SEQAN_FUNC MyersStateSimd():
        maxErrors(0),
        leftClip(0) {}
};

template <typename TWord>
struct MyersState
{
    // standard version
    // 1 pattern x 1 text
    //
    // patternBitmasks[0..sigma-1] is an array and
    // encode all bitmasks

    typedef TWord TSimdVector;
    typedef TWord TPatternBitmasks[5];
    TWord vp;
    TWord vn;
    TWord errors;
    TPatternBitmasks patternBitmasks;
    unsigned maxErrors;
    unsigned leftClip;

    SEQAN_FUNC MyersState():
        maxErrors(0),
        leftClip(0) {}
};


// non-SIMD wrappers

SEQAN_FUNC void clear(char &val)            { val = 0; }
SEQAN_FUNC void clear(signed char &val)     { val = 0; }
SEQAN_FUNC void clear(unsigned char &val)   { val = 0; }
    
template <typename TValue>
SEQAN_FUNC void
fill(TValue &val, TValue const &x)
{
    val = x;
}

template <typename TValue>
SEQAN_FUNC TValue
shiftRightLogical(TValue const &val, const int imm)
{
    return val >> imm;
}

template <typename TValue>
SEQAN_FUNC bool
testAllZeros(TValue const &val)
{
    return val == (TValue)0;
}

template <typename TValue>
SEQAN_FUNC bool
testAllOnes(TValue const &val)
{
    return val == (TValue)-1;
}

template <typename TVector, typename TPos>
SEQAN_FUNC typename Value<TVector>::Type
shuffleVector(TVector const &vector, TPos const &index)
{
    return vector[index];
}



template <typename TWord>
SEQAN_FUNC void
_myersPreInit(MyersState<TWord> & state) {}

template <typename TSimdVector>
SEQAN_FUNC void
_myersPreInit(MyersStateSimd<TSimdVector> & state)
{
    typedef typename Value<TSimdVector>::Type TWord;

    TWord orMask = (TWord)1 << (BitsPerValue<TWord>::VALUE - 1);
    for (int i = 0; i < 5; ++i)
    {
        clear(state.patternOr[i]);
        state.patternOr[i][i] = orMask;
    }
}

template <typename TWord, typename TValue, typename TShift, typename TSmallAlphabet>
SEQAN_FUNC void 
_myersAdjustBitmask(MyersState<TWord> &state, TValue const value, TShift, TSmallAlphabet)
{
    for (unsigned i = 0; i < 5 /*ValueSize<TValue>::VALUE*/; ++i)
        state.patternBitmasks[i] >>= 1;
   state.patternBitmasks[ordValue(value)] |= (TWord)1 << (BitsPerValue<TWord>::VALUE - 1);
}


template <typename TSimdVector>
SEQAN_FUNC void
transpose(TSimdVector matrix[])
{
    // we need a look-up table to reverse the lowest 4 bits
    // in order to place the permute the transposed rows
    static const unsigned char bitRev[] = {0,8,4,12,2,10,6,14,1,9,5,13,3,11,7,15};

    // transpose a 16x16 byte matrix
    //
    // matrix =
    // A0 A1 A2 ... Ae Af
    // B0 B1 B2 ... Be Bf
    // ...
    // P0 P1 P2 ... Pe Pf
    register TSimdVector tmp1[16];
    for (int i = 0; i < 8; ++i)
    {
        tmp1[i]   = _mm_unpacklo_epi8(matrix[2*i], matrix[2*i+1]);
        tmp1[i+8] = _mm_unpackhi_epi8(matrix[2*i], matrix[2*i+1]);
    }
    // tmp1[0]  = A0 B0 A1 B1 ... A7 B7
    // tmp1[1]  = C0 D0 C1 D1 ... C7 D7
    // ...
    // tmp1[7]  = O0 P0 O1 P1 ... O7 P7
    // tmp1[8]  = A8 B8 A9 B9 ... Af Bf
    // ...
    // tmp1[15] = O8 P8 O9 P9 ... Of Pf
    register TSimdVector tmp2[16];
    for (int i = 0; i < 8; ++i)
    {
        tmp2[i]   = _mm_unpacklo_epi16(tmp1[2*i], tmp1[2*i+1]);
        tmp2[i+8] = _mm_unpackhi_epi16(tmp1[2*i], tmp1[2*i+1]);
    }
    // tmp2[0]  = A0 B0 C0 D0 ... A3 B3 C3 D3
    // tmp2[1]  = E0 F0 G0 H0 ... E3 F3 G3 H3
    // ...
    // tmp2[3]  = M0 N0 O0 P0 ... M3 N3 O3 P3
    // tmp2[4]  = A8 B8 C8 D8 ... Ab Bb Cb Db
    // ...
    // tmp2[7]  = M8 N8 O8 P8 ... Mb Nb Ob Pb
    // tmp2[8]  = A4 B4 C4 D4 ... A7 B7 C7 D7
    // ..
    // tmp2[12] = Ac Bc Cc Dc ... Af Bf Cf Df
    // ...
    // tmp2[15] = Mc Nc Oc Pc ... Mf Nf Of Pf
    for (int i = 0; i < 8; ++i)
    {
        tmp1[i]   = _mm_unpacklo_epi32(tmp2[2*i], tmp2[2*i+1]);
        tmp1[i+8] = _mm_unpackhi_epi32(tmp2[2*i], tmp2[2*i+1]);
    }
    // tmp1[0]  = A0 B0 .... H0 A1 B1 .... H1
    // tmp1[1]  = I0 J0 .... P0 I1 J1 .... P1
    // ...
    // tmp1[4]  = A0 B0 .... H0 A1 B1 .... H1
    // tmp1[1]  = I0 J0 .... P0 I1 J1 .... P1
    for (int i = 0; i < 8; ++i)
    {
        matrix[bitRev[i]]   = _mm_unpacklo_epi64(tmp1[2*i], tmp1[2*i+1]);
        matrix[bitRev[i+8]] = _mm_unpackhi_epi64(tmp1[2*i], tmp1[2*i+1]);
    }
}


template <typename TSimdVector, typename TIterators, typename TOffset>
SEQAN_FUNC void
_gatherInterleavedImpl(
    TSimdVector *textChars,
    TIterators *iters,
    TOffset ofs,
    unsigned blockSize,
    True)
{
    for (unsigned i = 0; i < LENGTH<TSimdVector>::VALUE; ++i)
    {
        for (unsigned j = 0; j < blockSize; ++j)
            textChars[i][j] = ordValue(getValue(iters[i] + ofs + j));
//      textChars[i] = *reinterpret_cast<TSimdVector*>(textIters[i] + ofs);
    }
    transpose(textChars);
}

template <typename T>
SEQAN_FUNC int 
whatTypeIsIt(T x)
{
	return x;
}

template <typename TValue, typename TIterators, typename TOffset>
SEQAN_FUNC void
_gatherInterleavedImpl(
    TValue *textChars,
    TIterators *iters,
    TOffset ofs,
    unsigned blockSize,
    False)
{
    for (unsigned j = 0; j < blockSize; ++j)
       textChars[j] = ordValue(getValue(iters[0] + ofs + j));
}

template <typename TSimdVector/*, unsigned BLOCK_SIZE*/, typename TIterators, typename TOffset>
SEQAN_FUNC void
_gatherInterleaved(
    TSimdVector *textChars,
    TIterators *iters,
    TOffset ofs,
    unsigned blockSize)
{
    _gatherInterleavedImpl(
        textChars,
        iters,
        ofs,
        blockSize,
        typename Is<SimdVectorConcept<TSimdVector> >::Type());
}

template <
    typename TTextIterators,
    typename TSize,
//    typename TNeedle,
    typename TNeedle2,
    typename TState
//    typename TSpec,
//    typename TFinderCSP,
//    typename TPatternCSP,
//    typename TFindBeginPatternSpec
>
SEQAN_FUNC bool
simdMyersBandedDiagonal(
    TTextIterators *textIters,
    TSize columns,
	TNeedle2 const & needle, 
//    PatternState_<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > & state)
    TState & state)
{
    typedef FindPrefix  TSpec;
    typedef TNeedle2    TNeedle;

    typedef typename TState::TSimdVector        TSimdVector;
    typedef typename TState::TPatternBitmasks   TPatternBitmasks;

//	typedef Pattern<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > TPattern;
//	typedef typename TPattern::TWord TWord;
    typedef typename Value<TSimdVector>::Type TWord;
	typedef typename Iterator<TNeedle2 const, Standard>::Type TIter;
	typedef typename Value<TNeedle>::Type TValue;

//    _myersPreInit(state, typename MyersSmallAlphabet_<TValue>::Type());
    TPatternBitmasks patternBitmasks;
    
//    clear(patternBitmasks); //CUDA workaround

	typename Size<TNeedle>::Type const ndlLength = length(needle);
	
	// Initialize row 0 either with zeros or increasing numbers
	// This can be realized using the following DP pattern and
	// assuming character mismatches at rows -1, -2,...
	// Thus we initialize the bitmasks and VN with 0.
	// VP depends on global/local alignment
	//
	//  0  1  2  3  4   -2 -2 -2 -2 -2   (row -2)
	//  0  1  2  3  4   -1 -1 -1 -1 -1   (row -1)
	//  0  1  2  3  4    0  0  0  0  0   (row  0)
	//  1                1
	//	    global           local
	//
	//  VP = 100...      VP = 111...
	//

    register TSimdVector VP, VN;
    fill(VP, (TWord)((MyersUkkonenHP0_<TSpec>::VALUE == 1)? (TWord)1 << ((int)BitsPerValue<TWord>::VALUE-1): MaxValue<TWord>::VALUE)); // HP[0]==1 <-> global, HP[0]==0 <-> local
    clear(VN);
	
	// Errors are counted along the lowest diagonal and the
	// lowest row of the band.
	// The errors in the top-left corner are 0.
	//
	// 0 * * * *
	//   x * * * *
	//     x * * * *
	//       x x x x x
	//
	//       |-------|
	//     diagWidth + 1 = 5
	//
	// diagWidth = length(container(finder)) + state.leftClip + state.rightClip - length(needle)
	
	TIter ndlIter = begin(needle, Standard());
	TIter ndlEnd = end(needle, Standard());

	// The errors along the diagonal can only increase or stay the same.
	// There is only the last row of length diagWidth where errors can decrease.
	// If errors exceeds cutOff it cannot reach maxErrors again.
	

    register TSimdVector errors;
    register TSimdVector cutOff;
    register TSimdVector bit0;
    fill(bit0, (TWord)1);

    columns += state.leftClip;
	if (columns > ndlLength)
	{
		fill(cutOff, (TWord)(state.maxErrors + (columns - ndlLength))); // clipping case *0
        clear(errors);
	}
//    else  // CUDA workaround
//    {
//		fill(cutOff, (TWord)state.maxErrors);                           // clipping case *0
//        fill(errors, (TWord)(ndlLength - columns));
//		ndlEnd = ndlIter + columns;                                     // clipping case *1
//	}

//    const int BLOCK_SIZE = LENGTH<TSimdVector>::VALUE;
    const unsigned BLOCK_SIZE = 1;

    TSimdVector textChars[BLOCK_SIZE];
    for (unsigned pos = 0; ndlIter != ndlEnd; pos += BLOCK_SIZE)    //CUDA workaround
    {
		//////////////////////////////////////////////////////////////////
        // PART 1: go down the parallelogram
		//////////////////////////////////////////////////////////////////

        unsigned CURRENT_BLOCK_SIZE = BLOCK_SIZE;
        if (BLOCK_SIZE > 1)
            if (BLOCK_SIZE > (unsigned)(ndlEnd - ndlIter))
                CURRENT_BLOCK_SIZE = (unsigned)(ndlEnd - ndlIter);    //CUDA workaround

		/////////////////////////
        // DIAGONAL MYERS CORE
		
		// VP/VN --> D0  (original Myers)
        _gatherInterleaved(textChars, textIters, pos, CURRENT_BLOCK_SIZE);

//        if (testAllOnes(textChars) != 0) return false;

        for (unsigned ofs = 0; ofs < CURRENT_BLOCK_SIZE; ++ofs, ++ndlIter) //CUDA workaround
        {
            // adjust bitmask
//            patternBitmasks = shiftRightLogical(patternBitmasks, 1) | state.patternOr[ordValue(getValue(ndlIter))];
            _myersAdjustBitmask(state, 2, pos + ofs, typename MyersSmallAlphabet_<TValue>::Type());	// CUDA workaround: getValue(ndlIter)

//            register TWord X = _myersGetBitmask(state, ordValue(*finder), shift, typename MyersSmallAlphabet_<TValue>::Type()) | VN;
            register TSimdVector X = shuffleVector(patternBitmasks, textChars[ofs]) | VN;
            register TSimdVector D0 = ((VP + (X & VP)) ^ VP) | X;
            
            // adjust errors corresponding to rightmost bit of D0
            errors += shiftRightLogical(~D0, BitsPerValue<TWord>::VALUE - 1);

            // return false if all alignments early exceeded the error cutOff
            TSimdVector x = errors > cutOff;
            if (testAllOnes(x) != 0) return false;

            // D0 --> HP/HN  (original Myers)
            register TSimdVector HN = VP & D0;
            register TSimdVector HP = VN | ~(VP | D0);

            // moving register down corresponds to shifting HP/HN up (right shift)
            // HP/HN --> shift --> VP/VN (modified Myers)
            X = shiftRightLogical(D0, 1);
            VN = X & HP;
            VP = HN | ~(X | HP);
        }
    }
    state.vp = VP;
    state.vn = VN;
//    state.patternBitmasks = patternBitmasks;
    __builtin_memcpy(state.patternBitmasks, patternBitmasks, sizeof(TPatternBitmasks));
    return true;
}

template <
    typename TTextIterators,
    typename TSize,
//    typename TNeedle,
    typename TNeedle2,
    typename TSimdVector
//    typename TSpec,
//    typename TFinderCSP,
//    typename TPatternCSP,
//    typename TFindBeginPatternSpec
>
SEQAN_FUNC bool
simdMyersBandedHorizontal(
    TTextIterators textIters[],
    TSize textLength,
	TNeedle2 const & needle, 
//    PatternState_<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > & state)
    MyersStateSimd<TSimdVector> & state)
{
    typedef FindPrefix  TSpec;
    typedef TNeedle2    TNeedle;

//	typedef Pattern<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > TPattern;
//	typedef typename TPattern::TWord TWord;
    typedef typename Value<TSimdVector>::Type TWord;
	typedef typename Iterator<TNeedle2 const, Standard>::Type TIter;
	typedef typename Value<TNeedle>::Type TValue;

//    _myersPreInit(state, typename MyersSmallAlphabet_<TValue>::Type());

    register TSimdVector VP = state.vp;
    register TSimdVector VN = state.vn;
    register TSimdVector errors = state.errors;
    register TSimdVector patternBitmasks = state.patternBitmasks;
    register TSimdVector maxErrors;
    register TSimdVector bit0;

    clear(maxErrors);
    fill(maxErrors, state.maxErrors);
    fill(bit0, 1);

    TSimdVector maxPos;
    clear(maxPos);

    TSize pos = length(needle);
    textLength -= pos;


    TSimdVector textChars[LENGTH<TSimdVector>::VALUE];
//    for (unsigned i = 0; i < LENGTH<TSimdVector>::VALUE; ++i)
//        textChars[8] = *reinterpret_cast<TSimdVector*>(textIters[i] + pos);
//
    for (unsigned i = 0; i < LENGTH<TSimdVector>::VALUE; ++i)
    {
        for (unsigned ofs = 0; ofs < LENGTH<TSimdVector>::VALUE; ++ofs)
            textChars[i][ofs] = ordValue(getValue(textIters[i] + pos + ofs));
//            textChars[i] = *reinterpret_cast<TSimdVector*>(textIters[i] + pos);
    }
    transpose(textChars);

	for (unsigned pos = 0; pos < textLength; ++pos)
    {
		// PART 2: go right

		// normal Myers

//        register TWord X = _myersGetBitmask(state, ordValue(*finder), shift, typename MyersSmallAlphabet_<TValue>::Type()) | VN;
        register TSimdVector X = shuffleVector(patternBitmasks, textChars[pos]) | VN;
        register TSimdVector D0 = ((VP + (X & VP)) ^ VP) | X;
        register TSimdVector HN = VP & D0;
        register TSimdVector HP = VN | ~(VP | D0);

        X = (HP << 1) | bit0;
        VN = X & D0;
        VP = (HN << 1) | ~(X | D0);

        // adjust errors corresponding to rightmost bit of D0
        errors += (HP >> (BitsPerValue<TWord>::VALUE - 2)) & bit0;
        errors -= (HN >> (BitsPerValue<TWord>::VALUE - 2)) & bit0;

        // shift bitmasks and states

//        if (errors <= maxErrors)
//        {
//            state.vp = VP;
//            state.vn = VN;
//            state.errors = errors;
//            return true;
//        }
        register TSimdVector posVector;
        fill(posVector, pos);

        unsigned char maxMask = (errors <= maxErrors);
        maxErrors = blend(maxErrors, errors, maxMask);
        maxPos = blend(maxPos, pos, maxMask);
    }
    return false;
}

//parallelFor(pattern, Patterns)
//{
//    parallelFor(textPos, textPositions)
//    {
//        find(pattern,text[textPos], 3);
//    }
//}

//parallelFor(_1, Patterns)
//{
//    parallelFor(_2, textPositions)
//    {
//        find(_1,conv[_2], 3);


} // namespace seqan

#endif
