#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>
#include <algorithm>
#include <iostream>
#include <cstdlib>

#include "cuPrintf.cu"

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/find.h>

//#undef SEQAN_FUNC
//#define SEQAN_FUNC inline __device__ __host__

#include <seqan/sequence/adapt_thrust_vector.h>
#include <seqan/basic/basic_simd_vector.h>
#include <seqan/misc/misc_view.h>

#include "testset.h"
#include "myers.h"

using namespace seqan;


template <
    // typename TText,
	typename TPtrBegin,
	typename TPtrEnd,
    typename TPatterns,
    typename TLimits,
    typename TPos>
__global__ void verifyOnGPU(
    // View<TText, IteratorView> text,
	TPtrBegin textBegin, TPtrEnd textEnd,
    View<TPatterns, IteratorView> patterns,
    View<TLimits, IteratorView> limits,
    View<TPos, IteratorView> pos)
{	
	typedef View<TPtrBegin, IteratorView> TText;
	TText text(textBegin, textEnd);    
    
    typedef MyersState<unsigned char> TState;
    TState state;

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
	cuPrintf("index=%i", idx);
/*
    _myersPreInit(state);
    View<TPatterns> pattern(
        begin(patterns, Standard()) + limits[idx],
        begin(patterns, Standard()) + limits[idx+1]);
    
    typename Iterator<TText, Standard>::Type textIter = begin(text, Standard()) + pos[idx];

    simdMyersBandedDiagonal(
        &textIter, 
        length(pattern) + 3, 
        pattern, 
        state);
*/
}



void bench()
{
    typedef DnaString TText;
    typedef DnaString TPattern;
    typedef StringSet<TPattern, Owner<ConcatDirect<> > > TPatterns;
    
//    const unsigned verifications = 1u << 19;
//    const unsigned textLen = 1000000000;
    const unsigned verifications = 1u << 4;
    const unsigned textLen = 100000;
    const unsigned patternLen = 100;

    BenchParams params;
    params.patternAmount = verifications;       // simulate patterns
    params.maxErrors = 3;                       // currently unused
    params.findPositiveFrac=1.0;

    TText text;
    TPattern pattern;                           // search one pattern in multiple texts
    TPatterns patterns;                         // search multiple patterns
    String<unsigned> pos;


    ///////////////////////////////////////////////////////////////////////////////////
    // 0. DATA GENERATION
    ///////////////////////////////////////////////////////////////////////////////////

    {
        resize(text, textLen);                  // simulate a text of 1Gb
        resize(pattern, patternLen);            // simulate a text of 10Mb

        double start = cpuTime();
        generateRandomString(params, text);
        generateRandomString(params, pattern);
        generatePatterns(params, patterns, pos, text, patternLen);
//        generateMatches(params, pos, pattern, text, verifications);

        std::cout << "Data generation took\t" << cpuTime() - start << " seconds." << std::endl;        
    }


    ///////////////////////////////////////////////////////////////////////////////////
    // 1. BENCHMARK STANDARD CPU VERSION
    ///////////////////////////////////////////////////////////////////////////////////
    
    {
        typedef PatternState_<TPattern, Myers< AlignTextBanded<FindPrefix, NMatchesN_, NMatchesN_>, True, void> > TState;
        TState state;

        double start = cpuTime();
        for (int l=0;l<2;++l)
        {
            for (unsigned i = 0; i < verifications; ++i)
            {
        //    std::cout << pos[i]<< std::endl;
                Infix<TText>::Type textInfix = infix(text, _max((int)pos[i] - (int)params.maxErrors, 0), pos[i] + patternLen + params.maxErrors);
                Finder<Infix<TText>::Type> finder(textInfix);

//                _stateInit(finder, pattern, state);
                _stateInit(finder, patterns[i], state);
        //        _findMyersSmallPatternsBanded(finder, needle, state);
            }
        }
        std::cout << "Sequential search took\t" << cpuTime() - start << " seconds." << std::endl;


        typedef Iterator<TText, Standard>::Type TTextIter;
        String<TTextIter> textIter;

        resize(textIter, verifications);
        TTextIter textBegin = begin(text, Standard());
        for (unsigned i = 0; i < verifications; ++i)
            textIter[i] = textBegin + pos[i];
    }


    ///////////////////////////////////////////////////////////////////////////////////
    // 2. BENCHMARK SIMD CPU VERSION
    ///////////////////////////////////////////////////////////////////////////////////
    
    {

//        double start = cpuTime();
//        for (int l=0;l<2;++l)
//        {
//    //        typedef MyersStateSimd<unsigned char> TSimdState;
//            typedef MyersState<unsigned char> TSimdState;
//            TSimdState simdState;
//            _myersPreInit(simdState);
//            for (unsigned ofs = 0; ofs < verifications; ofs += 1)
//            {
//                // TODO prefetch next batch
//    //            if (ofs + LENGTH<SimdVector<unsigned char>::Type>::VALUE < verifications)
//    //                for (unsigned p = 0; p < LENGTH<SimdVector<unsigned char>::Type>::VALUE; ++p)
//    //                {
//    //                    _mm_prefetch(textIter[ofs + p], _MM_HINT_T0);
//    //                    _mm_prefetch(textIter[ofs + p] + 64, _MM_HINT_T0);
//    //                }
//
//                simdMyersBandedDiagonal(begin(textIter, Standard()) + ofs, patternLen + params.maxErrors, pattern, simdState);
//    //            simdMyersBandedHorizontal(begin(textIter, Standard()) + ofs, patternLen + params.maxErrors, pattern, simdState);
//            }
//        }
//        std::cout << "SIMD search took      \t" << cpuTime() - start << " seconds." << std::endl;
    }
    
    
    ///////////////////////////////////////////////////////////////////////////////////
    // 3. BENCHMARK GPU VERSION
    ///////////////////////////////////////////////////////////////////////////////////
    
    {

        // we need to use host_vectors as the space of SeqAn's iterators is not recognized
        thrust::host_vector<Dna> host_text(length(text));
        thrust::host_vector<unsigned char> host_patterns(length(patterns.concat));
        thrust::host_vector<unsigned> host_limits(length(patterns.limits));
        thrust::host_vector<unsigned> host_pos(length(pos));

        thrust::copy(begin(text, Standard()),               end(text, Standard()),              host_text.begin());
        thrust::copy(begin(patterns.concat, Standard()),    end(patterns.concat, Standard()),   host_patterns.begin());
        thrust::copy(begin(patterns.limits, Standard()),    end(patterns.limits, Standard()),   host_limits.begin());
        thrust::copy(begin(pos, Standard()),                end(pos, Standard()),               host_pos.begin());

        // copy text, patterns, and jobs to GPU
        thrust::device_vector<Dna> gpu_text(host_text);
        thrust::device_vector<unsigned char> gpu_patterns(host_patterns);
        thrust::device_vector<unsigned> gpu_limits(host_limits);
        thrust::device_vector<unsigned> gpu_pos(host_pos);

//        thrust::device_vector<unsigned char> gpu_text(length(text));
//        thrust::device_vector<unsigned char> gpu_patterns(length(patterns.concat));
//        thrust::device_vector<unsigned> gpu_limits(length(patterns.limits));
//        thrust::device_vector<unsigned> gpu_pos(length(pos));

//        thrust::copy(begin(text, Standard()),               end(text, Standard()),              gpu_text.begin());
//        thrust::copy(begin(patterns.concat, Standard()),    end(patterns.concat, Standard()),   gpu_patterns.begin());
//        thrust::copy(begin(patterns.limits, Standard()),    end(patterns.limits, Standard()),   gpu_limits.begin());
//        thrust::copy(begin(pos, Standard()),                end(pos, Standard()),               gpu_pos.begin());

        double start = cpuTime();
        for (int l=0;l<2;++l)
        {
            int block_size = 4;
            int n_blocks = (verifications + block_size - 1) / block_size;

			cudaPrintfInit();
            verifyOnGPU<<< n_blocks, block_size >>>(
				toRawView(gpu_text)._begin, 
				toRawView(gpu_text)._end, 
                toRawView(gpu_patterns), 
                toRawView(gpu_limits), 
                toRawView(gpu_pos));
            cudaDeviceSynchronize();
			cudaPrintfDisplay(stdout, true);
			cudaPrintfEnd();
        }
        std::cout << "GPU search took       \t" << cpuTime() - start << " seconds." << std::endl;
    }

	return;
}


int main(void)
{ 
  // generate random data serially
  thrust::host_vector<int> h_vec(100);
  std::generate(h_vec.begin(), h_vec.end(), rand);
  int x = thrust::reduce(h_vec.begin(), h_vec.end(), 0, thrust::plus<int>());
  std::cout << x << std::endl;

  // transfer to device and compute sum
  thrust::device_vector<int> d_vec = h_vec;
  x = thrust::reduce(d_vec.begin(), d_vec.end(), 0, thrust::plus<int>());
  std::cout << x << std::endl;

  bench();

  return 0;
}
