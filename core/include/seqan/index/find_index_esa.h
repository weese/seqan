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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_INDEX_ESA_FIND_H
#define SEQAN_HEADER_INDEX_ESA_FIND_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// ESA finders

	template < typename TText, typename TSpec >
	struct DefaultFinder< Index<TText, IndexEsa<TSpec> > > {
        typedef EsaFindMlr Type;	// standard suffix array finder is mlr-heuristic
    };

//////////////////////////////////////////////////////////////////////////////
// _findFirstIndex implementation

	template < typename TText, typename TSpec, typename TSpecFinder, typename TPattern >
	inline void 
	_findFirstIndex(
		Finder< Index<TText, TSpec>, TSpecFinder > &finder,
		TPattern const &pattern,
		EsaFindMlr const)
	{
		Index<TText, TSpec> &index = haystack(finder);
		indexRequire(index, EsaSA());
		finder.range = equalRangeSAIterator(indexText(index), indexSA(index), pattern);
	}

	template < typename TText, typename TSpec, typename TSpecFinder, typename TPattern >
	inline void 
	_findFirstIndex(
		Finder< Index<TText, TSpec>, TSpecFinder > &finder,
		TPattern const &pattern,
		EsaFindLcpe const)
	{
		Index<TText, TSpec> &index = haystack(finder);
		indexRequire(index, EsaSA());
		indexRequire(index, EsaLcpe());
		finder.range = equalRangeLcpeIterator(indexText(index), indexSA(index), indexLcpe(index), pattern);
	}

	template < typename TText, typename TSpec, typename TSpecFinder, typename TPattern >
	inline void
	_findFirstIndex(
		Finder< Index<TText, TSpec>, TSpecFinder > &finder,
		TPattern const &pattern,
		FinderSTree const)
	{
        typedef Index<TText, TSpec>                             TIndex;
		typedef typename Fibre<TIndex, FibreSA>::Type			TSA;
		typedef typename Iterator<TSA const, Standard>::Type	TIterator;
        
		TIndex &index = haystack(finder);
        typename Iterator<TIndex, TopDown<EmptyEdges> >::Type it(index);
        TIterator saIt = begin(indexSA(index), Standard());
        if (goDown(it, pattern))
        {
            Pair<typename Size<TIndex>::Type> rng = range(it);
			finder.range.i1 = saIt + rng.i1;
			finder.range.i2 = saIt + rng.i2;
        } else
        {
            finder.range.i1 = saIt;
            finder.range.i2 = saIt;
        }
	}


}
#endif

