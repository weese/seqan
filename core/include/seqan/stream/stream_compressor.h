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

#ifndef SEQAN_STREAM_STREAM_COMPRESSOR_H_
#define SEQAN_STREAM_STREAM_COMPRESSOR_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Classes
// ============================================================================

template <typename TAlgTag>
struct DefaultPageSize;

struct GZ {};
struct BGZF {};
struct BZ2 {};

struct CompressionContext<GZ>
{
    z_stream strm;
};

struct CompressionContext<BGZF>:
    CompressionContext<GZ>
{
    static const char header[] = {
        MagicHeader<BgzfFile>::VALUE[0], MagicHeader<BgzfFile>::VALUE[1], MagicHeader<BgzfFile>::VALUE[2],
        4, 0, 0, 0, 0, 0, -1, 6, 0, 'B', 'C', 2, 0, 0, 0
    };
    static const unsigned BLOCK_HEADER_LENGTH = sizeof(header); // == 18
    unsigned char headerPos;
};

template <>
struct DefaultPageSize<BGZF>
{
    const unsigned MAX_BLOCK_SIZE = 64 * 1024;
    const unsigned BLOCK_FOOTER_LENGTH = 8;
    const unsigned ZLIB_BLOCK_OVERHEAD = 5; // 5 bytes block overhead (see 3.2.4. at http://www.gzip.org/zlib/rfc-deflate.html)

    // Reduce the maximal input size, such that the compressed data 
    // always fits in one block even for level Z_NO_COMPRESSION.

    const unsigned VALUE = MAX_BLOCK_SIZE - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH - ZLIB_BLOCK_OVERHEAD;
};



template <typename TAlgTag>
Pager<TPager, Compress<TAlgTag> >
{
    TPager outPager;                // outbound pager
    PageTable<FixedSize> table;     // our page table

    Pager():
        table(DefaultPageSize<TAlgTag>::VALUE)
    {}

    Page & getPage (__int64 position)
    {
        Page *page;
        { 
            ScopedReadLock(table.lock);

            page = table[position];
            if (posInPage(position, page))                      // does the page exist yet?
                return page;
        }
        {
            ScopedWriteLock(table.lock);

            page = table[position];
            if (posInPage(position, page))                      // does the page exist yet?
                return page;

            page = new Page(table.rangeForPos(position));       // create new page
            reserve(page.data, table.pageSize);                 // allocate required memory
            table.insertPage(page);                             // insert page
            prevPage = prevPage(position);
        }
        return page;
    }
    
    void putPage (Page &page)
    {
        __int64 outPosition = 0;                                // compute start position in outbound pager
        if (page.range.begin != 0)
        {
            PageRange range = getPageRange(beginPosition(page.range) - 1);
            outPosition = endPosition(range);                   // wait for end position of the previous page
        }
        
        TCompressionContext ctx;
        initCompressionContext(ctx);

        Size<Page>::Type leftToCompress = length(page);
        while (leftToCompress != 0)
        {
            compress(toRange(
            Page &outPage = outPager.getPage(outPosition);

            auto handle = std::async(std::launch::async,
                                  parallel_sum<RAIter>, mid, end);
            compress
        }
    }
};

// ============================================================================
// Functions
// ============================================================================

compressInit(CompressionContext<GZ> &ctx)
{
    ctx.strm.zalloc = 0;
    ctx.strm.zfree = 0;
    int status = deflateInit2(&zs, Z_DEFAULT_COMPRESSION, Z_DEFLATED,
                              GZIP_WINDOW_BITS, Z_DEFAULT_MEM_LEVEL, Z_DEFAULT_STRATEGY);
    throw IOException("GZ deflateInit2() failed.");
}

compressInit(CompressionContext<BGZF> &ctx)
{
    compressInit(static_cast<CompressionContext<GZ> >(ctx));
    ctx.headerPos = 0;
}

template <typename TTarget, typename TSourceIterator>
inline typename Size<TTarget>::Type
compress(TTarget &target, TSourceIterator &source, CompressionContext<BGZF> &ctx)
{
    typedef typename Size<TTarget>::Type            TSize;
    typedef typename Chunk<TTarget>::Type           TTargetChunk;
    typedef typename Chunk<TSourceIterator>::Type   TSourceChunk;
    typedef typename Value<TSourceChunk>::Type      TSourceValue;

    TTargetChunk tChunk;
    TSourceChunk sChunk;

    if (ctx.headerPos < sizeof(ctx.header))
    {
        size_t headerLeft = sizeof(ctx.header) - ctx.headerPos;
        reserveChunk(target, headerLeft);

        tChunk = getChunk(target, Output());
        size_t size = std::min(headerLeft, length(tChunk));
        SEQAN_ASSERT_GT(size, 0u);

        std::copy(tChunk.begin, sChunk.begin, size);

        advanceChunk(target, size);
        ctx.headerPos += size;
        return size;
    }
    else
    {
        sChunk = getChunk(source, Input());
        tChunk = getChunk(target, Output());

        ctx.strm.next_in = static_cast<Bytef *>(sChunk.begin);
        ctx.strm.next_out = static_cast<Bytef *>(tChunk.begin);
        ctx.strm.avail_in = length(sChunk);
        ctx.strm.avail_out = length(tChunk);

        SEQAN_ASSERT_GT(ctx.strm.avail_out, 0u);

        int status = deflate(&zs, Z_NO_FLUSH);
        if (status != Z_OK)
            throw IOException("BGZF deflateInit2() failed.");

        source += length(sChunk) - ctx.strm.avail_in;
        size_t size = length(tChunk) - ctx.strm.avail_out;
        advanceChunk(target, size);
        return size;
    }

//    status = deflate(&zs, Z_FINISH);
//    bool rawDataTooBig = (status != Z_STREAM_END);
//
//    status = deflateEnd(&zs);
//    if (status != Z_OK)
//        throw IOException("BGZF deflateEnd() failed.");
//
//    if (!rawDataTooBig)
//    {
//        resize(page.raw, zs.total_out + BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH);
//        break;
//    }
}

#endif  // SEQAN_STREAM_STREAM_COMPRESSOR_H_

