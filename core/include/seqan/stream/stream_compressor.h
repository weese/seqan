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

template <typename TOutPager, typename TSpec>
struct Pager;

// ============================================================================
// Classes
// ============================================================================

template <typename TAlgTag>
struct Compress;

template <typename TAlgTag>
struct CompressionContext {};

template <typename TAlgTag>
struct DefaultPageSize;

#if SEQAN_HAS_ZLIB
#include <zlib.h>

template <>
struct CompressionContext<GZFile>
{
    z_stream strm;
};
#endif

template <>
struct CompressionContext<BgzfFile>:
    CompressionContext<GZFile>
{
    enum { BLOCK_HEADER_LENGTH = 18 };
    static const unsigned char header[BLOCK_HEADER_LENGTH];
    unsigned char headerPos;
};

const unsigned char CompressionContext<BgzfFile>::header[18] = {
    MagicHeader<BgzfFile>::VALUE[0], MagicHeader<BgzfFile>::VALUE[1], MagicHeader<BgzfFile>::VALUE[2],
    4, 0, 0, 0, 0, 0, -1, 6, 0, 'B', 'C', 2, 0, 0, 0
};

template <>
struct DefaultPageSize<BgzfFile>
{
    static const unsigned MAX_BLOCK_SIZE = 64 * 1024;
    static const unsigned BLOCK_FOOTER_LENGTH = 8;
    static const unsigned ZLIB_BLOCK_OVERHEAD = 5; // 5 bytes block overhead (see 3.2.4. at http://www.gzip.org/zlib/rfc-deflate.html)

    // Reduce the maximal input size, such that the compressed data 
    // always fits in one block even for level Z_NO_COMPRESSION.
    enum { BLOCK_HEADER_LENGTH = CompressionContext<BgzfFile>::BLOCK_HEADER_LENGTH };
    static const unsigned VALUE = MAX_BLOCK_SIZE - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH - ZLIB_BLOCK_OVERHEAD;
};


/*
template <typename TOutPager, typename TAlgTag>
Pager<TOutPager, Compress<TAlgTag> >
{
    TOutPager outPager;             // outbound pager
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
            Page &outPage = outPager.getPage(outPosition);

            auto handle = std::async(std::launch::async,
                                  parallel_sum<RAIter>, mid, end);
            compress
        }
    }
};
*/
// ============================================================================
// Functions
// ============================================================================

inline void
compressInit(CompressionContext<GZFile> &ctx)
{
    const int GZIP_WINDOW_BITS = -15; // no zlib header
    const int Z_DEFAULT_MEM_LEVEL = 8;

    ctx.strm.zalloc = NULL;
    ctx.strm.zfree = NULL;
    int status = deflateInit2(&ctx.strm, Z_DEFAULT_COMPRESSION, Z_DEFLATED,
                              GZIP_WINDOW_BITS, Z_DEFAULT_MEM_LEVEL, Z_DEFAULT_STRATEGY);
    if (status != Z_OK)
        throw IOException("GZFile deflateInit2() failed.");
}

inline void
compressInit(CompressionContext<BgzfFile> &ctx)
{
    compressInit(static_cast<CompressionContext<GZFile> &>(ctx));
    ctx.headerPos = 0;
}

template <typename TTarget, typename TSourceIterator>
inline typename Size<TTarget>::Type
compress(TTarget &target, TSourceIterator &source, CompressionContext<BgzfFile> &ctx)
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
        ctx.strm.avail_in = length(sChunk) * sizeof(TSourceValue);
        ctx.strm.avail_out = length(tChunk);

        SEQAN_ASSERT_GT(ctx.strm.avail_out, 0u);

        int status = deflate(&ctx.strm, Z_NO_FLUSH);
        if (status != Z_OK)
            throw IOException("BgzfFile deflateInit2() failed.");

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
//        throw IOException("BgzfFile deflateEnd() failed.");
//
//    if (!rawDataTooBig)
//    {
//        resize(page.raw, zs.total_out + BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH);
//        break;
//    }
}

template <typename TTarget, typename TSource>
inline typename Size<TTarget>::Type
compressAll(TTarget &target, TSource const &source, CompressionContext<BgzfFile> &ctx)
{
    typedef typename Value<TSource>::Type TSourceValue;

    const size_t BLOCK_HEADER_LENGTH = DefaultPageSize<BgzfFile>::BLOCK_HEADER_LENGTH;
    const size_t BLOCK_FOOTER_LENGTH = DefaultPageSize<BgzfFile>::BLOCK_FOOTER_LENGTH;

    // 1. COPY HEADER

    compressInit(ctx);

    // allocate output buffer
    resize(target, DefaultPageSize<BgzfFile>::MAX_BLOCK_SIZE, Exact());

    // copy header
    std::copy(&ctx.header[0], &ctx.header[BLOCK_HEADER_LENGTH], begin(target, Standard()));


    // 2. COMPRESS

    ctx.strm.next_in = static_cast<Bytef *>(static_cast<void *>(begin(source, Standard())));
    ctx.strm.next_out = static_cast<Bytef *>(static_cast<void const *>(begin(target, Standard())));
    ctx.strm.next_out += BLOCK_HEADER_LENGTH;
    ctx.strm.avail_in = length(source) * sizeof(TSourceValue);
    ctx.strm.avail_out = capacity(target) - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH;

    int status = deflate(&ctx.strm, Z_FINISH);
    if (status != Z_STREAM_END)
        throw IOException("Deflation failed. Compressed BgzfFile data is bigger than uncompressed data.");

    status = deflateEnd(&ctx.strm);
    if (status != Z_OK)
        throw IOException("BgzfFile deflateEnd() failed.");

    resize(target, capacity(target) - ctx.strm.avail_out);


    // 3. APPEND FOOTER

    // Set compressed length into buffer, compute CRC and write CRC into buffer.
    union {
        unsigned int i32;
        char raw[4];
    } tmp;

    tmp.i32 = crc32(0L, NULL, 0L);
    tmp.i32 = crc32(tmp.i32, static_cast<Bytef const *>(static_cast<void const *>(begin(source, Standard()))), length(source) * sizeof(TSourceValue));
    std::copy(&tmp.raw[0], &tmp.raw[4], begin(target, Standard()) + (length(target) - 8));

    tmp.i32 = length(source);
    std::copy(&tmp.raw[0], &tmp.raw[4], begin(target, Standard()) + (length(target) - 4));

    tmp.i32 = length(target) - 1;
    std::copy(&tmp.raw[0], &tmp.raw[4], begin(target, Standard()) + 16);

    return length(target);
}

}  // namespace seqan

#endif  // SEQAN_STREAM_STREAM_COMPRESSOR_H_

