/*
zipstream Library License:
--------------------------

The zlib/libpng License Copyright (c) 2003 Jonathan de Halleux.

This software is provided 'as-is', without any express or implied warranty. In no event will the authors be held liable for any damages arising from the use of this software.

Permission is granted to anyone to use this software for any purpose, including commercial applications, and to alter it and redistribute it freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.

2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.

3. This notice may not be removed or altered from any source distribution

Author: Jonathan de Halleux, dehalleux@pelikhan.com, 2003   (original zlib stream)
Author: David Weese, dave.weese@gmail.com, 2014             (extension to parallel block-wise compression in bgzf format)
*/

#ifndef BGZFSTREAM_HPP
#define BGZFSTREAM_HPP

#include <vector>
#include <iostream>
#include <algorithm>
#include <zlib.h>
#include "zutil.h"

namespace seqan {

/// default gzip buffer size,
/// change this to suite your needs
const size_t default_buffer_size = 4096;

const unsigned BGZF_MAX_BLOCK_SIZE = 64 * 1024;
const unsigned BGZF_BLOCK_HEADER_LENGTH = 18;
const unsigned BGZF_BLOCK_FOOTER_LENGTH = 8;
const unsigned ZLIB_BLOCK_OVERHEAD = 5; // 5 bytes block overhead (see 3.2.4. at http://www.gzip.org/zlib/rfc-deflate.html)

// Reduce the maximal input size, such that the compressed data 
// always fits in one block even for level Z_NO_COMPRESSION.
const unsigned BGZF_BLOCK_SIZE = BGZF_MAX_BLOCK_SIZE - BGZF_BLOCK_HEADER_LENGTH - BGZF_BLOCK_FOOTER_LENGTH - ZLIB_BLOCK_OVERHEAD;

/// Compression strategy, see bgzf doc.
enum EStrategy
{
	StrategyFiltered = 1,
	StrategyHuffmanOnly = 2,
	DefaultStrategy = 0
};


/** \brief A stream decorator that takes raw input and zips it to a ostream.

The class wraps up the inflate method of the bgzf library 1.1.4 http://www.gzip.org/bgzf/
*/
template<
	typename Elem, 
	typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = unsigned char,
    typename ByteAT = std::allocator<ByteT>
>	
class basic_bgzf_streambuf : public std::basic_streambuf<Elem, Tr> 
{
public:
	typedef std::basic_ostream<Elem, Tr>& ostream_reference;
    typedef ElemA char_allocator_type;
	typedef ByteT byte_type;
    typedef ByteAT byte_allocator_type;
	typedef byte_type* byte_buffer_type;
	typedef typename Tr::char_type char_type;
	typedef typename Tr::int_type int_type;

    struct Concatter
    {
        ostream_reference ostream;
        Mutex lock;
        unsigned waitForKey;
        unsigned nextKey;

        Concatter(ostream_reference ostream) :
            ostream(ostream),
            lock(false)
        {}
    };

    Concatter concatter;

    struct BgzfThreadContext
    {
        typedef std::vector<byte_type, byte_allocator_type> byte_vector_type;
        typedef std::vector<char_type, char_allocator_type> char_vector_type;

        struct BgzfCompressor
        {
            BgzfThreadContext *threadCtx;

            void operator()()
            {
                size_t outputLen = _compressBlock(
                    &threadCtx->outputBuffer[0], capacity(threadCtx->outputBuffer),
                    &threadCtx->buffer[0], threadCtx->size, threadCtx->ctx);
                while (true)
                {
                    ScopedLock<Mutex> scopedLock(threadCtx->concatter->lock);
                    if (threadCtx->concatter->waitForKey == threadCtx->waitForKey)
                    {
                        threadCtx->concatter->ostream.write(
                            (const char_type*) &(threadCtx->outputBuffer[0]),
                            outputLen);
                        threadCtx->concatter->waitForKey = threadCtx->waitForKey + 1;
                        break;
                    }
                }
            }
        };

        byte_vector_type                outputBuffer;
        char_vector_type                buffer;
        size_t                          size;
        CompressionContext<BgzfFile>    ctx;
        Thread<BgzfCompressor>          compress;

        Concatter *concatter;
        unsigned waitForKey;

        BgzfThreadContext():
            outputBuffer(BGZF_MAX_BLOCK_SIZE, 0),
            buffer(BGZF_BLOCK_SIZE / sizeof(char_type), 0)
        {}
    };

    BgzfThreadContext   *threadCtx;
    BgzfThreadContext   *ctx;
    size_t              thread;
    size_t              threads;

    /** Construct a zip stream
     * More info on the following parameters can be found in the bgzf documentation.
     */
    basic_bgzf_streambuf(ostream_reference ostream_,
                         size_t threads = 4) :
		concatter(ostream_),
        thread(0),
        threads(threads)
    {
        concatter.waitForKey = 0;
        concatter.nextKey = 0;
        threadCtx = new BgzfThreadContext[threads];
        for (unsigned i = 0; i < threads; ++i)
        {
            threadCtx[i].concatter = &concatter;
            threadCtx[i].compress.worker.threadCtx = &threadCtx[i];
        }
        ctx = &threadCtx[thread];
		this->setp(&(ctx->buffer[0]), &(ctx->buffer[ctx->buffer.size() - 1]));
    }
	
	~basic_bgzf_streambuf()
    {
		flush();
		concatter.ostream.flush();
        delete[] threadCtx;
    }

	int sync()
    {
		if (this->pptr() > this->pbase())
		{
			int c = overflow(EOF);
			if (c == EOF)
				return -1;
        }
        return 0;
    }

    bool compressBuffer(size_t size)
    {
        ctx->size = size;
        ctx->waitForKey = concatter.nextKey++;
        run(ctx->compress);

        thread = (thread + 1) % threads;
        ctx = &threadCtx[thread];
        waitFor(ctx->compress);
    }

    int_type overflow(int_type c)
    {
        int w = static_cast<int>(this->pptr() - this->pbase());
        if (c != EOF)
        {
            *this->pptr() = c;
            ++w;
        }
        if (compressBuffer(w))
        {
            this->setp(&(ctx->buffer[0]), &(ctx->buffer[ctx->buffer.size() - 1]));
            return c;
        }
        else
        {
            return EOF;
        }
    }

	/** flushes the zip buffer and output buffer.

	This method should be called at the end of the compression. Calling flush multiple times, will lower the
	compression ratio.
	*/
	std::streamsize flush()
    {
//        std::streamsize totalWrittenByteSize = compressBuffer();
//		concatter.ostream.flush();
//		return totalWrittenByteSize;
    }

	/// returns a reference to the output stream
	ostream_reference get_ostream() const	{ return concatter.ostream; };
	/// returns the latest bgzf error status
	int get_zerr() const					{ return m_err; };
private:
	size_t fill_input_buffer();

    int m_err;

};

/** \brief A stream decorator that takes compressed input and unzips it to a istream.

The class wraps up the deflate method of the bgzf library 1.1.4 http://www.gzip.org/bgzf/
*/
template<
	typename Elem, 
	typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = unsigned char,
    typename ByteAT = std::allocator<ByteT>
>	
class basic_unbgzf_streambuf : 
	public std::basic_streambuf<Elem, Tr> 
{
public:
	typedef std::basic_istream<Elem, Tr>& istream_reference;
    typedef ElemA char_allocator_type;
	typedef ByteT byte_type;
    typedef ByteAT byte_allocator_type;
	typedef byte_type* byte_buffer_type;
	typedef typename Tr::char_type char_type;
	typedef typename Tr::int_type int_type;
	typedef std::vector<byte_type, byte_allocator_type > byte_vector_type;
	typedef std::vector<char_type, char_allocator_type > char_vector_type;

     /** Construct a unzip stream
     * More info on the following parameters can be found in the bgzf documentation.
     */
	 basic_unbgzf_streambuf(istream_reference istream_);
	
	~basic_unbgzf_streambuf();

    int_type underflow();


	/// returns the compressed input istream
	istream_reference get_istream()	{	return m_istream;};
	/// returns the bgzf stream structure
	z_stream& get_bgzf_stream()		{	return m_bgzf_stream;};
	/// returns the latest bgzf error state
	int get_zerr() const					{	return m_err;};
	/// returns the crc of the uncompressed data so far 
	long get_crc() const					{	return m_crc;};
	/// returns the number of uncompressed bytes
	long get_out_size() const				{	return m_bgzf_stream.total_out;};
	/// returns the number of read compressed bytes
	long get_in_size() const				{	return m_bgzf_stream.total_in;};
private:
	void put_back_from_bgzf_stream();
	std::streamsize unbgzf_from_stream( char_type*, std::streamsize);

	size_t fill_input_buffer();

	istream_reference m_istream;
	z_stream m_bgzf_stream;
    int m_err;
	byte_vector_type m_input_buffer;
	char_vector_type m_buffer; 
	long m_crc;
};

/*! \brief Base class for zip ostreams

Contains a basic_bgzf_streambuf.
*/
template<
	typename Elem, 
	typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = unsigned char,
    typename ByteAT = std::allocator<ByteT>
>	
class basic_bgzf_ostreambase : virtual public std::basic_ios<Elem,Tr>
{
public:
	typedef std::basic_ostream<Elem, Tr>& ostream_reference;
	typedef basic_bgzf_streambuf<
        Elem,
        Tr,
        ElemA,
        ByteT,
        ByteAT
        > bgzf_streambuf_type;

    /** Construct a zip stream
     * More info on the following parameters can be found in the bgzf documentation.
     */
	basic_bgzf_ostreambase(ostream_reference ostream_)
		: m_buf(ostream_)
	{
		this->init(&m_buf );
	};
	
	/// returns the underlying zip ostream object
	bgzf_streambuf_type* rdbuf() { return &m_buf; };

	/// returns the bgzf error state
	int get_zerr() const					{	return m_buf.get_err();};
	/// returns the uncompressed data crc
	long get_crc() const					{	return m_buf.get_crc();};
	/// returns the compressed data size
	long get_out_size() const				{	return m_buf.get_out_size();};
	/// returns the uncompressed data size
	long get_in_size() const				{	return m_buf.get_in_size();};
private:
	bgzf_streambuf_type m_buf;
};

/*! \brief Base class for unzip istreams

Contains a basic_unbgzf_streambuf.
*/
template<
	typename Elem, 
	typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = unsigned char,
    typename ByteAT = std::allocator<ByteT>
>
class basic_bgzf_istreambase : virtual public std::basic_ios<Elem,Tr>
{
public:
	typedef std::basic_istream<Elem, Tr>& istream_reference;
	typedef basic_unbgzf_streambuf<
        Elem,
        Tr,
        ElemA,
        ByteT,
        ByteAT
        > unbgzf_streambuf_type;

	basic_bgzf_istreambase(istream_reference ostream_)
		: m_buf(ostream_)
	{
		this->init(&m_buf );
	};
	
	/// returns the underlying unzip istream object
	unbgzf_streambuf_type* rdbuf() { return &m_buf; };

	/// returns the bgzf error state
	int get_zerr() const					{	return m_buf.get_zerr();};
	/// returns the uncompressed data crc
	long get_crc() const					{	return m_buf.get_crc();};
	/// returns the uncompressed data size
	long get_out_size() const				{	return m_buf.get_out_size();};
	/// returns the compressed data size
	long get_in_size() const				{	return m_buf.get_in_size();};
private:
	unbgzf_streambuf_type m_buf;
};

/*! \brief A zipper ostream

This class is a ostream decorator that behaves 'almost' like any other ostream.

At construction, it takes any ostream that shall be used to output of the compressed data.

When finished, you need to call the special method zflush or call the destructor 
to flush all the intermidiate streams.

Example:
\code
// creating the target zip string, could be a fstream
ostringstream ostringstream_;
// creating the zip layer
bgzf_ostream zipper(ostringstream_);

	
// writing data	
zipper<<f<<" "<<d<<" "<<ui<<" "<<ul<<" "<<us<<" "<<c<<" "<<dum;
// zip ostream needs special flushing...
zipper.zflush();
\endcode
*/
template<
	typename Elem, 
	typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = unsigned char,
    typename ByteAT = std::allocator<ByteT>
>	
class basic_bgzf_ostream : 
	public basic_bgzf_ostreambase<Elem,Tr,ElemA,ByteT,ByteAT>, 
	public std::basic_ostream<Elem,Tr>
{
public:
	typedef basic_bgzf_ostreambase<
        Elem,Tr,ElemA,ByteT,ByteAT> bgzf_ostreambase_type;
	typedef std::basic_ostream<Elem,Tr> ostream_type;
    typedef ostream_type& ostream_reference;

	/** Constructs a zipper ostream decorator
	 *
	 * \param ostream_ ostream where the compressed output is written

	 When is_gbgzf_ is true, a gzip header and footer is automatically added.
	 */
	basic_bgzf_ostream(ostream_reference ostream_)
	: 
		bgzf_ostreambase_type(ostream_),
		ostream_type(this->rdbuf())
	{}

	/// flush inner buffer and zipper buffer
	basic_bgzf_ostream<Elem,Tr>& zflush()	
	{	
		this->flush(); this->rdbuf()->flush(); return *this;
	};

private:
    static void put_long(ostream_reference out_, unsigned long x_);
};

/*! \brief A zipper istream

This class is a istream decorator that behaves 'almost' like any other ostream.

At construction, it takes any istream that shall be used to input of the compressed data.

Simlpe example:
\code
// create a stream on zip string
istringstream istringstream_( ostringstream_.str());
// create unzipper istream
bgzf_istream unzipper( istringstream_);

// read and unzip
unzipper>>f_r>>d_r>>ui_r>>ul_r>>us_r>>c_r>>dum_r;
\endcode
*/
template<
	typename Elem, 
	typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = unsigned char,
    typename ByteAT = std::allocator<ByteT>
>
class basic_bgzf_istream : 
	public basic_bgzf_istreambase<Elem,Tr,ElemA,ByteT,ByteAT>, 
	public std::basic_istream<Elem,Tr>
{
public:
	typedef basic_bgzf_istreambase<
        Elem,Tr,ElemA,ByteT,ByteAT> bgzf_istreambase_type;
	typedef std::basic_istream<Elem,Tr> istream_type;
    typedef istream_type& istream_reference;
	typedef unsigned char byte_type;

	/** Construct a unzipper stream
	 *
	 * \param istream_ input buffer
	 */
	basic_bgzf_istream(istream_reference istream_)
	  : 
		bgzf_istreambase_type(istream_),
		istream_type(this->rdbuf()),
		m_is_gzip(false),
		m_gbgzf_data_size(0)
	{
 	      if (this->rdbuf()->get_zerr()==Z_OK)
			  check_header();
	};

	/// returns true if it is a gzip file
	bool is_gzip() const				{	return m_is_gzip;};
	/// return data size check
	bool check_data_size() const		{	return this->get_out_size() == m_gbgzf_data_size;};

	/// return the data size in the file 
	long get_gbgzf_data_size() const		{	return m_gbgzf_data_size;};
protected:
    static void read_long(istream_reference in_, unsigned long& x_);

	int check_header();
	bool m_is_gzip;
	unsigned long m_gbgzf_data_size;
};

/// A typedef for basic_bgzf_ostream<char>
typedef basic_bgzf_ostream<char> bgzf_ostream;
/// A typedef for basic_bgzf_ostream<wchar_t>
typedef basic_bgzf_ostream<wchar_t> bgzf_wostream;
/// A typedef for basic_bgzf_istream<char>
typedef basic_bgzf_istream<char> bgzf_istream;
/// A typedef for basic_bgzf_istream<wchart>
typedef basic_bgzf_istream<wchar_t> bgzf_wistream;

}  // namespace seqan

#include "bgzfstream.ipp"

#endif

