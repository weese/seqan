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
#ifndef BGZFSTREAM_IPP
#define BGZFSTREAM_IPP

#include "bgzfstream.hpp"
#include <sstream>

namespace seqan {

namespace detail{
	const int gz_magic[2] = {0x1f, 0x8b}; /* gzip magic header */

	/* gzip flag byte */
	const int gz_ascii_flag =  0x01; /* bit 0 set: file probably ascii text */
	const int gz_head_crc    = 0x02; /* bit 1 set: header CRC present */
	const int gz_extra_field = 0x04; /* bit 2 set: extra field present */
	const int gz_orig_name  =  0x08; /* bit 3 set: original file name present */
	const int gz_comment    =  0x10; /* bit 4 set: file comment present */
	const int gz_reserved   =  0xE0; /* bits 5..7: reserved */	

}

	template<
		typename Elem, 
		typename Tr,
		typename ElemA,
		typename ByteT,
		typename ByteAT
	>
    basic_unbgzf_streambuf<
		Elem,Tr,ElemA,ByteT,ByteAT
		>::basic_unbgzf_streambuf(istream_reference istream_)
	:  
		m_istream(istream_),
		m_input_buffer(BGZF_MAX_BLOCK_SIZE),
		m_buffer(BGZF_MAX_BLOCK_SIZE)
	{
		this->setg( &(m_buffer[0])+4,     // beginning of putback area
		    &(m_buffer[0])+4,     // read position
	        &(m_buffer[0])+4);    // end position    
	}

	template<
		typename Elem, 
		typename Tr,
		typename ElemA,
		typename ByteT,
		typename ByteAT
	>
    size_t basic_unbgzf_streambuf<
		Elem,Tr,ElemA,ByteT,ByteAT
		>::fill_input_buffer()
	{
		m_bgzf_stream.next_in=&(m_input_buffer[0]);
		m_istream.read( 
			(char_type*)(&(m_input_buffer[0])), 
			static_cast<std::streamsize>(m_input_buffer.size()/sizeof(char_type)) 
			); 
		return m_bgzf_stream.avail_in=m_istream.gcount()*sizeof(char_type);
	}


	template<
		typename Elem, 
		typename Tr,
		typename ElemA,
		typename ByteT,
		typename ByteAT
	>
    basic_unbgzf_streambuf<
		Elem,Tr,ElemA,ByteT,ByteAT
		>::~basic_unbgzf_streambuf()
	{
		inflateEnd(&m_bgzf_stream);
	}

	template<
		typename Elem, 
		typename Tr,
		typename ElemA,
		typename ByteT,
		typename ByteAT
	>
    typename basic_unbgzf_streambuf<
		Elem,Tr,ElemA,ByteT,ByteAT
		>::int_type 
		basic_unbgzf_streambuf<
			Elem,Tr,ElemA,ByteT,ByteAT
			>::underflow() 
	{ 
		if ( this->gptr() && ( this->gptr() < this->egptr()))
			return * reinterpret_cast<unsigned char *>( this->gptr());
     
       int n_putback = static_cast<int>(this->gptr() - this->eback());
       if ( n_putback > 4)
          n_putback = 4;
       memcpy( 
			&(m_buffer[0]) + (4 - n_putback), 
			this->gptr() - n_putback,
			n_putback*sizeof(char_type)
			);
  
	   int num = unbgzf_from_stream( 
		   &(m_buffer[0])+4, 
		   static_cast<std::streamsize>((m_buffer.size()-4)*sizeof(char_type))
		   );
        if (num <= 0) // ERROR or EOF
           return EOF;
    
        // reset buffer pointers
        this->setg( 
              &(m_buffer[0]) + (4 - n_putback),   // beginning of putback area
              &(m_buffer[0]) + 4,                 // read position
              &(m_buffer[0]) + 4 + num);          // end of buffer
    
         // return next character
         return* reinterpret_cast<unsigned char *>( this->gptr());
     }

	template<
		typename Elem, 
		typename Tr,
		typename ElemA,
		typename ByteT,
		typename ByteAT
	>
    std::streamsize basic_unbgzf_streambuf<
		Elem,Tr,ElemA,ByteT,ByteAT
		>::unbgzf_from_stream( 
			char_type* buffer_, 
			std::streamsize buffer_size_
			)
	{
		m_bgzf_stream.next_out=(byte_buffer_type)buffer_;
		m_bgzf_stream.avail_out=static_cast<uInt>(buffer_size_*sizeof(char_type));
		size_t count =m_bgzf_stream.avail_in;

		do
		{
			if (m_bgzf_stream.avail_in==0)
				count=fill_input_buffer();

			if (m_bgzf_stream.avail_in)
			{
				m_err = inflate( &m_bgzf_stream,  Z_SYNC_FLUSH );
			}
		} while (m_err==Z_OK && m_bgzf_stream.avail_out != 0 && count != 0);

		// updating crc
		m_crc = crc32( 
			m_crc, 
			(byte_buffer_type)buffer_,
			buffer_size_ - m_bgzf_stream.avail_out/sizeof(char_type)
			);	
		std::streamsize n_read = buffer_size_ - m_bgzf_stream.avail_out/sizeof(char_type);
		
		// check if it is the end
		if (m_err==Z_STREAM_END)
			put_back_from_bgzf_stream();				
		
		return n_read;
	}
	
	template<
		typename Elem, 
		typename Tr,
		typename ElemA,
		typename ByteT,
		typename ByteAT
	>
	void basic_unbgzf_streambuf<
		Elem,Tr,ElemA,ByteT,ByteAT
		>::put_back_from_bgzf_stream()
	{
		if (m_bgzf_stream.avail_in==0)
			return;

		m_istream.clear( std::ios::goodbit );
		m_istream.seekg(
			-static_cast<int>(m_bgzf_stream.avail_in),
			std::ios_base::cur
			);

		m_bgzf_stream.avail_in=0;
	}

    template<
        typename Elem, 
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    int _checkGZHeader(basic_unbgzf_streambuf<Elem, Tr, ElemA, ByteT, ByteAT> *buf)
    {
        int method; /* method byte */
        int flags;  /* flags byte */
        uInt len;
        int c;
        int err=0;
        z_stream& bgzf_stream = buf->get_bgzf_stream();

        std::basic_istream<Elem, Tr> &istream = buf->get_istream();
        bool m_is_gzip;

        /* Check the gzip magic header */
         for (len = 0; len < 2; len++) 
         {
            c = (int)istream.get();
            if (c != detail::gz_magic[len])
            {
                if (len != 0) 
                    istream.unget();
                if (c!= EOF) 
                {
                    istream.unget();
                }
            
                err = bgzf_stream.avail_in != 0 ? Z_OK : Z_STREAM_END;
                m_is_gzip = false;
                return err;
            }
        }

        m_is_gzip = true;
        method = (int)istream.get();
        flags = (int)istream.get();
        if (method != Z_DEFLATED || (flags & detail::gz_reserved) != 0)
        {
            err = Z_DATA_ERROR;
            return err;
        }

        /* Discard time, xflags and OS code: */
        for (len = 0; len < 6; len++) 
            istream.get();

        if ((flags & detail::gz_extra_field) != 0)
        { 
            /* skip the extra field */
            len  =  (uInt)istream.get();
            len += ((uInt)istream.get())<<8;
            /* len is garbage if EOF but the loop below will quit anyway */
            while (len-- != 0 && istream.get() != EOF) ;
        }
        if ((flags & detail::gz_orig_name) != 0)
        { 
            /* skip the original file name */
            while ((c = istream.get()) != 0 && c != EOF) ;
        }
        if ((flags & detail::gz_comment) != 0)
        {   
            /* skip the .gz file comment */
            while ((c = istream.get()) != 0 && c != EOF) ;
        }
        if ((flags & detail::gz_head_crc) != 0)
        {  /* skip the header crc */
            for (len = 0; len < 2; len++) 
                istream.get();
        }
        err = istream.eof() ? Z_DATA_ERROR : Z_OK;

        return err;
    }
    
        
    template<
		typename Elem, 
		typename Tr,
		typename ElemA,
		typename ByteT,
		typename ByteAT
	>
	int basic_bgzf_istream<
		Elem,Tr,ElemA,ByteT,ByteAT
		>::check_header()
	{
	    return _checkGZHeader(this->rdbuf());
	}

	template<
		typename Elem, 
		typename Tr
	>
    void _putBinaryLong(std::basic_ostream<Elem,Tr> & out_, unsigned long x_)
    {
		static const int size_ul = sizeof(unsigned long);
		static const int size_c = sizeof(typename Tr::char_type);
        static const int n_end = size_ul/size_c;
		out_.write(reinterpret_cast<typename Tr::char_type const*>(&x_), n_end);
    }
   
	template<
		typename Elem, 
		typename Tr,
		typename ElemA,
		typename ByteT,
		typename ByteAT
	>
    void basic_bgzf_istream<
			Elem,Tr,ElemA,ByteT,ByteAT
			>::read_long(
				istream_reference in_, 
			unsigned long& x_
			)
	{
		static const int size_ul = sizeof(unsigned long);
		static const int size_c = sizeof(typename Tr::char_type);
        static const int n_end = size_ul/size_c;
		in_.read(reinterpret_cast<char*>(&x_),n_end);
	}
    
    template<
        typename Elem, 
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    void _addGZHeader(basic_bgzf_streambuf<Elem, Tr, ElemA, ByteT, ByteAT> *buf)
    {
	    typename Tr::char_type zero=0;
	    
        buf->get_ostream()
			.put(static_cast<typename Tr::char_type>(detail::gz_magic[0]))
			.put(static_cast<typename Tr::char_type>(detail::gz_magic[1]))
			.put(static_cast<typename Tr::char_type>(Z_DEFLATED))
			.put(zero) //flags
			.put(zero).put(zero).put(zero).put(zero) // time
			.put(zero) //xflags
			.put(static_cast<typename Tr::char_type>(OS_CODE));
    }

    template<
        typename Elem, 
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    void _addGZFooter(basic_bgzf_streambuf<Elem, Tr, ElemA, ByteT, ByteAT> *buf)
	{
		_putBinaryLong( buf->get_ostream(), buf->get_crc() );
		_putBinaryLong( buf->get_ostream(), buf->get_in_size() );
	}

}  // namespace seqan

#endif