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
// Tests for the stream module.
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/stream.h>

#include "test_stream_tokenization.h"
#include "test_stream_virtual_stream.h"

using namespace seqan;

// ==========================================================================
// Types
// ========================================================================== 

// --------------------------------------------------------------------------
// Input Stream Types
// --------------------------------------------------------------------------

typedef
    TagList<std::fstream,
    TagList<std::ifstream
    > >
    InputStreamTypes;

SEQAN_TYPED_TEST_CASE(InputStreamTest, InputStreamTypes);

// --------------------------------------------------------------------------
// Output Stream Types
// --------------------------------------------------------------------------

typedef
    TagList<std::fstream,
    TagList<std::ofstream
    > >
    OutputStreamTypes;

SEQAN_TYPED_TEST_CASE(OutputStreamTest, OutputStreamTypes);

// --------------------------------------------------------------------------
// All Stream Types
// --------------------------------------------------------------------------

typedef
    TagList<std::fstream,
    TagList<std::ifstream,
    TagList<std::ofstream
    > > >
    AllStreamTypes;

SEQAN_TYPED_TEST_CASE(StreamTest, AllStreamTypes);

// ========================================================================== 
// Test Classes
// ========================================================================== 

// --------------------------------------------------------------------------
// Class StreamTest
// --------------------------------------------------------------------------

template <typename TStream>
class StreamTest : public Test
{
public:
    typedef TStream Type;
    
    TStream stream;
};

// --------------------------------------------------------------------------
// Class InputStreamTest
// --------------------------------------------------------------------------

template <typename TStream>
class InputStreamTest : public StreamTest<TStream> {};

// --------------------------------------------------------------------------
// Class OutputStreamTest
// --------------------------------------------------------------------------

template <typename TStream>
class OutputStreamTest : public StreamTest<TStream> {};

// ========================================================================== 
// Test Context
// ========================================================================== 

// --------------------------------------------------------------------------
// Class StreamTestContext
// --------------------------------------------------------------------------
// NOTE(esiragusa): The test context is shared by all tests of a given stream type.

template <typename TStream = void>
class StreamTestContext
{
public:
    CharString      outputFilename;
    CharString      inputFilename;
    CharString      content;
    unsigned        filesize;

    static bool     needToInitialize;

    // Get the singleton.
    static StreamTestContext const & get()
    {
        static StreamTestContext ctx;
        if (needToInitialize)
        {
            init(ctx);
            needToInitialize = false;
        }
        return ctx;
    }

protected:
    // The constructor is called on first call to the get().
    StreamTestContext() :
        outputFilename(SEQAN_TEMP_FILENAME()),
        inputFilename(SEQAN_TEMP_FILENAME()),
        content("This \t is a \r test \n")
    {
        // Testing for null chars is always good practice!
        appendValue(content, '\0');

        // Total file size is ~ 2.1 MB.
        filesize = (2 * Power<2, 20>::VALUE) - ((2 * Power<2, 20>::VALUE) % length(content));
    };

private:
    // Prevent copying the singleton.
    StreamTestContext(StreamTestContext const &);
    void operator=(StreamTestContext const &);
};

template <typename TStream>
bool StreamTestContext<TStream>::needToInitialize = true;

// --------------------------------------------------------------------------
// Function init(); Default
// --------------------------------------------------------------------------

template <typename TStream>
void init(StreamTestContext<TStream> & ctx)
{
    std::ofstream file;

    // Truncate file.
    file.open(toCString(ctx.inputFilename), std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
    file.exceptions(std::ios_base::failbit | std::ios_base::badbit);

    // Write content to file.
    for (unsigned i = 0; i < ctx.filesize / length(ctx.content); ++i)
        file.write(begin(ctx.content, Standard()), length(ctx.content));

    file.close();
}

// --------------------------------------------------------------------------
// Function init(); Stream<GZFile>
// --------------------------------------------------------------------------

//#ifdef SEQAN_HAS_ZLIB
//template <typename TStream>
//void init(StreamTestContext<Stream<GZFile> > & ctx)
//{
//    gzFile file = gzopen(toCString(ctx.inputFilename), "wb");
//    SEQAN_ASSERT(file != NULL);
//
//    for (unsigned i = 0; i < ctx.filesize / length(ctx.content); ++i)
//        gzwrite(file, begin(ctx.content, Standard()), length(ctx.content));
//
//    gzclose(file);
//}
//#endif

// ========================================================================== 
// Input Tests
// ========================================================================== 

// --------------------------------------------------------------------------
// Test open()
// --------------------------------------------------------------------------

SEQAN_TYPED_TEST(InputStreamTest, Open)
{
    typedef typename TestFixture::Type  TStream;

    {
        bool feature = HasStreamFeature<TStream, HasFilename>::VALUE;
        SEQAN_ASSERT(feature);
    }
    {
        bool feature = HasStreamFeature<TStream, IsInput>::VALUE;
        SEQAN_ASSERT(feature);
    }

    StreamTestContext<TStream> const & ctx = StreamTestContext<TStream>::get();

    SEQAN_ASSERT(open(this->stream, toCString(ctx.inputFilename), OPEN_RDONLY));
    streamGet(this->stream);

//    SEQAN_ASSERT_NOT(open(this->stream, toCString(ctx.inputFilename), OPEN_RDONLY));
}

// --------------------------------------------------------------------------
// Test streamGet()
// --------------------------------------------------------------------------

SEQAN_TYPED_TEST(InputStreamTest, Get)
{
    typedef typename TestFixture::Type  TStream;

    StreamTestContext<TStream> const & ctx = StreamTestContext<TStream>::get();
    open(this->stream, toCString(ctx.inputFilename));

    for (unsigned i = 0; i < ctx.filesize; i += length(ctx.content))
    {
        streamSeek(this->stream, i, SEEK_SET);
        SEQAN_ASSERT_EQ(streamGet(this->stream), front(ctx.content));
    }
}

// --------------------------------------------------------------------------
// Test streamEof()
// --------------------------------------------------------------------------

SEQAN_TYPED_TEST(InputStreamTest, Eof)
{
    typedef typename TestFixture::Type  TStream;

    StreamTestContext<TStream> const & ctx = StreamTestContext<TStream>::get();

    open(this->stream, toCString(ctx.inputFilename), OPEN_RDONLY);
    SEQAN_ASSERT_NOT(streamEof(this->stream));
    streamSeek(this->stream, 0u, SEEK_END);
    SEQAN_ASSERT(streamEof(this->stream));
}

// ========================================================================== 
// Output Tests
// ========================================================================== 

// --------------------------------------------------------------------------
// Test open()
// --------------------------------------------------------------------------

SEQAN_TYPED_TEST(OutputStreamTest, Open)
{
    typedef typename TestFixture::Type  TStream;

    {
        bool feature = HasStreamFeature<TStream, HasFilename>::VALUE;
        SEQAN_ASSERT(feature);
    }
    {
        bool feature = HasStreamFeature<TStream, IsOutput>::VALUE;
        SEQAN_ASSERT(feature);
    }

    StreamTestContext<TStream> const & ctx = StreamTestContext<TStream>::get();

    SEQAN_ASSERT(open(this->stream, toCString(ctx.outputFilename), OPEN_WRONLY));
    streamPut(this->stream, front(ctx.content));
    close(this->stream);
    SEQAN_ASSERT(open(this->stream, toCString(ctx.inputFilename), OPEN_WRONLY|OPEN_APPEND));
    SEQAN_ASSERT_EQ(streamTell(this->stream), ctx.filesize);
}

// --------------------------------------------------------------------------
// Test streamPut()
// --------------------------------------------------------------------------

SEQAN_TYPED_TEST(OutputStreamTest, Put)
{
    typedef typename TestFixture::Type  TStream;

    StreamTestContext<TStream> const & ctx = StreamTestContext<TStream>::get();
    open(this->stream, toCString(ctx.outputFilename), OPEN_WRONLY);

//    for (unsigned i = 0; i < ctx.filesize; i += length(ctx.content))
//    {
//        streamPut(this->stream, front(ctx.content));
//    }
//    SEQAN_ASSERT_EQ(streamTell(this->stream), 1u);
//    close(this->stream);
}

// ========================================================================== 
// Other Tests
// ========================================================================== 

// --------------------------------------------------------------------------
// Test streamTell()
// --------------------------------------------------------------------------

SEQAN_TYPED_TEST(StreamTest, Tell)
{
    typedef typename TestFixture::Type  TStream;

//    {
//        bool feature = HasStreamFeature<TStream, Seek<void> >::VALUE;
//        SEQAN_ASSERT(feature);
//    }
//    {
//        bool feature = HasStreamFeature<TStream, Tell>::VALUE;
//        SEQAN_ASSERT(feature);
//    }

    StreamTestContext<TStream> const & ctx = StreamTestContext<TStream>::get();
    open(this->stream, toCString(ctx.inputFilename));

    streamSeek(this->stream, 0u, SEEK_END);
    SEQAN_ASSERT_EQ(streamTell(this->stream), ctx.filesize);

    streamSeek(this->stream, 11u, SEEK_END);
    SEQAN_ASSERT_EQ(streamTell(this->stream), ctx.filesize + 11u);

    streamSeek(this->stream, 23u, SEEK_SET);
    SEQAN_ASSERT_EQ(streamTell(this->stream), 23u);

    streamSeek(this->stream, 23u, SEEK_CUR);
    SEQAN_ASSERT_EQ(streamTell(this->stream), 46u);
}

// ========================================================================== 
// Functions
// ========================================================================== 

int main(int argc, char const ** argv)
{
    TestSystem::init(argc, argv);
    return TestSystem::runAll();
}
