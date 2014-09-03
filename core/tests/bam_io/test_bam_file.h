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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef CORE_TESTS_BAM_IO_TEST_EASY_BAM_IO_H_
#define CORE_TESTS_BAM_IO_TEST_EASY_BAM_IO_H_

#include <sstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/bam_io.h>

// ---------------------------------------------------------------------------
// Read Header
// ---------------------------------------------------------------------------

void testBamIOBamFileReadHeader(char const * pathFragment)
{
    seqan::CharString filePath = SEQAN_PATH_TO_ROOT();
    append(filePath, pathFragment);

    seqan::BamFileIn bamIO(toCString(filePath));
    seqan::BamHeader header;

    readRecord(header, bamIO);

    SEQAN_ASSERT_EQ(length(header.records), 2u);
    SEQAN_ASSERT_EQ(header.records[0].type, seqan::BAM_HEADER_FIRST);
    SEQAN_ASSERT_EQ(header.records[0].tags[0].i1, "VN");
    SEQAN_ASSERT_EQ(header.records[0].tags[0].i2, "1.3");
    SEQAN_ASSERT_EQ(header.records[0].tags[1].i1, "SO");
    SEQAN_ASSERT_EQ(header.records[0].tags[1].i2, "coordinate");
    SEQAN_ASSERT_EQ(header.records[1].type, seqan::BAM_HEADER_REFERENCE);
    SEQAN_ASSERT_EQ(header.records[1].tags[0].i1, "SN");
    SEQAN_ASSERT_EQ(header.records[1].tags[0].i2, "REFERENCE");
    SEQAN_ASSERT_EQ(header.records[1].tags[1].i1, "LN");
    SEQAN_ASSERT_EQ(header.records[1].tags[1].i2, "10000");
    SEQAN_ASSERT_EQ(length(header.sequenceInfos), 1u);
    SEQAN_ASSERT_EQ(header.sequenceInfos[0].i1, "REFERENCE");
    SEQAN_ASSERT_EQ(header.sequenceInfos[0].i2, 10000);
}

SEQAN_DEFINE_TEST(test_bam_io_bam_file_sam_read_header)
{
    testBamIOBamFileReadHeader("/core/tests/bam_io/small.sam");
}

SEQAN_DEFINE_TEST(test_bam_io_bam_file_bam_read_header)
{
    testBamIOBamFileReadHeader("/core/tests/bam_io/small.bam");
}

// ---------------------------------------------------------------------------
// Read Records
// ---------------------------------------------------------------------------

void testBamIOBamFileReadRecords(char const * pathFragment)
{
    seqan::CharString filePath = SEQAN_PATH_TO_ROOT();
    append(filePath, pathFragment);

    seqan::BamFileIn bamIO(toCString(filePath));
    seqan::BamHeader header;

    readRecord(header, bamIO);

    seqan::String<seqan::BamAlignmentRecord> alignments;
    while (!atEnd(bamIO))
    {
        resize(alignments, length(alignments) + 1);
        readRecord(back(alignments), bamIO);
    }
    SEQAN_ASSERT_EQ(length(alignments), 3u);

    SEQAN_ASSERT_EQ(alignments[0].qName, "READ0");
    SEQAN_ASSERT_EQ(alignments[0].flag, 2u);
    SEQAN_ASSERT_EQ(alignments[0].rID, 0);
    SEQAN_ASSERT_EQ(alignments[0].beginPos, 0);
    SEQAN_ASSERT_EQ(alignments[0].mapQ, 8u);
    SEQAN_ASSERT_EQ(length(alignments[0].cigar), 3u);
    SEQAN_ASSERT_EQ(alignments[0].cigar[0].count, 5u);
    SEQAN_ASSERT_EQ(alignments[0].cigar[0].operation, 'M');
    SEQAN_ASSERT_EQ(alignments[0].cigar[1].count, 1u);
    SEQAN_ASSERT_EQ(alignments[0].cigar[1].operation, 'I');
    SEQAN_ASSERT_EQ(alignments[0].cigar[2].count, 4u);
    SEQAN_ASSERT_EQ(alignments[0].cigar[2].operation, 'M');
    SEQAN_ASSERT_EQ(alignments[0].rNextId, 0);
    SEQAN_ASSERT_EQ(alignments[0].pNext, 30);
    SEQAN_ASSERT_EQ(alignments[0].tLen, 40);
    SEQAN_ASSERT_EQ(alignments[0].seq, "AAAAAAAAAA");
    SEQAN_ASSERT_EQ(alignments[0].qual, "!!!!!!!!!!");
    SEQAN_ASSERT_EQ(length(alignments[0].tags), 0u);

    SEQAN_ASSERT_EQ(alignments[1].qName, "READ0");
    SEQAN_ASSERT_EQ(alignments[1].flag, 1u);
    SEQAN_ASSERT_EQ(alignments[1].rID, 0);
    SEQAN_ASSERT_EQ(alignments[1].beginPos, 1);
    SEQAN_ASSERT_EQ(alignments[1].mapQ, 8u);
    SEQAN_ASSERT_EQ(length(alignments[1].cigar), 3u);
    SEQAN_ASSERT_EQ(alignments[1].cigar[0].count, 5u);
    SEQAN_ASSERT_EQ(alignments[1].cigar[0].operation, 'M');
    SEQAN_ASSERT_EQ(alignments[1].cigar[1].count, 1u);
    SEQAN_ASSERT_EQ(alignments[1].cigar[1].operation, 'I');
    SEQAN_ASSERT_EQ(alignments[1].cigar[2].count, 4u);
    SEQAN_ASSERT_EQ(alignments[1].cigar[2].operation, 'M');
    SEQAN_ASSERT_EQ(alignments[1].rNextId, 0);
    SEQAN_ASSERT_EQ(alignments[1].pNext, 30);
    SEQAN_ASSERT_EQ(alignments[1].tLen, 40);
    SEQAN_ASSERT_EQ(alignments[1].seq, "AAAAAAAAAA");
    SEQAN_ASSERT_EQ(alignments[1].qual, "!!!!!!!!!!");
    SEQAN_ASSERT_EQ(length(alignments[1].tags), 0u);

    SEQAN_ASSERT_EQ(alignments[2].qName, "READ0");
    SEQAN_ASSERT_EQ(alignments[2].flag, 3u);
    SEQAN_ASSERT_EQ(alignments[2].rID, 0);
    SEQAN_ASSERT_EQ(alignments[2].beginPos, 2);
    SEQAN_ASSERT_EQ(alignments[2].mapQ, 8u);
    SEQAN_ASSERT_EQ(length(alignments[2].cigar), 3u);
    SEQAN_ASSERT_EQ(alignments[2].cigar[0].count, 5u);
    SEQAN_ASSERT_EQ(alignments[2].cigar[0].operation, 'M');
    SEQAN_ASSERT_EQ(alignments[2].cigar[1].count, 1u);
    SEQAN_ASSERT_EQ(alignments[2].cigar[1].operation, 'I');
    SEQAN_ASSERT_EQ(alignments[2].cigar[2].count, 4u);
    SEQAN_ASSERT_EQ(alignments[2].cigar[2].operation, 'M');
    SEQAN_ASSERT_EQ(alignments[2].rNextId, -1);
//    SEQAN_ASSERT_EQ(alignments[2].rNextId, alignments[2].INVALID_REFID);
    SEQAN_ASSERT_EQ(alignments[2].pNext, -1);
//    SEQAN_ASSERT_EQ(alignments[2].pNext, alignments[2].INVALID_POS);
    SEQAN_ASSERT_EQ(alignments[2].tLen, 0);
    SEQAN_ASSERT_EQ(alignments[2].seq, "AAAAAAAAAA");
    SEQAN_ASSERT_EQ(alignments[2].qual, "!!!!!!!!!!");
    SEQAN_ASSERT_EQ(length(alignments[2].tags), 0u);

    SEQAN_ASSERT_EQ(nameStore(context(bamIO))[0], "REFERENCE");
}

SEQAN_DEFINE_TEST(test_bam_io_bam_file_sam_read_records)
{
    testBamIOBamFileReadRecords("/core/tests/bam_io/small.sam");
}

SEQAN_DEFINE_TEST(test_bam_io_bam_file_bam_read_records)
{
    testBamIOBamFileReadRecords("/core/tests/bam_io/small.bam");
}

SEQAN_DEFINE_TEST(test_bam_io_bam_file_bam_read_ex1)
{
    seqan::CharString filePath = SEQAN_PATH_TO_ROOT();
    append(filePath, "/core/tests/bam_io/ex1.bam");

    seqan::BamFileIn bamIO(toCString(filePath));
    seqan::BamHeader header;

    readRecord(header, bamIO);

    SEQAN_ASSERT_EQ(nameStore(context(bamIO))[0], "seq1");
    SEQAN_ASSERT_EQ(nameStore(context(bamIO))[1], "seq2");
    SEQAN_ASSERT_EQ(context(bamIO).translateFile2GlobalRefId[0], 0u);
    SEQAN_ASSERT_EQ(context(bamIO).translateFile2GlobalRefId[1], 1u);

    seqan::BamAlignmentRecord record;
    seqan::String<int> counts;
    resize(counts, 2, 0);
    while (!atEnd(bamIO))
    {
        readRecord(record, bamIO);
        ++counts[record.rID];
//        seqan::CharString name = nameStore(context(bamIO))[record.rID];
//        std::cout << "Chrom: " << name << " (" << record.rID << ")" << std::endl;
    }
    SEQAN_ASSERT_EQ(counts[0], 1501);
    SEQAN_ASSERT_EQ(counts[1], 1806);
}

// ---------------------------------------------------------------------------
// Write Header
// ---------------------------------------------------------------------------

void testBamIOBamFileWriteHeader(char const * pathFragmentExpected)
{
    seqan::CharString filePath = SEQAN_PATH_TO_ROOT();
    append(filePath, pathFragmentExpected);

    seqan::CharString tmpPath = SEQAN_TEMP_FILENAME();
    if (seqan::endsWith(pathFragmentExpected, ".bam"))
        append(tmpPath, ".bam");
    else
        append(tmpPath, ".sam");

    // Initialize BamFile, build header.
    seqan::BamFileOut bamIO(toCString(tmpPath));

    seqan::BamHeader header;
    resize(header.sequenceInfos, 1);
    header.sequenceInfos[0].i1 = "REFERENCE";
    header.sequenceInfos[0].i2 = 10000;
    resize(header.records, 2);
    resize(header.records[0].tags, 2);
    header.records[0].type = seqan::BAM_HEADER_FIRST;
    header.records[0].tags[0].i1 = "VN";
    header.records[0].tags[0].i2 = "1.3";
    header.records[0].tags[1].i1 = "SO";
    header.records[0].tags[1].i2 = "coordinate";
    resize(header.records[1].tags, 2);
    header.records[1].type = seqan::BAM_HEADER_REFERENCE;
    header.records[1].tags[0].i1 = "SN";
    header.records[1].tags[0].i2 = "REFERENCE";
    header.records[1].tags[1].i1 = "LN";
    header.records[1].tags[1].i2 = "10000";
    writeRecord(bamIO, header);

    // Force writing of header on flush.
    close(bamIO);

    // Compare results.
    if (seqan::endsWith(pathFragmentExpected, ".bam"))
        SEQAN_ASSERT(seqan::_compareBinaryFiles(toCString(tmpPath), toCString(filePath)));
    else
        SEQAN_ASSERT(seqan::_compareTextFiles(toCString(tmpPath), toCString(filePath)));
}

SEQAN_DEFINE_TEST(test_bam_io_bam_file_sam_write_header)
{
    testBamIOBamFileWriteHeader("/core/tests/bam_io/header.sam");
}

SEQAN_DEFINE_TEST(test_bam_io_bam_file_bam_write_header)
{
    testBamIOBamFileWriteHeader("/core/tests/bam_io/header.bam");
}

// ---------------------------------------------------------------------------
// Write Records
// ---------------------------------------------------------------------------

void testBamIOBamFileWriteRecords(char const * pathFragmentExpected)
{
    seqan::CharString filePath = SEQAN_PATH_TO_ROOT();
    append(filePath, pathFragmentExpected);

    seqan::CharString tmpPath = SEQAN_TEMP_FILENAME();
    if (seqan::endsWith(pathFragmentExpected, ".bam"))
        append(tmpPath, ".bam");
    else
        append(tmpPath, ".sam");

    // Initialize BamFile, build header.
    seqan::BamFileOut bamIO(toCString(tmpPath));

    seqan::BamHeader header;
    resize(header.sequenceInfos, 1);
    header.sequenceInfos[0].i1 = "REFERENCE";
    header.sequenceInfos[0].i2 = 10000;
    resize(header.records, 2);
    resize(header.records[0].tags, 2);
    header.records[0].type = seqan::BAM_HEADER_FIRST;
    header.records[0].tags[0].i1 = "VN";
    header.records[0].tags[0].i2 = "1.3";
    header.records[0].tags[1].i1 = "SO";
    header.records[0].tags[1].i2 = "coordinate";
    resize(header.records[1].tags, 2);
    header.records[1].type = seqan::BAM_HEADER_REFERENCE;
    header.records[1].tags[0].i1 = "SN";
    header.records[1].tags[0].i2 = "REFERENCE";
    header.records[1].tags[1].i1 = "LN";
    header.records[1].tags[1].i2 = "10000";
    writeRecord(bamIO, header);

    // Construct first records.
    seqan::BamAlignmentRecord record;

    record.qName = "READ0";
    record.flag = 2;
    record.rID = 0;
    record.beginPos = 0;
    record.mapQ = 8;
    resize(record.cigar, 3);
    record.cigar[0].count = 5;
    record.cigar[0].operation = 'M';
    record.cigar[1].count = 1;
    record.cigar[1].operation = 'I';
    record.cigar[2].count = 4;
    record.cigar[2].operation = 'M';
    record.rNextId = 0;
    record.pNext = 30;
    record.tLen = 40;
    record.seq = "AAAAAAAAAA";
    record.qual = "!!!!!!!!!!";
    writeRecord(bamIO, record);

    record.qName = "READ0";
    record.flag = 1;
    record.rID = 0;
    record.beginPos = 1;
    record.mapQ = 8;
    resize(record.cigar, 3);
    record.cigar[0].count = 5;
    record.cigar[0].operation = 'M';
    record.cigar[1].count = 1;
    record.cigar[1].operation = 'I';
    record.cigar[2].count = 4;
    record.cigar[2].operation = 'M';
    record.rNextId = 0;
    record.pNext = 30;
    record.tLen = 40;
    record.seq = "AAAAAAAAAA";
    record.qual = "!!!!!!!!!!";
    writeRecord(bamIO, record);

    record.qName = "READ0";
    record.flag = 3;
    record.rID = 0;
    record.beginPos = 2;
    record.mapQ = 8;
    resize(record.cigar, 3);
    record.cigar[0].count = 5;
    record.cigar[0].operation = 'M';
    record.cigar[1].count = 1;
    record.cigar[1].operation = 'I';
    record.cigar[2].count = 4;
    record.cigar[2].operation = 'M';
    record.rNextId = seqan::BamAlignmentRecord::INVALID_REFID;
    record.pNext = seqan::BamAlignmentRecord::INVALID_POS;
    record.tLen = seqan::BamAlignmentRecord::INVALID_LEN;
    record.seq = "AAAAAAAAAA";
    record.qual = "!!!!!!!!!!";
    writeRecord(bamIO, record);

    // Force writing of everything.
    close(bamIO);

    // Compare results.
    if (seqan::endsWith(pathFragmentExpected, ".bam"))
        SEQAN_ASSERT(seqan::_compareBinaryFiles(toCString(tmpPath), toCString(filePath)));
    else
        SEQAN_ASSERT(seqan::_compareTextFiles(toCString(tmpPath), toCString(filePath)));
}

SEQAN_DEFINE_TEST(test_bam_io_bam_file_sam_write_records)
{
    testBamIOBamFileWriteRecords("/core/tests/bam_io/small.sam");
}

SEQAN_DEFINE_TEST(test_bam_io_bam_file_bam_write_records)
{
    testBamIOBamFileWriteRecords("/core/tests/bam_io/small.bam");
}

// ---------------------------------------------------------------------------
// File Size / Byte Position In File
// ---------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_bam_io_bam_file_sam_file_size)
{
    seqan::CharString filePath = SEQAN_PATH_TO_ROOT();
    append(filePath, "/core/tests/bam_io/small.sam");

    seqan::BamFileIn bamFile(toCString(filePath));

    SEQAN_ASSERT_EQ(position(bamFile), 0u);

    seqan::BamHeader header;
    readRecord(header, bamFile);

    SEQAN_ASSERT_EQ(position(bamFile), 51u);
//    SEQAN_ASSERT_EQ(fileSize(bamFile), 226u);

    seqan::BamAlignmentRecord record;
    readRecord(record, bamFile);

    SEQAN_ASSERT_EQ(position(bamFile), 110u);

    readRecord(record, bamFile);

    SEQAN_ASSERT_EQ(position(bamFile), 169u);

    readRecord(record, bamFile);

    SEQAN_ASSERT_EQ(position(bamFile), 226u);
}

SEQAN_DEFINE_TEST(test_bam_io_bam_file_bam_file_size)
{
    seqan::CharString filePath = SEQAN_PATH_TO_ROOT();
    append(filePath, "/core/tests/bam_io/small.bam");

    seqan::BamFileIn bamFile(toCString(filePath));

    SEQAN_ASSERT_EQ(position(bamFile), 0u);

    seqan::BamHeader header;
    readRecord(header, bamFile);

//    SEQAN_ASSERT_EQ(fileSize(bamFile), 181u);
    SEQAN_ASSERT_EQ(position(bamFile), 0x0051u);

    seqan::BamAlignmentRecord record;
    readRecord(record, bamFile);

    SEQAN_ASSERT_EQ(position(bamFile), 0x0096u);  // [block begin in bam file] << 16 + [local offset]

    readRecord(record, bamFile);

    SEQAN_ASSERT_EQ(position(bamFile), 0x00dbu);

    readRecord(record, bamFile);

    SEQAN_ASSERT_EQ(position(bamFile), 0x0120u);
}

#endif  // #ifndef CORE_TESTS_BAM_IO_TEST_EASY_BAM_IO_H_