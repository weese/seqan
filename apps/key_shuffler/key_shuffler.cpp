// ==========================================================================
//                                  key_shuffler
// ==========================================================================
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Shuffle numeric idententifiers in tsv files
// ==========================================================================

#include <sstream>

#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <seqan/modifier.h>
#include <seqan/sequence.h>
#include <seqan/tsv_io.h>

#include <cryptopp/rc5.h>
#include <cryptopp/simeck.h>
#include <cryptopp/osrng.h>
#include <cryptopp/secblock.h>
#include <cryptopp/modes.h>

// --------------------------------------------------------------------------
// Class key_shufflerOptions
// --------------------------------------------------------------------------

struct key_shufflerOptions
{
    // Verbosity level.  0 - quiet, 1 - normal, 2 - verbose, 3 - very verbose.
    int verbosity;

    // Path to tsv input file.
    seqan::CharString inTsvPath;

    // Path to tsv output file.
    seqan::CharString outPath;

    // Path to key file.
    seqan::CharString keyPath;

    seqan::StringSet<seqan::CharString> bigintCols, intCols;

    key_shufflerOptions() :
        verbosity(1)
    {
    }
};

// --------------------------------------------------------------------------
// Function parseArgs()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseArgs(key_shufflerOptions & options,
          int argc,
          char ** argv)
{
    seqan::ArgumentParser parser("key_shuffler");
    setShortDescription(parser, "Shuffle numeric identifiers in tsv files..");
    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);
    setCategory(parser, "Utilities");

    addUsageLine(parser, "[\\fIOPTIONS\\fP] [\\fB-o\\fP \\fIOUT.tsv\\fP] \\fIIN.tsv\\fP");
    addDescription(parser, "A tool to permute numeric identifiers and prevent mapping from user ids back to persons\"");
    addDescription(parser, "while maintaining the structure of dataset.");

    // The only argument is the input file.
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE, "IN"));

    // Only tsv files are allowed as input.
    setValidValues(parser, 0, seqan::TsvFileIn::getFileExtensions());

    // TODO(holtgrew): I want a custom help text!
    // addOption(parser, seqan::ArgParseOption("h", "help", "This helpful screen."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Verbose, log to STDERR."));
    hideOption(parser, "verbose");
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Very verbose, log to STDERR."));
    hideOption(parser, "very-verbose");

    addSection(parser, "Output Options");
    addOption(parser, seqan::ArgParseOption("o", "out-path",
                                            "Path to the resulting file.  If omitted, result is printed to stdout in FastQ format.",
                                            seqan::ArgParseOption::OUTPUT_FILE, "TSV"));
    setValidValues(parser, "out-path", seqan::TsvFileOut::getFileExtensions());

    addSection(parser, "Shuffle Options");
    addOption(parser, seqan::ArgParseOption("b", "bigint", "Shuffle big integer column.",
                                            seqan::ArgParseOption::STRING, "COL", true));
    addOption(parser, seqan::ArgParseOption("i", "int", "Shuffle integer column.",
                                            seqan::ArgParseOption::STRING, "COL", true));
    addOption(parser, seqan::ArgParseOption("k", "key-path",
                                            "Path to the key file.",
                                            seqan::ArgParseOption::OUTPUT_FILE, "KEY"));
    setDefaultValue(parser, "key-path", "key.bin");

    seqan::ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    getArgumentValue(options.inTsvPath, parser, 0);

    seqan::CharString tmp;
    getOptionValue(tmp, parser, "out-path");

    if (isSet(parser, "out-path"))
        getOptionValue(options.outPath, parser, "out-path");

    getOptionValue(options.keyPath, parser, "key-path");

    options.intCols = getOptionValues(parser, "int");
    options.bigintCols = getOptionValues(parser, "bigint");

    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;

    return res;
}

// ---------------------------------------------------------------------------
// Function main()
// ---------------------------------------------------------------------------

int main(int argc, char ** argv)
{
    double startTime = 0;

    // Parse command line.
    key_shufflerOptions options;
    seqan::ArgumentParser::ParseResult res = parseArgs(options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;  // 1 on errors, 0 otherwise

    // -----------------------------------------------------------------------
    // Show options.
    // -----------------------------------------------------------------------
    if (options.verbosity >= 2)
    {
        std::cerr << "____OPTIONS___________________________________________________________________\n"
                  << "\n"
                  << "VERBOSITY    " << options.verbosity << "\n"
                  << "IN           " << options.inTsvPath << "\n"
                  << "OUT          " << options.outPath << "\n";
    }

    // -----------------------------------------------------------------------
    // Open Files.
    // -----------------------------------------------------------------------
    seqan::TsvFileIn inFile;
    seqan::TsvFileOut outFile;

    bool openRes = false;
    if (!empty(options.inTsvPath))
        openRes = open(inFile, toCString(options.inTsvPath));
    else
        openRes = open(inFile, std::cin);
    if (!openRes)
    {
        std::cerr << "ERROR: Problem opening input file.\n";
        return 1;
    }

    if (!empty(options.outPath))
        openRes = open(outFile, toCString(options.outPath));
    else
        openRes = open(outFile, std::cout, seqan::Tsv());
    if (!openRes)
    {
        std::cerr << "ERROR: Problem opening output file.\n";
        return 1;
    }

    seqan::String<unsigned char, seqan::MMap<> > keyFile;
    if (!open(keyFile, toCString(options.keyPath), seqan::OPEN_RDWR | seqan::OPEN_CREATE | seqan::OPEN_APPEND))// | seqan::OPEN_CREATE))
    {
        std::cout << "Could not open for writing\n";
        return 1;
    }

    // -----------------------------------------------------------------------
    // Read and Write Filtered.
    // -----------------------------------------------------------------------
    startTime = seqan::sysTime();

    // Copy header.
    seqan::TsvHeader header;
    readHeader(header, inFile);
    writeHeader(outFile, header);

    seqan::CharString shuffleType;

    resize(shuffleType, length(header), '.');
    for (unsigned i = 0; i < length(header); ++i) {
        for (unsigned j = 0; j < length(options.intCols); ++j) {
            if (header[i] == options.intCols[j]) {
                shuffleType[i] = 'i';
                break;
            }
        }
        for (unsigned j = 0; j < length(options.bigintCols); ++j) {
            if (header[i] == options.bigintCols[j]) {
                shuffleType[i] = 'b';
                break;
            }
        }
    }

    if (options.verbosity >= 2)
    {
        std::cerr << "SHUFFLE TYPE " << shuffleType << "\n";
    }


    CryptoPP::AutoSeededRandomPool prng;
    CryptoPP::SecByteBlock key(128);
    if (empty(keyFile)) {
        if (options.verbosity >= 1) {
            std::cerr << "Create new key" << std::endl;
        }
        // create new key and write to file
        prng.GenerateBlock(key, key.size());
        resize(keyFile, key.size(), seqan::Exact());
        std::copy(key.begin(), key.end(), begin(keyFile, seqan::Standard()));
        shrinkToFit(keyFile);
    } else
    {
        // read key from file
        key.Assign(begin(keyFile, seqan::Standard()), length(keyFile));
    }

    CryptoPP::ECB_Mode<CryptoPP::SIMECK32>::Encryption e32;
    e32.SetKey(key, CryptoPP::SIMECK32::MAX_KEYLENGTH);
    CryptoPP::StreamTransformationFilter encryptor32(e32, NULL);

    if (CryptoPP::SIMECK32::MAX_KEYLENGTH > key.size()) {
        std::cerr << "SIMECK32 key too short!" << std::endl;
        return 1;
    }

    CryptoPP::ECB_Mode<CryptoPP::RC5>::Encryption e64;
    e64.SetKey(key, key.size());
    CryptoPP::StreamTransformationFilter encryptor64(e64, NULL);

    seqan::TsvRecord record;
    while (!atEnd(inFile))
    {
        readRecord(record, inFile);

        // fix trailing header bug
        if (length(record) == length(header) * 2 - 1) {
            unsigned len = length(header);
            if (endsWith(record[len - 1], header[0])) {
                resize(record[len - 1], length(record[len - 1]) - length(header[0]));
                resize(record, length(header));
            }
        }

        uint32_t i;
        uint64_t b;
        for (unsigned j = 0; j < length(header) && j < length(record); ++j) {
            switch (shuffleType[j]) {
                case 'i':
                    i = seqan::lexicalCast<uint32_t>(record[j]);

                    encryptor32.Initialize();
                    encryptor32.PutWord32(i);
                    encryptor32.MessageEnd();
                    // if (encryptor32.GetWord32(i) != sizeof(i)) {
                    //     std::cerr << "Somethings wrong here (32)!" << std::endl;
                    //     return 1;
                    // }

                    clear(record[j]);
                    appendNumber(record[j], i);
                    break;
                case 'b':
                    b = seqan::lexicalCast<uint64_t>(record[j]);

                    encryptor64.Initialize();
                    encryptor64.PutWord32(b >> 32);
                    encryptor64.PutWord32(b & 0xffffffff);
                    encryptor64.MessageEnd();
                    uint32_t lo, hi;

                    if (encryptor64.GetWord32(hi) != sizeof(hi) || encryptor64.GetWord32(lo) != sizeof(lo)) {
                        std::cerr << "Somethings wrong here (64)!" << std::endl;
                        return 1;
                    }

                    // b = ((uint64_t)hi << 32) | lo;
                    clear(record[j]);
                    appendNumber(record[j], b);
                    break;
            }
        }

        writeRecord(outFile, record);
    }

    if (options.verbosity >= 2)
        std::cerr << "Took " << (seqan::sysTime() - startTime) << " s\n";

    return 0;
}
