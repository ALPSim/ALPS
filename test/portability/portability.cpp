// SPDX-License-Identifier: MIT
#include <alps/utility/temporary_filename.hpp>
#include <alps/utility/os.hpp>
#include <alps/parapack/filelock.h>
#include <alps/ngs/sleep.hpp>
#include <alps/osiris/xdrcore.h>
#include <alps/osiris/xdrdump.h>
#include <alps/parser/xmlstream.h>
#include <alps/parser/parser.h>
#include <alps/random/parallel/lcg64.hpp>
#include <alps/random/rngfactory.h>
#include <boost/filesystem/operations.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <memory>
#include <limits>
#include <stdexcept>
#include <sstream>
#include <iostream>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
struct temporary_file {
    std::string path = alps::temporary_filename(
        (alps::temp_directory_path() / "alps-portability-").string());
    ~temporary_file() { boost::system::error_code ec; boost::filesystem::remove(path, ec); }
};
}

int main() try {
    alps::lcg64a parallel_rng, same_seed;
    for (int i = 0; i < 100; ++i)
        require(parallel_rng() == same_seed(), "parallel RNG seeding must be reproducible");
    std::unique_ptr<alps::buffered_rng_base> rng(alps::rng_factory.create("mt19937"));
    rng->seed(42);
    const auto first_value = (*rng)();
    rng->seed(42);
    require(first_value == (*rng)(), "factory RNG must restart after seeding");
    rng->seed();
    const auto default_value = (*rng)();
    rng->seed();
    require(default_value == (*rng)(), "default reseeding must discard buffered values");
    alps::pseudo_des seed_source(7), same_seed_source(7);
    rng->seed(seed_source);
    const auto generated_seed_value = (*rng)();
    rng->seed(same_seed_source);
    require(generated_seed_value == (*rng)(), "generator reseeding must discard buffered values");
    alps::buffered_rng<boost::mt19937> sequence_rng;
    std::array<uint32_t, boost::mt19937::state_size> seed_words{};
    for (std::size_t i = 0; i < seed_words.size(); ++i) seed_words[i] = uint32_t(i + 1);
    sequence_rng.seed(seed_words.begin(), seed_words.end());
    const auto sequence_value = sequence_rng();
    sequence_rng.seed(seed_words.begin(), seed_words.end());
    require(sequence_value == sequence_rng(), "sequence reseeding must discard buffered values");
    std::istringstream text("value;");
    require(alps::read_until(text, ';') == "value", "SDK parser API must be usable by consumers");
    require(alps::precision(std::numeric_limits<double>::signaling_NaN(), 6) == "nan",
            "XML must use portable NaN spelling");
    require(alps::precision(-std::numeric_limits<double>::infinity(), 6) == "-inf",
            "XML must preserve infinity sign");
    temporary_file first, second;
    require(first.path != second.path, "temporary names must be distinct");
    require(boost::filesystem::is_regular_file(first.path), "temporary_filename must reserve a file");
    alps::filelock owner(first.path, true, 0);
    alps::filelock competitor(first.path);
    require(competitor.locked(), "a second owner must see the lock");
    bool refused = false;
    try { competitor.lock(0); } catch (const std::logic_error&) { refused = true; }
    require(refused, "exclusive lock must reject a second owner");
    owner.release();
    competitor.lock(0);
    competitor.release();
    require(!competitor.locked(), "release must remove the lock");

    auto start = std::chrono::steady_clock::now();
    alps::sleep(20000000);
    require(std::chrono::steady_clock::now() - start >= std::chrono::milliseconds(5),
            "portable sleep must wait");

    std::unique_ptr<std::FILE, decltype(&std::fclose)> file(
        std::fopen(first.path.c_str(), "w+b"), &std::fclose);
    require(bool(file), "open XDR fixture");
    XDR stream;
    xdrstdio_create(&stream, file.get(), XDR_ENCODE);
    int32_t signed_value = -2;
    uint32_t unsigned_value = UINT32_C(0xfedcba98);
    float single = 1.0f;
    double real = -2.5;
    require(xdr_int32_t(&stream, &signed_value) && xdr_uint32_t(&stream, &unsigned_value)
            && xdr_float(&stream, &single) && xdr_double(&stream, &real), "encode XDR fixture");
    xdr_destroy(&stream);
    const std::array<unsigned char, 20> expected = {{
        0xff,0xff,0xff,0xfe, 0xfe,0xdc,0xba,0x98,
        0x3f,0x80,0x00,0x00, 0xc0,0x04,0x00,0x00,0x00,0x00,0x00,0x00}};
    std::array<unsigned char, 20> actual{};
    std::rewind(file.get());
    require(std::fread(actual.data(), 1, actual.size(), file.get()) == actual.size(), "read XDR bytes");
    require(actual == expected, "XDR must use standard big-endian wire bytes");
    std::rewind(file.get());
    xdrstdio_create(&stream, file.get(), XDR_DECODE);
    signed_value = 0; unsigned_value = 0; single = 0; real = 0;
    require(xdr_int32_t(&stream, &signed_value) && xdr_uint32_t(&stream, &unsigned_value)
            && xdr_float(&stream, &single) && xdr_double(&stream, &real), "decode XDR fixture");
    require(signed_value == -2 && unsigned_value == UINT32_C(0xfedcba98)
            && single == 1.0f && real == -2.5, "XDR round trip");
    xdr_destroy(&stream);

    // A 32-bit negative XDR long must sign-extend on LP64 and LLP64 alike.
    std::rewind(file.get());
    xdrstdio_create(&stream, file.get(), XDR_DECODE);
    long wide = 0;
    require(xdr_long(&stream, &wide) && wide == -2, "XDR long sign extension");
    xdr_destroy(&stream);

    // The high bit must survive on both LP64 and Windows LLP64. This also
    // catches undefined signed left shifts in the old 64-bit decoder.
    file.reset();
    const long long minimum = std::numeric_limits<long long>::min();
    const unsigned long long maximum = std::numeric_limits<unsigned long long>::max();
    {
        alps::OXDRFileDump dump(first.path);
        dump << minimum << maximum;
    }
    {
        alps::IXDRFileDump dump(first.path);
        long long decoded_min = 0;
        unsigned long long decoded_max = 0;
        dump >> decoded_min >> decoded_max;
        require(decoded_min == minimum && decoded_max == maximum, "XDR 64-bit limits");
    }
} catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
}
