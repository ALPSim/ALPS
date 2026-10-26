#ifdef _WIN32
#include <windows.h>
#endif
#include <boost/spirit/include/classic_core.hpp>
#include <alps/parapack/integer_range.h>
#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/params.hpp>
#include <alps/utilities/encode.hpp>
#include <alps/accumulators.hpp>
#include <alps/parser/xslt_path.h>
#include <alps/osiris/xdrdump.h>
#include <alps/numeric/functional.hpp>
#include <boost/filesystem/operations.hpp>
#include <vector>

template<class Accumulator> bool check_accumulator() {
    alps::accumulators::accumulator_set measurements;
    measurements << Accumulator("x");
    for (double value : {1.0, 2.0, 3.0})
        measurements["x"] << value;
    const alps::accumulators::result_set results(measurements);
    return results["x"].count() == 3 && results["x"].mean<double>() == 2.0;
}

int main(int argc, char** argv) {
    if (argc == 2) {
        for (const char* resource : {"lattices.xml", "models.xml", "ALPS.xsl"}) {
            if (!boost::filesystem::equivalent(alps::search_xml_library_path(resource),
                                               boost::filesystem::path(argv[1]) / resource))
                return 1;
        }
        return 0;
    }
    const std::string archive_name = "a/path with spaces";
    if (alps::hdf5_name_decode(alps::hdf5_name_encode(archive_name)) != archive_name)
        return 1;
    alps::params parameters;
    parameters["count"] = 3;
    const alps::integer_range<int> range("[2:5]");
    const std::vector<double> expected{1.0, 2.0, 3.0};
    {
        alps::OXDRFileDump output("sdk-contract.xdr");
        output << 42;
    }
    int checkpoint_value = 0;
    {
        alps::IXDRFileDump input("sdk-contract.xdr");
        input >> checkpoint_value;
    }
    {
        alps::hdf5::archive output("sdk-contract.h5", "w");
        output["/values"] << expected;
    }
    std::vector<double> actual;
    alps::hdf5::archive input("sdk-contract.h5", "r");
    input["/values"] >> actual;
    return actual == expected && checkpoint_value == 42 && int(parameters["count"]) == 3
        && (range.min)() == 2 && (range.max)() == 5
        && check_accumulator<alps::accumulators::MeanAccumulator<double>>()
        && check_accumulator<alps::accumulators::NoBinningAccumulator<double>>()
        && check_accumulator<alps::accumulators::LogBinningAccumulator<double>>()
        && check_accumulator<alps::accumulators::FullBinningAccumulator<double>>() ? 0 : 1;
}
