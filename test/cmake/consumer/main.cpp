#ifdef _WIN32
#include <windows.h>
#endif
#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/ngs/params.hpp>
#include <alps/ngs/accumulator/feature/binning_analysis.hpp>
#include <alps/ngs/accumulator/feature/max_num_binning.hpp>
#include <alps/parser/xslt_path.h>
#include <boost/filesystem/operations.hpp>
#include <vector>

using namespace alps::accumulator;
using Count = impl::Accumulator<double, count_tag, impl::AccumulatorBase<double>>;
using Mean = impl::Accumulator<double, mean_tag, Count>;
using Error = impl::Accumulator<double, error_tag, Mean>;
using Binning = impl::Accumulator<double, binning_analysis_tag, Error>;
using MaxBinning = impl::Accumulator<double, max_num_binning_tag, Error>;

template<class Accumulator> bool check_accumulator() {
    Accumulator accumulator;
    accumulator(1.0);
    accumulator(2.0);
    accumulator(3.0);
    return accumulator.count() == 3 && accumulator.mean() == 2.0;
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
    alps::params parameters;
    parameters["count"] = 3;
    const std::vector<double> expected{1.0, 2.0, 3.0};
    {
        alps::hdf5::archive output("sdk-contract.h5", "w");
        output["/values"] << expected;
    }
    std::vector<double> actual;
    alps::hdf5::archive input("sdk-contract.h5", "r");
    input["/values"] >> actual;
    return actual == expected && int(parameters["count"]) == 3
        && check_accumulator<Mean>() && check_accumulator<Error>()
        && check_accumulator<Binning>() && check_accumulator<MaxBinning>() ? 0 : 1;
}
