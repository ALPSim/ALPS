// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
// The deprecated include spelling must stay self-contained in the SDK and keep
// writing observables to the historical results path.
#include <alps/python/save_observable_to_hdf5.hpp>

#include <boost/filesystem.hpp>

#include <stdexcept>
#include <string>

struct observable {
    std::string representation() const { return "Energy"; }
    void save(alps::hdf5::archive & ar) const { ar["mean/value"] << 1.5; }
};

int main() {
    const std::string filename = "legacy-save-observable.h5";
    boost::filesystem::remove(filename);
    alps::python::save_observable_to_hdf5(observable(), filename);
    double value = 0;
    {
        alps::hdf5::archive ar(filename, "r");
        ar["/simulation/results/Energy/mean/value"] >> value;
    }
    boost::filesystem::remove(filename);
    if (value != 1.5) throw std::runtime_error("observable was not written to /simulation/results");
}
