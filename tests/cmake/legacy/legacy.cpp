// Uses libalps and the hdf5 component, as an ALPS 2.x consumer would.
#include <alps/parameter.h>
#include <alps/hdf5/archive.hpp>
#include <iostream>

int main() {
    alps::Parameters parameters;
    parameters["L"] = 4;
    if (!parameters.defined("L"))
        return 1;
    int written = 4;
    int read = 0;
    {
        alps::hdf5::archive archive("legacy.h5", "w");
        archive["/L"] << written;
    }
    {
        alps::hdf5::archive archive("legacy.h5", "r");
        archive["/L"] >> read;
    }
    std::cout << "L = " << read << std::endl;
    return read == written ? 0 : 1;
}
