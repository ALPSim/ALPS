// SPDX-License-Identifier: MIT
#include <alps/gf/gf.hpp>
#include <alps/hdf5/archive.hpp>
#include <complex>
#include <cstdio>
#include <stdexcept>

int main() {
    typedef alps::gf::matsubara_mesh<alps::gf::mesh::POSITIVE_ONLY> matsubara_mesh;
    alps::gf::omega_gf gf(matsubara_mesh(10.0, 8));
    gf.initialize();
    alps::gf::matsubara_index omega;
    omega = 3;
    gf(omega) = std::complex<double>(3, 4);
    {
        alps::hdf5::archive archive("gf_contract.h5", "w");
        gf.save(archive, "/gf");
    }
    alps::gf::omega_gf loaded;
    {
        alps::hdf5::archive archive("gf_contract.h5", "r");
        loaded.load(archive, "/gf");
    }
    std::remove("gf_contract.h5");
    if (loaded != gf || loaded(omega) != std::complex<double>(3, 4))
        throw std::runtime_error("GF archive round trip failed");
}
