// Copyright (C) 2005-2026 by Emanuel Gull
// ALPS Project: https://alps.comp-phys.org/
// SPDX-License-Identifier: MIT

#include<alps/hdf5.hpp>
#include"./vector.h"
#include"./matrix.h"
#ifndef BLAS_DUMP_HEADER
#define BLAS_DUMP_HEADER
namespace blas{
alps::hdf5::archive & serialize(alps::hdf5::archive & ar, std::string const & p, blas::matrix const &M);
alps::hdf5::archive & serialize(alps::hdf5::archive & ar, std::string const & p, blas::matrix &M);
alps::hdf5::archive & serialize(alps::hdf5::archive & ar, std::string const & p, blas::vector const &M);
alps::hdf5::archive & serialize(alps::hdf5::archive & ar, std::string const & p, blas::vector &M);
alps::hdf5::archive & serialize(alps::hdf5::archive & ar, std::string const & p, blas::general_matrix<double> const &M);
alps::hdf5::archive & serialize(alps::hdf5::archive & ar, std::string const & p, blas::general_matrix<double> &M);
alps::hdf5::archive & serialize(alps::hdf5::archive & ar, std::string const & p, blas::general_matrix<std::complex<double> > const &M);
alps::hdf5::archive & serialize(alps::hdf5::archive & ar, std::string const & p, blas::general_matrix<std::complex<double> > &M);
}
#endif
