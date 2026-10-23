/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1998-2018 ALPS Collaboration
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

// Adapted from CQMP/Maxent src/maxent_kernel.hpp, revision
// e727650a1f0a427e95efebc409d907b4c84011ec, on 2026-10-14.
// Uses ALPS's HDF5 archive and retains ALPS input/output conventions.
#pragma once
#include<vector>
#include <alps/params.hpp>
#include<iostream>
#include"maxent_matrix_def.hpp"
///enum that enumerates if we're in time or in frequency
enum dataspace_type{
  time_dataspace,
  frequency_dataspace,
  legendre_dataspace
};
///we have a range of different kernels that all need different input parameters. This enum enumerates them.
enum kernel_type{
  time_fermionic_kernel,
  time_bosonic_kernel,
  time_boris_kernel,
  time_fermionic_legendre_kernel,
  time_bosonic_legendre_kernel, //TODO determine boris->legendre
  legendre_fermionic_kernel,
  legendre_bosonic_kernel,
  frequency_fermionic_ph_kernel,
  frequency_bosonic_ph_kernel,
  frequency_anomalous_ph_kernel,
  frequency_fermionic_kernel,
  frequency_bosonic_kernel,
  frequency_anomalous_kernel
};
class kernel{
public:
  kernel(alps::params &p, const vector_type& freq, vector_type &inputGrid);

  ///getter function for the kernel matrix
  const matrix_type &operator()()const{return K_;}
  ///getter function for kernel type
  kernel_type getKernelType() const{return ktype_;}
  ///getter function for dataspace type
  dataspace_type getDataspaceType() const{return dtype_;}
	///getter function for tau_points
	const vector_type getTauPoints() const{return tau_points_;}
    
private:
  ///figure out which kernel is to be used
  void set_kernel_type(const std::string &dataspace_name, const std::string &kernel_name,
                       bool ph_symmetry);
  ///set up kernel with the legendre transform
  void setup_legendre_kernel(const alps::params &p, const vector_type& freq, const int lmax);
  ///number; of Matsubara points
  int ndat_;
  ///number of real frequency points
  int nfreq_;
  ///temperature
  double T_;
  ///type of the kernel (see enums above)
  kernel_type ktype_;
  ///type of the dataspace
  dataspace_type dtype_;
  ///kernel matrix
  matrix_type K_;
	///array of tau points
	vector_type tau_points_;
};
