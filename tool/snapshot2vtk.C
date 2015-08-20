/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Library / snapshot2vtk tool
*
* Copyright (C) 2012-2015 by Synge Todo <wistaria@comp-phys.org>
*
* This software is part of the ALPS libraries, published under the ALPS
* Library License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Library License along with
* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

// ALPS snapshot format (xdr):
//
//   lattice_dimension(int)  state_dimension(int)  number_of_sites(int)
//   coordinate_of_site1(vector of double)  state_of_site1(vector of double)
//   coordinate_of_site2(vector of double)  state_of_site2(vector of double)
//   coordinate_of_site3(vector of double)  state_of_site3(vector of double)
//   ...

#include <alps/osiris.h>
#include <alps/utility/vectorio.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/regex.hpp>
#include <complex>
#include <iostream>
#include <string>
#include <vector>

double angle(double x, double y) {
  double g = arg(std::complex<double>(x, y)) / (2 * M_PI);
  return (g >= 0) ? g : (g + 1);
}

int main(int argc, char **argv) {
  if (argc == 1) {
    std::cerr << "Usage: " << argv[0] << " snapshot files\n";
    return 127;
  }

  for (int i = 1; i < argc; ++i) {
    std::string basename = regex_replace(std::string(argv[i]), boost::regex("\\.xdr$"), "");
    boost::filesystem::path xdrfile =
      complete(boost::filesystem::path(basename + ".xdr")).normalize();
    boost::filesystem::path vtkfile =
      complete(boost::filesystem::path(basename + ".vtk")).normalize();
    if (!exists(xdrfile)) {
      std::cerr << "Warning: " << xdrfile << " not found. Skipped.\n";
      continue;
    }
    if (exists(vtkfile)) {
      std::cerr << "Warning: " << vtkfile << " exists. Skipped.\n";
      continue;
    }
      
    alps::IXDRFileDump xdr(xdrfile);
    int dim_lattice, dim_state, num_sites;
    xdr >> dim_lattice >> dim_state >> num_sites;
    if (dim_lattice > 3) {
      std::cerr << "Warning: lattice dimension is too large (>3). Skipped\n";
      continue;
    }
    if (dim_state > 3) {
      std::cerr << "Warning: state dimension is too large (>3). Skipped\n";
      continue;
    }
    std::vector<std::vector<double> > coords;
    std::vector<std::vector<double> > states;
    std::vector<double> coord;
    std::vector<double> state;
    for (int s = 0; s < num_sites; ++s) {
      xdr >> coord >> state;
      coord.resize(3, 0.0);
      state.resize(3, 0.0);
      coords.push_back(coord);
      states.push_back(state);
    }

    // shift center of mass of the lattice to the origin
    std::vector<double> center(3, 0.0);
    for (int s = 0; s < num_sites; ++s) {
      for (int j = 0; j < 3; ++j) center[j] += coords[s][j];
    }
    for (int j = 0; j < 3; ++j) center[j] /= num_sites;
    for (int s = 0; s < num_sites; ++s) {
      for (int j = 0; j < 3; ++j) coords[s][j] -= center[j];
    }

    boost::filesystem::ofstream vtk(vtkfile);
    vtk << "# vtk DataFile Version 3.0\n"
        << "snapshot\n"
        << "ASCII\n"
        << "DATASET UNSTRUCTURED_GRID\n"
        << "POINTS " << num_sites << " float\n";
    for (int s = 0; s < num_sites; ++s) vtk << alps::write_vector(coords[s]) << std::endl;
    if (dim_state == 1) {
      vtk << "POINT_DATA " << num_sites << "\n"
          << "SCALARS spins float\n"
          << "LOOKUP_TABLE default\n";
      for (int s = 0; s < num_sites; ++s) vtk << states[s][0] << std::endl;
    } else {
      vtk << "POINT_DATA " << num_sites << "\n"
          << "VECTORS spins float\n";
      for (int s = 0; s < num_sites; ++s) vtk << alps::write_vector(states[s]) << std::endl;
      vtk << "SCALARS angle float\n"
          << "LOOKUP_TABLE default\n";
      for (int s = 0; s < num_sites; ++s) vtk << angle(states[s][0], states[s][1]) << std::endl;
    }
  }
};
