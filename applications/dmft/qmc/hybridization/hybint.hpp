/****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2012 by Emanuel Gull <gull@pks.mpg.de>,
 *                   
 *  based on an earlier version by Philipp Werner and Emanuel Gull
 *
 *
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 * 
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
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

#ifndef HYB_INT_HPP
#define HYB_INT_HPP
#include <vector>
#include <fstream>
#include <alps/ngs.hpp>
#include "../green_function.h"

//the interaction matrix keeps track of the impurity density-density interactions (other interactions are not possible in this code). Two general methods: either specify U (and optionally J and U'), or define a matrix and write it into a file, from where it is read in.
class interaction_matrix{
public:
  //constructor
  interaction_matrix(const alps::params &p);
  //element access routines
  double &operator()(std::size_t flavor_i, std::size_t flavor_j){
    return val_[flavor_i*n_orbitals_+flavor_j];
  }
  
  const double &operator() (std::size_t flavor_i, std::size_t flavor_j)const {
    return val_[flavor_i*n_orbitals_+flavor_j];
  }
  const int &n_orbitals() const{return n_orbitals_;}
  void apply_shift(const double shift);
private:
  //get the usual density-density form of the interaction
  void assemble(const double U, const double Uprime, const double J);
  std::vector<double> val_; //interaction matrix values
  int n_orbitals_;
};

//This class handles the 'mu' term (chemical potential term). In the simplest case, 'mu' is just a constant for all orbitals.
//The values can be different through a magnetic or christal field, or due to a double counting which is orbitally dependent.
class chemical_potential{
public:
  chemical_potential(const alps::params &p){
    extern int global_mpi_rank;
    val_.resize(p["N_ORBITALS"], p["MU"]|0.);
    if(p.defined("MU_VECTOR")){
      if(p.defined("MU") && !global_mpi_rank){ std::cout << "Warning::parameter MU_VECTOR defined, ignoring parameter MU" << std::flush << std::endl; };
      std::string mufilename=p["MU_VECTOR"].cast<std::string>();
      if(p.defined("MU_IN_HDF5") && p["MU_IN_HDF5"].cast<bool>()){//attempt to read from h5 archive
        alps::hdf5::archive ar(mufilename, alps::hdf5::archive::READ);
        ar>>alps::make_pvp("/MUvector",val_);
      }
      else{//read from text file
        std::ifstream mu_file(mufilename.c_str());
        if(!mu_file.good()) throw std::runtime_error("Problem reading in MU_VECTOR.");
        std::size_t i=0;
        double MU_i;
        for(i=0; i<n_orbitals(); ++i){
          mu_file>>MU_i;
          val_[i]=MU_i;
          if(!mu_file.good()) throw std::runtime_error("Problem reading in MU_VECTOR.");
        }
      }
    }
  }
  std::size_t n_orbitals(void)const { return val_.size(); }
  const double &operator[] (std::size_t flavor)const {
    return val_[flavor];
  }
  void apply_shift(const double shift){
    for(std::size_t i=0; i<n_orbitals(); ++i)
      val_[i]+=shift; //apply shift
  }
private:
  std::vector<double> val_;
};
std::ostream &operator<<(std::ostream &os, const interaction_matrix &U);
std::ostream &operator<<(std::ostream &os, const chemical_potential &mu);

#endif
