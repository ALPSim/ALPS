 /*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Sebastian Fuchs <fuchs@theorie.physik.uni-goettingen.de>
 *                              Matthias Troyer <troyer@comp-phys.org>
 *
 *
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/


//Data structure for repulsion for density density bands.
#ifndef U_MATRIX_H
#define U_MATRIX_H
#include "types.h"
#include "alps/parameter.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

/// The density-density interaction U(i,j) n_i n_j / 2 between flavors i, j.
///
/// It is read from the file U_MATRIX, which holds the dense FLAVORS x FLAVORS
/// matrix row by row, as the hybridization and ct-hyb-segment solvers read
/// it. Without U_MATRIX it is built from U, and for more than two flavors
/// from U, J and U' (default U - 2J).
class U_matrix{
public:
  U_matrix(const alps::Parameters &parms) :
    ns_(parms.value_or_default("SITES", 1)),
    nf_(parms.value_or_default("FLAVORS", 2)),
    val_(nf_*nf_, 0.), //default: non-interacting.
    n_nonzero_(0)
  {
    if(parms.defined("U_MATRIX")){
      read(parms["U_MATRIX"]);
    } else if (nf_==1) {
      //special case: only 1 orbital
      operator()(0,0)=static_cast<double>(require(parms, "U"));
    }else if (nf_==2) {
      //your ordinary two site problem
      double U=static_cast<double>(require(parms, "U"));
      operator()(0,0)=0; operator()(1,1)=0;
      operator()(0,1)=U; operator()(1,0)=U;
    } else {
      double U=static_cast<double>(require(parms, "U"));
      double J=static_cast<double>(require(parms, "J"));
      double Uprime = parms.value_or_default("U'", U-2*J);
      assemble(U, Uprime, J);
    }
    for (unsigned i=0; i<nf_*nf_; ++i)
      if (val_[i]!=0) 
        n_nonzero_++;
  }
  

  void assemble(const double U, const double Uprime, const double J){
    //this implements the U matrix for the special case of n_flavor/2 degenerate bands
    if (ns_!=1 || nf_%2!=0)
      throw std::invalid_argument("U_matrix: U, J and U' build the interaction for SITES = 1 and an even "
                                  "number of FLAVORS; give U_MATRIX otherwise.");
    for(spin_t i=0;i<nf_;i+=2){
      operator()(i  , i  ) = 0; //Pauli
      operator()(i+1, i+1) = 0; //Pauli
      operator()(i  , i+1) = U; //Hubbard repulsion same band
      operator()(i+1, i  ) = U; //Hubbard repulsion same band
      for(spin_t j=0; j<nf_; j+=2){
        if(j==i) 
          continue;
        operator()(i  ,j  ) = Uprime-J; //Hubbard repulsion interband same spin 
        operator()(i+1,j+1) = Uprime-J; //Hubbard repulsion interband same spin
        operator()(i  ,j+1) = Uprime; //Hubbard repulsion interband opposite spin (this used to be '+J', the rest of the world uses '-J' -> changed to be consistent).
        operator()(i+1,j  ) = Uprime; //Hubbard repulsion interband opposite spin
      }
    }
  } 
  
  double &operator()(spin_t flavor_i, spin_t flavor_j){
    return val_[flavor_i*nf_+flavor_j];
    }
  
  const double &operator() (spin_t flavor_i, spin_t flavor_j)const {
    return val_[flavor_i*nf_+flavor_j];
  }
  
  spin_t nf()const {return nf_;}
  spin_t ns()const {return ns_;}

  /// The shift of the orbital energy of flavor f that puts half filling at
  /// MU = 0: sum_j U(f,j) / 2. It differs between flavors when the orbitals
  /// are not equivalent.
  double mu_shift(spin_t f) const {
    double shift=0;
    for (spin_t j=0; j<nf_; ++j)
      shift += operator()(f,j);
    return shift/2;
  }
  /// Whether every flavor has the same mu_shift.
  bool uniform_mu_shift() const {
    for (spin_t f=1; f<nf_; ++f)
      if (mu_shift(f)!=mu_shift(0))
        return false;
    return true;
  }

  inline int n_nonzero() const{return n_nonzero_;}

private:
  alps::StringValue require(const alps::Parameters &parms, const std::string &name) const {
    if (!parms.defined(name))
      throw std::invalid_argument("U_matrix: FLAVORS = " + std::to_string(nf_)
                                  + " needs " + name + ", or the interaction matrix as U_MATRIX.");
    return parms[name];
  }

  /// Reads the dense matrix, row by row, and checks it.
  void read(const std::string &file) {
    std::ifstream u_file(file.c_str());
    if (!u_file.is_open())
      throw std::runtime_error("U_MATRIX: cannot open the file " + file + ".");
    std::vector<double> entries;
    std::string word;
    while (u_file >> word) {
      std::istringstream in(word);
      double x;
      if (!(in >> x) || !in.eof())
        throw std::runtime_error("U_MATRIX " + file + ": entry " + std::to_string(entries.size()+1)
                                 + " (\"" + word + "\") is not a number.");
      entries.push_back(x);
    }
    const std::size_t n = std::size_t(nf_)*nf_;
    if (entries.size() != n) {
      std::ostringstream msg;
      msg << "U_MATRIX " << file << ": has " << entries.size() << " entries; the dense matrix for FLAVORS = "
          << nf_ << " has " << n << ", row by row.";
      if (looks_like_index_triples(entries))
        msg << " The file looks like lines \"i j U_ij\", which dmft read before; write the dense matrix "
               "instead, as the hybridization and ct-hyb-segment solvers read it.";
      throw std::runtime_error(msg.str());
    }
    val_ = entries;
    for (spin_t i=0; i<nf_; ++i) {
      if (operator()(i,i) != 0)
        throw std::runtime_error("U_MATRIX " + file + ": the diagonal entry U(" + std::to_string(i) + ","
                                 + std::to_string(i) + ") is not 0; a flavor does not interact with itself.");
      for (spin_t j=0; j<i; ++j)
        if (operator()(i,j) != operator()(j,i))
          throw std::runtime_error("U_MATRIX " + file + ": U(" + std::to_string(i) + "," + std::to_string(j)
                                   + ") and U(" + std::to_string(j) + "," + std::to_string(i)
                                   + ") differ; the interaction matrix is symmetric.");
    }
  }

  /// Whether the entries read as lines "i j U_ij" with flavor indices.
  bool looks_like_index_triples(std::vector<double> const &entries) const {
    if (entries.empty() || entries.size() % 3 != 0)
      return false;
    for (std::size_t k=0; k<entries.size(); k+=3)
      for (std::size_t c=0; c<2; ++c) {
        const double index = entries[k+c];
        if (index != std::floor(index) || index < 0 || index >= nf_)
          return false;
      }
    return true;
  }

  site_t ns_;
  spin_t nf_;
  std::vector<double> val_;
  int n_nonzero_;
};

std::ostream &operator<<(std::ostream &os, const U_matrix &U);
#endif
//U_MATRIX_H
