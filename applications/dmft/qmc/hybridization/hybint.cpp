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

#include"hybint.hpp"


interaction_matrix::interaction_matrix(const alps::params &p){
  extern int global_mpi_rank;
  n_orbitals_=p["N_ORBITALS"];
  val_.resize(n_orbitals_*n_orbitals_,0.);
  //if the parameter U_MATRIX is defined: read in the U_MATRIX from file
  if(p.defined("U_MATRIX")){
    if(p.defined("U") && !global_mpi_rank){ std::cout << "Warning::parameter U_MATRIX defined, ignoring parameter U" << std::flush << std::endl; };
    std::string ufilename=p["U_MATRIX"].cast<std::string>();
    if(p.defined("UMATRIX_IN_HDF5") && p["UMATRIX_IN_HDF5"].cast<bool>()){//attempt to read from h5 archive
      alps::hdf5::archive ar(ufilename, alps::hdf5::archive::READ);
      ar>>alps::make_pvp("/Umatrix",val_);
    }
    else{//read from text file
      std::ifstream u_file(ufilename.c_str());
      if(!u_file.good()) throw std::runtime_error("problem reading in U_MATRIX.");
      double U_ij;
      for(int i=0; i<n_orbitals(); ++i)
        for(int j=0; j<n_orbitals(); ++j){
          u_file>>U_ij;
          operator()(i,j)=U_ij;
          if(!u_file.good()) throw std::runtime_error("problem reading in U_MATRIX.");
        }
    }
  }else{
    if(!p.defined("U")) throw std::invalid_argument("please specify either U (and, optionally, U' and J) or a file U_MATRIX in your parameter file");
    double U=(double)(p["U"]);
    double J=(double)(p["J"]|0.);
    double Uprime = (p["U'"]|(U-2*J));
    assemble(U, Uprime, J);
  }
}

void interaction_matrix::apply_shift(const double shift){
  for(int i=0; i<n_orbitals(); ++i)
    for(int j=0; j<n_orbitals(); ++j)
      if(i!=j) operator()(i,j)+=shift; //apply shift
}

void interaction_matrix::assemble(const double U, const double Uprime, const double J){
  if(Uprime==U && J==0){
     for(int i=0;i<n_orbitals_;++i){
       for(int j=0;j<n_orbitals_;++j){
         operator()(i,j)=(i==j)?0:U;
       }
     }
   }else{
  if(n_orbitals_%2!=0){
    std::cerr<<"n_orbitals is: "<<n_orbitals_<<std::endl;
    throw std::logic_error("extend assemble or write interaction matrix to file for odd # orbitals");
  }
  for(int i=0;i<n_orbitals_;i+=2){
    operator()(i  , i  ) = 0; //Pauli
    operator()(i+1, i+1) = 0; //Pauli
    operator()(i  , i+1) = U; //Hubbard repulsion same band
    operator()(i+1, i  ) = U; //Hubbard repulsion same band
    for(int j=0; j<n_orbitals_; j+=2){
      if(j==i) 
        continue;
      operator()(i  ,j  ) = Uprime-J; //Hubbard repulsion interband same spin 
      operator()(i+1,j+1) = Uprime-J; //Hubbard repulsion interband same spin
      operator()(i  ,j+1) = Uprime; //Hubbard repulsion interband opposite spin (this used to be '+J', the rest of the world uses '-J' -> changed to be consistent).
      operator()(i+1,j  ) = Uprime; //Hubbard repulsion interband opposite spin
    }
  }
}
} 

std::ostream &operator<<(std::ostream &os, const interaction_matrix &U){
  os<<"(effective) U matrix with "<<U.n_orbitals()<<" orbitals: "<<std::endl;
  for(int i=0;i<U.n_orbitals();++i){
    for(int j=0;j<U.n_orbitals();++j){
      os<<U(i,j)<<" ";
    }
    os<<std::endl;
  }
  return os;
}
std::ostream &operator<<(std::ostream &os, const chemical_potential &mu){
  os<<"(effective) chemical potential with "<<mu.n_orbitals()<<" orbitals: "<<std::endl;
  for(std::size_t i=0;i<mu.n_orbitals();++i){
    os<<mu[i]<<" ";
  }
  os<<std::endl;
  return os;
}

