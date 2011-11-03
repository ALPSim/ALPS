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


//Data structure for repulsion for density density bands.
#ifndef U_MATRIX_H
#define U_MATRIX_H
#include "types.h"
#include "alps/parameter.h"

class U_matrix{
public:
  U_matrix(const alps::Parameters &parms) :
    ns_(parms.value_or_default("SITES", 1)),
    nf_(parms.value_or_default("FLAVORS", 2)),
      n_nonzero_(0), mu_shift_(0)
  {
    val_ = new double[nf_*nf_];
    for(unsigned i=0; i<nf_*nf_; ++i) 
      val_[i]=0; //default: non-interacting.
    if(parms.defined("U_MATRIX")){
      std::string ufilename(parms["U_MATRIX"]);
      std::ifstream u_file(ufilename.c_str());
      assert(u_file.is_open());
      int i;
      int j;
      double U_ij;
      while(u_file>>i>>j>>U_ij){
        operator()(i,j)=U_ij;
      }
    } else if (nf_==1) {
      //special case: only 1 orbital
      operator()(0,0)=(double)(parms["U"]);
    }else if (nf_==2) {
      //your ordinary two site problem
      assert(parms.defined("U"));
      double U=(double)(parms["U"]);
      operator()(0,0)=0; operator()(1,1)=0;
      operator()(0,1)=U; operator()(1,0)=U;
    } else {
      assert(parms.defined("U") && parms.defined("J"));     
      double U=(double)(parms["U"]);
      double J=(double)(parms["J"]);
      double Uprime = parms.value_or_default("U'", U-2*J);
      assemble(U, Uprime, J);
    }
    for (unsigned i=0; i<nf_*nf_; ++i)
      if (val_[i]!=0) 
        n_nonzero_++;
    for (unsigned i=0; i<nf_; ++i) 
      mu_shift_ += operator()(i,0);
    mu_shift_ /= 2;
  }
  

  void assemble(const double U, const double Uprime, const double J){
    //this implements the U matrix for the special case of n_flavor/2 degenerate bands
    assert(ns_==1);
    assert(nf_%2==0);
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
  
  ~U_matrix(){
    delete[] val_;
  }
  
  double &operator()(spin_t flavor_i, spin_t flavor_j){
    return val_[flavor_i*nf_+flavor_j];
    }
  
  const double &operator() (spin_t flavor_i, spin_t flavor_j)const {
    return val_[flavor_i*nf_+flavor_j];
  }
  
  spin_t nf()const {return nf_;}
  spin_t ns()const {return ns_;}
  double mu_shift() const { return mu_shift_; }

  inline int n_nonzero() const{return n_nonzero_;}

private:
  double *val_;
  site_t ns_;
  spin_t nf_;
  int n_nonzero_;
  double mu_shift_;
};

std::ostream &operator<<(std::ostream &os, const U_matrix &U);
//U_MATRIX_H
#endif 
