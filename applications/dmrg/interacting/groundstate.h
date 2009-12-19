/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2003-2004 by Salvatore R. Manmana <Salva@theo3.physik.uni-stuttgart.de>,
*                            Reinhard M. Noack <Reinhard.Noack@physik.uni-marburg.de>
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

#ifndef GROUND_H
#define GROUND_H

#include "headers.h"
#include "system.h"
#include "superblock.h"

// vectorspace using wavefunctions in matrix-representation as needed for DMRG:
class DMRGVecspace {
  
 public:
  
  typedef WaveFunction vector_type; 
  typedef double scalar_type;
  typedef int size_type;
  
  DMRGVecspace(size_type n, size_type m):n_(n), m_(m){}
  
  inline size_type vec_dimension() const {
    return n_ * m_;
  }

  void project(vector_type&) const {}
  
  vector_type new_vector() const {
    vector_type psi_new(n_,m_);
    psi_new.clear();
    return psi_new;
  }
  
 private:
  size_type n_; 
  size_type m_;
};


namespace ietl{
  
  // Ian's proposal:
  template <typename Gen>
    inline void generate(WaveFunction& psi, Gen& g)
    {
      typedef WaveFunction::iterator1 iterator1;
      iterator1 iEnd = psi.end1();
      for(iterator1 i = psi.begin1(); i != iEnd; ++i)  // prefer ++i to i++ 
        {
          std::generate(i.begin(), i.end(), g);
        }
    }
  
  // rewrite using iterators:
  inline double dot(const WaveFunction& psi1 , const WaveFunction& psi2) {
    
    double prod = 0.0;
    
    for(unsigned int i = 0; i < psi1.size1(); i++)
      for(unsigned int j = 0; j < psi1.size2(); j++)
        prod += psi1(i,j) * psi2(i,j);
    
    return prod;
    
  }
    
  inline double two_norm(const WaveFunction& psi) {
    
    double norm = sqrt(ietl::dot(psi,psi));
    return norm;
    
  }
  
  inline void scale(WaveFunction& psi, double a){
    psi = a * psi;
  }
  
  inline void add(WaveFunction& psi1, WaveFunction& psi2){
    psi2 = psi1 + psi2;
  }
  
  inline void copy(WaveFunction& psi1, WaveFunction& psi2){
    psi2 = psi1;
  }
  
}// matches namespace ietl

namespace ietl{
  inline void mult(const Superblock& S, 
                   const WaveFunction& psi_in, WaveFunction& psi) { 
    S.multiply(psi_in, psi);
  }
}

#include <ietl/lanczos.h>

double GetGroundState(Superblock S, WaveFunction& psi_new){
  
  typedef boost::lagged_fibonacci607 Gen;  
  Gen mygen;
  
  DMRGVecspace vec(psi_new.size1(),psi_new.size2());
  
  ietl::lanczos<Superblock,DMRGVecspace> lanczos(S,vec); 
  
  // Creation of an iteration object:  
  int max_iter = 10 * psi_new.size1() * psi_new.size2(); 

  double rel_tol = 500*std::numeric_limits<double>::epsilon(), 
    abs_tol = std::pow(std::numeric_limits<double>::epsilon(),2./3);  
  
  int n_lowest_eigenval = 1; 

  std::vector<double> eigen, err; std::vector<int> multiplicity;  
  
  ietl::lanczos_iteration_nlowest<double> 
    iter(max_iter, n_lowest_eigenval, rel_tol, abs_tol);
  
  try{
    lanczos.calculate_eigenvalues(iter,mygen);
    eigen = lanczos.eigenvalues(); err = lanczos.errors(); 
    multiplicity = lanczos.multiplicities();
    std::cerr<<"number of iterations: "<<iter.iterations()<<"\n";
    cerr << "ground state energy = " << eigen[0] << endl;
  }
  catch (std::runtime_error& e) {std::cerr << e.what() << "\n";} 
  
  // call of eigenvectors function follows:   
  std::vector<double>::iterator start = eigen.begin(); std::vector<double>::iterator end = eigen.begin()+1;
  std::vector<WaveFunction> eigenvectors; ietl::Info<double> info; 
  
  try {
    lanczos.eigenvectors(start,end,std::back_inserter(eigenvectors),info,mygen,max_iter); 
  }
  catch (std::runtime_error& e) {std::cerr << e.what() << "\n";}  
  
  psi_new = *eigenvectors.begin(); // fetch the ground state
  
  std::cerr << "Information about the eigenvectors computations:\n";
  for(int i = 0; i < info.size(); i++) {
    std::cerr << "m1(" << i+1 << "): " << info.m1(i) << ", m2(" << i+1 << "): "
              << info.m2(i) << ", ma(" << i+1 << "): " << info.ma(i) << " eigenvalue("
              << i+1 << "): " << info.eigenvalue(i) << " residual(" << i+1 << "): "
              << info.residual(i) << " error_info(" << i+1 << "): "
              << info.error_info(i) << "\n\n";
  
    if(info.error_info(i) != 0){
      cerr << "Could not calculate ground state vector for V1!" << endl;
      return 1;
    }
  }
  
  double energy = *eigen.begin();               
  return energy;
  
} 

#endif
