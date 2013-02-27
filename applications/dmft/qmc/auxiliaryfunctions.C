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

#include "types.h"
#include "green_function.h"
#include "U_matrix.h"
#include <fstream>
#include <sstream>
#include <math.h>



std::ostream &operator<<(std::ostream &os, const multiple_vector_type &v){
  os<<std::setprecision(20);
  for(unsigned int i=0;i<v.first.size();++i)
    os<<i<<"\t"<<v.first[i]<<"\t"<<v.second[i]<<std::endl;
  return os;
}

std::ostream &operator<<(std::ostream &os, const multiple_complex_vector_type &v){
  os<<std::setprecision(20);
  for(unsigned int i=0;i<v.first.size();++i)
    os<<i<<"\t"<<v.first[i].imag()<<"\t"<<v.second[i].imag()<<std::endl;
  return os;
}

std::ostream &operator<<(std::ostream &os, const double_vector &v){
  os<<std::setprecision(20);
  for(unsigned int i=0;i<v.size();++i)
    os<<i<<"\t"<<v[i]<<std::endl;
  return os;
}

std::ostream &operator<<(std::ostream &os, const complex_vector &v){
  os<<std::setprecision(20);
  for(unsigned int i=0;i<v.size();++i)
    os<<i<<"\t"<<v[i].real()<<"\t"<<v[i].imag()<<std::endl;
  return os;
}




void print_imag_green_matsubara(std::ostream &os, const multiple_complex_vector_type &v, const double beta){
  os<<std::setprecision(20);
  for(unsigned int i=0;i<v.first.size();++i)
    os<<(2*i+1)*M_PI/beta<<"\t"<<v.first[i].imag()<<"\t"<<v.second[i].imag()<<std::endl;
}




void print_imag_green_matsubara(std::ostream &os, const matsubara_green_function_t &v, const double beta, const shape_t shape)
{
  os<<std::setprecision(20);
  for(unsigned int i=0;i<v.nfreq();++i){
    os<<(2*i+1)*M_PI/beta<<"\t";
    switch (shape) {
      case diagonal:
        for (unsigned int k=0; k<v.nsite(); ++k) 
          for(unsigned int f=0;f<v.nflavor();++f)
            os<<v(i,k,k,f).imag()<<"\t";
        break;
      case blockdiagonal:
        for (unsigned int k=0; k<v.nsite(); k+=2)
          for (unsigned int p1=0; p1<=1; p1++)
            for (unsigned int p2=0; p2<=1; p2++)
              for(unsigned int f=0;f<v.nflavor();++f)
                os<<v(i,k+p1,k+p2,f).imag()<<"\t";
        break;
      case nondiagonal:
        for (unsigned int s1=0; s1<v.nsite(); ++s1)
          for (unsigned int s2=0; s2<v.nsite(); ++s2) 
            for(unsigned int f=0;f<v.nflavor();++f)
              os<<v(i,s1,s2,f).imag()<<"\t";
    }
    os<<std::endl;
  }
}




void print_selfenergy_matsubara(std::ostream &os, const matsubara_green_function_t &g_bare, 
                                const matsubara_green_function_t g_dressed, const double beta, const shape_t shape)
{
  os<<std::setprecision(20);
  switch (shape) {
    case diagonal:
      for(unsigned int i=0;i<g_bare.nfreq();++i){
        os<<(2*i+1)*M_PI/beta<<"\t";
        for (unsigned int k=0; k<g_bare.nsite(); ++k) 
          for(unsigned int f=0;f<g_bare.nflavor();++f)
            os<<(1./g_bare(i,k,k,f) - 1./g_dressed(i,k,k,f)).real()<<"\t"<<(1./g_bare(i,k,k,f) - 1./g_dressed(i,k,k,f)).imag()<<"\t";
        os<<std::endl;
      }
      break;
    case blockdiagonal:
    {
      throw std::invalid_argument("please use cluster framework for block diag GF");
    }
    case nondiagonal:
    {
      throw std::invalid_argument("please use cluster framework for full nondiag GF");
    }
  }
}




void print_quasiparticle_estimate(std::ostream &os, const multiple_complex_vector_type &g_bare, 
                                  const multiple_complex_vector_type &g_dressed, const double beta){
  os<<"Quasiparticle weight estimate Zeta for Z, using sigma(pi T) and sigma(3 Pi T):"<<std::endl;
  double sigma_1=1./g_bare.first[0].imag() - 1./g_dressed.first[0].imag();
  double sigma_2=1./g_bare.first[1].imag() - 1./g_dressed.first[1].imag();
  os<<"1: "<<"\t"<<1./(1.-sigma_1*beta/M_PI)<<"\t"<<1./(1.-(sigma_2)*beta/M_PI/3.)<<std::endl;
}




void print_quasiparticle_estimate(std::ostream &os, const matsubara_green_function_t &g_bare, 
                                  const matsubara_green_function_t &g_dressed, const double beta)
{
  os<<"Quasiparticle weight estimate Zeta for Z, using sigma(pi T) and sigma(3 Pi T):"<<std::endl;
  for(unsigned int k=0; k<g_bare.nsite(); ++k) {
    for(unsigned int f=0;f<g_bare.nflavor();++f){
      double sigma_1=1./g_bare(0,k,k,f).imag() - 1./g_dressed(0,k,k,f).imag();
      double sigma_2=1./g_bare(1,k,k,f).imag() - 1./g_dressed(1,k,k,f).imag();
      os<<k<<":"<<"\t"<<1./(1.-sigma_1*beta/M_PI)<<"\t"<<1./(1.-(sigma_2)*beta/M_PI/3.)<<std::endl;
    }
  }
}




void print_real_green_matsubara(std::ostream &os, const multiple_complex_vector_type &v, const double beta){
  os<<std::setprecision(20);
  for(unsigned int i=0;i<v.first.size();++i){
    os<<(2*i+1)*M_PI/beta<<"\t"<<v.first[i].real()<<"\t"<<v.second[i].real()<<std::endl;
  }
}




void print_real_green_matsubara(std::ostream &os, const matsubara_green_function_t &v, const double beta, const shape_t shape){
  os<<std::setprecision(20);
  for(unsigned int i=0;i<v.nfreq();++i){
    os<<(2*i+1)*M_PI/beta<<"\t";
    switch (shape) {
      case diagonal:
        for (unsigned int k=0; k<v.nsite(); ++k) 
          for(unsigned int f=0;f<v.nflavor();++f)
            os<<v(i,k,k,f).real()<<"\t";
        break;
      case blockdiagonal:
        for (unsigned int k=0; k<v.nsite(); k+=2)
          for (unsigned int p1=0; p1<=1; p1++)
            for (unsigned int p2=0; p2<=1; p2++)
              for(unsigned int f=0;f<v.nflavor();++f)
                os<<v(i,k+p1,k+p2,f).real()<<"\t";
        break;
      case nondiagonal:
        for (unsigned int s1=0; s1<v.nsite(); ++s1)
          for (unsigned int s2=0; s2<v.nsite(); ++s2) 
            for(unsigned int f=0;f<v.nflavor();++f)
              os<<v(i,s1,s2,f).real()<<"\t";
    }
    os<<std::endl;
  }
}




void print_green_itime(std::ostream &os, const multiple_vector_type &v, const double beta){
  os<<std::setprecision(20);
  for(unsigned int i=0;i<v.first.size();++i){
    os<<i*beta/(v.first.size()-1)<<"\t"<<v.first[i]<<"\t"<<v.second[i]<<std::endl;
  }
}




void print_green_itime(std::ostream &os, const itime_green_function_t &v, const double beta, const shape_t shape){
  os<<std::setprecision(20);
  for(unsigned int i=0;i<v.ntime();++i){
    os<<beta*i/(v.ntime()-1)<<"\t";
    switch (shape) {
      case diagonal:
        for (unsigned int k=0; k<v.nsite(); ++k) 
          for(unsigned int f=0;f<v.nflavor();++f)
            os<<v(i,k,k,f)<<"\t";
        break;
      case blockdiagonal:
        for (unsigned int k=0; k<v.nsite(); k+=2)
          for (unsigned int p1=0; p1<=1; p1++)
            for (unsigned int p2=0; p2<=1; p2++)
              for(unsigned int f=0;f<v.nflavor();++f)
                os<<v(i,k+p1,k+p2,f)<<"\t";
        break;
      case nondiagonal:
        for (unsigned int s1=0; s1<v.nsite(); ++s1)
          for (unsigned int s2=0; s2<v.nsite(); ++s2) 
            for(unsigned int f=0;f<v.nflavor();++f)
              os<<v(i,s1,s2,f)<<"\t";
    }
    os<<std::endl;
  }
}

void print_all_green_functions(std::string const &basename, const int iteration_ctr, const matsubara_green_function_t &G0_omega, 
                               const matsubara_green_function_t &G_omega, const itime_green_function_t &G0_tau, 
                               const itime_green_function_t &G_tau, const double beta, const shape_t shape, 
                               const std::string suffix)
{
  std::ostringstream G0omega_name, G0omegareal_name, G0tau_name, Gomega_name, 
  Gomegareal_name, Gtau_name, Gtau_name_2, Gomega_name_2, selfenergy_name;
  G0omega_name<<"G0_omega"<<suffix<<"_"<<iteration_ctr;//<<"_"<<process_id;
  G0omegareal_name<<"G0_omegareal"<<suffix<<"_"<<iteration_ctr;//<<"_"<<process_id;
  G0tau_name<<"G0_tau"<<suffix<<"_"<<iteration_ctr;//<<"_"<<process_id;
  Gomega_name<<"G_omega"<<suffix<<"_"<<iteration_ctr;//<<"_"<<process_id;
  Gomegareal_name<<"G_omegareal"<<suffix<<"_"<<iteration_ctr;//<<"_"<<process_id;
  Gtau_name<<"G_tau"<<suffix<<"_"<<iteration_ctr;//<<"_"<<process_id;
  Gomega_name_2<<"G_omega"<<suffix;//<<"_"<<process_id;
  Gtau_name_2<<"G_tau"<<suffix;//<<"_"<<process_id;
  selfenergy_name<<"selfenergy"<<suffix<<"_"<<iteration_ctr;//<<"_"<<process_id;
  std::ofstream G0_omega_file(G0omega_name.str().c_str());
  std::ofstream G0_omegareal_file(G0omegareal_name.str().c_str());
  std::ofstream G0_tau_file(G0tau_name.str().c_str());
  std::ofstream G_omega_file(Gomega_name.str().c_str());
  std::ofstream G_omegareal_file(Gomegareal_name.str().c_str());
  std::ofstream G_tau_file(Gtau_name.str().c_str());
  std::ofstream G_omega_file_2(Gomega_name_2.str().c_str());
  std::ofstream G_tau_file_2(Gtau_name_2.str().c_str());
  std::ofstream selfenergy_file(selfenergy_name.str().c_str());
  assert(G0_omega_file.is_open() && G0_tau_file.is_open() && G_omega_file.is_open() && G_tau_file.is_open());
  print_imag_green_matsubara(G0_omega_file, G0_omega, beta, shape);
  print_real_green_matsubara(G0_omegareal_file,G0_omega,beta, shape);
  print_green_itime(G0_tau_file,G0_tau, beta, shape);
  print_imag_green_matsubara(G_omega_file,G_omega,beta, shape);
  print_real_green_matsubara(G_omegareal_file,G_omega,beta, shape);
  print_green_itime(G_tau_file,G_tau,beta, shape);
  print_imag_green_matsubara(G_omega_file_2,G_omega,beta, shape);
  print_green_itime(G_tau_file_2,G_tau,beta, shape);
  print_selfenergy_matsubara(selfenergy_file, G0_omega, G_omega, beta, shape);
  if (shape==diagonal)
    print_quasiparticle_estimate(std::cout, G_omega, G0_omega, beta);
  
  alps::hdf5::archive ar(basename+".h5", "a");
  //writeout into hf5 file, using /simulation/iteration/ path
  std::stringstream basepath; basepath<<"/simulation/iteration/"<<iteration_ctr<<"/results/";
  G_tau.write_hdf5(ar,basepath.str()+"G_tau");
  G_omega.write_hdf5(ar,basepath.str()+"G_omega");
  G0_tau.write_hdf5(ar,basepath.str()+"G0_tau");
  G0_omega.write_hdf5(ar,basepath.str()+"G0_omega");
  
  //writeout into hf5 file, using /simulation/results/ path
  std::stringstream basepath2; basepath2<<"/simulation/results/";
  G_tau.write_hdf5(ar,basepath2.str()+"G_tau");
  G_omega.write_hdf5(ar,basepath2.str()+"G_omega");
  G0_tau.write_hdf5(ar,basepath2.str()+"G0_tau");
  G0_omega.write_hdf5(ar,basepath2.str()+"G0_omega");
}


void print_tau_green_functions(std::string const &basename, const int iteration_ctr, const itime_green_function_t &G0_tau, const itime_green_function_t &G_tau, const double beta,
                               const shape_t shape, const std::string suffix){
  std::ostringstream G0tau_name, Gtau_name, Gtau_name_2;
  G0tau_name<<"G0_tau_"<<iteration_ctr;//<<"_"<<process_id;
  Gtau_name<<"G_tau_"<<iteration_ctr;//<<"_"<<process_id;
  Gtau_name_2<<"G_tau";//<<"_"<<process_id;
  std::ofstream G0_tau_file(G0tau_name.str().c_str());
  std::ofstream G_tau_file(Gtau_name.str().c_str());
  std::ofstream G_tau_file_2(Gtau_name_2.str().c_str());
  print_green_itime(G0_tau_file,G0_tau,beta,shape);
  print_green_itime(G_tau_file,G_tau,beta,shape);
  print_green_itime(G_tau_file_2,G_tau,beta,shape);
  
  alps::hdf5::archive ar(basename+".h5", "a");
  //writeout into hf5 file, using /simulation/iteration/ path
  std::stringstream basepath; basepath<<"/simulation/iteration/"<<iteration_ctr<<"/results/";
  G_tau.write_hdf5(ar,basepath.str()+"G_tau");
  G0_tau.write_hdf5(ar,basepath.str()+"G0_tau");
  
  //writeout into hf5 file, using /simulation/results/ path
  std::stringstream basepath2; basepath2<<"/simulation/results/";
  G_tau.write_hdf5(ar,basepath2.str()+"G_tau");
  G0_tau.write_hdf5(ar,basepath2.str()+"G0_tau");
}

void print_dressed_tau_green_functions(std::string const &basename, const int iteration_ctr, const itime_green_function_t &G_tau, const double beta, 
                                       const shape_t shape, const std::string suffix){
  std::ostringstream Gtau_name, Gtau_name_2;
  Gtau_name<<"G_tau"<<suffix<<"_"<<iteration_ctr;//<<"_"<<process_id;
  Gtau_name_2<<"G_tau"<<suffix;//<<"_"<<process_id;
  std::ofstream G_tau_file(Gtau_name.str().c_str());
  std::ofstream G_tau_file_2(Gtau_name_2.str().c_str());
  print_green_itime(G_tau_file,G_tau,beta,shape);
  print_green_itime(G_tau_file_2,G_tau,beta,shape);
  
  alps::hdf5::archive ar(basename+".h5", "a");
  //writeout into hf5 file, using /simulation/iteration/ path
  std::stringstream basepath; basepath<<"/simulation/iteration/"<<iteration_ctr<<"/results/";
  G_tau.write_hdf5(ar,basepath.str()+"G_tau");
  
  //writeout into hf5 file, using /simulation/results/ path
  std::stringstream basepath2; basepath2<<"/simulation/results/";
  G_tau.write_hdf5(ar,basepath2.str()+"G_tau");
}

std::ostream &operator<<(std::ostream &os, const green_function<double> &v){
  os<<std::setprecision(20);
  for(unsigned int o=0;o<v.nfreq();++o){
    os<<o<<"\t";
    for(unsigned int i=0;i<v.nsite();++i)
      for(unsigned int j=0;j<v.nsite();++j)
        for(unsigned int z=0;z<v.nflavor();++z)
          os<<v(o,i,j,z)<<"\t";
    os<<std::endl;
    
  }
  return os;
}
std::istream &operator>>(std::istream &is, green_function<double> &v){
  double index;
  for(unsigned int o=0;o<v.nfreq();++o){
    is>>index;
    for(unsigned int i=0;i<v.nsite();++i)
      for(unsigned int j=0;j<v.nsite();++j)
        for(unsigned int z=0;z<v.nflavor();++z)
          is>>v(o,i,j,z);
  }
  return is;
}

/*std::istream &operator>>(std::istream &is, green_function<std::complex<double> > &v){
 double index;
 for(unsigned int o=0;o<v.nfreq();++o){
 is>>index;
 for(unsigned int i=0;i<v.nsite();++i)
 for(unsigned int j=0;j<v.nsite();++j)
 for(unsigned int z=0;z<v.nflavor();++z) {
 double re,im;
 is >>re >> im;
 v(o,i,j,z)=std::complex<double>(re,im);
 }
 }
 return is;
 }*/

std::ostream &operator<<(std::ostream &os, const green_function<std::complex<double> > &v){
  os<<std::setprecision(20);
  for(unsigned int o=0;o<v.nfreq();++o){
    os<<o<<"\t";
    for(unsigned int i=0;i<v.nsite();++i)
      for(unsigned int j=0;j<v.nsite();++j)
        for(unsigned int z=0;z<v.nflavor();++z)
          os<<(v(o,i,j,z).real())<<"\t"<<(v(o,i,j,z).imag())<<"\t";
    os<<std::endl;    
  }
  return os;
}

inline double G_times_minusG(const int i, const int N, const multiple_vector_type &G_tau){
  return G_tau.first[i]*(-1 -G_tau.first[N-i]);
}
std::ostream &operator<<(std::ostream &os, const U_matrix &U){
  os<<"[ ";
  for(spin_t flavor_i=0;flavor_i<U.nf();++flavor_i){
    os<<"[ ";
    for(spin_t flavor_j=0;flavor_j<U.nf();++flavor_j){
      os<<U(flavor_i,flavor_j)<<"\t";
    }
    os<<"] "<<std::endl;
  }
  os<<"]"<<std::endl;
  return os;
}
/*alps::hdf5::archive& alps::hdf5::save(alps::hdf5::archive& ar, const std::string& basepath, const matsubara_green_function_t &gf){
 std::cout<<"trying to store gf to: "<<basepath<<std::endl;
 }
 alps::hdf5::archive& alps::hdf5::save(alps::hdf5::archive& ar, const std::string& basepath, const itime_green_function_t&gf){
 std::cout<<"trying to store gf to: "<<basepath<<std::endl;
 }*/
