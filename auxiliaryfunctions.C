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
#include "nambu_green_function.h"
#include <fstream>
#include <sstream>
#include <math.h>
#include "alps/numeric/matrix.hpp"

#include "alps/hdf5/archive.hpp"
#include "alps/hdf5/vector.hpp"
#include "alps/hdf5/pointer.hpp"
#include "alps/hdf5/complex.hpp"

#include "blas_classes/blasheader.h"
#include "blas_classes/general_matrix.h"


void matrix_inversion(alps::numeric::matrix<std::complex<double> >&A, alps::numeric::matrix<std::complex<double> >&A_inv);
void print_imag_green_matsubara(std::ostream &os, const matsubara_green_function_t &v, const double beta)
{
  os<<std::setprecision(14);
  if(!v.get_header().empty()){
	  os<<v.get_header()<<std::endl;
  }
  for(unsigned int i=0;i<v.nfreq();++i){
    os<<(2*i+1)*M_PI/beta<<" ";
    if (v.shape()==diagonal) {
      for (unsigned int k=0; k<v.nsite(); ++k)
        for(unsigned int f=0;f<v.nflavor();++f)
          os<<v(i,k,k,f).imag()<<" ";
    }else if(v.shape()==blockdiagonal){
      for (unsigned int k=0; k<v.nsite(); k+=2)
        for (unsigned int p1=0; p1<=1; p1++)
          for (unsigned int p2=0; p2<=1; p2++)
            for(unsigned int f=0;f<v.nflavor();++f)
              os<<v(i,k+p1,k+p2,f).imag()<<" ";
    }else if(v.shape()==nondiagonal){
      for (unsigned int s1=0; s1<v.nsite(); ++s1)
        for (unsigned int s2=0; s2<v.nsite(); ++s2)
          for(unsigned int f=0;f<v.nflavor();++f)
            os<<v(i,s1,s2,f).imag()<<" ";
    }
    os<<std::endl;
  }
}

void matrix_inversion(alps::numeric::matrix<std::complex<double> >&A, alps::numeric::matrix<std::complex<double> >&A_inv){
  fortran_int_t size=(fortran_int_t)(A.num_rows());
  if(size==0){
    return;
  }
  //build the unit matrix
  A_inv=alps::numeric::matrix<std::complex<double> >::identity_matrix(size);
  //destructive in-place solver - will overwrite unit matrix
  if(size==1){
    A_inv(0,0)=1./A(0,0);
  }else{
    int ipiv[size];
    int info;
#ifdef MY_BLAS_DEFS
    using namespace lapack;
#endif
   //boost::numeric::bindings::lapack::
    zgesv_(&size, &size, &(A(0,0)),&size, ipiv, &(A_inv(0,0)), &size, &info);
    if(info!=0) std::cerr<<"gesv went wrong!!"<<std::endl;
  }
}



void print_selfenergy_matsubara(std::ostream &os, const matsubara_green_function_t &g_bare,
                                const matsubara_green_function_t g_dressed, const double beta)
{
  os<<std::setprecision(14);
  bool causal = true;
    if(g_dressed.shape()==diagonal){
     for(unsigned int i=0;i<g_bare.nfreq();++i){
      os<<(2*i+1)*M_PI/beta<<"\t";
      for (unsigned int k=0; k<g_bare.nsite(); ++k){
        for(unsigned int f=0;f<g_bare.nflavor();++f){
          os<<(1./g_bare(i,k,k,f) - 1./g_dressed(i,k,k,f)).real()<<"\t"<<(1./g_bare(i,k,k,f) - 1./g_dressed(i,k,k,f)).imag()<<"\t";
          if( i < g_bare.nfreq()/2. && (1./g_bare(i,k,k,f) - 1./g_dressed(i,k,k,f)).imag() > 0){
        	  causal = false;
          }
        }
      }
      os<<std::endl;
     }
    }else if(g_dressed.shape()==blockdiagonal){
      
      alps::numeric::matrix<std::complex<double> > Gup(2,2), G0up(2,2), Gdn(2,2), G0dn(2,2);
      alps::numeric::matrix<std::complex<double> > Ginvup(2,2), G0invup(2,2), Ginvdn(2,2), G0invdn(2,2);
      alps::numeric::matrix<std::complex<double> > Sigmaup(2,2);
      alps::numeric::matrix<std::complex<double> > Sigmadn(2,2);
      for(unsigned int i=0;i<g_bare.nfreq();++i){
       if(g_bare.nflavor()!=2) throw std::logic_error("not understood how to handle blockdiag with >2 flavors.");
       os<<(2*i+1)*M_PI/beta<<" ";
       for (unsigned int k=0; k<g_bare.nsite(); k+=2){
        int q=k+1;
        
        //initialize G
        Gup(0,0)=g_dressed(i,k,k,0); Gup(0,1)=g_dressed(i,k,q,0); Gup(1,0)=g_dressed(i,q,k,0); Gup(1,1)=g_dressed(i,q,q,0);
        Gdn(0,0)=g_dressed(i,k,k,1); Gdn(0,1)=g_dressed(i,k,q,1); Gdn(1,0)=g_dressed(i,q,k,1); Gdn(1,1)=g_dressed(i,q,q,1);
        
        //initialize G0
        G0up(0,0)=g_bare(i,k,k,0); G0up(0,1)=g_bare(i,k,q,0); G0up(1,0)=g_bare(i,q,k,0); G0up(1,1)=g_bare(i,q,q,0);
        G0dn(0,0)=g_bare(i,k,k,1); G0dn(0,1)=g_bare(i,k,q,1); G0dn(1,0)=g_bare(i,q,k,1); G0dn(1,1)=g_bare(i,q,q,1);
        
        //perform inversion
        matrix_inversion(Gup, Ginvup);
        matrix_inversion(G0up, G0invup);
        matrix_inversion(Gdn, Ginvdn);
        matrix_inversion(G0dn, G0invdn);
        
        //get sigma
        for(unsigned int i=0;i<2;++i){
          for(unsigned int j=0;j<2;++j){
            Sigmaup(i,j)=G0invup(i,j)-Ginvup(i,j);
            Sigmadn(i,j)=G0invdn(i,j)-Ginvdn(i,j);
          }
        }
        for(unsigned int i=0;i<2;++i){
          for(unsigned int j=0;j<2;++j){
            os<<Sigmaup(i,j).real()<<" "<<Sigmaup(i,j).imag()<<"  "<<Sigmadn(i,j).real()<<" "<<Sigmadn(i,j).imag()<<" ";
          }
        }
      }
      os<<std::endl;
     }
    }else if(g_dressed.shape()==nondiagonal){
      int nsite=g_dressed.nsite();
      alps::numeric::matrix<std::complex<double> > Gup(nsite,nsite), G0up(nsite,nsite), Gdn(nsite,nsite), G0dn(nsite,nsite);
      alps::numeric::matrix<std::complex<double> > Ginvup(nsite,nsite), G0invup(nsite,nsite), Ginvdn(nsite,nsite), G0invdn(nsite,nsite);
      alps::numeric::matrix<std::complex<double> > Sigmaup(nsite,nsite);
      alps::numeric::matrix<std::complex<double> > Sigmadn(nsite,nsite);
      for(unsigned int w=0;w<g_dressed.nfreq();++w){
        os<<(2*w+1)*M_PI/beta<<" ";
        for(int i=0;i< nsite;++i){
          for(int j=0;j< nsite;++j){
            Gup(i  , j)=g_dressed(w,i,j,0);
            Gdn(i  , j)=g_dressed(w,i,j,1);
            G0up(i  , j)=g_bare(w,i,j,0);
            G0dn(i  , j)=g_bare(w,i,j,1);
          }
        }
        matrix_inversion(Gup, Ginvup);
        matrix_inversion(G0up, G0invup);
        matrix_inversion(Gdn, Ginvdn);
        matrix_inversion(G0dn, G0invdn);
        
        for(int i=0;i<nsite;++i){
          for(int j=0;j<nsite;++j){
            Sigmaup(i,j)=G0invup(i,j)-Ginvup(i,j);
            Sigmadn(i,j)=G0invdn(i,j)-Ginvdn(i,j);
          }
        }
        for(int i=0;i<nsite;++i){
          for(int j=0;j<nsite;++j){
            os<<Sigmaup(i,j).real()<<" "<<Sigmaup(i,j).imag()<<"  "<<Sigmadn(i,j).real()<<" "<<Sigmadn(i,j).imag()<<" ";
          }
        }
        os<<std::endl;
      }
    }

    if(!causal){
    	std::cout<<"Warning! Elements of the imaginary part of the self-energy are positive!"<<std::endl;
    }
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







void print_real_green_matsubara(std::ostream &os, const matsubara_green_function_t &v, const double beta){
  os<<std::setprecision(14);
  if(!v.get_header().empty()){
	  os<<v.get_header()<<std::endl;
  }
  for(unsigned int i=0;i<v.nfreq();++i){
    os<<(2*i+1)*M_PI/beta<<"\t";
    if(v.shape()==diagonal){
      for (unsigned int k=0; k<v.nsite(); ++k)
        for(unsigned int f=0;f<v.nflavor();++f)
          os<<v(i,k,k,f).real()<<"\t";
    }else if(v.shape()==blockdiagonal){
      for (unsigned int k=0; k<v.nsite(); k+=2)
        for (unsigned int p1=0; p1<=1; p1++)
          for (unsigned int p2=0; p2<=1; p2++)
            for(unsigned int f=0;f<v.nflavor();++f)
              os<<v(i,k+p1,k+p2,f).real()<<"\t";
    }else if(v.shape()==nondiagonal){
      for (unsigned int s1=0; s1<v.nsite(); ++s1)
        for (unsigned int s2=0; s2<v.nsite(); ++s2)
          for(unsigned int f=0;f<v.nflavor();++f)
            os<<v(i,s1,s2,f).real()<<"\t";
    }
    os<<std::endl;
  }
}








void print_green_itime(std::ostream &os, const itime_green_function_t &v, const double beta){
  os<<std::setprecision(14);
  if(!v.get_header().empty()){
	  os<<v.get_header()<<std::endl;
  }
  for(unsigned int i=0;i<v.ntime();++i){
    os<<beta*i/(v.ntime()-1)<<" ";
    if (v.shape()==diagonal){
      for (unsigned int k=0; k<v.nsite(); ++k)
        for(unsigned int f=0;f<v.nflavor();++f)
          os<<v(i,k,k,f)<<" ";
    }else if(v.shape()==blockdiagonal){
      for (unsigned int k=0; k<v.nsite(); k+=2)
        for (unsigned int p1=0; p1<=1; p1++)
          for (unsigned int p2=0; p2<=1; p2++)
            for(unsigned int f=0;f<v.nflavor();++f)
              os<<v(i,k+p1,k+p2,f)<<" ";
    }else if(v.shape()==nondiagonal){
      for (unsigned int s1=0; s1<v.nsite(); ++s1)
        for (unsigned int s2=0; s2<v.nsite(); ++s2)
          for(unsigned int f=0;f<v.nflavor();++f)
            os<<v(i,s1,s2,f)<<" ";
    }
    os<<std::endl;
  }
}

void print_all_green_functions(std::string const &basename, const int iteration_ctr, const matsubara_green_function_t &G0_omega,
                               const matsubara_green_function_t &G_omega, const itime_green_function_t &G0_tau,
                               const itime_green_function_t &G_tau, const double beta,
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
  
  if(!G0_omega_file.is_open() || !G0_tau_file.is_open() || !G_omega_file.is_open() || !G_tau_file.is_open()) throw std::runtime_error("error opening gf files.");
  print_imag_green_matsubara(G0_omega_file, G0_omega, beta);
  print_real_green_matsubara(G0_omegareal_file,G0_omega,beta);
  print_green_itime(G0_tau_file,G0_tau, beta);
  print_imag_green_matsubara(G_omega_file,G_omega,beta);
  print_real_green_matsubara(G_omegareal_file,G_omega,beta);
  print_green_itime(G_tau_file,G_tau,beta);
  print_imag_green_matsubara(G_omega_file_2,G_omega,beta);
  print_green_itime(G_tau_file_2,G_tau,beta);
  print_selfenergy_matsubara(selfenergy_file, G0_omega, G_omega, beta);
  if (G_omega.shape()==diagonal)
    print_quasiparticle_estimate(std::cout, G_omega, G0_omega, beta);
  
}


void print_tau_green_functions(const int iteration_ctr, const itime_green_function_t &G0_tau, const itime_green_function_t &G_tau, const double beta){
  std::ostringstream G0tau_name, Gtau_name, Gtau_name_2;
  G0tau_name<<"G0_tau_"<<iteration_ctr;//<<"_"<<process_id;
  Gtau_name<<"G_tau_"<<iteration_ctr;//<<"_"<<process_id;
  Gtau_name_2<<"G_tau";//<<"_"<<process_id;
  std::ofstream G0_tau_file(G0tau_name.str().c_str());
  std::ofstream G_tau_file(Gtau_name.str().c_str());
  std::ofstream G_tau_file_2(Gtau_name_2.str().c_str());
  print_green_itime(G0_tau_file,G0_tau,beta);
  print_green_itime(G_tau_file,G_tau,beta);
  print_green_itime(G_tau_file_2,G_tau,beta);
}



void print_dressed_tau_green_functions(const int iteration_ctr, const itime_green_function_t &G_tau, const double beta,
                                       const std::string suffix){
  std::ostringstream Gtau_name, Gtau_name_2;
  Gtau_name<<"G_tau"<<suffix<<"_"<<iteration_ctr;//<<"_"<<process_id;
  Gtau_name_2<<"G_tau"<<suffix;//<<"_"<<process_id;
  std::ofstream G_tau_file(Gtau_name.str().c_str());
  std::ofstream G_tau_file_2(Gtau_name_2.str().c_str());
  print_green_itime(G_tau_file,G_tau,beta);
  print_green_itime(G_tau_file_2,G_tau,beta);
}

std::ostream &operator<<(std::ostream &os, const std::pair<green_function<double>,double> &v_and_beta){
  os<<std::setprecision(14);
  double dtau=v_and_beta.second/(v_and_beta.first.ntime()-1);
  if(!(v_and_beta.first.get_header().empty())){
	  os<<v_and_beta.first.get_header()<<std::endl;
  }

  for(unsigned int o=0;o<v_and_beta.first.ntime();++o){
    os<<o*dtau<<"\t";
    for(unsigned int i=0;i<v_and_beta.first.nsite();++i)
      for(unsigned int j=0;j<v_and_beta.first.nsite();++j)
        for(unsigned int z=0;z<v_and_beta.first.nflavor();++z)
          os<<v_and_beta.first(o,i,j,z)<<"\t";
    os<<std::endl;
  }
  return os;
}
std::ostream &operator<<(std::ostream &os, const std::pair<green_function<std::complex<double> >, double> &v_and_beta){
  os<<std::setprecision(14);
  if(!(v_and_beta.first.get_header().empty())){
	  os<<v_and_beta.first.get_header()<<std::endl;
  }
  for(unsigned int o=0;o<v_and_beta.first.ntime();++o){
    os<<(2.*o+1)*M_PI/v_and_beta.second<<"\t";
    if(v_and_beta.first.shape()==nondiagonal){
      for(unsigned int i=0;i<v_and_beta.first.nsite();++i){
        for(unsigned int j=0;j<v_and_beta.first.nsite();++j){
          for(unsigned int z1=0;z1<v_and_beta.first.nflavor();++z1){
            os<<v_and_beta.first(o,i,j,z1).real()<<" "<<v_and_beta.first(o,i,j,z1).imag()<<" ";
          }
        }
      }
    }else if(v_and_beta.first.shape()==blockdiagonal){
      for (unsigned int k=0; k<v_and_beta.first.nsite(); k+=2){
        for (unsigned int p1=0; p1<=1; p1++){
          for (unsigned int p2=0; p2<=1; p2++){
            for(unsigned int f=0;f<v_and_beta.first.nflavor();++f){
              os<<v_and_beta.first(o,k+p1  ,k+p2  ,f).real()<<" "<<v_and_beta.first(o,k+p1  ,k+p2  ,f).imag()<<" ";
            }
          }
        }
      }
    }else if(v_and_beta.first.shape()==diagonal){
      for(unsigned int i=0;i<v_and_beta.first.nsite();++i){
        for(unsigned int z1=0;z1<v_and_beta.first.nflavor();++z1){
          os<<v_and_beta.first(o,i,i,z1).real()<<" "<<v_and_beta.first(o,i,i,z1).imag()<<" ";
        }
      }
    }else{
      throw std::logic_error("shape not understood in write.");
    }
    os<<std::endl;
  }
  return os;
}
std::istream &operator>>(std::istream &is, green_function<std::complex<double> > &gf){
  double ignored;
  for(unsigned int o=0;o<gf.ntime();++o){
    is>>ignored>>std::ws;
    if(gf.shape()==nondiagonal){
      for(unsigned int i=0;i<gf.nsite();++i){
        for(unsigned int j=0;j<gf.nsite();++j){
          for(unsigned int z1=0;z1<gf.nflavor();++z1){
            is>>gf(o,i,j,z1).real()>>gf(o,i,j,z1).imag();
          }
        }
      }
    }else if(gf.shape()==blockdiagonal){
      for (unsigned int k=0; k<gf.nsite(); k+=2){
        for (unsigned int p1=0; p1<=1; p1++){
          for (unsigned int p2=0; p2<=1; p2++){
            for(unsigned int f=0;f<gf.nflavor();++f){
              is>>gf(o,k+p1  ,k+p2  ,f).real()>>gf(o,k+p1  ,k+p2  ,f).imag();
            }
          }
        }
      }
    }else if(gf.shape()==diagonal){
      for(unsigned int i=0;i<gf.nsite();++i){
        for(unsigned int z1=0;z1<gf.nflavor();++z1){
          is>>gf(o,i,i,z1).real()>>gf(o,i,i,z1).imag();
        }
      }
    }else{
      throw std::logic_error("shape not understood in read.");
    }
    is>>std::ws;
  }
  return is;
}

std::ostream &operator<<(std::ostream &os, const std::pair<nambu_green_function<double>, double> &v_and_beta){
  os<<std::setprecision(14);
  double dtau=v_and_beta.second/(v_and_beta.first.ntime()-1);
  for(unsigned int o=0;o<v_and_beta.first.ntime();++o){
    os<<o*dtau<<"\t";
    for(unsigned int i=0;i<v_and_beta.first.nsite();++i)
      for(unsigned int j=0;j<v_and_beta.first.nsite();++j)
        for(unsigned int z1=0;z1<v_and_beta.first.nflavor();++z1)
          for(unsigned int z2=0;z2<v_and_beta.first.nflavor();++z2)
            os<<v_and_beta.first(o,i,j,z1,z2)<<"\t";
    os<<std::endl;
  }
  return os;
}
std::ostream &operator<<(std::ostream &os, const std::pair<nambu_green_function<std::complex<double> >, double> &v_and_beta){
  os<<std::setprecision(14);
  for(unsigned int o=0;o<v_and_beta.first.ntime();++o){
    os<<(2.*o+1)*M_PI/v_and_beta.second<<"\t";
    if(v_and_beta.first.shape()==diagonal){
      for(unsigned int i=0;i<v_and_beta.first.nsite();++i)
        for(unsigned int z1=0;z1<v_and_beta.first.nflavor();++z1)
          for(unsigned int z2=0;z2<v_and_beta.first.nflavor();++z2)
            os<<v_and_beta.first(o,i,i,z1,z2).real()<<" "<<v_and_beta.first(o,i,i,z1,z2).imag()<<" ";
    }else{
      for(unsigned int i=0;i<v_and_beta.first.nsite();++i)
        for(unsigned int j=0;j<v_and_beta.first.nsite();++j)
          for(unsigned int z1=0;z1<v_and_beta.first.nflavor();++z1)
            for(unsigned int z2=0;z2<v_and_beta.first.nflavor();++z2)
              os<<v_and_beta.first(o,i,j,z1,z2).real()<<" "<<v_and_beta.first(o,i,j,z1,z2).imag()<<" ";
    }
    os<<std::endl;
  }
  return os;
}

std::istream &operator>>(std::istream &is, green_function<double> &v){
  double index;
  for(unsigned int o=0;o<v.nfreq();++o){
    is>>index;
    for(unsigned int i=0;i<v.nsite();++i)
      if(v.shape()==nondiagonal){
        for(unsigned int j=0;j<v.nsite();++j)
          for(unsigned int z=0;z<v.nflavor();++z)
            is>>v(o,i,j,z);
      }else{
        for(unsigned int z=0;z<v.nflavor();++z)
          is>>v(o,i,i,z);
      }
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
  os<<std::setprecision(14);
  for(unsigned int o=0;o<v.nfreq();++o){
    os<<o<<"\t";
    for(unsigned int i=0;i<v.nsite();++i)
      if(v.shape()==nondiagonal){
        for(unsigned int j=0;j<v.nsite();++j)
          for(unsigned int z=0;z<v.nflavor();++z)
            os<<(v(o,i,j,z).real())<<"\t"<<(v(o,i,j,z).imag())<<"\t";
      }else{
        for(unsigned int z=0;z<v.nflavor();++z)
          os<<(v(o,i,i,z).real())<<"\t"<<(v(o,i,i,z).imag())<<"\t";
        
      }
    os<<std::endl;
  }
  return os;
}


//////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////THE NAMBU PART/////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
void print_green_itime(std::ostream &os, const nambu_itime_green_function_t &v, const double beta, const shape_t shape){
  os<<std::setprecision(14);
  for(unsigned int i=0;i<v.ntime();++i){
    os<<beta*i/(v.ntime()-1)<<" ";
    switch (shape) {
      case diagonal:
        for (unsigned int k=0; k<v.nsite(); ++k)
          for(unsigned int f1=0;f1<v.nflavor();++f1)
            for(unsigned int f2=0;f2<v.nflavor();++f2)
              os<<v(i,k,k,f1,f2)<<" ";
        break;
      case blockdiagonal:
        for (unsigned int k=0; k<v.nsite(); k+=2)
          for (unsigned int p1=0; p1<=1; p1++)
            for (unsigned int p2=0; p2<=1; p2++)
              for(unsigned int f1=0;f1<v.nflavor();++f1)
                for(unsigned int f2=0;f2<v.nflavor();++f2)
                  os<<v(i,k+p1,k+p2,f1,f2)<<" ";
        break;
      case nondiagonal:
        for (unsigned int s1=0; s1<v.nsite(); ++s1)
          for (unsigned int s2=0; s2<v.nsite(); ++s2)
            for(unsigned int f1=0;f1<v.nflavor();++f1)
              for(unsigned int f2=0;f2<v.nflavor();++f2)
                os<<v(i,s1,s2,f1,f2)<<" ";
    }
    os<<std::endl;
  }
}

void print_imag_green_matsubara(std::ostream &os, const nambu_matsubara_green_function_t &v, const double beta, const shape_t shape)
{
  os<<std::setprecision(14);
  for(unsigned int i=0;i<v.nfreq();++i){
    os<<(2*i+1)*M_PI/beta<<" ";
    switch (shape) {
      case diagonal:
        for (unsigned int k=0; k<v.nsite(); ++k)
          for(unsigned int f1=0;f1<v.nflavor();++f1)
            for(unsigned int f2=0;f2<v.nflavor();++f2)
              os<<v(i,k,k,f1,f2).imag()<<" ";
        break;
      case blockdiagonal:
        for (unsigned int k=0; k<v.nsite(); k+=2)
          for (unsigned int p1=0; p1<=1; p1++)
            for (unsigned int p2=0; p2<=1; p2++)
              for(unsigned int f1=0;f1<v.nflavor();++f1)
                for(unsigned int f2=0;f2<v.nflavor();++f2)
                  os<<v(i,k+p1,k+p2,f1,f2).imag()<<" ";
        break;
      case nondiagonal:
        for (unsigned int s1=0; s1<v.nsite(); ++s1)
          for (unsigned int s2=0; s2<v.nsite(); ++s2)
            for(unsigned int f1=0;f1<v.nflavor();++f1)
              for(unsigned int f2=0;f2<v.nflavor();++f2)
                os<<v(i,s1,s2,f1,f2).imag()<<" ";
    }
    os<<std::endl;
  }
}
void print_real_green_matsubara(std::ostream &os, const nambu_matsubara_green_function_t &v, const double beta, const shape_t shape){
  os<<std::setprecision(14);
  for(unsigned int i=0;i<v.nfreq();++i){
    os<<(2*i+1)*M_PI/beta<<" ";
    switch (shape) {
      case diagonal:
        for (unsigned int k=0; k<v.nsite(); ++k)
          for(unsigned int f1=0;f1<v.nflavor();++f1)
            for(unsigned int f2=0;f2<v.nflavor();++f2)
              os<<v(i,k,k,f1,f2).real()<<" ";
        break;
      case blockdiagonal:
        for (unsigned int k=0; k<v.nsite(); k+=2)
          for (unsigned int p1=0; p1<=1; p1++)
            for (unsigned int p2=0; p2<=1; p2++)
              for(unsigned int f1=0;f1<v.nflavor();++f1)
                for(unsigned int f2=0;f2<v.nflavor();++f2)
                  os<<v(i,k+p1,k+p2,f1,f2).real()<<" ";
        break;
      case nondiagonal:
        for (unsigned int s1=0; s1<v.nsite(); ++s1)
          for (unsigned int s2=0; s2<v.nsite(); ++s2)
            for(unsigned int f1=0;f1<v.nflavor();++f1)
              for(unsigned int f2=0;f2<v.nflavor();++f2)
                os<<v(i,s1,s2,f1,f2).real()<<" ";
    }
    os<<std::endl;
  }
}
void print_selfenergy_matsubara(std::ostream &os, const nambu_matsubara_green_function_t &g_bare,
                                const nambu_matsubara_green_function_t g_dressed, const double beta, const shape_t shape)
{
  os<<std::setprecision(14);
  switch (shape) {
    case diagonal:
      for(unsigned int i=0;i<g_bare.nfreq();++i){
        os<<(2*i+1)*M_PI/beta<<"\t";
        for (unsigned int k=0; k<g_bare.nsite(); ++k){
          std::complex<double> g0_block[4]={g_bare(i,k,k,0,0),g_bare(i,k,k,0,1),g_bare(i,k,k,1,0),g_bare(i,k,k,1,1)};
          std::complex<double> det_g0=g0_block[0]*g0_block[3]-g0_block[1]*g0_block[2];
          std::complex<double>  g_block[4]={g_dressed(i,k,k,0,0),g_dressed(i,k,k,0,1),g_dressed(i,k,k,1,0),g_dressed(i,k,k,1,1)};
          std::complex<double> det_g=g_block[0]*g_block[3]-g_block[1]*g_block[2];
          std::complex<double> sigma[4];
          sigma[0]= g0_block[3]/det_g0-g_block[3]/det_g;
          sigma[1]=-g0_block[1]/det_g0+g_block[1]/det_g;
          sigma[2]=-g0_block[2]/det_g0+g_block[2]/det_g;
          sigma[3]= g0_block[0]/det_g0-g_block[0]/det_g;
          os<<sigma[0].real()<<" "<<sigma[0].imag()<<" "<<sigma[1].real()<<" "<<sigma[1].imag()<<" ";
          os<<sigma[2].real()<<" "<<sigma[2].imag()<<" "<<sigma[3].real()<<" "<<sigma[3].imag()<<" ";
        }
        os<<std::endl;
      }
      break;
    case blockdiagonal:
    {
      throw std::logic_error("please use cluster framework for block diag sigma");
    }
    case nondiagonal:
    {
      int nsite=g_dressed.nsite();
      for(unsigned int w=0;w<g_dressed.nfreq();++w){
        os<<(2*w+1)*M_PI/beta<<"\t";
        alps::numeric::matrix<std::complex<double> > G(nsite*2, nsite*2);
        alps::numeric::matrix<std::complex<double> > G_inv(nsite*2, nsite*2);
        alps::numeric::matrix<std::complex<double> > G0(nsite*2, nsite*2);
        alps::numeric::matrix<std::complex<double> > G0_inv(nsite*2, nsite*2);
        for(unsigned int i=0;i<(unsigned int) nsite;++i){
          for(unsigned int j=0;j<(unsigned int) nsite;++j){
            G(2*i  , 2*j  )=g_dressed(w,i,j,0,0);
            G(2*i  , 2*j+1)=g_dressed(w,i,j,0,1);
            G(2*i+1, 2*j  )=g_dressed(w,i,j,1,0);
            G(2*i+1, 2*j+1)=g_dressed(w,i,j,1,1);
            G0(2*i  , 2*j  )=g_bare(w,i,j,0,0);
            G0(2*i  , 2*j+1)=g_bare(w,i,j,0,1);
            G0(2*i+1, 2*j  )=g_bare(w,i,j,1,0);
            G0(2*i+1, 2*j+1)=g_bare(w,i,j,1,1);
          }
        }
        matrix_inversion(G, G_inv);
        matrix_inversion(G0, G0_inv);
        for(unsigned int i=0;i<(unsigned int)nsite;++i){
          for(unsigned int j=0;j<(unsigned int)nsite;++j){
            std::complex<double> sigma00=G0_inv(2*i  , 2*i  )-G_inv(2*i  ,2*i  );
            std::complex<double> sigma01=G0_inv(2*i  , 2*i+1)-G_inv(2*i  ,2*i+1);
            std::complex<double> sigma10=G0_inv(2*i+1, 2*i  )-G_inv(2*i+1,2*i  );
            std::complex<double> sigma11=G0_inv(2*i+1, 2*i+1)-G_inv(2*i+1,2*i+1);
            os<<sigma00.real()<<" "<<sigma00.imag()<<" "<<sigma01.real()<<" "<<sigma01.imag()<<" ";
            os<<sigma10.real()<<" "<<sigma10.imag()<<" "<<sigma11.real()<<" "<<sigma11.imag()<<" ";
          }
        }
        os<<std::endl;
      }
      //throw std::invalid_argument("please use cluster framework for full nondiag GF");
    }
  }
}
void print_all_nambu_green_functions(std::string const &basename, const int iteration_ctr, const nambu_matsubara_green_function_t &G0_omega,
                                     const nambu_matsubara_green_function_t &G_omega, const nambu_itime_green_function_t &G0_tau,
                                     const nambu_itime_green_function_t &G_tau, const double beta, const shape_t shape,
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
}

template<>  void green_function<double>::write_hdf5(alps::hdf5::archive &ar, const std::string &path) const{
    ar<<alps::make_pvp(path+"/nt",nt_);
    ar<<alps::make_pvp(path+"/ns",ns_);
    ar<<alps::make_pvp(path+"/nf",nf_);
    ar<<alps::make_pvp(path+"/shape",((int)shape_));
    std::stringstream subpath; subpath<<path<<"/values/mean";
    ar<<alps::make_pvp(subpath.str(), val_, nt_*(shape_==diagonal?ns_:ns_*ns_)*nf_);
  }
template<>  void green_function<double>::read_hdf5(alps::hdf5::archive &ar, const std::string &path) {
    unsigned int nt, ns, nf,shape;
    clear();
    ar>>alps::make_pvp(path+"/nt",nt);
    ar>>alps::make_pvp(path+"/ns",ns);
    ar>>alps::make_pvp(path+"/nf",nf);
    ar>>alps::make_pvp(path+"/shape",shape);
    if(nt!=nt_ || ns!=ns_ || nf!=nf_){ std::cerr<<path<<" nt: "<<nt_<<" new: "<<nt<<" ns: "<<ns_<<" "<<ns<<" nf: "<<nf_<<" "<<nf<<" dimensions do not match."<<std::endl; throw std::runtime_error("Green's function read in: dimensions do not match."); }
    if(shape!=shape_){ std::cerr<<path<<" shape: "<<shape_<<" new: "<<shape<<" shapes do not match."<<std::endl; throw std::runtime_error("Green's function read in: shapes do not match."); }
    std::stringstream subpath; subpath<<path<<"/values/mean";
    ar>>alps::make_pvp(subpath.str(), val_, nt_*(shape_==diagonal?ns_:ns_*ns_)*nf_);
  }

template<>  void green_function<std::complex<double> >::write_hdf5(alps::hdf5::archive &ar, const std::string &path) const{
    ar<<alps::make_pvp(path+"/nt",nt_);
    ar<<alps::make_pvp(path+"/ns",ns_);
    ar<<alps::make_pvp(path+"/nf",nf_);
    ar<<alps::make_pvp(path+"/shape",((int)shape_));
    std::stringstream subpath; subpath<<path<<"/values/mean";
    ar<<alps::make_pvp(subpath.str(), val_, nt_*(shape_==diagonal?ns_:ns_*ns_)*nf_);
  }
template<>  void green_function<std::complex<double> >::read_hdf5(alps::hdf5::archive &ar, const std::string &path) {
    unsigned int nt, ns, nf,shape;
    clear();
    ar>>alps::make_pvp(path+"/nt",nt);
    ar>>alps::make_pvp(path+"/ns",ns);
    ar>>alps::make_pvp(path+"/nf",nf);
    ar>>alps::make_pvp(path+"/shape",shape);
    if(nt!=nt_ || ns!=ns_ || nf!=nf_){ std::cerr<<path<<" nt: "<<nt_<<" new: "<<nt<<" ns: "<<ns_<<" "<<ns<<" nf: "<<nf_<<" "<<nf<<" dimensions do not match."<<std::endl; throw std::runtime_error("Green's function read in: dimensions do not match."); }
    if(shape!=shape_){ std::cerr<<path<<" shape: "<<shape_<<" new: "<<shape<<" shapes do not match."<<std::endl; throw std::runtime_error("Green's function read in: shapes do not match."); }
    std::stringstream subpath; subpath<<path<<"/values/mean";
    ar>>alps::make_pvp(subpath.str(), val_, nt_*(shape_==diagonal?ns_:ns_*ns_)*nf_);
  }

template<>  void nambu_green_function<std::complex<double> >::write_hdf5(alps::hdf5::archive &ar, const std::string &path) const{
    ar<<alps::make_pvp(path+"/nt",nt_);
    ar<<alps::make_pvp(path+"/ns",ns_);
    ar<<alps::make_pvp(path+"/nf",nf_);
    std::stringstream subpath; subpath<<path<<"/values/mean";
    //std::cout<<path<<" nt: "<<nt_<<" new: "<<nt_<<" ns: "<<ns_<<" "<<ns_<<" nf: "<<nf_<<" "<<nf_<<" total size should be: "<<nt_*ns_*ns_*nf_*nf_<<std::endl;
    ar<<alps::make_pvp(subpath.str(), val_, nt_*ns_*ns_*nf_*nf_);
  }
template<>  void nambu_green_function<std::complex<double> >::read_hdf5(alps::hdf5::archive &ar, const std::string &path) {
    unsigned int nt, ns, nf;
    clear();
    ar>>alps::make_pvp(path+"/nt",nt);
    ar>>alps::make_pvp(path+"/ns",ns);
    ar>>alps::make_pvp(path+"/nf",nf);
    if(nt!=nt_ || ns!=ns_ || nf!=nf_){ std::cerr<<path<<" nt: "<<nt_<<" new: "<<nt<<" ns: "<<ns_<<" "<<ns<<" nf: "<<nf_<<" "<<nf<<" dimensions do not match."<<std::endl; throw std::runtime_error("Green's function read in: dimensions do not match."); }
    //std::cout<<path<<" nt: "<<nt_<<" new: "<<nt<<" ns: "<<ns_<<" "<<ns<<" nf: "<<nf_<<" "<<nf<<" total size should be: "<<nt_*ns_*ns_*nf_*nf_<<std::endl;
    std::stringstream subpath; subpath<<path<<"/values/mean";
    ar>>alps::make_pvp(subpath.str(), val_, nt_*ns_*ns_*nf_*nf_);
  }
template<>  void nambu_green_function<double>::write_hdf5(alps::hdf5::archive &ar, const std::string &path) const{
    ar<<alps::make_pvp(path+"/nt",nt_);
    ar<<alps::make_pvp(path+"/ns",ns_);
    ar<<alps::make_pvp(path+"/nf",nf_);
    std::stringstream subpath; subpath<<path<<"/values/mean";
    //std::cout<<path<<" nt: "<<nt_<<" new: "<<nt_<<" ns: "<<ns_<<" "<<ns_<<" nf: "<<nf_<<" "<<nf_<<" total size should be: "<<nt_*ns_*ns_*nf_*nf_<<std::endl;
    ar<<alps::make_pvp(subpath.str(), val_, nt_*ns_*ns_*nf_*nf_);
  }
  template<> void nambu_green_function<double>::read_hdf5(alps::hdf5::archive &ar, const std::string &path) {
    unsigned int nt, ns, nf;
    clear();
    ar>>alps::make_pvp(path+"/nt",nt);
    ar>>alps::make_pvp(path+"/ns",ns);
    ar>>alps::make_pvp(path+"/nf",nf);
    if(nt!=nt_ || ns!=ns_ || nf!=nf_){ std::cerr<<path<<" nt: "<<nt_<<" new: "<<nt<<" ns: "<<ns_<<" "<<ns<<" nf: "<<nf_<<" "<<nf<<" dimensions do not match."<<std::endl; throw std::runtime_error("Green's function read in: dimensions do not match."); }
    //std::cout<<path<<" nt: "<<nt_<<" new: "<<nt<<" ns: "<<ns_<<" "<<ns<<" nf: "<<nf_<<" "<<nf<<" total size should be: "<<nt_*ns_*ns_*nf_*nf_<<std::endl;
    std::stringstream subpath; subpath<<path<<"/values/mean";
    ar>>alps::make_pvp(subpath.str(), val_, nt_*ns_*ns_*nf_*nf_);
  }
