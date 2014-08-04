/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Matthias Troyer <troyer@comp-phys.org>
 *                              Sebastian Fuchs <fuchs@theorie.physik.uni-goettingen.de>
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

#ifndef GREEN_FUNCTION_H
#define GREEN_FUNCTION_H
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <cstdlib>
#include <string.h>
#ifdef USE_MPI
#include <mpi.h>
#endif

namespace alps{namespace hdf5{class archive; } }

#include "types.h"


//Matsubara GF: use T=std::complex<double>
//Imaginary time: use T=double
template <typename T> class green_function{
  
public:
  //construction and destruction, assignement and copy constructor
  ///constructor: how many time slices, how many sites, how many flavors
  green_function(unsigned int ntime, unsigned int nsite, unsigned int nflavor, shape_t shape):nt_(ntime), ns_(nsite), nf_(nflavor),
  ntnsns_(ntime*nsite*nsite), ntns_(ntime*nsite), shape_(shape), header_("") {
    val_=(shape_==diagonal?new T[nt_*ns_*nf_]:new T[nt_*ns_*ns_*nf_]);
    clear();
  }
  ///destructor
  ~green_function(){
    delete [] val_;
  }
  ///copy constructor
  green_function(const green_function &g):nt_(g.nt_), ns_(g.ns_), nf_(g.nf_), ntnsns_(g.ntnsns_), ntns_(g.ntns_), shape_(g.shape_), header_(g.header_){
    val_=(shape_==diagonal?new T[nt_*ns_*nf_]:new T[nt_*ns_*ns_*nf_]);
    operator=(g);
  }
  ///operator= (assignement operator)
  const green_function &operator=(const green_function &g){
    if(ns_==!g.ns_)
      throw std::logic_error("in operator= : number of sites is different.");
    if(nf_!=g.nf_)
      throw std::logic_error("in operator= : number of flavors is different.");
    if(nt_!=g.nt_)
      throw std::logic_error("in operator= : number of frequencies or times is different.");
    if(shape_!=g.shape_ && ns_>1)
      throw std::logic_error("in operator= : shape is different.");
    memcpy(val_, g(), sizeof(T)*nt_*(shape_==diagonal?ns_:ns_*ns_)*nf_);
    return *this;
  }
  void clear(){ memset(val_, 0, (shape_==diagonal?ns_:ns_*ns_)*nt_*nf_*sizeof(T)); }
  //access of vectors and elements
  ///specialization for only one site: access element with given time and flavor
  inline T &operator()(unsigned int t, unsigned int flavor){return val_[t+nt_*flavor];}
  ///specialization for only one site: return const reference to element with given time and flavor
  inline const T &operator()(unsigned int t, unsigned int flavor)const{return val_[t+nt_*flavor];}
  
  ///return an entire vector of times for a given flavor
  inline T *operator()(unsigned int flavor){return val_+(shape_==diagonal?ntns_:ntnsns_)*flavor;}
  
  ///access element with given time, site 1, site 2, and flavor
  inline T &operator()(unsigned int t, unsigned int site1, unsigned int site2, unsigned int flavor){return val_[t+nt_*site1+(shape_==diagonal?ntns_*flavor:ntns_*site2+ntnsns_*flavor)];}
  ///access element with given time, site 1, site 2, and flavor (const reference)
  inline const T &operator()(unsigned int t, unsigned int site1, unsigned int site2, unsigned int flavor)const{return val_[t+nt_*site1+(shape_==diagonal?ntns_*flavor:ntns_*site2+ntnsns_*flavor)];}
  ///return an entire vector of imaginary time values for a given site 1, site2, flavor
  inline T *operator()(unsigned int site1, unsigned int site2, unsigned int flavor){return val_+nt_*site1+(shape_==diagonal?ntns_*flavor:ntns_*site2+ntnsns_*flavor);}
  
  ///get all values at once
  inline const T *operator()() const {return val_;}
  
  //size information
  ///how many flavors do we have? (flavors are usually spins, GF of different flavors are zero)
  inline const unsigned int &nflavor()const{return nf_;}
  ///return # of sites
  inline const unsigned int &nsite()const{return ns_;}
  ///return # of imaginary time values
  inline const unsigned int &ntime()const{return nt_;}
  ///return shape
  inline const shape_t &shape()const{return shape_;}
  ///return # of matsubara frequencies. Exactly equivalent to ntime().
  ///In the case of a Matsubara GF 'ntime' sounds odd -> define 'nfreq' instead.
  inline const unsigned int &nfreq()const{return nt_;} //nfreq is an alias to ntime - more intuitive use for Matsubara GF
  void read(const char *filename);
  void write(const char *filename) const;
  void write_hdf5(alps::hdf5::archive &ar, const std::string &path) const;
  void read_hdf5(alps::hdf5::archive &ar, const std::string &path) ;
  double diff(const green_function<T> &gf2){
    double maxdiff=0.;
    for(int i=0; i<ntnsns_*nf_;++i){
      double d=std::abs(val_[i]-gf2.val_[i]);
      maxdiff=d>maxdiff?d:maxdiff;
    }
    return maxdiff;
  }
  const std::string get_header()const{return header_;}
  void set_header(std::string header){header_=header.insert(0,"#");}

#ifdef USE_MPI
  void broadcast(){
    if(shape_==nondiagonal || shape_==blockdiagonal)
      MPI_Bcast( val_, ntnsns_*nf_*sizeof(T)/sizeof(double), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    else
      MPI_Bcast( val_, ntns_*nf_*sizeof(T)/sizeof(double), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
#endif
  
private:
  //const values
  const unsigned int nt_; ///imag time points
  const unsigned int ns_; ///number of sites
  const unsigned int nf_; ///number of flavors
  const unsigned int ntnsns_; ///nt*ns*ns
  const unsigned int ntns_; ///nt*ns
  // the actual values and errors.
  T *val_;
  const shape_t shape_;
  //Header string for Greens function files, labels the columns of the file with momentum or coordinates
  std::string header_;
};
typedef green_function<std::complex<double> > matsubara_green_function_t;
typedef green_function<double> itime_green_function_t;

std::ostream &operator<<(std::ostream &os, const std::pair<green_function<double> ,double>&v_and_beta);
std::ostream &operator<<(std::ostream &os, const std::pair<green_function<std::complex<double> > ,double>&v_and_beta);

std::istream &operator>>(std::istream &is, green_function<double> &v);
std::istream &operator>>(std::istream &is, green_function<std::complex<double> > &v);

///compute kinetic energy

template<typename T> void green_function<T>::read(const char *filename){
  std::ifstream in_file(filename);
  if(!in_file.is_open()){ throw(std::invalid_argument("input file could not be opened!")); }
  double ignored=0;
  for(unsigned int i=0;i<nt_;++i){
    in_file>>ignored; //read first entry, which could be # matsubara frequencies, or tau-point, or N/beta*tau, or...
    if(shape()==nondiagonal){
      for(unsigned int s0=0; s0<ns_; ++s0) {
        for(unsigned int s1=0; s1<ns_; ++s1){
          for(unsigned int f=0; f<nf_; ++f) {
            in_file>>operator()(i, s0, s1, f)>>std::ws; //read the actual value
          }
        }
      }
    }else if (shape()==diagonal){
      for(unsigned int s0=0; s0<ns_; ++s0) {
        for(unsigned int f=0; f<nf_; ++f) {
          in_file>>operator()(i, s0, s0, f)>>std::ws; //read the actual value
        }
      }
    }else if (shape()==blockdiagonal){
      for(unsigned int s0=0; s0<ns_; s0+=2) {
        for(unsigned int f=0; f<nf_; ++f) {
          in_file>>operator()(i, s0, s0, f) >>operator()(i, s0, s0+1, f)>>operator()(i, s0+1, s0, f)>>operator()(i, s0, s0+1, f)>>std::ws;
        }
      }
    }
  }
}

template<typename T> void green_function<T>::write(const char *filename) const{
  std::ofstream out_file(filename);
  out_file<<std::setprecision(14);
  if(!out_file.is_open()){ std::cerr<<"output file "<<filename<<" could not be opened!"<<std::endl; exit(1);}
  for(unsigned int i=0;i<nt_;++i){
    out_file << i << " ";
    if(shape_==nondiagonal){
      for(unsigned int s0=0; s0<ns_; ++s0) {
        for(unsigned int s1=0; s1<ns_; ++s1){
          for(unsigned int f=0; f<nf_; ++f) {
            out_file<<operator()(i, s0, s1, f) << " ";
          }
        }
      }
    }else if(shape_==blockdiagonal){
      for(unsigned int s0=0; s0<ns_; s0+=2) {
        for(unsigned int f=0; f<nf_; ++f) {
          out_file<<operator()(i, s0, s0, f) << " "<<operator()(i, s0, s0+1, f)<<" "<<operator()(i, s0+1, s0, f)<<" "<<operator()(i, s0, s0+1, f)<<" ";
        }
      }
    }
    else if(shape_==diagonal){
      for(unsigned int s0=0; s0<ns_; ++s0) {
        for(unsigned int f=0; f<nf_; ++f) {
          out_file<<operator()(i, s0, s0, f) << " ";
        }
      }
    }else throw std::logic_error("don't know shape in write.");
    out_file << std::endl;
  }
}



void print_all_green_functions(std::string const &basename, const int iteration_ctr, const matsubara_green_function_t &G0_omega,
                               const matsubara_green_function_t &G_omega, const itime_green_function_t &G0_tau,
                               const itime_green_function_t &G_tau, const double beta,
                               const std::string suffix="");
void print_real_green_matsubara(std::ostream &os, const matsubara_green_function_t &v, const double beta);
void print_imag_green_matsubara(std::ostream &os, const matsubara_green_function_t &v, const double beta);
void print_dressed_tau_green_functions(const int iteration_ctr, const itime_green_function_t &G_tau, const double beta,
                                       const std::string suffix="");
void print_tau_green_functions(const int iteration_ctr, const itime_green_function_t &G0_tau, const itime_green_function_t &G_tau, const double beta);

#endif
