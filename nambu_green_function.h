#ifndef nambu_green_function_H
#define nambu_green_function_H
#include "types.h"
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <iostream>
#ifdef USE_MPI
#include <mpi.h>
#endif

namespace alps{namespace hdf5{class archive; } }

//Matsubara GF: use T=std::complex<double>
//Imaginary time: use T=double
template <typename T> class nambu_green_function{
  
public:
	//construction and destruction, assignement and copy constructor
	///constructor: how many time slices, how many sites, how many flavors
	nambu_green_function(unsigned int ntime, unsigned int nsite, unsigned int nflavor, shape_t shape):nt_(ntime), ns_(nsite), nf_(nflavor),
   ntns_(ntime*nsite), ntnsns_(ntime*nsite*nsite), ntnsnsnf_(ntime*nsite*nsite*nflavor), shape_(shape){
		val_=new T[nt_*ns_*ns_*nf_*nf_];
	}
	///destructor
	~nambu_green_function(){
		delete [] val_;
	}
	///copy constructor
	nambu_green_function(const nambu_green_function &g):nt_(g.nt_), ns_(g.ns_), nf_(g.nf_), ntns_(g.ntns_), ntnsns_(g.ntnsns_), ntnsnsnf_(g.ntnsnsnf_), shape_(g.shape_){
		val_=new T[nt_*ns_*ns_*nf_*nf_];
		operator=(g);
	}
	///operator= (assignement operator)
	const nambu_green_function &operator=(const nambu_green_function &g){
		memcpy(val_, g(), sizeof(T)*nt_*ns_*ns_*nf_*nf_);
    return *this;
  }
  void clear(){ memset(val_, 0, ns_*ns_*nt_*nf_*nf_*sizeof(T)); }
	//access of vectors and elements
	///return an entire vector of times for a given flavor
	///access element with given time, site 1, site 2, and flavor
  inline T &operator()(unsigned int t, unsigned int site1, unsigned int site2, unsigned int flavor1, unsigned int flavor2){
    return val_[t+nt_*site1+ntns_*site2+ntnsns_*flavor1+ntnsnsnf_*flavor2];
  }
	///access element with given time, site 1, site 2, and flavor (const reference)
  inline const T &operator()(unsigned int t, unsigned int site1, unsigned int site2, unsigned int flavor1, unsigned int flavor2)const{
    return val_[t+nt_*site1+ntns_*site2+ntnsns_*flavor1+ntnsnsnf_*flavor2];
  }
	///return an entire vector of imaginary time values for a given site 1,
  //site2, flavor1, flavor2
  inline T *operator()(unsigned int site1, unsigned int site2, unsigned int flavor1, unsigned int flavor2){
    return val_+nt_*site1+ntns_*site2+ntnsns_*flavor1+ntnsnsnf_*flavor2;
  }
	///get all values at once
	inline const T *operator()() const {return val_;}
	///get all errors at once
	inline const T *error() const {return err_;}
  
	//size information
	///how many flavors do we have? (flavors are usually spins, GF of different flavors are zero)
	inline const unsigned int &nflavor()const{return nf_;}
	///return # of sites
	inline const unsigned int &nsite()const{return ns_;}
	///return # of imaginary time values
	inline const unsigned int &ntime()const{return nt_;}
	///return # of matsubara frequencies. Exactly equivalent to ntime().
	///In the case of a Matsubara GF 'ntime' sounds odd -> define 'nfreq' instead.
	inline const unsigned int &nfreq()const{return nt_;} //nfreq is an alias to ntime - more intuitive use for Matsubara GF
  
  void read(const char *filename);
  void write(const char *filename);
  void write_hdf5(alps::hdf5::archive &ar, const std::string &path) const;
  void read_hdf5(alps::hdf5::archive &ar, const std::string &path);
#ifdef USE_MPI
    void broadcast(){
      MPI_Bcast( val_, ntnsnsnf_*nf_*sizeof(T)/sizeof(double), MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast( err_, ntnsnsnf_*nf_*sizeof(T)/sizeof(double), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
#endif

  const shape_t &shape() const{return shape_;}
private:
	//const values
	const unsigned int nt_; ///imag time points
	const unsigned int ns_; ///number of sites
	const unsigned int nf_; ///number of flavors
	const unsigned int ntns_; ///nt*ns
	const unsigned int ntnsns_; ///nt*ns*ns
	const unsigned int ntnsnsnf_; ///nt*ns*ns*nf
  const shape_t shape_;
	// the actual values and errors.
  T *val_;
  T *err_;
  
};
typedef nambu_green_function<std::complex<double> > nambu_matsubara_green_function_t;
typedef nambu_green_function<double> nambu_itime_green_function_t;
///write out imag time Green function
std::ostream &operator<<(std::ostream &os, const std::pair<nambu_green_function<double>, double> &v_and_beta);
std::ostream &operator<<(std::ostream &os, const std::pair<nambu_green_function<std::complex<double> >, double> &v_and_beta);
///read in imag time Green function
std::istream &operator>>(std::istream &is, const nambu_green_function<double> &v);
///write out Matsubara Green function
std::ostream &operator<<(std::ostream &os, const nambu_green_function<std::complex<double> > &v);
///read in Matsubara Green function
std::istream &operator>>(std::istream &is, const nambu_green_function<std::complex<double> > &v);

template<typename T> void nambu_green_function<T>::read(const char *filename){
  std::ifstream in_file(filename);
  if(!in_file.is_open()){ std::cerr<<"input file "<<filename<<" could not be opened!"<<std::endl; exit(1);}
  double ignored;
  for(unsigned int i=0;i<nt_;++i){
    in_file>>ignored; //read first entry, which could be # matsubara frequencies, or tau-point, or N/beta*tau, or...
    if(shape_==diagonal){
      for(unsigned int s0=0; s0<ns_; ++s0) {
        for(unsigned int f1=0; f1<nf_; ++f1){
          for(unsigned int f2=0; f2<nf_; ++f2) {
            in_file>>operator()(i, s0, s0, f1,f2)>>std::ws; 
          }
        }
      }
      in_file>>std::ws;
    }else{
      for(unsigned int s0=0; s0<ns_; ++s0) {
        for(unsigned int s1=0; s1<ns_; ++s1){
          for(unsigned int f1=0; f1<nf_; ++f1) {
            for(unsigned int f2=0; f2<nf_; ++f2) {
              in_file>>operator()(i, s0, s1, f1, f2); //read the actual value
            }
          }
        }
      }
    }
  }
}

template<typename T> void nambu_green_function<T>::write(const char *filename){
  std::ofstream out_file(filename);
  out_file<<std::setprecision(14);
  if(!out_file.is_open()){ std::cerr<<"output file "<<filename<<" could not be opened!"<<std::endl; exit(1);}
  for(unsigned int i=0;i<nt_;++i){
    out_file << i << " ";
    if(shape_==diagonal){
      for(unsigned int s0=0; s0<ns_; ++s0) {
        for(unsigned int f1=0; f1<nf_; ++f1){
          for(unsigned int f2=0; f2<nf_; ++f2) {
            out_file<<operator()(i, s0, s0, f1,f2)<<" ";
          }
        }
      }
    }else{
      for(unsigned int s0=0; s0<ns_; ++s0) {
        for(unsigned int s1=0; s1<ns_; ++s1){
          for(unsigned int f1=0; f1<nf_; ++f1) {
            for(unsigned int f2=0; f2<nf_; ++f2) {
              out_file<<operator()(i, s0, s1, f1,f2) << " "; 
            }
          }
        }
      }
    }
    out_file << std::endl;
  }
}

void print_all_nambu_green_functions(std::string const &basename, const int iteration_ctr, const nambu_matsubara_green_function_t &G0_omega,
                                     const nambu_matsubara_green_function_t &G_omega, const nambu_itime_green_function_t &G0_tau, 
                                     const nambu_itime_green_function_t &G_tau, const double beta, const shape_t shape=diagonal,
                                     const std::string suffix="");
void print_real_green_matsubara(std::ostream &os, const nambu_matsubara_green_function_t &v, const double beta, const shape_t shape=diagonal);
void print_imag_green_matsubara(std::ostream &os, const nambu_matsubara_green_function_t &v, const double beta, const shape_t shape=diagonal);
#endif
