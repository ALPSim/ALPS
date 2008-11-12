#ifndef GREEN_FUNCTION_H
#define GREEN_FUNCTION_H
#include "types.h"
#include <fstream>
//Matsubara GF: use T=std::complex<double>
//Imaginary time: use T=double
template <typename T> class green_function{

	public:
	//construction and destruction, assignement and copy constructor
	///constructor: how many time slices, how many sites, how many flavors
	green_function(unsigned int ntime, unsigned int nsite, unsigned int nflavor):nt_(ntime), ns_(nsite), nf_(nflavor),
		 ntnsns_(ntime*nsite*nsite), ntns_(ntime*nsite){
		val_=new T[nt_*ns_*ns_*nf_];
		err_=new T[nt_*ns_*ns_*nf_];
	}
	///specialization: constructor for problems with only one site
	green_function(unsigned int ntime, unsigned int nflavor):nt_(ntime), ns_(1), nf_(nflavor),
		ntnsns_(ntime), ntns_(ntime){
		val_=new T[nt_*nf_];
		err_=new T[nt_*nf_];
	}
	///destructor
	~green_function(){
		delete [] val_;
		delete [] err_;
	}
	///copy constructor
	green_function(const green_function &g):nt_(g.nt_), ns_(g.ns_), nf_(g.nf_), ntnsns_(g.ntnsns_), ntns_(g.ntns_){
		val_=new T[nt_*ns_*ns_*nf_];
		err_=new T[nt_*ns_*ns_*nf_];
		operator=(g);
	}
	///operator= (assignement operator)
	const green_function &operator=(const green_function &g){
		memcpy(val_, g(), sizeof(T)*nt_*ns_*ns_*nf_);
		memcpy(err_, g.error(), sizeof(T)*nt_*ns_*ns_*nf_);
                return *this;
        }
  void clear(){ memset(val_, 0, ns_*ns_*nt_*nf_*sizeof(T)); }
	//access of vectors and elements
	///specialization for only one site: access element with given time and flavor
	inline T &operator()(unsigned int t, unsigned int flavor){return val_[t+nt_*flavor];}
	///specialization for only one site: return const reference to element with given time and flavor
	inline const T &operator()(unsigned int t, unsigned int flavor)const{return val_[t+nt_*flavor];}

	///return an entire vector of times for a given flavor
	inline T *operator()(unsigned int flavor){return val_+ntnsns_*flavor;}
	
	//error access
	inline T &error(unsigned int t, unsigned int flavor){return err_[t+nt_*flavor];}
	inline const T &error(unsigned int t, unsigned int flavor)const{return err_[t+nt_*flavor];}
	inline T *errors(unsigned int flavor){return err_+nt_*flavor;}
	///access element with given time, site 1, site 2, and flavor
        inline T &operator()(unsigned int t, unsigned int site1, unsigned int site2, unsigned int flavor){return val_[t+nt_*site1+ntns_*site2+ntnsns_*flavor];}
	///access element with given time, site 1, site 2, and flavor (const reference)
        inline const T &operator()(unsigned int t, unsigned int site1, unsigned int site2, unsigned int flavor)const{return val_[t+nt_*site1+ntns_*site2+ntnsns_*flavor];}
	///return an entire vector of imaginary time values for a given site 1, site2, flavor
        inline T *operator()(unsigned int site1, unsigned int site2, unsigned int flavor){return val_+nt_*site1+ntns_*site2+ntnsns_*flavor;}

        inline T &error(unsigned int t, unsigned int site1, unsigned int site2, unsigned int flavor){return err_[t+nt_*site1+ntns_*site2+ntnsns_*flavor];}
        inline const T &error(unsigned int t, unsigned int site1, unsigned int site2, unsigned int flavor)const{return err_[t+nt_*site1+ntns_*site2+ntnsns_*flavor];}
        inline T *errors(unsigned int site1, unsigned int site2, unsigned int flavor){return err_+nt_*site1+ntns_*site2+ntnsns_*flavor;}

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
        std::pair<std::vector<T>,std::vector<T> > to_multiple_vector() const;
        void from_multiple_vector(const std::pair<std::vector<T>,std::vector<T> > &mv);
        ///mirror a GF at zero.
        /*void mirror_real(){
            for(unsigned int i=0;i<ntnsns_*nf_;++i){
              val_[i].real()*=-1;
            }
        }
        void mirror_imag(){
            for(unsigned int i=0;i<ntnsns_*nf_;++i){
              val_[i].imag()*=-1;
            }
        }*/
              

        ///hobby implementation of an iterator with begin, end, and operator++
        //that goes through all entries of the GF
        /*class iterator{
          public:
            iterator(){}
            iterator(unsigned int i){pos_=i;}
            static iterator begin() {return iterator(0);}
            static iterator end()   {return iterator(nf_*ns_*nt_);}
            void operator++(){pos_++;}
            T &operator&()            {return val_[pos_];
            const T &operator&()const {return val_[pos_];
            iterator &operator=(const iterator &it){pos_=it.pos_;}
          private:
            unsigned int pos_;
        };*/
	private:
	//const values
	const unsigned int nt_; ///imag time points
	const unsigned int ns_; ///number of sites
	const unsigned int nf_; ///number of flavors
	const unsigned int ntnsns_; ///nt*ns*ns
	const unsigned int ntns_; ///nt*ns
	// the actual values and errors.
        T *val_;
        T *err_;

};
typedef green_function<std::complex<double> > matsubara_green_function_t;
typedef green_function<double> itime_green_function_t;
///write out imag time Green function
std::ostream &operator<<(std::ostream &os, const green_function<double> &v);
///read in imag time Green function
std::istream &operator>>(std::istream &is, green_function<double> &v);
///write out Matsubara Green function
std::ostream &operator<<(std::ostream &os, const green_function<std::complex<double> > &v);
///read in Matsubara Green function
std::istream &operator>>(std::istream &is, green_function<std::complex<double> > &v);

///compute kinetic energy
double kinetic_energy(const multiple_vector_type &G_tau, const double &beta, const double &t);

template<typename T> void green_function<T>::read(const char *filename){
  std::ifstream in_file(filename);
  assert(in_file.is_open()); //this is not good enough when people use -DNDEBUG.
  if(!in_file.is_open()){ std::cerr<<"input file "<<filename<<" could not be opened!"<<std::endl; exit(1);}
  double ignored=0;
  for(unsigned int i=0;i<nt_;++i){
    in_file>>ignored; //read first entry, which could be # matsubara frequencies, or tau-point, or N/beta*tau, or...
    for(unsigned int s0=0; s0<ns_; ++s0) {
      for(unsigned int s1=0; s1<ns_; ++s1){
	for(unsigned int f=0; f<nf_; ++f) {
	  in_file>>operator()(i, s0, s1, f)>>std::ws; //read the actual value
	}
      }
    }
  }
}

template<typename T> void green_function<T>::write(const char *filename){
  std::ofstream out_file(filename);
  assert(out_file.is_open()); //this is not good enough when people use -DNDEBUG.
  if(!out_file.is_open()){ std::cerr<<"output file "<<filename<<" could not be opened!"<<std::endl; exit(1);}
  for(unsigned int i=0;i<nt_;++i){
    out_file << i << " ";
    for(unsigned int s0=0; s0<ns_; ++s0) {
      for(unsigned int s1=0; s1<ns_; ++s1){
	for(unsigned int f=0; f<nf_; ++f) {
	  out_file<<operator()(i, s0, s1, f) << " "; 
	}
      }
    }
    out_file << std::endl;
  }
}

///for the transition period from multiple vectors to this data structure only.
template<typename T> std::pair<std::vector<T>,std::vector<T> > green_function<T>::to_multiple_vector() const{
  assert(ns_==1 && nf_<=2);
  std::pair<std::vector<T>,std::vector<T> > mv;
  mv.first.resize(nt_);
  mv.second.resize(nt_);
  for(unsigned int i=0;i<nt_;++i){
    mv.first[i]=operator()(i,0,0,0);
    mv.second[i]=nf_==1?operator()(i,0,0,0):operator()(i,0,0,1);
  }
  return mv;
}
template<typename T> void green_function<T>::from_multiple_vector(const std::pair<std::vector<T>,std::vector<T> >&mv){
  assert(ns_==1 && nf_<=2);
  for(unsigned int i=0;i<nt_;++i){
    operator()(i,0,0,0)=mv.first[i];
    if(nf_==2)
      operator()(i,0,0,1)=mv.second[i];
  }
}

void print_all_green_functions(const int iteration_ctr, const matsubara_green_function_t &G0_omega,
			       const matsubara_green_function_t &G_omega, const itime_green_function_t &G0_tau, 
			       const itime_green_function_t &G_tau, const double beta, const bool diagonal=true,
			       const std::string suffix="");
void print_real_green_matsubara(std::ostream &os, const matsubara_green_function_t &v, const double beta, const bool diagonal=true);
void print_imag_green_matsubara(std::ostream &os, const matsubara_green_function_t &v, const double beta, const bool diagonal=true);
void print_dressed_tau_green_functions(const int iteration_ctr, const itime_green_function_t &G_tau, const double beta, 
				       const bool diagonal=true, const std::string suffix="");
#endif
