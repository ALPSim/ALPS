#include "types.h"
//proposal: use this class for the storage of Green's functions.
//Matsubara GF: use T=std::complex<double>
//Imaginary time: use T=double
template <typename T> class green_function{

	public:
	//construction and destruction, assignement and copy constructor
	///constructor: how many time slices, how many sites, how many flavors
	green_function(uint ntime, uint nsite, uint nflavor):nt_(ntime), ns_(nsite), nf_(nflavor),
		 ntnsns_(ntime*nsite*nsite), ntns_(ntime*nsite){
		val_=new T[nt_*ns_*ns_*nf_];
		err_=new T[nt_*ns_*ns_*nf_];
	}
	///specialization: constructor for problems with only one site
	green_function(uint ntime, uint nflavor):nt_(ntime), ns_(1), nf_(nflavor),
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
		green_function(g.ntime(), g.nsite(), g.nflavor());
		operator=(g);
	}
	///operator= (assignement operator)
	const green_function &operator=(const green_function &g){
		memcpy(val_, g(), sizeof(T)*nt_*ns_*ns_*nf_);
		memcpy(err_, g.error(), sizeof(T)*nt_*ns_*ns_*nf_);
                return *this;
        }
	//access of vectors and elements
	///specialization for only one site: access element with given time and flavor
	inline T &operator()(uint t, uint flavor){return val_[t+nt_*flavor];}
	///specialization for only one site: return const reference to element with given time and flavor
	inline const T &operator()(uint t, uint flavor)const{return val_[t+nt_*flavor];}
	///specialization for only one site: return an entire vector of times for a given flavor
	inline T *operator()(uint flavor){return val_+nt_*flavor;}
	
	//error access
	inline T &error(uint t, uint flavor){return err_[t+nt_*flavor];}
	inline const T &error(uint t, uint flavor)const{return err_[t+nt_*flavor];}
	inline T *errors(uint flavor){return err_+nt_*flavor;}
	///access element with given time, site 1, site 2, and flavor
        inline T &operator()(uint t, uint site1, uint site2, uint flavor){return val_[t+nt_*site1+ntns_*site2+ntnsns_*flavor];}
	///access element with given time, site 1, site 2, and flavor (const reference)
        inline const T &operator()(uint t, uint site1, uint site2, uint flavor)const{return val_[t+nt_*site1+ntns_*site2+ntnsns_*flavor];}
	///return an entire vector of imaginary time values for a given site 1, site2, flavor
        inline T *operator()(uint site1, uint site2, uint flavor){return val_+nt_*site1+ntns_*site2+ntnsns_*flavor;}

        inline T &error(uint t, uint site1, uint site2, uint flavor){return err_[t+nt_*site1+ntns_*site2+ntnsns_*flavor];}
        inline const T &error(uint t, uint site1, uint site2, uint flavor)const{return err_[t+nt_*site1+ntns_*site2+ntnsns_*flavor];}
        inline T *errors(uint site1, uint site2, uint flavor){return err_+nt_*site1+ntns_*site2+ntnsns_*flavor;}

	///get all values at once
	inline const T *operator()() const {return val_;}
	///get all errors at once
	inline const T *error() const {return err_;}

	//size information
	///how many flavors do we have? (flavors are usually spins, GF of different flavors are zero)
	inline const uint &nflavor()const{return nf_;}
	///return # of sites
	inline const uint &nsite()const{return ns_;}
	///return # of imaginary time values
	inline const uint &ntime()const{return nt_;}
	///return # of matsubara frequencies. Exactly equivalent to ntime().
	///In the case of a Matsubara GF 'ntime' sounds odd -> define 'nfreq' instead.
	inline const uint &nfreq()const{return nt_;} //nfreq is an alias to ntime - more intuitive use for Matsubara GF
	private:
	//const values
	const uint nt_;
	const uint ns_;
	const uint nf_;
	const uint ntnsns_;
	const uint ntns_;

	// the actual values and errors.
        T *val_;
        T *err_;

};
typedef green_function<std::complex<double> > matsubara_green_function_t;
typedef green_function<double> itime_green_function_t;
std::ostream &operator<<(std::ostream &os, const green_function<double> &v);
std::ostream &operator<<(std::ostream &os, const green_function<std::complex<double> > &v);


double kinetic_energy(const multiple_vector_type &G_tau, const double &beta, const double &t);
