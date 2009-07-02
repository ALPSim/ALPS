#include<iostream>
#include"green_function.h"
#include"fouriertransform.h"
#include"cluster.h"
#include<alps/parameter.h>

typedef class moments{
  public:
    moments(int flavors, std::vector<double> &epsav, std::vector<double> &epssqav):
      c0_(flavors,0.),
      c1_(flavors,0.),
      c2_(flavors,0.),
      c3_(flavors,0.),
      flavors_(flavors){
      for(int i=0;i<epsav.size();++i){
        c0_[i]=0;
        c1_[i]=epssqav[i]-epsav[i]*epsav[i];
        c2_[i]=0;
        c3_[i]=0;
      }
    }
    double ftau(const double &tau, const double &beta, const int &index, const bool pos=false) const{ //pos in case G>0
      return (pos?1.:-1.)*0.5*c1_[index]+ (c2_[index]*0.25)*(-beta+2.*tau) + (c3_[index]*0.25)*(beta*tau-tau*tau);
    }
    std::complex<double> fomega(const std::complex<double> &iw, const int &index) const{
      std::complex<double> iwsq=iw*iw;
      return c1_[index]/iw + c2_[index]/(iwsq) + c3_[index]/(iw*iwsq);
    }
    const std::vector<double> &c0()const{return c0_;}
    const std::vector<double> &c1()const{return c1_;}
    const std::vector<double> &c2()const{return c2_;}
    const std::vector<double> &c3()const{return c3_;}
    const int& flavors() const {return flavors_;}
  private:
    std::vector<double> c0_,c1_,c2_,c3_;
    int flavors_;
} moments;

void fourier_to_itime(const matsubara_green_function_t hyb_matsubara, itime_green_function_t &hyb_itime, const moments &m, double beta, const int N_omega_cutoff);
int main(int argc, char **argv){
  if(argc!=3) throw std::invalid_argument("wrong number of arguments! use ./program_name parms input_g0");

  std::ifstream parm_file(argv[1]);
  alps::Parameters parms(parm_file);
 
  boost::shared_ptr<FourierTransformer> fourier_ptr;
  SimpleDCATransformer clusterhandler(parms);
  FourierTransformer::generate_transformer(parms, fourier_ptr, clusterhandler.epsav(), clusterhandler.epssqav());

  std::cout<<"moments are: "<<std::endl;
  std::vector<double> epsav=clusterhandler.epsav();
  std::vector<double> epssqav=clusterhandler.epssqav();
  for(int i=0;i<epsav.size();++i){
    std::cout<<i<<" "<<epsav[i]<<" "<<epssqav[i]<<std::endl;
  }

  std::ifstream G0_file(argv[2]);
  int n_tau  =static_cast<int>(parms["N"]);
  int n_sites=static_cast<int>(parms.value_or_default("SITES",1));
  double beta=static_cast<double>(parms["BETA"]);
  int n_matsubara=static_cast<int>(parms["NMATSUBARA"]);
  double mu=static_cast<double>(parms["MU"]);
  int N_omega_cutoff=50;

  matsubara_green_function_t g0(n_matsubara, n_sites, 2);
  matsubara_green_function_t hyb(n_matsubara, n_sites, 2);
  itime_green_function_t hyb_itime(n_tau, n_sites, 2);
  g0.read(argv[2]);

  //std::cout<<g0<<std::endl;

  for(int freq=0;freq<n_matsubara;++freq){
    std::complex<double> iomegan=std::complex<double>(0.,(2.*freq+1)*M_PI/beta);
    for(int spin=0;spin<2;++spin){
      for(int site=0;site<n_sites;++site){
        hyb(freq,site,site,spin)=iomegan+mu -epsav[site]-1./g0(freq,site,site,spin);
      }
    }
  }
  moments hybridization_moments(n_sites, epsav, epssqav);
  fourier_to_itime(hyb, hyb_itime, hybridization_moments, beta, N_omega_cutoff);
 std::cout<<hyb_itime<<std::endl;


}

void fourier_to_itime(const matsubara_green_function_t hyb_matsubara, itime_green_function_t &hyb_itime, const moments &m, double beta, const int N_omega_cutoff){
//  std::cout<<"hyb matsubara is: "<<hyb_matsubara<<std::endl;
  matsubara_green_function_t hyb_matsubara_nomodel(hyb_matsubara);
  hyb_itime.clear();
  double t_res=(double)hyb_itime.ntime();
  for(uint i=0;i<hyb_matsubara.nfreq();++i){
    std::complex<double> iomegan(0., (2.*i+1)*M_PI/beta);
    for(uint f=0;f<hyb_matsubara.nflavor();++f){
      for(uint s=0;s<hyb_matsubara.nsite();++s){
        hyb_matsubara_nomodel(i,s,s,f)-=m.fomega(iomegan,s);
      }
    }
  }
  for(uint f=0;f<hyb_matsubara.nflavor();++f){
    for(uint s=0;s<hyb_matsubara.nsite();++s){
      for(uint i=0;i<hyb_itime.ntime();++i){
        for(uint n=0;n<hyb_matsubara.nfreq();++n){
          //double nw=1;
          double fermi=1.;//(std::exp((n-N_omega_cutoff)/nw)+1);
          hyb_itime(i,s,s,f)+=2./beta*(std::complex<double>(cos((2*n+1)*M_PI*i/t_res), -sin((2*n+1)*M_PI*i/t_res))*hyb_matsubara_nomodel(n,s,s,f)).real()*fermi;
        }
      double tau=i*beta/t_res;
      hyb_itime(i,s,s,f)+=m.ftau(tau, beta,s);
      }
    }
  }
  //std::cout<<"hyb itime is: "<<hyb_itime<<std::endl;
}
