#include<iostream>
#include"types.h"
#include"alps/parameter.h"
#include"green_function.h"
#include"cluster.h"
#include <boost/program_options.hpp>
#include "alps/numeric/matrix.hpp"

double diff(matsubara_green_function_t &gf,matsubara_green_function_t &gf2);
void read_parameters(int &iteration, std::string &pfilename, std::string &input_hdf5, std::string &output_hdf5, bool &nonint, int argc, char **argv);
void invert_2x2(const alps::numeric::matrix<std::complex<double> > &M, alps::numeric::matrix<std::complex<double> > &Minv, const std::complex<double> &alpha, const std::complex<double> &beta);
void invert_2x2(const alps::numeric::matrix<std::complex<double> > &M, alps::numeric::matrix<std::complex<double> > &Minv);
void evaluate_AFM_RBZ(const matsubara_green_function_t &selfenergy, matsubara_green_function_t &gloc, matsubara_green_function_t &g0, const std::vector<std::vector<std::vector<double> > > &epsilon_k_kx_ky, int N, double beta, double mu, double h);
void evaluate_AFM_2site(const matsubara_green_function_t &selfenergy, matsubara_green_function_t &gloc, matsubara_green_function_t &g0, const std::vector<std::vector<std::vector<double> > > &epsilon_k_kx_ky, int N, double beta, double mu, double h);
int main(int argc, char **argv){
  
  std::string pfilename, input_hdf5, output_hdf5;
  double mu, beta, h, t;
  int nfreq;
  bool nonint;
  int ns=4, N=20;
  int iteration=0;
  
  read_parameters(iteration, pfilename, input_hdf5, output_hdf5, nonint, argc, argv);
  std::ifstream parm_file(pfilename.c_str());
  if(!parm_file.is_open()) throw std::runtime_error("parm file not open.");
  alps::Parameters parms(parm_file);
  
  mu=parms["MU"];
  beta=parms["BETA"];
  nfreq=parms["NMATSUBARA"];
  t=parms["t"];
  h=nonint?parms["H_INIT"]:parms["H"];

  boost::shared_ptr<DCATransformer> clusterhandler;
  DCATransformer::generate_transformer(parms, clusterhandler);
  
  matsubara_green_function_t selfenergy(nfreq,ns, 2, blockdiagonal); selfenergy.clear();
  matsubara_green_function_t g0(nfreq,ns, 2, blockdiagonal); g0.clear();
  matsubara_green_function_t gloc(nfreq,ns, 2, blockdiagonal); gloc.clear();
  matsubara_green_function_t g0_rs(nfreq,ns, 2, nondiagonal); g0_rs.clear();
  matsubara_green_function_t gloc_rs(nfreq,ns, 2, nondiagonal); gloc_rs.clear();
  
  alps::numeric::matrix<std::complex<double> > G0inv_2site_up(2,2),G0inv_2site_dn(2,2),G0inv_kq_up(2,2),G0inv_kq_dn(2,2),Ginv_kq_up(2,2),Ginv_kq_dn(2,2);
  alps::numeric::matrix<std::complex<double> > G_kq_up(2,2),G_kq_dn(2,2),G0_kq_up(2,2),G0_kq_dn(2,2),G_kq_up_tmp(2,2),G_kq_dn_tmp(2,2);
  
  if(!nonint){
    std::cout<<"computing sigma"<<std::endl;
    alps::hdf5::archive ar(input_hdf5, alps::hdf5::archive::READ);
    g0_rs.read_hdf5(ar, "/G0_omega_rs");
    //g0.read_hdf5(ar, "/G0_omega_ks");
    g0=clusterhandler->transform_into_k_space(g0_rs);
    {
      std::stringstream g0_file_name; g0_file_name<<"green0_function_lastit_"<<iteration;
      std::ofstream g0_file(g0_file_name.str().c_str());
      g0_file<<std::make_pair(g0, beta);
    }
    gloc_rs.read_hdf5(ar, "/G_omega");
    gloc=clusterhandler->transform_into_k_space(gloc_rs);
    {
      std::stringstream g_file_name; g_file_name<<"green_function_lastit_"<<iteration;
      std::ofstream g_file(g_file_name.str().c_str());
      g_file<<std::make_pair(gloc, beta);
    }
    for (int i=0;i<nfreq;++i){
      for(int k=0;k<ns;k+=2){
        G0_kq_up(0,0)=g0(i,k  ,k  ,0);
        G0_kq_up(0,1)=g0(i,k  ,k+1,0);
        G0_kq_up(1,0)=g0(i,k+1,k  ,0);
        G0_kq_up(1,1)=g0(i,k+1,k+1,0);
        
        G0_kq_dn(0,0)=g0(i,k  ,k  ,1);
        G0_kq_dn(0,1)=g0(i,k  ,k+1,1);
        G0_kq_dn(1,0)=g0(i,k+1,k  ,1);
        G0_kq_dn(1,1)=g0(i,k+1,k+1,1);
        
        G_kq_up(0,0)=gloc(i,k  ,k  ,0);
        G_kq_up(0,1)=gloc(i,k  ,k+1,0);
        G_kq_up(1,0)=gloc(i,k+1,k  ,0);
        G_kq_up(1,1)=gloc(i,k+1,k+1,0);
        
        G_kq_dn(0,0)=gloc(i,k  ,k  ,1);
        G_kq_dn(0,1)=gloc(i,k  ,k+1,1);
        G_kq_dn(1,0)=gloc(i,k+1,k  ,1);
        G_kq_dn(1,1)=gloc(i,k+1,k+1,1);
        
        invert_2x2(G0_kq_up,G0inv_kq_up);
        invert_2x2(G_kq_up ,Ginv_kq_up);
        invert_2x2(G0_kq_dn,G0inv_kq_dn);
        invert_2x2(G_kq_dn ,Ginv_kq_dn);
        selfenergy(i,k  ,k  ,0)=G0inv_kq_up(0,0)-Ginv_kq_up(0,0);
        selfenergy(i,k  ,k+1,0)=G0inv_kq_up(0,1)-Ginv_kq_up(0,1);
        selfenergy(i,k+1,k  ,0)=G0inv_kq_up(1,0)-Ginv_kq_up(1,0);
        selfenergy(i,k+1,k+1,0)=G0inv_kq_up(1,1)-Ginv_kq_up(1,1);

        selfenergy(i,k  ,k  ,1)=G0inv_kq_dn(0,0)-Ginv_kq_dn(0,0);
        selfenergy(i,k  ,k+1,1)=G0inv_kq_dn(0,1)-Ginv_kq_dn(0,1);
        selfenergy(i,k+1,k  ,1)=G0inv_kq_dn(1,0)-Ginv_kq_dn(1,0);
        selfenergy(i,k+1,k+1,1)=G0inv_kq_dn(1,1)-Ginv_kq_dn(1,1);
}
    }
    std::stringstream sigma_file_name; sigma_file_name<<"selfenergy_"<<iteration;
    std::ofstream sigma_file(sigma_file_name.str().c_str());
    sigma_file<<std::make_pair(selfenergy, beta);
  }
  
  
  std::vector<std::vector<std::vector<double> > > epsilon_k_kx_ky(ns, std::vector<std::vector<double> >(N, std::vector<double>(N, 0.)));
  std::vector<std::vector<std::vector<bool> > > inside_RBZ(ns, std::vector<std::vector<bool> >(N, std::vector<bool>(N, false)));
  std::vector<int > inside_RBZ_count(ns, 0);
  std::vector<double> KX(ns), KY(ns);
  int inside_RBZ_tot=0.;
  KX[0]=0; KX[1]=M_PI; KX[2]=0   ; KX[3]=M_PI;
  KY[1]=0; KY[1]=M_PI; KY[2]=M_PI; KY[3]=0   ;
  for(int k=0;k<ns;k++){
    for(int kx=0;kx<N;++kx){
      for(int ky=0;ky<N;++ky){
        double kxcoord=KX[k]+(kx-N/2.)*M_PI/N; if(kxcoord>M_PI) kxcoord -= 2*M_PI; if(kxcoord<-M_PI) kxcoord += 2*M_PI;
        double kycoord=KY[k]+(ky-N/2.)*M_PI/N; if(kycoord>M_PI) kycoord -= 2*M_PI; if(kycoord<-M_PI) kycoord += 2*M_PI;
        epsilon_k_kx_ky[k][kx][ky]=-2*t*(cos(kxcoord)+cos(kycoord));
        if((kxcoord+kycoord <= M_PI) && (kxcoord+kycoord > -M_PI) && (kxcoord-kycoord <=M_PI) && (kxcoord-kycoord>-M_PI)){
          inside_RBZ[k][kx][ky]=true;
          inside_RBZ_count[k]++;
          inside_RBZ_tot++;
        }
      }
    }
  }
  
  std::cout<<"beta: "<<beta<<std::endl;
  if(inside_RBZ_tot != N*N*ns/2) throw std::logic_error("problem with BZ tiling: not half the points inside RBZ!");

  //evaluate_AFM_RBZ(selfenergy, gloc, g0, epsilon_k_kx_ky, N,  beta, mu, h);
  evaluate_AFM_2site(selfenergy, gloc, g0, epsilon_k_kx_ky, N,  beta, mu, h);
  alps::hdf5::archive ar(output_hdf5, alps::hdf5::archive::WRITE);
  parms["INFILE"]=output_hdf5;
  parms["READWRITE_REAL_SPACE"]=1;
  g0=clusterhandler->transform_into_real_space(g0);
  ar<<alps::make_pvp("/parameters", parms);
  g0.write_hdf5(ar, "/G0");
  
  {
    std::stringstream g_file_name; g_file_name<<"green_function_"<<iteration;
    std::ofstream g_file(g_file_name.str().c_str());
    g_file<<std::make_pair(gloc, beta);
  }
  {
    std::stringstream g0_file_name; g0_file_name<<"green0_function_"<<iteration;
    std::ofstream g0_file(g0_file_name.str().c_str());
    g0_file<<std::make_pair(g0, beta);
  }
}

void read_parameters(int &iteration, std::string &pfilename, std::string &input_hdf5, std::string &output_hdf5, bool &nonint, int argc, char **argv){
  namespace po = boost::program_options;
  po::options_description desc("Allowed options");
  nonint=false;
  desc.add_options()
  ("help", "show this help")
  ("pfile", po::value<std::string>(&pfilename)->default_value("ctaux.param"), "parameter file")
  ("input_hdf5", po::value<std::string>(&input_hdf5)->default_value("sc_in.h5"), "input hdf5 file")
  ("output_hdf5", po::value<std::string>(&output_hdf5)->default_value("sc_out.h5"), "output hdf5 file")
  ("nonint", "non-interacting if specified.")
  ("iteration", po::value<int>(&iteration)->default_value(0), "iteration");
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout<<desc;
    exit(1);
  }
  if (vm.count("nonint")) {
    nonint=true;
  }
  
  if(!vm.count("pfile")) throw std::runtime_error("you need to specify the parameter file as pfile.");
}
void invert_2x2(const alps::numeric::matrix<std::complex<double> > &M, alps::numeric::matrix<std::complex<double> > &Minv){
  invert_2x2(M, Minv, 1., 0.);
}
///Minv = \alpha M^{-1}+\beta Minv
void invert_2x2(const alps::numeric::matrix<std::complex<double> > &M, alps::numeric::matrix<std::complex<double> > &Minv, const std::complex<double> &alpha, const std::complex<double> &beta){
  std::complex<double> det=M(0,0)*M(1,1)-M(0,1)*M(1,0);
  Minv(0,0)=alpha*( M(1,1))/det+beta*Minv(0,0);
  Minv(0,1)=alpha*(-M(0,1))/det+beta*Minv(0,1);
  Minv(1,0)=alpha*(-M(1,0))/det+beta*Minv(1,0);
  Minv(1,1)=alpha*( M(0,0))/det+beta*Minv(1,1);
}

void evaluate_AFM_RBZ(const matsubara_green_function_t &selfenergy, matsubara_green_function_t &gloc, matsubara_green_function_t &g0, const std::vector<std::vector<std::vector<double> > > &epsilon_k_kx_ky, int N, double beta, double mu, double h){

  alps::numeric::matrix<std::complex<double> > G0inv_2site_up(2,2),G0inv_2site_dn(2,2),G0inv_kq_up(2,2),G0inv_kq_dn(2,2),Ginv_kq_up(2,2),Ginv_kq_dn(2,2);
  alps::numeric::matrix<std::complex<double> > G_kq_up(2,2),G_kq_dn(2,2),G0_kq_up(2,2),G0_kq_dn(2,2),G_kq_up_tmp(2,2),G_kq_dn_tmp(2,2);
  for (int i=0;i<g0.nfreq();++i){
    std::complex<double> iomegan(0., (2.*i+1)*M_PI/beta);
    for(int k=0;k<g0.nsite();k+=2){
      G_kq_up*=0.;
      G_kq_dn*=0.;
      for(int kx=0;kx<N;++kx){
        for(int ky=0;ky<N;++ky){
          G0inv_kq_up(0,0)=iomegan+mu-epsilon_k_kx_ky[k][kx][ky];
          G0inv_kq_up(0,1)=h;
          G0inv_kq_up(1,0)=h;
          G0inv_kq_up(1,1)=iomegan+mu-epsilon_k_kx_ky[k+1][kx][ky];
          
          G0inv_kq_dn(0,0)=iomegan+mu-epsilon_k_kx_ky[k][kx][ky];
          G0inv_kq_dn(0,1)=-h;
          G0inv_kq_dn(1,0)=-h;
          G0inv_kq_dn(1,1)=iomegan+mu-epsilon_k_kx_ky[k+1][kx][ky];
          
          Ginv_kq_up(0,0) = (G0inv_kq_up(0,0)-selfenergy(i,k  ,k  ,0));
          Ginv_kq_up(0,1) = (G0inv_kq_up(0,1)-selfenergy(i,k  ,k+1,0));
          Ginv_kq_up(1,0) = (G0inv_kq_up(1,0)-selfenergy(i,k+1,k  ,0));
          Ginv_kq_up(1,1) = (G0inv_kq_up(1,1)-selfenergy(i,k+1,k+1,0));
          
          Ginv_kq_dn(0,0) = (G0inv_kq_dn(0,0)-selfenergy(i,k  ,k  ,1));
          Ginv_kq_dn(0,1) = (G0inv_kq_dn(0,1)-selfenergy(i,k  ,k+1,1));
          Ginv_kq_dn(1,0) = (G0inv_kq_dn(1,0)-selfenergy(i,k+1,k  ,1));
          Ginv_kq_dn(1,1) = (G0inv_kq_dn(1,1)-selfenergy(i,k+1,k+1,1));
          
          invert_2x2(Ginv_kq_up, G_kq_up, 1./(N*N), 1.); //invert and add G
          invert_2x2(Ginv_kq_dn, G_kq_dn, 1./(N*N), 1.);
        }
      }
      invert_2x2(G_kq_up,Ginv_kq_up); //compute the inverse local G
      invert_2x2(G_kq_dn,Ginv_kq_dn); //compute the inverse local G
      
      gloc(i,k  ,k  ,0)=G_kq_up(0,0);
      gloc(i,k  ,k+1,0)=G_kq_up(0,1);
      gloc(i,k+1,k  ,0)=G_kq_up(1,0);
      gloc(i,k+1,k+1,0)=G_kq_up(1,1);
      
      gloc(i,k  ,k  ,1)=G_kq_dn(0,0);
      gloc(i,k  ,k+1,1)=G_kq_dn(0,1);
      gloc(i,k+1,k  ,1)=G_kq_dn(1,0);
      gloc(i,k+1,k+1,1)=G_kq_dn(1,1);
      
      
      G0inv_kq_up(0,0)=selfenergy(i,k  ,k  ,0)+Ginv_kq_up(0,0);
      G0inv_kq_up(0,1)=selfenergy(i,k  ,k+1,0)+Ginv_kq_up(0,1);
      G0inv_kq_up(1,0)=selfenergy(i,k+1,k  ,0)+Ginv_kq_up(1,0);
      G0inv_kq_up(1,1)=selfenergy(i,k+1,k+1,0)+Ginv_kq_up(1,1);
      
      G0inv_kq_dn(0,0)=selfenergy(i,k  ,k  ,1)+Ginv_kq_dn(0,0);
      G0inv_kq_dn(0,1)=selfenergy(i,k  ,k+1,1)+Ginv_kq_dn(0,1);
      G0inv_kq_dn(1,0)=selfenergy(i,k+1,k  ,1)+Ginv_kq_dn(1,0);
      G0inv_kq_dn(1,1)=selfenergy(i,k+1,k+1,1)+Ginv_kq_dn(1,1);
      
      invert_2x2(G0inv_kq_up,G0_kq_up); //compute the local G0
      invert_2x2(G0inv_kq_dn,G0_kq_dn); //compute the local G0
      g0(i,k  ,k  ,0)=G0_kq_up(0,0);
      g0(i,k  ,k+1,0)=G0_kq_up(0,1);
      g0(i,k+1,k  ,0)=G0_kq_up(1,0);
      g0(i,k+1,k+1,0)=G0_kq_up(1,1);
      
      g0(i,k  ,k  ,1)=G0_kq_dn(0,0);
      g0(i,k  ,k+1,1)=G0_kq_dn(0,1);
      g0(i,k+1,k  ,1)=G0_kq_dn(1,0);
      g0(i,k+1,k+1,1)=G0_kq_dn(1,1);
    }
  }
}

void evaluate_AFM_2site(const matsubara_green_function_t &selfenergy, matsubara_green_function_t &gloc, matsubara_green_function_t &g0, const std::vector<std::vector<std::vector<double> > > &epsilon_k_kx_ky, int N, double beta, double mu, double h){
  
  alps::numeric::matrix<std::complex<double> > G0inv_2site_up(2,2),G0inv_2site_dn(2,2),G0inv_kq_up(2,2),G0inv_kq_dn(2,2),Ginv_kq_up(2,2),Ginv_kq_dn(2,2);
  alps::numeric::matrix<std::complex<double> > G_kq_up(2,2),G_kq_dn(2,2),G0_kq_up(2,2),G0_kq_dn(2,2),G_kq_up_tmp(2,2),G_kq_dn_tmp(2,2);
  alps::numeric::matrix<std::complex<double> > U(2,2); U(0,0)=1.; U(1,0)=1.; U(0,1)=1.; U(1,1)=-1.; U/=sqrt(2.);
  for (int i=0;i<g0.nfreq();++i){
    std::complex<double> iomegan(0., (2.*i+1)*M_PI/beta);
    for(int k=0;k<g0.nsite();k+=2){
      G_kq_up*=0.;
      G_kq_dn*=0.;
      for(int kx=0;kx<N;++kx){
        for(int ky=0;ky<N;++ky){
          G0inv_2site_up(0,0)=iomegan+mu+h;
          G0inv_2site_up(0,1)=epsilon_k_kx_ky[k][kx][ky];
          G0inv_2site_up(1,0)=epsilon_k_kx_ky[k][kx][ky];
          G0inv_2site_up(1,1)=iomegan+mu-h;
          
          G0inv_2site_dn(0,0)=iomegan+mu-h;
          G0inv_2site_dn(0,1)=epsilon_k_kx_ky[k][kx][ky];
          G0inv_2site_dn(1,0)=epsilon_k_kx_ky[k][kx][ky];
          G0inv_2site_dn(1,1)=iomegan+mu+h;
          
          G0inv_kq_up=U*G0inv_2site_up*U;
          G0inv_kq_dn=U*G0inv_2site_dn*U;
          
          int K=k; int Q=k+1;

          Ginv_kq_up(0,0) = (G0inv_kq_up(0,0)-selfenergy(i,K  ,K  ,0));
          Ginv_kq_up(0,1) = (G0inv_kq_up(0,1)-selfenergy(i,K  ,Q  ,0));
          Ginv_kq_up(1,0) = (G0inv_kq_up(1,0)-selfenergy(i,Q  ,K  ,0));
          Ginv_kq_up(1,1) = (G0inv_kq_up(1,1)-selfenergy(i,Q  ,Q  ,0));
          
          Ginv_kq_dn(0,0) = (G0inv_kq_dn(0,0)-selfenergy(i,K  ,K  ,1));
          Ginv_kq_dn(0,1) = (G0inv_kq_dn(0,1)-selfenergy(i,K  ,Q  ,1));
          Ginv_kq_dn(1,0) = (G0inv_kq_dn(1,0)-selfenergy(i,Q  ,K  ,1));
          Ginv_kq_dn(1,1) = (G0inv_kq_dn(1,1)-selfenergy(i,Q  ,Q  ,1));
          
          invert_2x2(Ginv_kq_up, G_kq_up, 1./(N*N), 1.); //invert and add G
          invert_2x2(Ginv_kq_dn, G_kq_dn, 1./(N*N), 1.);
        }
      }
      invert_2x2(G_kq_up,Ginv_kq_up); //compute the inverse local G
      invert_2x2(G_kq_dn,Ginv_kq_dn); //compute the inverse local G
      
      gloc(i,k  ,k  ,0)=G_kq_up(0,0);
      gloc(i,k  ,k+1,0)=G_kq_up(0,1);
      gloc(i,k+1,k  ,0)=G_kq_up(1,0);
      gloc(i,k+1,k+1,0)=G_kq_up(1,1);
      
      gloc(i,k  ,k  ,1)=G_kq_dn(0,0);
      gloc(i,k  ,k+1,1)=G_kq_dn(0,1);
      gloc(i,k+1,k  ,1)=G_kq_dn(1,0);
      gloc(i,k+1,k+1,1)=G_kq_dn(1,1);
      
      
      G0inv_kq_up(0,0)=selfenergy(i,k  ,k  ,0)+Ginv_kq_up(0,0);
      G0inv_kq_up(0,1)=selfenergy(i,k  ,k+1,0)+Ginv_kq_up(0,1);
      G0inv_kq_up(1,0)=selfenergy(i,k+1,k  ,0)+Ginv_kq_up(1,0);
      G0inv_kq_up(1,1)=selfenergy(i,k+1,k+1,0)+Ginv_kq_up(1,1);
      
      G0inv_kq_dn(0,0)=selfenergy(i,k  ,k  ,1)+Ginv_kq_dn(0,0);
      G0inv_kq_dn(0,1)=selfenergy(i,k  ,k+1,1)+Ginv_kq_dn(0,1);
      G0inv_kq_dn(1,0)=selfenergy(i,k+1,k  ,1)+Ginv_kq_dn(1,0);
      G0inv_kq_dn(1,1)=selfenergy(i,k+1,k+1,1)+Ginv_kq_dn(1,1);
      
      invert_2x2(G0inv_kq_up,G0_kq_up); //compute the local G0
      invert_2x2(G0inv_kq_dn,G0_kq_dn); //compute the local G0
      
      //std::cout<<"G0 kq up: "<<G0_kq_up<<std::endl;

      g0(i,k  ,k  ,0)=G0_kq_up(0,0);
      g0(i,k  ,k+1,0)=G0_kq_up(0,1);
      g0(i,k+1,k  ,0)=G0_kq_up(1,0);
      g0(i,k+1,k+1,0)=G0_kq_up(1,1);
      
      g0(i,k  ,k  ,1)=G0_kq_dn(0,0);
      g0(i,k  ,k+1,1)=G0_kq_dn(0,1);
      g0(i,k+1,k  ,1)=G0_kq_dn(1,0);
      g0(i,k+1,k+1,1)=G0_kq_dn(1,1);
    }
  }
}
