#include <iostream>
#include <fstream>
#include "types.h"
#include "fouriertransform.h"
//call this: four N beta mu G_tau and read G_omega from stdout. G_tau has to
//include beta, G_tau[N]=beta, i.e. G_tau has N+1 points.
int main(int argc, char** argv){
  if(argc!=6){
    std::cerr<<"use: four N beta mu t G_tau"<<std::endl;
    exit(1);
  }
  int N=atoi(argv[1]);
  double beta=atof(argv[2]);
  double mu=atof(argv[3]);
  double t=atof(argv[4]);
  std::cerr<<"using N: "<<N<<" beta ="<<beta<<" mu= "<<mu<<" and t= "<<t<<std::endl;
  BetheFourierTransformer trans(beta, mu, t);
  std::ifstream G_tau_file(argv[5]);
  assert(G_tau_file.is_open());
  itime_green_function_t G_tau(N+1,1,1);
  matsubara_green_function_t G_omega(N,1,1);
  for (int i=0; i<N+1; i++) {
    double ignored;
    G_tau_file >> ignored;
    G_tau_file >> G_tau(i,0);
  }
  std::cerr<<"using G_tau: "<<G_tau<<std::endl;
  trans.forward_ft(G_tau, G_omega); 
  for (int i=0; i<N; i++) {
    std::cout.precision(20);
    std::cout<<(2*i+1)*M_PI/beta<<" "<<G_omega(i,0).imag()<<std::endl;
  }
  //std::cout<<G_omega<<std::endl;
}
