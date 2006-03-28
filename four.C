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
  std::ifstream G_tau_file(argv[5]);
  assert(G_tau_file.is_open());
  vector_type G_tau;
  complex_vector_type G_omega;
  G_tau.resize(N+1);
  G_omega.resize(N);
  for (int i=0; i<N+1; i++) {
    double ignored;
    G_tau_file >> ignored;
    G_tau_file >> G_tau[i];
  }
  std::cerr<<"using G_tau: "<<G_tau<<std::endl;
  forward_FT(beta, -mu, mu*mu+t*t, G_tau, G_omega); 
  std::cout<<G_omega<<std::endl;
}
