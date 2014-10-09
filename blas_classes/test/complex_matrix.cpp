#include <sys/types.h>
#include <iostream>
#include "../general_matrix.h"


int main(int argc, char **argv){
  std::cout<<"general matrix test: "<<std::endl;

  blas::general_matrix<std::complex<double> > M(2);
  M(0,0)=1.;
  M(0,1)=2.;
  M(1,0)=3.;
  M(1,1)=4.;
  std::cout<<"starting with matrix: "<<M<<std::endl;
  blas::general_matrix<std::complex<double> > M2(2);
  M2(0,0)=7.;
  M2(0,1)=2.;
  M2(1,0)=3.;
  M2(1,1)=4.;
  std::cout<<"starting with matrix M2: "<<M2<<std::endl;
  
  std::cout<<"difference: "<<M2-M<<std::endl;
  std::cout<<"sum: "<<M2+M<<std::endl;

  blas::general_matrix<std::complex<double> > M3(M2);
  M2.invert();
  std::cout<<M3*M2<<std::endl;

  /*blas::general_matrix<std::complex<double> > M4(2), M5;
  M4(0,0)=1.;
  M4(0,1)=2.;
  M4(1,0)=3.;
  M4(1,1)=4.;
  M5=M4;
  M4.invert();
  std::cout<<M4*M5<<std::endl;*/
} 
