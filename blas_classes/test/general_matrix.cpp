#include <sys/types.h>
#include <general_matrix.h>


int main(int argc, char **argv){
  std::cout<<"general matrix test: "<<std::endl;

  blas::general_matrix<double> M(2);
  M(0,0)=1;
  M(0,1)=2;
  M(1,0)=3;
  M(1,1)=4;
  std::cout<<"starting with matrix: "<<M<<std::endl;
  blas::general_matrix<double> M2(2);
  M2(0,0)=1;
  M2(0,1)=2;
  M2(1,0)=3;
  M2(1,1)=4;
  std::cout<<"starting with matrix M2: "<<M2<<std::endl;
  
  std::cout<<"difference: "<<M2-M<<std::endl;
  std::cout<<"sum: "<<M2+M<<std::endl;

  blas::general_matrix<double> M3(M2);
  M2.invert();
  std::cout<<M3*M2<<std::endl;

  blas::general_matrix<std::complex<double> > M4(2), M5;
  M4(0,0)=1.;
  M4(0,1)=2.;
  M4(1,0)=3.;
  M4(1,1)=4.;
  M5=M4;
  M4.invert();
  std::cout<<M4*M5<<std::endl;
  blas::general_matrix<std::complex<double> > M6(2,3), M7(3,1);
  M6(0,0)=1.; M6(0,1)=2.; 
  M6(0,1)=3.; M6(1,1)=5.; 
  M6(0,2)=6.; M6(2,1)=7.;
  M7(0,0)=3.; M7(1,0)=2.; M7(2,0)=3.;
  std::cout<<M6*M7<<std::endl;
} 
