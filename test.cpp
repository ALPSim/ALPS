#include <alps/config.h> // needed to set up correct bindings
#include <boost/numeric/bindings/ublas.hpp>
#include <boost/numeric/bindings/lapack/driver/gesv.hpp>

#include "fouriertransform.h"
#include "cluster.h"
#include <valarray>
#include <alps/parameter.h>
#include <gtest/gtest.h>
typedef boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> dense_matrix;

TEST(Fouriertransform, nambuCDMFT){

  int n_matsubara=1024;
  int n_flavor=2;
  int n_site=4;
  double beta =10;
  double mu=1.5;
  double h=0;
  std::vector<std::vector<double> > eps_f(n_flavor, std::vector<double>(n_site*n_site,0.00));
  std::vector<std::vector<double> > epssq_f(n_flavor, std::vector<double>(n_site*n_site,0.0001));
  nambu_matsubara_green_function_t G_omega(n_matsubara, n_site, 2, nondiagonal);
  //compare it with existing correct result
  G_omega.read("G0_omega_input_r");
  nambu_itime_green_function_t G_tau_standard(n_matsubara+1, n_site, 2, nondiagonal);
  G_tau_standard.read("G0_tau_input");
  
  ClusterG0FourierTransformer testFourier(beta,mu,h,n_flavor,n_site,eps_f,epssq_f);
  nambu_itime_green_function_t G_tau(n_matsubara+1, n_site, 2, nondiagonal);
  testFourier.frequency_to_time_ft(G_tau,G_omega);
  double error=10e-6;
  for(int i=0;i<n_matsubara;i++){
    for(int s1=0;s1<n_site;s1++){
      for(int s2=0;s2<n_site;s2++){
        EXPECT_LT(abs(G_tau(i,s1,s2,0,0)-G_tau_standard(i,s1,s2,0,0)),error);
      }
    }
  }
  //check the value at tau=0. G(0)=2/beta*sum(G.real)
  dense_matrix G_0(n_site,n_site);
  G_0.clear();
  for(int s1=0;s1<n_site;s1++){
    for(int s2=0;s2<n_site;s2++){
      for(int freq=0;freq<n_matsubara;freq++){
        G_0(s1,s2)+=G_omega(freq,s1,s2,0,0).real()*2/beta;
      }
      EXPECT_LT(abs(G_tau(0,s1,s2,0,0)-G_0(s1,s2)),error);
    }
  }
  
  //check the value at tau=beta/2. G(beta/2)=-2/beta*sum(G.Imag*(-1)^n)
  G_0.clear();
  for(int s1=0;s1<n_site;s1++){
    for(int s2=0;s2<n_site;s2++){
      for(int freq=0;freq<n_matsubara;freq++){
        double sign=-2.*(freq%2)+1;
        G_0(s1,s2)+=-G_omega(freq,s1,s2,0,0).imag()*2/beta*sign;
      }
      EXPECT_LT(abs(G_tau(n_matsubara/2,s1,s2,0,0)-G_0(s1,s2)),error);
    }
  }

}


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
