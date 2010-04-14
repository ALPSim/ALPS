/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Philipp Werner <werner@itp.phys.ethz.ch>
 *                              Emanuel Gull <gull@phys.columbia.edu>
 *                              Matthias Troyer <troyer@comp-phys.org>
 *
 *
 * THIS SOFTWARE NEEDS AN APPROPRIATE LICENSE BLOCK HERE
 *****************************************************************************/
#include "types.h"
#include <alps/model.h>
#include <alps/lattice.h>
#include "impurity.h"
#include "moves.h"
#include "fouriertransform.h"
#include "xml.h"
#include <alps/alea.h>
//#include "matrix_builder.h"
using namespace std;
using namespace alps;

inline int compute_perm(double *times){
  //very naive sorting for finding out sign of permutation.
  int sign=1;
  if(times[0]>times[1]){std::swap(times[0], times[1]); sign*=-1; } //permute once - minus sign
  if(times[0]>times[2]){std::swap(times[0], times[2]); sign*= 1; } //permute twice - plus sign
  if(times[0]>times[3]){std::swap(times[0], times[3]); sign*=-1; } //permute thrice - minus sign
  if(times[1]>times[2]){std::swap(times[1], times[2]); sign*=-1; }
  if(times[1]>times[3]){std::swap(times[1], times[3]); sign*=+1; }
  if(times[2]>times[3]){std::swap(times[2], times[3]); sign*=-1; }
  //std::cout<<"sorted times: "<<times[0]<<" "<<times[1]<<" "<<times[2]<<" "<<times[3]<<" "<<sign<<std::endl;
  return sign;
}
void HybridizationRun::measure_fourpoint()
{
  
  // increment sweep count
  
  static double BETA ( static_cast<double>(parms["BETA"]));
  static int N4point ( static_cast<int>(parms["N4point"]));
  static int FLAVORS ( static_cast<int>(parms["FLAVORS"]));
  //static int N_order ( static_cast<int>(parms["N_ORDER"]));
  static int N_meas ( static_cast<int>(parms["N_MEAS"]));
  //static int N_shift (static_cast<int>(parms["N_SHIFT"]));
  
  /*************FOR BETA != 1**************/
  BETA=1;
  
  times full_segment(0,BETA);
  
  segment_container_t::const_iterator it1, it2;
  
  int N41=N4point+1;
  int N42=N41*N41;
  int N43=N42*N41;
  int N44=N42*N42;
  std::valarray<double> fourpoint_array;
  std::valarray<double> fourpoint_ud_array;
  
  fourpoint_array.resize(FLAVORS*N44);
  fourpoint_ud_array.resize(N44);
  for(unsigned i=0;i<fourpoint_array.size();++i) fourpoint_array[i]=0;
  for(unsigned i=0;i<fourpoint_ud_array.size();++i) fourpoint_ud_array[i]=0;
  
  //double sign_meas=0, s=1;
  G_meas = 0;
  
  assert(FLAVORS==2); //we specifically program the up-down function for two flavors.
  //std::cout<<"M: "<<M[0]<<std::endl;
  
  for (int i=0; i<N_meas; i++) {
    for (int j=0; j<FLAVORS; j++) {
      if (segments[j].size()>0) {
        //double N_div_beta=N/BETA;
        double N4_div_beta=N4point/BETA;
        double *start_times=new double[segments[j].size()];
        double *end_times=new double[segments[j].size()];
        int i=0;
        for(it1=segments[j].begin();it1!=segments[j].end();++it1){
          start_times[i]=it1->t_start();
          end_times[i]=it1->t_end();
          ++i;
        }
        for (int i=0; i<(int)M[j].size1(); i++) {
          for (int k=0; k<(int)M[j].size1(); k++) {
            if (M[j](k,i)!=0.) {
              double argument = end_times[i]-start_times[k];
              double bubble_sign=1.;
              if (argument < 0) {
                bubble_sign = -1.;
                argument += BETA;
              }
              
              //Fourpoint, up and down sector:
              if(j==0){
                double *start_times_1=new double[segments[1].size()];
                double *end_times_1=new double[segments[1].size()];
                int z=0;
                for(it1=segments[1].begin();it1!=segments[1].end();++it1){
                  start_times_1[z]=it1->t_start();
                  end_times_1[z]=it1->t_end();
                  ++z;
                }
                for (int p=0; p<(int)M[1].size1(); p++) {
                  for (int q=0; q<(int)M[1].size1(); q++) {
                    if (M[1](p,q)!=0.) {
                      double argument_1 = end_times_1[p]-start_times_1[q];
                      double bubble_sign_1=1.;
                      if (argument_1 < 0) {
                        bubble_sign_1 = -1.;
                        argument += BETA;
                      }
                      double prefactor=1;
                      int index_1=(int)(end_times[i]*N4_div_beta+0.5);
                      int index_2=(int)(end_times_1[p]*N4_div_beta+0.5);
                      int index_3=(int)(start_times[k]*N4_div_beta+0.5);
                      int index_4=(int)(start_times_1[q]*N4_div_beta+0.5);
                      if(index_1==0 || index_1==N4point) prefactor *=2; //these bins are just half the size of the others
                      if(index_2==0 || index_2==N4point) prefactor *=2;
                      if(index_3==0 || index_3==N4point) prefactor *=2;
                      if(index_4==0 || index_4==N4point) prefactor *=2;
                      //if two entries are in the same bin order them:
                      if(index_1==index_3 && end_times[i] > start_times[k]) prefactor*=-1;
                      if(index_2==index_4 && end_times_1[p] > start_times_1[q]) prefactor*=-1;
                      //compute sign of permutation -- this is imaginary time
                      //ordering
                      double times[4]={end_times[i], end_times_1[p], start_times[k], start_times_1[q]};
                      prefactor*=compute_perm(times); //the 'bubble sign' for the fourpoint function - correct ordering of four operators.
                      if(end_times[i]>start_times[k]) prefactor *=-1; //get rid of the time ordering for the green function
                      if(end_times_1[p]>start_times_1[q]) prefactor *=-1;
                      fourpoint_ud_array[index_1*N43+index_2*N42+index_3*N41+index_4]+=prefactor*M[0](k, i)*M[1](q, p);
                    }
                  }
                }
                delete[] start_times_1;
                delete[] end_times_1;
              }
            } //if M != 0
          } //for k
        } //for i
        delete [] start_times;
        delete [] end_times;
      }// if there are segments.
      
    }
    
    
  }
  int N4pow4=N4point*N4point*N4point*N4point;
  measurements.get<vec_obs_t>("Greens_4point_ud") << (fourpoint_ud_array*(N4pow4/(double)N_meas));//*sign;
  
}

void HybridizationSimFrequency::write_fourpoint() const{
  std::ofstream fpudfile("fourpoint_ud.dat");
  int N4point=parms["N4point"];
  assert(fpudfile.is_open());
  alps::RealVectorObsevaluator g4udobseval=get_measurements()["Greens_4point_ud"];
  std::valarray<double> fourpoint_ud_mean=g4udobseval.mean();
  std::valarray<double> fourpoint_ud_err=g4udobseval.error();
  int N4=N4point+1; 
  for(int i=0;i<N4point+1;++i){ //c_up
    for(int j=0;j<N4point+1;++j){ //c_down
      for(int k=0;k<N4point+1;++k){ //c_up_dagger
        for(int l=0;l<N4point+1;++l){ //c_down_dagger
          int z=N4*N4*N4*i+N4*N4*j+N4*k+l;
          fpudfile<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<fourpoint_ud_mean[z]<<
          " "<<fourpoint_ud_err[z]<</*" "<<fourpoint_tau[z]<<*/std::endl;
        }
      }
    }
  }
}
