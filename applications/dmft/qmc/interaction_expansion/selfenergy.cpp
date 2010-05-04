/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Sebastian Fuchs <fuchs@theorie.physik.uni-goettingen.de>
 *                              Matthias Troyer <troyer@comp-phys.org>
 *
 *
 * THIS SOFTWARE NEEDS AN APPROPRIATE LICENSE BLOCK HERE
 *****************************************************************************/

#include "interaction_expansion.hpp"
#include <complex>
#include <alps/alea.h>
#include <alps/alea/simpleobseval.h>
#include <alps/scheduler/montecarlo.h>
#include <alps/osiris/dump.h>
#include <alps/osiris/std/vector.h>


typedef alps::SignedObservable<alps::RealVectorObservable> signed_vec_obs_t;
typedef alps::RealVectorObservable vec_obs_t;
typedef alps::SimpleObservable<double,alps::DetailedBinning<double> > simple_obs_t;
typedef const alps::SimpleObservable<double,alps::DetailedBinning<double> > const_simple_obs_t;


#ifdef SSE
#include<emmintrin.h>
class twocomplex{
public:
  inline twocomplex(){};
  inline twocomplex(const std::complex<double>&p, const std::complex<double> &q){
    r=_mm_loadl_pd(r, &(p.real())); //load two complex numbers.
    r=_mm_loadh_pd(r, &(q.real()));
    i=_mm_loadl_pd(i, &(p.imag()));
    i=_mm_loadh_pd(i, &(q.imag()));
  }
  inline void store(std::complex<double> &p, std::complex<double> &q){
    _mm_store_sd(&(p.real()), r);
    _mm_store_sd(&(p.imag()), i);
    _mm_storeh_pd(&(q.real()), r);
    _mm_storeh_pd(&(q.imag()), i);
  }
  __m128d r;
  __m128d i;
};

inline twocomplex fastcmult(const twocomplex &a, const twocomplex &b)
{
  twocomplex c;
  c.r = _mm_sub_pd(_mm_mul_pd(a.r,b.r), _mm_mul_pd(a.i,b.i));
  c.i = _mm_add_pd(_mm_mul_pd(a.r,b.i), _mm_mul_pd(a.i,b.r));
  return c;
}
#endif



void InteractionExpansionRun::compute_W_matsubara()
{
  static std::vector<std::vector<std::valarray<std::complex<double> > > >Wk(n_zone); 
  for(int z=0;z<n_zone;++z){
    Wk[z].resize(n_site);
    for(int j=0;j<n_site;++j){
      Wk[z][j].resize(n_matsubara);
      memset(&(Wk[z][j][0]), 0, sizeof(std::complex<double>)*(n_matsubara));
    }
  }
  measure_Wk(Wk, n_matsubara);
  measure_densities();
}



void InteractionExpansionRun::measure_Wk(std::vector<std::vector<std::valarray<std::complex<double> > > >& Wk, 
                                 const unsigned int nfreq) 
{
  for (int z=0; z<n_zone; ++z) {                     
    for (int k=0; k<n_site; k++) {                   
      for(int p=0;p<M[z].size();++p){         
        site_t site_p=M[z].creators()[p].s();
        M[z].creators()[p].compute_exp(n_matsubara, +1);
        for(int q=0;q<M[z].size();++q){       
          site_t site_q=M[z].creators()[q].s();
          M[z].annihilators()[q].compute_exp(n_matsubara, -1);
          std::complex<double>* Wk_z_k1_k2 = &Wk[z][k][0];
          const std::complex<double>* exparray_creators = M[z].creators()[p].exp_iomegat();
          const std::complex<double>* exparray_annihilators = M[z].annihilators()[q].exp_iomegat();
          std::complex<double> tmp = M[z](p,q);
#ifndef SSE
#pragma ivdep
          for(int o=0; o<nfreq; ++o){      
            *Wk_z_k1_k2++ += (*exparray_creators++)*(*exparray_annihilators++)*tmp;
          }
#else
#pragma ivdep
          for(int o=0;o<nfreq;o+=2) {
            twocomplex exp_c(*exparray_creators++,    *exparray_creators++); //load it all into xmm registers
            twocomplex exp_a(*exparray_annihilators++,*exparray_annihilators++);
            twocomplex tmp2(tmp,tmp);
            twocomplex product=fastcmult(fastcmult(exp_c,exp_a),tmp2);
            std::complex<double> Wk1, Wk2;
            product.store(Wk1, Wk2);
            *Wk_z_k1_k2++ += Wk1;
            *Wk_z_k1_k2++ += Wk2;
          }
#endif
        }
      }
    }
  }
  for(int flavor=0;flavor<n_zone;++flavor){
    for (int k=0; k<n_site; k++) {                   
      std::stringstream Wk_real_name, Wk_imag_name;
      Wk_real_name  << "Wk_real_"  << flavor << "_" << k << "_" << k;
      Wk_imag_name  << "Wk_imag_"  << flavor << "_" << k << "_" << k;
      std::valarray<double> Wk_real(nfreq);
      std::valarray<double> Wk_imag(nfreq);
      for (int w=0; w<nfreq; ++w) {
        Wk_real[w] = Wk[flavor][k][w].real();
        Wk_imag[w] = Wk[flavor][k][w].imag();
      }
      measurements.get<alps::SignedObservable<vec_obs_t> >(Wk_real_name.str().c_str()) <<Wk_real*sign;
      measurements.get<alps::SignedObservable<vec_obs_t> >(Wk_imag_name.str().c_str()) <<Wk_imag*sign;
    }
  }
}



void InteractionExpansionRun::measure_densities()
{
  std::vector< std::vector<double> > dens(n_zone);
  for(int z=0;z<n_zone;++z){
    dens[z].resize(n_site);
    memset(&(dens[z][0]), 0., sizeof(double)*(n_site));
  }
  double tau = beta*random_01();
  for (int z=0; z<n_zone; ++z) {                 
    double g0_tauj[M[z].size()];
    double M_g0_tauj[M[z].size()];
    double g0_taui[M[z].size()];
    for (int s=0;s<n_site;++s) {             
      for (int j=0;j<M[z].size();++j) 
	g0_tauj[j] = green0_spline(M[z].creators()[j].t()-tau, z, M[z].creators()[j].s(), s);
      for (int i=0;i<M[z].size();++i) 
        g0_taui[i] = green0_spline(tau-M[z].annihilators()[i].t(),z, s, M[z].annihilators()[i].s());
      if (M[z].size()>0)
	M[z].right_multiply(g0_tauj, M_g0_tauj);
      dens[z][s] += green0_spline(0,z,s,s);
      for (int j=0;j<M[z].size();++j) 
	dens[z][s] -= g0_taui[j]*M_g0_tauj[j]; 
    }
  }
  std::valarray<double> densities(0., n_zone);
  for (int z=0; z<n_zone; ++z) {                  
    std::valarray<double> densmeas(n_site);
    for (int i=0; i<n_site; ++i) {
      densities[z] += dens[z][i];
      densmeas[i] = 1+dens[z][i];
    }
    measurements.get<alps::SignedObservable<vec_obs_t> >("densities_"+boost::lexical_cast<std::string>(z)) << densmeas*sign;
    densities[z] /= n_site;
    densities[z] = 1 + densities[z];
  }
  measurements.get<alps::SignedObservable<vec_obs_t> >("densities") << densities*sign;
  double density_correlation = 0.;
  for (int i=0; i<n_site; ++i) {
    density_correlation += (1+dens[0][i])*(1+dens[1][i]);
  }
  density_correlation /= n_site;
  measurements["density_correlation"] << density_correlation*sign;
  std::valarray<double> ninj(n_site*n_site*4);
  for (int i=0; i<n_site; ++i) {
    for (int j=0; j<n_site; ++j) {
      ninj[i*n_site+j] = (1+dens[0][i])*(1+dens[0][j]);
      ninj[i*n_site+j+1] = (1+dens[0][i])*(1+dens[1][j]);
      ninj[i*n_site+j+2] = (1+dens[1][i])*(1+dens[0][j]);
      ninj[i*n_site+j+3] = (1+dens[1][i])*(1+dens[1][j]);
    }
  }
  measurements.get<alps::SignedObservable<vec_obs_t> >("n_i n_j") << ninj*sign;
}



void InteractionExpansionRun::compute_W_itime()
{
  static std::vector<std::vector<std::vector<std::valarray<double> > > >W_z_i_j(n_zone); 
  //first index: zone. Second index: momentum. Third index: self energy tau point.
  std::vector<std::vector<double> >density(n_zone);
  for(int z=0;z<n_zone;++z){
    W_z_i_j[z].resize(n_site);
    density[z].resize(n_site);
    for(int i=0;i<n_site;++i){
      W_z_i_j[z][i].resize(n_site);
      for(int j=0;j<n_site;++j){
        W_z_i_j[z][i][j].resize(n_self+1);
        memset(&(W_z_i_j[z][i][j][0]), 0, sizeof(double)*(n_self+1));
      }
    }
  }
  int ntaupoints=10; //# of tau points at which we measure.
  double tau_2[ntaupoints];
  for(int i=0; i<ntaupoints;++i) 
    tau_2[i]=beta*random_01();
  for(int z=0;z<n_zone;++z){                  //loop over zone
    double g0_tauj[M[z].size()*ntaupoints];
    double M_g0_tauj[M[z].size()*ntaupoints];
    double g0_taui[M[z].size()];
    for(int s2=0;s2<n_site;++s2){             //site loop - second site.
      for(int j=0;j<ntaupoints;++j) { 
	for(int i=0;i<M[z].size();++i){ //G0_{s_p s_2}(tau_p - tau_2) where we set t2=0.
	  g0_tauj[i*ntaupoints+j]=green0_spline(M[z].creators()[i].t()-tau_2[j], z, M[z].creators()[i].s(), s2);
	}
      }
      if (M[z].size()>0)
	M[z].matrix_right_multiply(g0_tauj, M_g0_tauj, ntaupoints);
      for(int j=0;j<ntaupoints;++j) {
	for(int p=0;p<M[z].size();++p){       //operator one
	  double sgn=1;
	  double delta_tau=M[z].creators()[p].t()-tau_2[j];
	  if(delta_tau<0){ 
	    sgn=-1; 
	    delta_tau+=beta; 
	  }
	  int bin=(int)(delta_tau/beta*n_self+0.5);
	  site_t site_p=M[z].creators()[p].s();
	  W_z_i_j[z][site_p][s2][bin]+=M_g0_tauj[p*ntaupoints+j]*sgn;
	}
      }
      for(int i=0;i<M[z].size();++i){
        g0_taui[i]=green0_spline(tau_2[0]-M[z].annihilators()[i].t(),z, s2, M[z].annihilators()[i].s());
      }
      density[z][s2]=green0_spline(0,z,s2,s2);
      for(int j=0;j<M[z].size();++j){
	density[z][s2]-= g0_taui[j]*M_g0_tauj[j*ntaupoints]  ; 
      }
    }
  }
  if(is_thermalized()){
    for(int flavor=0;flavor<n_zone;++flavor){
      for(int i=0;i<n_site;++i){
        for(int j=0;j<n_site;++j){
          std::stringstream W_name;
          W_name  <<"W_"  <<flavor<<"_"<<i<<"_"<<j;
          measurements.get<alps::SignedObservable<vec_obs_t> >(W_name  .str().c_str()) 
	    << W_z_i_j[flavor][i][j]*(sign/ntaupoints);
        }
        std::stringstream density_name;
        density_name<<"density_"<<flavor<<"_"<<i;
        measurements[density_name.str().c_str()]<<density[flavor][i]*sign;
        if(n_zone==2){ //then we know how to compute Sz^2
          std::stringstream sz_name, sz2_name, sz0_szj_name;
          sz_name<<"Sz_"<<i; sz2_name<<"Sz2_"<<i; sz0_szj_name<<"Sz0_Sz"<<i;
          measurements[sz_name.str().c_str()]<<(density[0][i]-density[1][i])*sign;
          measurements[sz2_name.str().c_str()]<<(density[0][i]-density[1][i])*(density[0][i]-density[1][i])*sign;
          measurements[sz0_szj_name.str().c_str()]<<(density[0][0]-density[1][0])*(density[0][i]-density[1][i])*sign;
        }
      }
    }
  }
}




void InteractionExpansionSim::evaluate_selfenergy_measurement_matsubara(const alps::ObservableSet &gathered_measurements, 
                                                                matsubara_green_function_t &green_matsubara_measured,
                                                                const matsubara_green_function_t &bare_green_matsubara, 
                                                                std::vector<double>& densities,
                                                                const double &beta, const int n_site, 
                                                                const int n_zone, const int n_matsubara) const
{
  double max_error = 0.;
  std::cout<<"evaluating self energy measurement."<<std::endl;
  matsubara_green_function_t Wk(n_matsubara, n_site, n_zone);
  Wk.clear();
  matsubara_green_function_t reduced_bare_green_matsubara(n_matsubara, n_site, n_zone);
  reduced_bare_green_matsubara.clear();
  for(int z=0;z<n_zone;++z){
    for (int k=0; k<n_site; k++) {                   
      std::stringstream Wk_real_name, Wk_imag_name;
      Wk_real_name  <<"Wk_real_"  <<z<<"_"<<k << "_" << k;
      Wk_imag_name  <<"Wk_imag_"  <<z<<"_"<<k << "_" << k;
      alps::RealVectorObsevaluator Weval_real=gathered_measurements[Wk_real_name.str().c_str()];
      alps::RealVectorObsevaluator Weval_imag=gathered_measurements[Wk_imag_name.str().c_str()];
      Weval_real /= beta*n_site;
      Weval_imag /= beta*n_site;
      std::valarray<double> mean_real = Weval_real.mean();
      std::valarray<double> mean_imag = Weval_imag.mean();
      for(int w=0;w<n_matsubara;++w)
        Wk(w, k, k, z) = std::complex<double>(mean_real[w], mean_imag[w]);
      for(int w=0;w<n_matsubara;++w)
        reduced_bare_green_matsubara(w, k, k, z) = bare_green_matsubara(w, k, k, z);
      std::valarray<double> error_real = Weval_real.error();
      std::valarray<double> error_imag = Weval_imag.error();
      for (int e=0; e<error_real.size(); ++e) {
        double ereal = error_real[e];
        double eimag = error_imag[e];
        double error = (ereal >= eimag) ? ereal : eimag;
        max_error = (error > max_error) ? error : max_error;
      }
    }
  }
  std::cout << "Maximal error in Wk: " << max_error << std::endl;
  green_matsubara_measured.clear();
  for(int z=0;z<n_zone;++z)
    for (int k=0; k<n_site; k++)                    
      for(int w=0;w<n_matsubara;++w)
        green_matsubara_measured(w,k,k, z) = bare_green_matsubara(w,k,k,z) 
          - bare_green_matsubara(w,k,k,z) * bare_green_matsubara(w,k,k,z) * Wk(w,k,k,z);
  std::valarray<double> dens = alps::RealVectorObsevaluator(gathered_measurements["densities"]).mean();
  for (int z=0; z<n_zone; ++z) 
    densities[z] = dens[z];
}




void InteractionExpansionSim::evaluate_selfenergy_measurement_itime_rs(const alps::ObservableSet &gathered_measurements, 
                                                               itime_green_function_t &green_result,
                                                               const itime_green_function_t &green0, 
                                                               const double &beta, const int n_site, 
                                                               const int n_zone, const int n_tau, const int n_self) const
{
  std::cout<<"evaluating self energy measurement: itime, real space."<<std::endl;
  clock_t time0=clock();
  std::vector<std::vector<std::vector<std::valarray<double> > > >W_z_i_j(n_zone); 
  //first index: zone. Second index: momentum. Third index: self energy tau point.
  double max_error = 0.;
  for(int z=0;z<n_zone;++z){
    W_z_i_j[z].resize(n_site);
    for(int i=0;i<n_site;++i){
      W_z_i_j[z][i].resize(n_site);
    }
    for(int i=0;i<n_site;++i){
      for(int j=0;j<n_site;++j){
        std::stringstream W_name;
        W_name<<"W_"<<z<<"_"<<i<<"_"<<j;
        alps::RealVectorObsevaluator W_eval  =gathered_measurements[W_name.  str().c_str()];
        W_z_i_j[z][i][j].resize(n_self+1);
        std::valarray<double> tmp=(W_eval.mean());
	std::valarray<double> errorvec = W_eval.error();
        for(int k=0;k<n_self+1;++k){
          W_z_i_j[z][i][j][k]=tmp[k];
	  double error = errorvec[k];
	  max_error = error>max_error ? error : max_error;
        }
      }
    }
    std::cout << "Maximum error in W: " << max_error << std::endl;
    for(int i=0;i<n_site;++i){
      for(int j=0;j<n_site;++j){
        for(int k=0;k<n_tau;++k){
          green_result(k,i,j,z)=green0(k,i,j,z);
          double timek=(k/(double)n_tau)*beta;
          for(int l=0;l<=n_self;++l){
            double timel;
            if(l==0){
              timel=(0.5/(double)n_self)*beta; //middle position of our first bin.
            }else if(l==n_self){
              timel=((n_self-0.5)/(double)n_self)*beta; //middle position of our last bin.
            }else{
              timel=(l/(double)n_self)*beta; //middle position of our remaining bins.
            }
            for(int p=0;p<n_site;++p){
              //this will make a tiny error if the bins are equal.
              green_result(k,i,j,z)-=green0_spline(green0, timek-timel,i,p,z, n_tau, beta)*W_z_i_j[z][p][j][l];
            }
          }
        }
        green_result(n_tau,i,j,z)=(i==j?-1:0)-green_result(0,i,j,z);
      }
    }
  }
  if(n_zone==2){  
    alps::RealObsevaluator Sz_obs = gathered_measurements["Sz_0"];
    for(int i=1;i<n_site;i++) {
      std::stringstream Sz_name;
      Sz_name << "Sz_" << i;
      alps::RealObsevaluator Sz_i_obs = gathered_measurements[Sz_name.str().c_str()];
      Sz_i_obs *= ((i%2==0) ? 1 : -1);
      Sz_obs += Sz_i_obs;
    }
    Sz_obs /= n_site;
    std::ofstream szstream("staggered_sz", std::ios::app);
    szstream << Sz_obs.mean() << "\t" << Sz_obs.error() << std::endl;
  }
  clock_t time1=clock();
  std::cout<<"evaluate of SE measurement took: "<<(time1-time0)/(double)CLOCKS_PER_SEC<<std::endl;
}



template<class X, class Y> inline Y linear_interpolate(const X x0, const X x1, const Y y0, const Y y1, const X x)
{
  return y0 + (x-x0)/(x1-x0)*(y1-y0);
}





double InteractionExpansionSim::green0_spline(const itime_green_function_t &green0, const itime_t delta_t, 
                                              const int s1, const int s2, const spin_t z, int n_tau, double beta) const
{
  double temperature=1./beta;
  double almost_zero=1.e-12;
  double n_tau_inv=1./n_tau;
  if(delta_t*delta_t < almost_zero){
    return green0(0,s1,s2,z);
  } else if(delta_t>0){
    int time_index_1 = (int)(delta_t*n_tau*temperature);
    int time_index_2 = time_index_1+1;
    return linear_interpolate((double)time_index_1*beta*n_tau_inv, (double)time_index_2*beta*n_tau_inv,green0(time_index_1,s1,s2,z),
      green0(time_index_2,s1,s2,z),delta_t);
  } else{
    int time_index_1 = (int)((beta+delta_t)*n_tau*temperature);
    int time_index_2 = time_index_1+1;
    return -linear_interpolate((double)time_index_1*beta*n_tau_inv, (double)time_index_2*beta*n_tau_inv, green0(time_index_1,s1,s2,z),
      green0(time_index_2, s1,s2,z), delta_t+beta);
  }
}

