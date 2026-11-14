/****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2012 by Emanuel Gull <gull@pks.mpg.de>,
 *
 *  based on an earlier version by Philipp Werner and Emanuel Gull
 *
 *
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 *
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#include "hybmatrix.hpp"
#include "nfft3.h"

//compute the hybridization weight change when an operator pair is inserted
void hybmatrix::measure_Gw(std::vector<double> &Gwr, std::vector<double> &Gwi , std::vector<double> &Fwr, std::vector<double> &Fwi , const std::map<double,double> &F_prefactor, double sign) const{
  static std::vector<double> cdagger_times(size()); cdagger_times.resize(size());
  static std::vector<double> c_times(size()); c_times.resize(size());
  static std::vector<std::complex<double> > cdagger_exp(size()); cdagger_exp.resize(size());
  static std::vector<std::complex<double> > c_exp(size()); c_exp.resize(size());

  //create map of creator and annihilator times
  for (hyb_map_t::const_iterator it= c_index_map_.begin(); it != c_index_map_.end(); ++it) {
    c_times[it->second] = it->first;
  }
  for (hyb_map_t::const_iterator it= cdagger_index_map_.begin(); it != cdagger_index_map_.end(); ++it) {
    cdagger_times[it->second] = it->first;
  }
  int size2=size()*size();
  int n_omega_meas=Gwr.size();
  //nfft calculation
  static std::vector<nfft_plan> plan_vector(size()+1);
  static std::vector<bool> plan_initialized(size()+1, false);
  if(size()>=(int)(plan_vector.size())){
    plan_vector.resize(size()+1);
    plan_initialized.resize(size()+1, false);
    //std::cout<<"resizing up to: "<<size()<<std::endl;
  }

  // init a one dimensional plan
  {
    std::size_t nfft_flags_later = PRE_PHI_HUT| PRE_PSI| FFT_OUT_OF_PLACE;
    std::size_t fftw_flags= FFTW_MEASURE| FFTW_DESTROY_INPUT;
    static int n[1];
    static int m;
    if(!plan_initialized[size()]){
      //simple initialization
      nfft_init_1d(&(plan_vector[size()]),2*n_omega_meas, size2);
      //make sure we don't release memory and the plan, reset by plan w/o MALLOC flags set
      plan_vector[size()].fftw_flags=fftw_flags;
      plan_vector[size()].flags=nfft_flags_later;
      m=plan_vector[size()].m;
      n[0]=plan_vector[size()].n[0];
      plan_initialized[size()]=true;
      //std::cout<<"initializing: "<<size()<<std::endl;
    }else{
      int N[1];
      N[0]=2*n_omega_meas;
      int M_total=size2;
      nfft_init_guru(&(plan_vector[size()]), 1, N, M_total, n, m,  nfft_flags_later, fftw_flags);
    }
  }

  /** init the nonequidistant times*/
  for(int i=0; i<size();++i){
    for(int j=0; j<size();++j){
      (plan_vector[size()]).x[i*size()+j]=(c_times[i]-cdagger_times[j])/(2.*beta_);
    }
  }

  /** precompute psi, the entries of the matrix B */
  if((plan_vector[size()]).flags & PRE_ONE_PSI){ //(PRE_LIN_PSI| PRE_FG_PSI| PRE_PSI| PRE_FULL_PSI)
    nfft_precompute_one_psi(&(plan_vector[size()]));
  }

  //factorized exp factor: exp(c_times[p]*M_PI/beta)
  //                       exp(2*M_PI*n_omega_meas*c_times[p]/(2 beta))
  //consolidate:           exp(c_times[p]*(1+n_omega_meas)*M_PI/beta)
  for(int i=0;i<size();++i){ c_exp      [i]=std::exp(std::complex<double>(0., (c_times[i]*(1+n_omega_meas)*M_PI/beta_))); }
  for(int i=0;i<size();++i){ cdagger_exp[i]=std::exp(std::complex<double>(0., -(cdagger_times[i]*(1+n_omega_meas)*M_PI/beta_))); }

  //these are O(N^2) exps. Previously we had O(N*Nmatsubara) exps. Can we cut
  //some of these and replace them by O(N) exps? -> probably yes!
  static std::vector<std::complex<double> > coeff_f; coeff_f.resize(size2);
  for(int p=0;p<size();++p){
    double f_pref=(F_prefactor.find(c_times[p]))->second;
    for(int q=0;q<size();++q){
      //std::complex<double> times_exp=std::exp(std::complex<double>(0., (c_times[p]-cdagger_times[q])*M_PI/beta_+2*M_PI*n_omega_meas*nfft.x[p*size()+q]));
      std::complex<double> times_exp=c_exp[p]*cdagger_exp[q];
      std::complex<double> M_ji = operator() (q, p) * sign;
      std::complex<double> coeff=M_ji*times_exp;
      (plan_vector[size()]).f[p*size()+q][0]=coeff.real();
      (plan_vector[size()]).f[p*size()+q][1]=coeff.imag();
      coeff_f[p*size()+q]=coeff*f_pref;
    }
  }

  // approx. adjoint and show the result
  nfft_adjoint(&(plan_vector[size()]));
  for(int wn=0;wn<n_omega_meas;wn++){
    Gwr[wn]+=-(plan_vector[size()]).f_hat[2*wn][0]/beta_;
    Gwi[wn]+=-(plan_vector[size()]).f_hat[2*wn][1]/beta_;
  }

  for(int p=0;p<size();++p){
    for(int q=0;q<size();++q){
      (plan_vector[size()]).f[p*size()+q][0]=coeff_f[p*size()+q].real();
      (plan_vector[size()]).f[p*size()+q][1]=coeff_f[p*size()+q].imag();
    }
  }
  nfft_adjoint(&(plan_vector[size()]));
  for(int wn=0;wn<n_omega_meas;wn++){
    Fwr[wn]+=-(plan_vector[size()]).f_hat[2*wn][0]/beta_;
    Fwi[wn]+=-(plan_vector[size()]).f_hat[2*wn][1]/beta_;
  }
  // finalise the one dimensional plan
  nfft_finalize(&plan_vector[size()]); //store plan for later!
  //std::cout<<"exiting measurement."<<std::endl;

}

void hybmatrix::measure_G2w(std::vector<std::complex<double> > &G2w, std::vector<std::complex<double> >&F2w, int N_w2, int N_w_aux, const std::map<double,double> &F_prefactor) const{
  static std::vector<double> cdagger_times(size()); cdagger_times.resize(size());
  static std::vector<double> c_times(size()); c_times.resize(size());
  static std::vector<std::complex<double> > cdagger_exp(size()); cdagger_exp.resize(size());
  static std::vector<std::complex<double> > c_exp(size()); c_exp.resize(size());

  //create map of creator and annihilator times
  for (hyb_map_t::const_iterator it= c_index_map_.begin(); it != c_index_map_.end(); ++it) {
    c_times[it->second] = it->first;
  }
  for (hyb_map_t::const_iterator it= cdagger_index_map_.begin(); it != cdagger_index_map_.end(); ++it) {
    cdagger_times[it->second] = it->first;
  }

  int n_omega_meas=N_w_aux;//this is the effective number of fermionic frequencies needed

  //nfft calculation
  nfft_plan nfft_c;
  nfft_plan nfft_cdagger;

  // init a one dimensional plan
  nfft_init_1d(&nfft_c,n_omega_meas*2, size());
  nfft_init_1d(&nfft_cdagger,n_omega_meas*2, size());

  // init the nonequidistant times
  for(int i=0; i<size();++i){
    nfft_c      .x[i]= c_times[i]      /beta_-0.5;
    nfft_cdagger.x[i]=-cdagger_times[i]/beta_+0.5;
  }

  // precompute psi, the entries of the matrix B
  if(nfft_c.flags & PRE_ONE_PSI){
    nfft_precompute_one_psi(&nfft_c      );
    nfft_precompute_one_psi(&nfft_cdagger);
  }

  for(int i=0;i<size();++i){ c_exp      [i]=std::exp(std::complex<double>(0., (c_times[i]*M_PI/beta_))); }
  for(int i=0;i<size();++i){ cdagger_exp[i]=std::exp(std::complex<double>(0., -(cdagger_times[i]*M_PI/beta_))); }

  std::vector<std::complex<double> > M_tilde_im(size()*n_omega_meas, 0.);
  std::vector<std::complex<double> > coeff_G_ij(size()*size(), 0.);
  std::vector<std::complex<double> > coeff_F_ij(size()*size(), 0.);
  std::vector<std::complex<double> > coeff_im(size()*n_omega_meas, 0.);

  std::vector<int> freq(n_omega_meas); for(int i=0;i<n_omega_meas;++i) freq[i]=i-n_omega_meas;

  //first preproc step
  for(int i=0;i<size();++i){
    double f_pref=(F_prefactor.find(c_times[i]))->second;
    for(int j=0;j<size();++j){
      coeff_G_ij[i*size()+j]=operator() (j, i)*cdagger_exp[j];
      coeff_F_ij[i*size()+j]=operator() (j, i)*cdagger_exp[j]*f_pref;
    }
  }
  if(measure_g2w_){
    //first ft step as nfft, for G
    for(int i=0;i<size();++i){
      for(int j=0;j<size();++j){
        nfft_cdagger.f[j][0]=coeff_G_ij[i*size()+j].real();
        nfft_cdagger.f[j][1]=coeff_G_ij[i*size()+j].imag();
      }
      nfft_adjoint(&nfft_cdagger);
      for(int m=0;m<n_omega_meas;++m){
        std::complex<double> prefactor=std::exp(std::complex<double>(0.,-M_PI*freq[m]));
        M_tilde_im[i*n_omega_meas+m]=std::complex<double>(nfft_cdagger.f_hat[m+n_omega_meas-N_w2/2][0], nfft_cdagger.f_hat[m+n_omega_meas-N_w2/2][1])*prefactor;
      }
    }

    //second preproc step for G
    for(int i=0;i<size();++i){
      for(int m=0;m<n_omega_meas;++m){
        coeff_im[i*n_omega_meas+m]=M_tilde_im[i*n_omega_meas+m]*c_exp[i];
      }
    }

    //second ft step as nfft for G
    for(int m=0;m<n_omega_meas;++m){
      for(int i=0;i<size();++i){
        nfft_c.f[i][0]=coeff_im[i*n_omega_meas+m].real();
        nfft_c.f[i][1]=coeff_im[i*n_omega_meas+m].imag();
      }
      nfft_adjoint(&nfft_c);
      for(int n=0;n<n_omega_meas;++n){
        std::complex<double> prefactor=std::exp(std::complex<double>(0., M_PI*freq[n]));
        G2w[n*n_omega_meas+m]=std::complex<double>(nfft_c.f_hat[n+n_omega_meas-N_w2/2][0], nfft_c.f_hat[n+n_omega_meas-N_w2/2][1])*prefactor;
      }
    }
  }
  if(measure_h2w_){
    //first ft step as nfft, for F
    for(int i=0;i<size();++i){
      for(int j=0;j<size();++j){
        nfft_cdagger.f[j][0]=coeff_F_ij[i*size()+j].real();
        nfft_cdagger.f[j][1]=coeff_F_ij[i*size()+j].imag();
      }
      nfft_adjoint(&nfft_cdagger);
      for(int m=0;m<n_omega_meas;++m){
        std::complex<double> prefactor=std::exp(std::complex<double>(0.,-M_PI*freq[m]));
        M_tilde_im[i*n_omega_meas+m]=std::complex<double>(nfft_cdagger.f_hat[m+n_omega_meas-N_w2/2][0], nfft_cdagger.f_hat[m+n_omega_meas-N_w2/2][1])*prefactor;
      }
    }
    //second preproc step for F
    for(int i=0;i<size();++i){
      for(int m=0;m<n_omega_meas;++m){
        coeff_im[i*n_omega_meas+m]=M_tilde_im[i*n_omega_meas+m]*c_exp[i];
      }
    }

    //second ft step as nfft for F
    for(int m=0;m<n_omega_meas;++m){
      for(int i=0;i<size();++i){
        nfft_c.f[i][0]=coeff_im[i*n_omega_meas+m].real();
        nfft_c.f[i][1]=coeff_im[i*n_omega_meas+m].imag();
      }
      nfft_adjoint(&nfft_c);
      for(int n=0;n<n_omega_meas;++n){
        std::complex<double> prefactor=std::exp(std::complex<double>(0., M_PI*freq[n]));
        F2w[n*n_omega_meas+m]=std::complex<double>(nfft_c.f_hat[n+n_omega_meas-N_w2/2][0], nfft_c.f_hat[n+n_omega_meas-N_w2/2][1])*prefactor;
      }
    }
  }
  nfft_finalize(&nfft_c);
  nfft_finalize(&nfft_cdagger);
}


//correct, for debugging.
/*void hybmatrix::measure_G2w_nfft(std::vector<std::complex<double> > &G2w, std::vector<std::complex<double> >&F2w, int N_w2, int N_w_aux, const std::map<double,double> &F_prefactor) const{
 std::cout<<"entering measure G2w."<<std::endl;
 if(N_w2 !=N_w_aux) throw std::logic_error("understand why aux and w2 can be different!");
 if(N_w2 %2!=0) throw std::logic_error("understand why  w2 can be odd!");
 clock_t t_nfft_start=clock();
 std::vector<double> cdagger_times(size());
 std::vector<double> c_times(size());

 //create map of creator and annihilator times
 for (hyb_map_t::const_iterator it= c_index_map_.begin(); it != c_index_map_.end(); ++it) {
 c_times[it->second] = it->first;
 }
 for (hyb_map_t::const_iterator it= cdagger_index_map_.begin(); it != cdagger_index_map_.end(); ++it) {
 cdagger_times[it->second] = it->first;
 }

 int size2=size()*size();
 int n_omega_meas=N_w2;

 //std::cout<<"N_w2: "<<N_w2<<" aux: "<<N_w_aux<<std::endl;
 //nfft calculation
 nfft_plan nfft_c;
 nfft_plan nfft_cdagger;

 // init a one dimensional plan
 nfft_init_1d(&nfft_c,n_omega_meas, size());
 nfft_init_1d(&nfft_cdagger,n_omega_meas, size());

 // init the nonequidistant times
 for(int i=0; i<size();++i){
 nfft_c      .x[i]= c_times[i]      /beta_-0.5;
 nfft_cdagger.x[i]=-cdagger_times[i]/beta_+0.5;
 }

 // precompute psi, the entries of the matrix B
 if(nfft_c.flags & PRE_ONE_PSI){
 nfft_precompute_one_psi(&nfft_c      );
 nfft_precompute_one_psi(&nfft_cdagger);
 }

 //for(int i=0;i<size();++i){ c_exp      [i]=std::exp(std::complex<double>(0., (c_times[i]*M_PI/beta_))); }
 //for(int i=0;i<size();++i){ cdagger_exp[i]=std::exp(std::complex<double>(0., -(cdagger_times[i]*M_PI/beta_))); }


 std::vector<std::complex<double> > M_tilde_im(size()*n_omega_meas, 0.);
 std::vector<std::complex<double> > coeff_ij(size()*size(), 0.);
 std::vector<std::complex<double> > coeff_im(size()*n_omega_meas, 0.);
 std::vector<std::complex<double> > M_tilde_nm(n_omega_meas*n_omega_meas, 0.);
 std::vector<int> freq(n_omega_meas); for(int i=0;i<n_omega_meas;++i) freq[i]=i-n_omega_meas/2;

 //first preproc step
 for(int i=0;i<size();++i){
 for(int j=0;j<size();++j){
 coeff_ij[i*size()+j]=operator() (j, i)*std::exp(std::complex<double>(0, -M_PI/beta_*cdagger_times[j]));
 }
 }

 //first ft step
 for(int i=0;i<size();++i){
 for(int m=0;m<n_omega_meas;++m){
 for(int j=0;j<size();++j){
 M_tilde_im[i*n_omega_meas+m]+=coeff_ij[i*size()+j]*std::exp(std::complex<double>(0, -(2.*M_PI*freq[m])/beta_*cdagger_times[j]));
 //M_tilde_im[i*n_omega_meas+m]+=coeff_ij[i*size()+j]*std::exp(std::complex<double>(0, (2.*M_PI*freq[m])*(-cdagger_times[j]/beta_+0.5)));
 }
 }
 }
 //first ft step as nfft
 for(int i=0;i<size();++i){
 for(int j=0;j<size();++j){
 nfft_cdagger.f[j][0]=coeff_ij[i*size()+j].real();
 nfft_cdagger.f[j][1]=coeff_ij[i*size()+j].imag();
 }
 nfft_adjoint(&nfft_cdagger);
 //std::cout<<i<<std::endl;
 for(int m=0;m<n_omega_meas;++m){
 std::complex<double> prefactor=std::exp(std::complex<double>(0.,-M_PI*freq[m]));
 std::complex<double> val=std::complex<double>(nfft_cdagger.f_hat[m][0], nfft_cdagger.f_hat[m][1])*prefactor;
 //std::cout<<m<<" "<<M_tilde_im[i*n_omega_meas+m]<<" "<<val<<std::endl;
 }
 }
 //second preproc step
 for(int i=0;i<size();++i){
 for(int m=0;m<n_omega_meas;++m){
 coeff_im[i*n_omega_meas+m]=M_tilde_im[i*n_omega_meas+m]*std::exp(std::complex<double>(0, M_PI/beta_*c_times[i]));
 }
 }

 //second ft step
 for(int n=0;n<n_omega_meas;++n){
 for(int m=0;m<n_omega_meas;++m){
 for(int i=0;i<size();++i){
 M_tilde_nm[n*n_omega_meas+m]+=coeff_im[i*n_omega_meas+m]*std::exp(std::complex<double>(0, (2.*M_PI*freq[n])/beta_*c_times[i]));
 }
 }
 }
 //second ft step as nfft
 //first ft step as nfft
 for(int m=0;m<n_omega_meas;++m){
 for(int i=0;i<size();++i){
 nfft_c.f[i][0]=coeff_im[i*n_omega_meas+m].real();
 nfft_c.f[i][1]=coeff_im[i*n_omega_meas+m].imag();
 }
 nfft_adjoint(&nfft_c);
 //std::cout<<i<<std::endl;
 for(int n=0;n<n_omega_meas;++n){
 std::complex<double> prefactor=std::exp(std::complex<double>(0., M_PI*freq[n]));
 std::complex<double> val=std::complex<double>(nfft_c.f_hat[n][0], nfft_c.f_hat[n][1])*prefactor;
 //std::cout<<n<<" "<<M_tilde_nm[n*n_omega_meas+m]<<" "<<val<<std::endl;
 M_tilde_nm[n*n_omega_meas+m]=val;
 }
 }


 std::cout<<"result is: "<<M_tilde_nm[0*N_w_aux+0]<<" "<<M_tilde_nm[0*N_w_aux+1]<<" "<<M_tilde_nm[1*N_w_aux+0]<<" "<<M_tilde_nm[1*N_w_aux+1]<<" "<<std::endl;

 }*/
