/****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2012 by Emanuel Gull <gull@pks.mpg.de>,
 *                       Hartmut Hafermann <hafermann@cpht.polytechnique.fr>
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

  for(int i=0;i<size();++i){ c_exp      [i]=std::exp(std::complex<double>(0,  M_PI*c_times      [i]/beta_)); }
  for(int i=0;i<size();++i){ cdagger_exp[i]=std::exp(std::complex<double>(0, -M_PI*cdagger_times[i]/beta_)); }
  //measures the Fourier transform of G(tau-tau'):=-<T c(tau) c^dagger(tau')>
  for (int i = 0; i < size(); i++) {
    //note: strictly F_prefactor should be evaluated for cdagger_times for F(tau)
    //here we evaluate it instead for the annihilator times
    //this works only as long as F is diagonal (i.e. for diagonal hybridization)
    //see comments on get_F_prefactor in hyblocal.cpp
    double f_pref=(F_prefactor.find(c_times[i]))->second;
    for (int j = 0; j < size(); j++) {
      std::complex<double> M_ji = operator() (j, i) * sign;
      std::complex<double> exp=c_exp[i]*cdagger_exp[j];
      std::complex<double> dexp=exp*exp;
      for(std::size_t wn=0; wn<Gwr.size(); wn++){
        //std::complex<double> meas = -M_ji*std::exp(std::complex<double>(0,(2.*wn+1)*M_PI/beta_*(c_times[i]-cdagger_times[j])))/beta_;
        std::complex<double> meas = -M_ji*exp/beta_;
        Gwr[wn] += meas.real();
        Gwi[wn] += meas.imag();
        std::complex<double> fmeas = f_pref*meas;
        Fwr[wn] += fmeas.real();
        Fwi[wn] += fmeas.imag();
        exp*=dexp;
      }
    }
  }
}


void hybmatrix::measure_G2w(std::vector<std::complex<double> > &G2w, std::vector<std::complex<double> >&F2w, int N_w2, int N_w_aux, const std::map<double,double> &F_prefactor) const{
  static std::vector<double> cdagger_times(size()); cdagger_times.resize(size());
  static std::vector<double> c_times(size()); c_times.resize(size());
  static std::vector<std::complex<double> > c_exp_ini(size()); c_exp_ini.resize(size());
  static std::vector<std::complex<double> > c_exp_inc(size()); c_exp_inc.resize(size());
  static std::vector<std::complex<double> > cdagger_exp_ini(size()); cdagger_exp_ini.resize(size());
  static std::vector<std::complex<double> > cdagger_exp_inc(size()); cdagger_exp_inc.resize(size());

  //create map of creator and annihilator times
  for (hyb_map_t::const_iterator it= c_index_map_.begin(); it != c_index_map_.end(); ++it) {
    c_times[it->second] = it->first;
  }
  for (hyb_map_t::const_iterator it= cdagger_index_map_.begin(); it != cdagger_index_map_.end(); ++it) {
    cdagger_times[it->second] = it->first;
  }

  memset(&(G2w[0]),0, G2w.size()*sizeof(std::complex<double>));
  memset(&(F2w[0]),0, F2w.size()*sizeof(std::complex<double>));

  double w_ini = (2*(-N_w2/2)+1)*M_PI/beta_;
  double w_inc = 2*M_PI/beta_;
  for(int i=0;i<size();++i){ c_exp_ini      [i]=std::exp(std::complex<double>(0,  w_ini*c_times      [i])); }
  for(int i=0;i<size();++i){ c_exp_inc      [i]=std::exp(std::complex<double>(0,  w_inc*c_times      [i])); }
  for(int i=0;i<size();++i){ cdagger_exp_ini[i]=std::exp(std::complex<double>(0, -w_ini*cdagger_times[i])); }
  for(int i=0;i<size();++i){ cdagger_exp_inc[i]=std::exp(std::complex<double>(0, -w_inc*cdagger_times[i])); }

  //measures the Fourier transform of G(tau,tau'):=-<T c(tau) c^dagger(tau')>
  double f_pref=0.0;
  for (int i = 0; i < size(); i++) {
    if(measure_h2w_) f_pref=(F_prefactor.find(c_times[i]))->second;
    std::complex<double> exp1_ini=c_exp_ini[i];
    std::complex<double> exp1_inc=c_exp_inc[i];
    for (int j = 0; j < size(); j++) {
      std::complex<double> exp2_ini=cdagger_exp_ini[j];
      std::complex<double> exp2_inc=cdagger_exp_inc[j];

      std::complex<double> exp1=exp1_ini;
      std::complex<double> exp2=exp2_ini;

      std::complex<double> M_ji = operator() (j, i);

      for(int w1n=0; w1n<N_w_aux; w1n++){
        for(int w2n=0; w2n<N_w_aux; w2n++){
          if(measure_g2w_ || measure_h2w_){
            std::complex<double> meas = M_ji*exp1*exp2;
            G2w[w1n*N_w_aux+w2n] += meas;
          }
          if(measure_h2w_){
            std::complex<double> fmeas = f_pref*M_ji*exp1*exp2;
            F2w[w1n*N_w_aux+w2n] += fmeas;
          }
          exp2*=exp2_inc;//update exp2
        }
        exp2=exp2_ini;//reset exp2
        exp1*=exp1_inc;//update exp1
      }

    }
  }
}



