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

#include"hyb.hpp"
#include"alps/numeric/vector_functions.hpp"

using namespace alps::numeric;

#ifdef ALPS_NGS_USE_NEW_ALEA
  using alps::accumulator::max_bin_number;
  #define NUM_BINS_CONSTRUCTOR_ARG max_bin_number=NUM_BINS
#else
  #define NUM_BINS_CONSTRUCTOR_ARG NUM_BINS
#endif

typedef alps::accumulator::RealVectorObservable vec_obs_t;
typedef alps::accumulator::RealObservable obs_t;

void hybridization::create_measurements(){//called once in the constructor

  //basic measurements for all orbitals
  //  std::cerr << "NUM_BINS = " << NUM_BINS << std::endl;
  measurements<< vec_obs_t("order_histogram_total",NUM_BINS_CONSTRUCTOR_ARG);
  measurements<< vec_obs_t("sector_statistics",NUM_BINS_CONSTRUCTOR_ARG);
  measurements<< obs_t("Sign",NUM_BINS_CONSTRUCTOR_ARG);
  
  g2wr_names.resize(n_orbitals); g2wi_names.resize(n_orbitals);
  h2wr_names.resize(n_orbitals); h2wi_names.resize(n_orbitals);
  g2wr.resize(n_orbitals); g2wi.resize(n_orbitals);
  h2wr.resize(n_orbitals); h2wi.resize(n_orbitals);

  G2w.resize(n_orbitals, std::vector<std::complex<double> >()); //resized conditionally
  F2w.resize(n_orbitals, std::vector<std::complex<double> >()); //below

  if(MEASURE_g2w){
    g2wr.resize(N_w2*N_w2*N_W, 0.);
    g2wi.resize(N_w2*N_w2*N_W, 0.);
  }
  if(MEASURE_h2w){
    h2wr.resize(N_w2*N_w2*N_W, 0.);
    h2wi.resize(N_w2*N_w2*N_W, 0.);
  }

  nnt_names.resize(n_orbitals);
  nnw_re_names.resize(n_orbitals);
  nn_names.resize(n_orbitals);
  nnt.resize(n_orbitals);
  nnw_re.resize(n_orbitals);
  nn.resize(n_orbitals);

  //additional measurements, per orbital
  for(std::size_t i=0;i<n_orbitals;++i){
    //g in tau
    std::stringstream g_name; g_name<<"g_"<<i; g_names.push_back(g_name.str());
    //f in tau
    std::stringstream f_name; f_name<<"f_"<<i; f_names.push_back(f_name.str());

    //density
    std::stringstream density_name; density_name<<"density_"<<i; density_names.push_back(density_name.str());

    //histogram
    std::stringstream order_name; order_name<<"order_"<<i; order_names.push_back(order_name.str());
    std::stringstream order_histogram_name; order_histogram_name<<"order_histogram_"<<i; order_histogram_names.push_back(order_histogram_name.str());

    //g in matsubara, f (for improved estimator) in matsubara
    std::stringstream gwr_name; gwr_name<<"gw_re_"<<i; gwr_names.push_back(gwr_name.str());
    std::stringstream gwi_name; gwi_name<<"gw_im_"<<i; gwi_names.push_back(gwi_name.str());
    std::stringstream fwr_name; fwr_name<<"fw_re_"<<i; fwr_names.push_back(fwr_name.str());
    std::stringstream fwi_name; fwi_name<<"fw_im_"<<i; fwi_names.push_back(fwi_name.str());

    //Legendre coefficients for g and f
    std::stringstream gl_name; gl_name<<"gl_"<<i; gl_names.push_back(gl_name.str());
    std::stringstream fl_name; fl_name<<"fl_"<<i; fl_names.push_back(fl_name.str());

    //nnt and nnw
    nnt[i].resize(i+1);
    nnw_re[i].resize(i+1);
    nn[i].resize(i);

    if(MEASURE_g2w || MEASURE_h2w) G2w[i].resize(N_w_aux*N_w_aux);
    if(MEASURE_h2w) F2w[i].resize(N_w_aux*N_w_aux);

    //initialize measurements for observable names
    measurements << vec_obs_t(g_name.str(),NUM_BINS_CONSTRUCTOR_ARG);
    measurements << vec_obs_t(f_name.str(),NUM_BINS_CONSTRUCTOR_ARG);
    measurements << obs_t(density_name.str(),NUM_BINS_CONSTRUCTOR_ARG);

    measurements << vec_obs_t(order_histogram_name.str(),NUM_BINS_CONSTRUCTOR_ARG);
    measurements << obs_t(order_name.str(),NUM_BINS_CONSTRUCTOR_ARG);

    measurements << vec_obs_t(gwr_name.str(),NUM_BINS_CONSTRUCTOR_ARG);
    measurements << vec_obs_t(gwi_name.str(),NUM_BINS_CONSTRUCTOR_ARG);
    measurements << vec_obs_t(fwr_name.str(),NUM_BINS_CONSTRUCTOR_ARG);
    measurements << vec_obs_t(fwi_name.str(),NUM_BINS_CONSTRUCTOR_ARG);

    measurements << vec_obs_t(gl_name.str(),NUM_BINS_CONSTRUCTOR_ARG);
    measurements << vec_obs_t(fl_name.str(),NUM_BINS_CONSTRUCTOR_ARG);

    if(MEASURE_nn){
      for(std::size_t j=0;j<i;++j){//j<i not j<=i
        std::stringstream nn_name; nn_name<<"nn_"<<i<<"_"<<j; nn_names[i].push_back(nn_name.str());
        nn[i][j]=0.;
        measurements << obs_t(nn_name.str(),NUM_BINS_CONSTRUCTOR_ARG);
      }
    }
    for(std::size_t j=0;j<=i;++j){//two-particle quantities
      if(MEASURE_nnt){
        std::stringstream nnt_name; nnt_name<<"nnt_"<<i<<"_"<<j; nnt_names[i].push_back(nnt_name.str());
        nnt[i][j].resize(N_nn+1, 0.);
        measurements << vec_obs_t(nnt_name.str(),NUM_BINS_CONSTRUCTOR_ARG);
      }
      if(MEASURE_nnw){
        std::stringstream nnw_re_name; nnw_re_name<<"nnw_re_"<<i<<"_"<<j; nnw_re_names[i].push_back(nnw_re_name.str());
        nnw_re[i][j].resize(N_W, 0.);
        measurements << vec_obs_t(nnw_re_name.str(),NUM_BINS_CONSTRUCTOR_ARG);
      }
      if(MEASURE_g2w){    //the two-particle Green's function is large and error is usually not needed -> declare as simple observable
        std::stringstream g2wr_name; g2wr_name<<"g2w_re_"<<i<<"_"<<j; g2wr_names[i].push_back(g2wr_name.str());
        std::stringstream g2wi_name; g2wi_name<<"g2w_im_"<<i<<"_"<<j; g2wi_names[i].push_back(g2wi_name.str());
#ifdef ALPS_NGS_USE_NEW_ALEA
        measurements << alps::accumulator::SimpleRealVectorObservable(g2wr_name.str());
        measurements << alps::accumulator::SimpleRealVectorObservable(g2wi_name.str());
#else
        measurements << alps::ngs::SimpleRealVectorObservable(g2wr_name.str());
        measurements << alps::ngs::SimpleRealVectorObservable(g2wi_name.str());
#endif
      }
      if(MEASURE_h2w){
        std::stringstream h2wr_name; h2wr_name<<"h2w_re_"<<i<<"_"<<j; h2wr_names[i].push_back(h2wr_name.str());
        std::stringstream h2wi_name; h2wi_name<<"h2w_im_"<<i<<"_"<<j; h2wi_names[i].push_back(h2wi_name.str());
#ifdef ALPS_NGS_USE_NEW_ALEA
        measurements << alps::accumulator::SimpleRealVectorObservable(h2wr_name.str());
        measurements << alps::accumulator::SimpleRealVectorObservable(h2wi_name.str());
#else
        measurements << alps::ngs::SimpleRealVectorObservable(h2wr_name.str());
        measurements << alps::ngs::SimpleRealVectorObservable(h2wi_name.str());
#endif
      }
    }
  }
  measurements.reset(true);

  //initialize measurement vectors
  sgn=0.;
  order_histogram.resize(n_orbitals, std::vector<double> (N_hist_orders, 0.));
  orders.resize(n_orbitals, 0.);
  order_histogram_total.resize(N_hist_orders, 0.);
  densities.resize(n_orbitals, 0.);
  sector_statistics.resize(1<<n_orbitals, 0.);

  if(MEASURE_time){
    G.resize(n_orbitals, std::vector<double>(N_t+1, 0.));
    F.resize(n_orbitals, std::vector<double>(N_t+1, 0.));
  }
  if(MEASURE_freq){
    Gwr.resize(n_orbitals, std::vector<double>(N_w, 0.));
    Gwi.resize(n_orbitals, std::vector<double>(N_w, 0.));
    Fwr.resize(n_orbitals, std::vector<double>(N_w, 0.));
    Fwi.resize(n_orbitals, std::vector<double>(N_w, 0.));
  }
  if(MEASURE_legendre){
    Gl.resize(n_orbitals, std::vector<double>(N_l, 0.));
    Fl.resize(n_orbitals, std::vector<double>(N_l, 0.));
  }
  if(MEASURE_nnt) n_vectors.resize(n_orbitals, std::vector<double>(N_nn+1, 0.));

  F_prefactor.resize(n_orbitals);

}// create measurements


void hybridization::measure(){
  if(!is_thermalized()) return;

  accumulate_order();

  accumulate_G();

  if(!MEASURE_time && (MEASURE_freq || MEASURE_legendre || MEASURE_h2w))//F_prefactor is computed in update() if time measurement is turned on
    local_config.get_F_prefactor(F_prefactor);//compute segment overlaps in local config

  measure_Gw(F_prefactor);
  accumulate_Gw();

  measure_Gl(F_prefactor);
  accumulate_Gl();

  measure_sector_statistics();
  accumulate_sector_statistics();

  //measure 2-particle quantities
  measure_nn();
  accumulate_nn();

  measure_nnt();
  accumulate_nnt();

  measure_nnw();
  accumulate_nnw();

  if(MEASURE_g2w  || MEASURE_h2w) measure_G2w(F_prefactor); //accumulated during measurement to save memory

  sweep_count = sweeps;
}

void hybridization::measure_order(){
  //compute the order and store it in the vectors for the histograms
  sgn+=sign;
  local_config.measure_density(densities, sign);
  for(std::size_t i=0;i<n_orbitals;++i){
    double order=local_config.order(i);
    orders[i]+=order;
    if(order<order_histogram[i].size()){
      order_histogram[i][order]++;
      order_histogram_total[order]++;
    }
  }
}

void hybridization::accumulate_order(){
  measurements["order_histogram_total"]<<(order_histogram_total/N_meas);
  memset(&(order_histogram_total[0]), 0, sizeof(double)*order_histogram_total.size());
  measurements["Sign"]<<(sgn/N_meas); sgn=0.;
  for(std::size_t i=0;i<n_orbitals;++i){
    measurements[order_names[i]]<<(orders[i]/N_meas);
    measurements[density_names[i]]<<(densities[i]/N_meas);
    measurements[order_histogram_names[i]]<<(order_histogram[i]/N_meas);
    orders[i]=0.;
    densities[i]=0;
    memset(&(order_histogram[i][0]), 0, sizeof(double)*order_histogram[i].size());
  }
}

//measure the Green's function
void hybridization::measure_G(std::vector<std::map<double,double> > &F_prefactor){
  if(!MEASURE_time) return;
  //delegate the actual measurement to the hybridization configuration
  hyb_config.measure_G(G, F, F_prefactor, sign);
}

void hybridization::accumulate_G(){
  if(!MEASURE_time) return;
  for(std::size_t i=0;i<n_orbitals;++i){
    measurements[g_names[i]]<<(N_t*G[i]/(beta*beta*N_meas));
    measurements[f_names[i]]<<(N_t*F[i]/(beta*beta*N_meas));
    memset(&(G[i][0]), 0, sizeof(double)*G[i].size());
    memset(&(F[i][0]), 0, sizeof(double)*F[i].size());
  }
}

void hybridization::measure_sector_statistics(){
  if(!MEASURE_sector_statistics) return;
  local_config.measure_sector_statistics(sector_statistics, sign);
}

void hybridization::accumulate_sector_statistics(){
  if(!MEASURE_sector_statistics) return;
  measurements["sector_statistics"]<<sector_statistics;
  memset(&(sector_statistics[0]),0, sector_statistics.size()*sizeof(double));
}

void hybridization::measure_nn(){
  if(!MEASURE_nn) return;
  for(std::size_t i=0;i<n_orbitals;++i)
    for(std::size_t j=0;j<i;++j){//i==j would simply yield the density, which we measure separately
      nn[i][j]+=local_config.measure_nn(i,j)*sign;
    }
}

void hybridization::accumulate_nn(){
  if(!MEASURE_nn) return;
  for(std::size_t i=0;i<n_orbitals;++i)
    for(std::size_t j=0;j<i;++j){//i==j would simply yield the density, which we measure separately
      measurements[nn_names[i][j]]<<nn[i][j];
      nn[i][j]=0;
    }
}

void hybridization::measure_nnt(){
  if(!MEASURE_nnt) return;
  local_config.get_density_vectors(n_vectors);

  for(std::size_t i=0;i<n_orbitals;++i){
    for(std::size_t j=0;j<=i;++j){
      if(n_vectors[j][0]>0){//uses time translational invariance, i.e. n_j fixed at tau=0; looping over all entries of n_j is costly and measurements are correlated
//        for(int n=0;n<=N_nn;++n) nnt[n]=n_vectors[i][0]*n_vectors[j][n]*sign;
        for(std::size_t n=0; n<=N_nn; ++n) nnt[i][j][n]+=n_vectors[i][n]*sign;
      }
    }
  }
}

void hybridization::accumulate_nnt(){
  if(!MEASURE_nnt) return;
  for(std::size_t i=0;i<n_orbitals;++i){
    for(std::size_t j=0;j<=i;++j){
      measurements[nnt_names[i][j]]<<nnt[i][j];
      memset(&(nnt[i][j][0]),0, nnt[i][j].size()*sizeof(double));
    }
  }
}

void hybridization::measure_nnw(){
  if(!MEASURE_nnw) return;
  for(std::size_t i=0;i<n_orbitals;++i){
    for(std::size_t j=0;j<=i;++j){
      if(local_config.density(j,0.0)>0.){
        local_config.measure_nnw(i,nnw_re[i][j], sign);
      }
    }
  }
}

void hybridization::accumulate_nnw(){
  if(!MEASURE_nnw) return;
  for(std::size_t i=0;i<n_orbitals;++i){
    for(std::size_t j=0;j<=i;++j){
      measurements[nnw_re_names[i][j]]<<nnw_re[i][j];
      memset(&(nnw_re[i][j][0]),0, nnw_re[i][j].size()*sizeof(double));
    }
  }
}

void hybridization::measure_Gw(std::vector<std::map<double,double> > &F_prefactor){
  if(!MEASURE_freq) return;

  //compute frequency quantities
  hyb_config.measure_Gw(Gwr,Gwi, Fwr, Fwi, F_prefactor, sign);
}

void hybridization::accumulate_Gw(){
  if(!MEASURE_freq) return;
    for(std::size_t i=0;i<n_orbitals;++i){
      measurements[gwr_names[i]]<<Gwr[i];
      measurements[gwi_names[i]]<<Gwi[i];
      measurements[fwr_names[i]]<<Fwr[i];
      measurements[fwi_names[i]]<<Fwi[i];
      memset(&(Gwr[i][0]),0, Gwr[i].size()*sizeof(double));
      memset(&(Gwi[i][0]),0, Gwr[i].size()*sizeof(double));
      memset(&(Fwr[i][0]),0, Gwr[i].size()*sizeof(double));
      memset(&(Fwi[i][0]),0, Gwr[i].size()*sizeof(double));
    }
}

void hybridization::measure_Gl(std::vector<std::map<double,double> > &F_prefactor){
  if(!MEASURE_legendre) return;
  //compute legendre quantities
  hyb_config.measure_Gl(Gl, Fl, F_prefactor, sign);
}

void hybridization::accumulate_Gl(){
  if(!MEASURE_legendre) return;
    for(std::size_t i=0;i<n_orbitals;++i){
      measurements[gl_names[i]]<<Gl[i];
      measurements[fl_names[i]]<<Fl[i];
      memset(&(Gl[i][0]),0, Gl[i].size()*sizeof(double));
      memset(&(Fl[i][0]),0, Fl[i].size()*sizeof(double));
    }
}

void hybridization::measure_G2w(std::vector<std::map<double,double> > &F_prefactor){

  //compute two-frequency quantities
  hyb_config.measure_G2w(G2w, F2w, N_w2, N_w_aux, F_prefactor);

  for(std::size_t i=0;i<n_orbitals;++i){
    for(std::size_t j=0;j<=i;++j){//we measure only for j<=i since results for ij and ji are exactly the same (no gain through averaging)
      for(std::size_t w2n=0; w2n<N_w2;++w2n)
        for(std::size_t w3n=0; w3n<N_w2;++w3n)
          for(std::size_t Wn=0; Wn<N_W;++Wn){
            int w1n=w2n+Wn; int w4n=w3n+Wn;
            int index=Wn*N_w2*N_w2 + w2n*N_w2 + w3n;
            if(MEASURE_g2w){
              std::complex<double> meas =G2w[i][w1n*N_w_aux+w2n]*G2w[j][w3n*N_w_aux+w4n]; // M12M34
              if(i==j)             meas-=G2w[i][w1n*N_w_aux+w4n]*G2w[i][w3n*N_w_aux+w2n]; //-M14M32
              meas/=beta; meas*=sign;
              g2wr[index] += meas.real();
              g2wi[index] += meas.imag();
            }
            if(MEASURE_h2w){
              //std::cout<<i<<" "<<j<<" "<<w2n<<" "<<w3n<<" "<<Wn<<" "<<std::endl;
              //std::cout<<"G: "<<G2w[j][w3n*N_w_aux+w4n]<<" size: "<< G2w[j].size()<<" index: "<<w3n*N_w_aux+w4n<<std::endl;
              //std::cout<<"F: "<<F2w[i][w1n*N_w_aux+w2n]<<" size: "<< F2w[i].size()<<" index: "<<w1n*N_w_aux+w2n<<std::endl;
              //std::cout<<"done G and F"<<std::endl;
              std::complex<double> meas_h =F2w[i][w1n*N_w_aux+w2n]*G2w[j][w3n*N_w_aux+w4n]; // n1M12M34
              if(i==j)             meas_h-=F2w[i][w1n*N_w_aux+w4n]*G2w[i][w3n*N_w_aux+w2n]; //-n1M14M32
              meas_h/=beta; meas_h*=sign;
              if(i >= F2w.size() || j >= G2w.size()) throw std::logic_error("size is too large!");
              if(w1n*N_w_aux+w2n>= F2w[i].size() || w3n*N_w_aux+w4n>=G2w[j].size()) throw std::logic_error("size 2 is too large");
              //if(isinf(meas_h.real()) || isinf(meas_h.imag()) || isnan(meas_h.real()) || isnan(meas_h.imag())) throw std::runtime_error("inf or nan in measuring of h");
              h2wr[index] += meas_h.real();
              h2wi[index] += meas_h.imag();
              //std::cout<<i<<" "<<j<<" "<<w2n<<" "<<w3n<<" "<<Wn<<" done."<<std::endl;
            }
          }//Wn
      if(MEASURE_g2w){
        measurements[g2wr_names[i][j]]<<g2wr;
        measurements[g2wi_names[i][j]]<<g2wi;
        memset(&(g2wr[0]),0, g2wr.size()*sizeof(double));
        memset(&(g2wi[0]),0, g2wi.size()*sizeof(double));
      }
      if(MEASURE_h2w){
        measurements[h2wr_names[i][j]]<<h2wr;
        measurements[h2wi_names[i][j]]<<h2wi;
        memset(&(h2wr[0]),0, h2wr.size()*sizeof(double));
        memset(&(h2wi[0]),0, h2wi.size()*sizeof(double));
      }
    }//j
  }//i
}

/*
void hybridization::accumulate_G2w(){
  for(std::size_t i=0;i<n_orbitals;++i){
    for(std::size_t j=0;j<=i;++j){//we measure only for j<=i since results for ij and ji are exactly the same (no gain through averaging)
      if(MEASURE_g2w){
        measurements[g2wr_names[i][j]]<<g2wr[i][j];
        measurements[g2wi_names[i][j]]<<g2wi[i][j];
        memset(&(g2wr[i][j][0]),0, g2wr[i][j].size()*sizeof(double));
        memset(&(g2wi[i][j][0]),0, g2wi[i][j].size()*sizeof(double));
      }
      if(MEASURE_h2w){
        measurements[h2wr_names[i][j]]<<h2wr[i][j];
        measurements[h2wi_names[i][j]]<<h2wi[i][j];
        memset(&(h2wr[i][j][0]),0, h2wr[i][j].size()*sizeof(double));
        memset(&(h2wi[i][j][0]),0, h2wi[i][j].size()*sizeof(double));
      }
    }
  }
}
*/

