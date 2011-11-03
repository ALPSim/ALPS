/************************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2011 by Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Emanuel Gull <gull@phys.columbia.edu>,
 *                              Hartmut Hafermann <hafermann@cpht.polytechnique.fr>
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
 ************************************************************************************/

#include "impurity.h"
#include "update.h"
//#include "exponential.h"

/*
implements the measurement of Green's function in the Legendre orthogonal
polynomial basis, PRB 84, 075145 (2011) and the improved estimators for
the self-energy and vertex function, arXiv:1108.1936.
*/

using namespace std;

void hybridization::create_measurements(){
  //basic measurements
  measurements << alps::ngs::RealVectorObservable("n")
               << alps::ngs::RealVectorObservable("order")
               << alps::ngs::RealVectorObservable("Greens")
               << alps::ngs::RealObservable("sign")
               << alps::ngs::RealVectorObservable("matrix_size");

  //additional measurements
  if(MEASURE_gw){ 
    measurements << alps::ngs::RealVectorObservable("gw_re");
    measurements << alps::ngs::RealVectorObservable("gw_im");
  }
  if(MEASURE_fw){ 
    measurements << alps::ngs::RealVectorObservable("fw_re");
    measurements << alps::ngs::RealVectorObservable("fw_im");
  }
  if(MEASURE_gl) 
    measurements << alps::ngs::RealVectorObservable("gl");
  if(MEASURE_fl) 
    measurements << alps::ngs::RealVectorObservable("fl");
  
  if(MEASURE_g2w){ 
    measurements << alps::ngs::RealVectorObservable("g2w_re");
    measurements << alps::ngs::RealVectorObservable("g2w_im");
  }
  if(MEASURE_hw){ 
    measurements << alps::ngs::RealVectorObservable("hw_re");
    measurements << alps::ngs::RealVectorObservable("hw_im");
  }
  if(MEASURE_nn) 
    measurements << alps::ngs::RealVectorObservable("nn");
  if(MEASURE_nnt) 
    measurements << alps::ngs::RealVectorObservable("nnt");

return;
}


void hybridization::resize_measurement_vectors(int crank){

  //basic measurements
  n_meas.resize(FLAVORS, 0.);
  order_meas.resize(N_order*FLAVORS, 0.);
  n_vectors.resize(FLAVORS);
  matrix_size.resize(FLAVORS, 0.);

  //additional measurements
  if(MEASURE_gw){ gw_meas_re.resize(FLAVORS*N_w, 0.); gw_meas_im.resize(FLAVORS*N_w); if(!crank) cout << "measuring g(w)" << endl;}
  if(MEASURE_fw){ fw_meas_re.resize(FLAVORS*N_w, 0.); fw_meas_im.resize(FLAVORS*N_w); if(!crank) cout << "measuring f(w)" << endl;}
  if(MEASURE_gl){ gl_meas.resize(FLAVORS*N_l, 0.); if(!crank) cout << "measuring g(l)" << endl;}
  if(MEASURE_fl){ fl_meas.resize(FLAVORS*N_l, 0.); if(!crank) cout << "measuring f(l)" << endl;}

//  if(MEASURE_nn || MEASURE_nnt) n_vectors.resize(FLAVORS,std::vector<double>(N_nn+1));
  if(MEASURE_nn || MEASURE_nnt){ n_vectors.resize(FLAVORS);
    for(int f=0; f<FLAVORS; ++f) n_vectors[f].resize(N_nn+1);
  }

//  if(MEASURE_nn) nn_meas.resize(FLAVORS*FLAVORS, 0.);
  if(MEASURE_nn){ nn_meas.resize((FLAVORS*(FLAVORS+1))/2, 0.); if(!crank) cout << "measuring <nn>" << endl; }
  if(MEASURE_nnt){
    nnt_meas.resize(((FLAVORS*(FLAVORS+1))/2)*(N_nn+1));//store only for f1>=f2
    if(!crank) cout << "measuring <n(tau)n(0)>" << endl;
  }

  N_w_aux=N_w2+N_W-1;
  if(MEASURE_g2w){ g2w_meas_re.resize(FLAVORS*FLAVORS*N_w2*N_w2*N_W, 0.);
                   g2w_meas_im.resize(FLAVORS*FLAVORS*N_w2*N_w2*N_W, 0.);
                   mw_meas.resize(FLAVORS*N_w_aux*N_w_aux); 
                   if(!crank) cout << "measuring g2(w)" << endl;
                 }
  if(MEASURE_hw){
                  hw_meas_re.resize(FLAVORS*FLAVORS*N_w2*N_w2*N_W, 0.);
                  hw_meas_im.resize(FLAVORS*FLAVORS*N_w2*N_w2*N_W, 0.);
                  nmw_meas.resize(FLAVORS*N_w_aux*N_w_aux); 
                  if(!crank) cout << "measuring h(w)" << endl; 
                }

return;
}


void hybridization::set_measurement_vectors(){

//  int fine_grid=1000*BETA;
//  static TABULATED_FERMIONIC_IMAGINARY_EXPONENTIALS EXP(BETA,N_w,fine_grid); //elements can be accessed as exp(wn,tau)
//  static TABULATED_BOSONIC_IMAGINARY_EXPONENTIALS EXPW(BETA,N_W,fine_grid); //elements can be accessed as exp(Wn,tau)

  double s=1;
  static const std::complex<double> I(0.,1.);

  if(MEASURE_g2w){//measure <T c c^* c c^*> (a1 c1 a2 c2)
    memset(&( mw_meas[0]),0,  mw_meas.size()*sizeof(std::complex<double>));
    memset(&(nmw_meas[0]),0, nmw_meas.size()*sizeof(std::complex<double>));
  }
  for(int f=0; f<FLAVORS; f++){
    s *= sign[f];
    segment_container_t::iterator itc1, ita1;
    static times full_segment(0,BETA);
    n_meas[f] += compute_overlap(full_segment, segments[f], full_line[f], BETA)/BETA;//measure density
    if(segments[f].size()<(std::size_t)N_order)
      order_meas[f*N_order+segments[f].size()]++; // measure perturbation order
    matrix_size[f]+=M[f].size(); //measure Matrix size
    double tau_1, tau_2;
    if(segments[f].size()>0){//at least one segment in the configuration
      for(int a1=0; a1<(int)M[f].size1(); a1++){
        (a1==0 ? ita1 = segments[f].begin() : ita1++); //1st annihilator
        tau_1=ita1->t_end(); double pref(0.0);
        if(MEASURE_fw || MEASURE_fl || MEASURE_hw)
          for(int f1=0; f1<FLAVORS; f1++)
            pref += 0.5*(u(f1,f)+u(f,f1))*get_occupation(segments[f1],full_line[f1],tau_1,BETA);

        for(int c1=0; c1<(int)M[f].size1(); c1++){
          (c1==0 ? itc1 = segments[f].begin() : itc1++); //1st creator
          tau_2=itc1->t_start();
          double m_matrix = M[f](c1,a1);
          double arg = tau_1-tau_2; //measures -<T c(tau_1) c^*(tau_2)>
          //measurement of G(tau)
          double argument=arg; double bubble_sign=1;
          if(argument > 0) bubble_sign = 1;
          else{
            bubble_sign = -1;
            argument += BETA;
          }
          //all values in the interval [0,BETA/(2N) ) go to bin zero
          //all value int the interval (BETA-BETA/(2N),BETA] go to bin N
          int index = (int)(argument/BETA*N+0.5);
          G_meas[f*(N+1)+index] += m_matrix*bubble_sign;//measures +<T c c^*>; need to make definitions consistent
          if(MEASURE_gl){
            double x=2.0*argument/BETA-1.0;
            double pl_2=1; double pl_1=x; double legendre_p;
            for(int l=0; l<N_l; l++){
              if(l==0) legendre_p=1;
              else if(l==1) legendre_p=x;
              else{ 
                legendre_p=((2*l-1)*x*pl_1-(l-1)*pl_2)/static_cast<double>(l);//l
                pl_2=pl_1; //l-2
                pl_1=legendre_p; //l-1
               }
              gl_meas[f*N_l+l] -= m_matrix*legendre_p*bubble_sign;
              if(MEASURE_fl){
                fl_meas[f*N_l+l] -= m_matrix*pref*legendre_p*bubble_sign;
              }
            }
          }//measure gl
          if(MEASURE_gw){
            std::complex<double> exp_=exp( -1.0*I*M_PI*arg/BETA);
            std::complex<double> dexp=exp( 2.0*I*M_PI*arg/BETA);
            for(int wn=0; wn<N_w; wn++){
              exp_*=dexp;
              std::complex<double> meas = -m_matrix*exp_;//note the -
//              std::complex<double> meas = -m_matrix*EXP(wn,arg);//note the -
              gw_meas_re[f*N_w+wn] += meas.real();
              gw_meas_im[f*N_w+wn] += meas.imag();
              if(MEASURE_fw){
                std::complex<double> fmeas(0.0);
//                fmeas -= m_matrix*pref*EXP(wn,arg);
                fmeas -= m_matrix*pref*exp_;
                fw_meas_re[f*N_w+wn] += fmeas.real();
                fw_meas_im[f*N_w+wn] += fmeas.imag();
              }
            }
          }//measure_gw
          if(MEASURE_g2w){//measure <T c c^* c c^*> (a1 c1 a2 c2); factorized measurement
            static const int Nwh=N_w2/2;
            std::complex<double> dexp1=exp( 2.0*I*M_PI*tau_1/BETA);
            std::complex<double> dexp2=exp(-2.0*I*M_PI*tau_2/BETA);
            std::complex<double> expinit1=exp(I*( (2*(-Nwh-1)+1)*M_PI*tau_1/BETA));
            std::complex<double> expinit2=exp(-I*( (2*(-Nwh-1)+1)*M_PI*tau_2/BETA));
            std::complex<double> exp1=expinit1;//fast update of the exponential
            for(int w1n=0; w1n<N_w_aux; w1n++){
              exp1*=dexp1;
              std::complex<double> exp2=expinit2;
              //std::complex<double> exp1=EXP(w1n-Nwh,tau_1);//exp(I w1 t)
              for(int w2n=0; w2n<N_w_aux; w2n++){
                exp2*=dexp2;
                //std::complex<double> exp2=EXP(w2n-Nwh,-tau_2);//exp(-I w2 t)
                mw_meas[f*N_w_aux*N_w_aux+w1n*N_w_aux+w2n]+=m_matrix*exp1*exp2;
                if(MEASURE_hw) nmw_meas[f*N_w_aux*N_w_aux+w1n*N_w_aux+w2n]+=pref*m_matrix*exp1*exp2;//prefactor has opposite sign for h
              }
            }
          }//MEASURE_g2w
        }
      }
    }//segments.size()>0
  }//f

  if(MEASURE_g2w){
    for(int f1=0; f1<FLAVORS; f1++)
      for(int f2=0; f2<FLAVORS; f2++)//measure updn and dnup for averaging
        for(int w2n=0; w2n<N_w2; w2n++)
          for(int w3n=0; w3n<N_w2; w3n++)
            for(int Wn=0; Wn<N_W; Wn++){
              int w1n=w2n+Wn; int w4n=w3n+Wn;
//              int index=(FLAVORS*f1+f2) *N_w2*N_w2*N_W + w2n *N_w2*N_W + w3n *N_W + Wn;
              int index=Wn*FLAVORS*FLAVORS*N_w2*N_w2 + (FLAVORS*f1+f2) *N_w2*N_w2 + w2n*N_w2 + w3n;
              std::complex<double> meas =mw_meas[f1*N_w_aux*N_w_aux+w1n*N_w_aux+w2n]*mw_meas[f2*N_w_aux*N_w_aux+w3n*N_w_aux+w4n];
              if(f1==f2)   meas-=mw_meas[f1*N_w_aux*N_w_aux+w1n*N_w_aux+w4n]*mw_meas[f1*N_w_aux*N_w_aux+w3n*N_w_aux+w2n];
              g2w_meas_re[index] += meas.real();
              g2w_meas_im[index] += meas.imag();
              if(MEASURE_hw){
                std::complex<double> meas_h =nmw_meas[f1*N_w_aux*N_w_aux+w1n*N_w_aux+w2n]*mw_meas[f2*N_w_aux*N_w_aux+w3n*N_w_aux+w4n];
                if(f1==f2)  meas_h-=nmw_meas[f1*N_w_aux*N_w_aux+w1n*N_w_aux+w4n]*mw_meas[f1*N_w_aux*N_w_aux+w3n*N_w_aux+w2n];
                hw_meas_re[index] += meas_h.real();
                hw_meas_im[index] += meas_h.imag();
              }
            }//Wn
  }//end::MEASURE_g2w

  sign_meas += s;

  if(MEASURE_nn || MEASURE_nnt){  
    segment_container_t::iterator it;
    for(int flavor=0; flavor<FLAVORS; ++flavor){
      for(int i=0; i<N_nn+1; i++) n_vectors[flavor][i]=1;//reset
      if(segments[flavor].size()==0){
        if(full_line[flavor]==0){
        for(std::size_t i=0; i<n_vectors[flavor].size(); ++i)
          n_vectors[flavor][i]=0;
        }
      }
      else{
        it=segments[flavor].end(); it--;
        if(it->t_end()<it->t_start()) n_vectors[flavor][0]=1;//last segment winds around the circle //n(0)=1
        else n_vectors[flavor][0]=0;//n(0)=0
        // mark segment start and end points
        int index;
        for(it=segments[flavor].begin(); it!=segments[flavor].end(); it++){
          index = (int)(it->t_start()/BETA*N_nn+1);
          n_vectors[flavor][index] *= -1;
          index = (int)(it->t_end()/BETA*N_nn+1);
          n_vectors[flavor][index] *= -1;
        }
        // fill vector with occupation number
        for(std::size_t i=1; i<n_vectors[flavor].size(); i++){
          if(n_vectors[flavor][i]==-1) n_vectors[flavor][i]=1-n_vectors[flavor][i-1];//segment starts or ends -> occupation number changes
          else n_vectors[flavor][i]=n_vectors[flavor][i-1];//on the same segment -> occupation number identical to that of previous segment
        }
      }
    }

  //this gives the same result, but is somewhat slower
  /*
  for(int flavor=0; flavor<FLAVORS; ++flavor){
    for(int i=0; i<N_nn+1; i++){ 
      double tau = i*BETA/static_cast<double>(N_nn);
      n_vectors[flavor][i]=get_occupation(segments[flavor],full_line[flavor],tau,BETA);
  }
  }
*/
  }//end MEASURE_nn || MEASURE_nnt

  if(MEASURE_nn){
  // compute n(0)n(0)
    int pos=0;
    for(int flavor1=0; flavor1<FLAVORS; ++flavor1){
      for(int flavor2=0; flavor2<=flavor1; ++flavor2){
        for(int i=0;i<N_nn+1;++i){
          nn_meas[pos] += n_vectors[flavor1][i]*n_vectors[flavor2][i]; 
        }
        pos++;
      }
    }
  }//end::MEASURE_nn

  if(MEASURE_nnt){
    // compute n(\tau)n(0)
    int pos=0;
    for(int flavor1=0; flavor1<FLAVORS; ++flavor1){
      for(int flavor2=0; flavor2<=flavor1; ++flavor2){
      //for(int i=0; i<N_nn+1; ++i){// COMPUTING ALL PAIRS i,j MAY BE TOO SLOW
        int i=0;
        for(int index=0; index<N_nn+1; ++index){
          int j=i+index;
          if(j>N_nn) j -= N_nn; //no sign change, this correlator is bosonic
          nnt_meas[pos+index] += n_vectors[flavor1][j]*n_vectors[flavor2][i];
        }
      //}
        pos += (N_nn+1);
      }
    }
  }//end::MEASURE_nnt

return;
}//set_measurements


void hybridization::do_measurements(){
  if(is_thermalized()){
  //basic measurements
    order_meas /= N_meas;
    measurements["order"] << order_meas;
    
    sign_meas /= N_meas;
    measurements["sign"] << sign_meas;
    //if(sign_meas != 1.) throw std::runtime_error("negative sign encountered. The current code is not able to deal with this (in the segment representation such a situation should not arise for diagonal hybridizations). Do you know what you're doing?!?");
    
    n_meas /= N_meas;
    measurements["n"]<< n_meas;
    
    matrix_size /= N_meas;
    measurements["matrix_size"]<< matrix_size;

    G_meas *= (1.*N)/N_meas/(BETA*BETA);
    measurements["Greens"] << (G_meas); //*sign;

    //additional measurements
    if(MEASURE_gw){
      gw_meas_re *= 1./(N_meas*BETA);
      gw_meas_im *= 1./(N_meas*BETA);
      measurements["gw_re"] << gw_meas_re;
      measurements["gw_im"] << gw_meas_im;
    }
    if(MEASURE_fw){
      fw_meas_re *= 1./(N_meas*BETA);
      fw_meas_im *= 1./(N_meas*BETA);
      measurements["fw_re"] << fw_meas_re;
      measurements["fw_im"] << fw_meas_im;
    }
    if(MEASURE_gl){
      gl_meas *=1./(N_meas*BETA);
      measurements["gl"] << gl_meas;
    }
    if(MEASURE_fl){
      fl_meas *=1./(N_meas*BETA);
      measurements["fl"] << fl_meas;
    }
    if(MEASURE_nn){
      nn_meas*=1./(N_meas*(N_nn+1));//also divide zero-time correlator bin N_nn+1, since we measure at N_nn+1 times
      measurements["nn"] << nn_meas;
    }
    if(MEASURE_nnt){
      nnt_meas*=1./(N_meas);
      measurements["nnt"] << nnt_meas;
    }
    if(MEASURE_g2w){
      g2w_meas_re *= 1./(N_meas*BETA);
      g2w_meas_im *= 1./(N_meas*BETA);
      measurements["g2w_re"] << g2w_meas_re;
      measurements["g2w_im"] << g2w_meas_im;
    }
    if(MEASURE_hw){
      hw_meas_re *= 1./(N_meas*BETA);
      hw_meas_im *= 1./(N_meas*BETA);
      measurements["hw_re"] << hw_meas_re;
      measurements["hw_im"] << hw_meas_im;
    }
}//if thermalized

    //reset
    reset_measurement_vectors();
return;
}


void hybridization::reset_measurement_vectors(){
  //basic measurements
  memset(&(n_meas[0]),0, n_meas.size()*sizeof(double));
  memset(&(order_meas[0]),0, order_meas.size()*sizeof(double));
  memset(&(matrix_size[0]),0, matrix_size.size()*sizeof(double));
  memset(&(G_meas[0]),0, G_meas.size()*sizeof(double));
  sign_meas=0;

  //additional measurements
  if(MEASURE_gw){ memset(&(gw_meas_re[0]),0, gw_meas_re.size()*sizeof(double));
		  memset(&(gw_meas_im[0]),0, gw_meas_im.size()*sizeof(double));
  }
  if(MEASURE_fw){ memset(&(fw_meas_re[0]),0, fw_meas_re.size()*sizeof(double));
		  memset(&(fw_meas_im[0]),0, fw_meas_im.size()*sizeof(double));
  }
  if(MEASURE_gl) memset(&(gl_meas[0]),0, gl_meas.size()*sizeof(double));
  if(MEASURE_fl) memset(&(fl_meas[0]),0, fl_meas.size()*sizeof(double));

  if(MEASURE_nn) memset(&(nn_meas[0]), 0, nn_meas.size()*sizeof(double));
  if(MEASURE_nnt) memset(&(nnt_meas[0]), 0, nnt_meas.size()*sizeof(double));

  if(MEASURE_g2w){ memset(&(g2w_meas_re[0]),0, g2w_meas_re.size()*sizeof(double));
		   memset(&(g2w_meas_im[0]),0, g2w_meas_im.size()*sizeof(double));
  }
  if(MEASURE_hw){ memset(&(hw_meas_re[0]),0, hw_meas_re.size()*sizeof(double));
		  memset(&(hw_meas_im[0]),0, hw_meas_im.size()*sizeof(double));
  }

return;
}


