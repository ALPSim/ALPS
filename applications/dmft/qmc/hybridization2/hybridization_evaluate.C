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

#include <iostream>
#include <fstream>
#include <complex>
#include <cmath>
#include <alps/alea.h>
#include <alps/hdf5.hpp>
#include <alps/hdf5/pointer.hpp>
#include <alps/hdf5/complex.hpp>
#include <alps/hdf5/vector.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/constants/constants.hpp>

using namespace boost::math;

const std::complex<double> i_c(0.0,1.0);
const double pi = constants::pi<double>();

std::complex<double> t(int n, int l){
return (sqrt(2*l+1)/sqrt(2*n+1)) * exp(i_c*(n+0.5)*pi) * pow(i_c,l) * cyl_bessel_j(l+0.5,(n+0.5)*pi);
}

template <typename T>
T sqr(T x){
  return x*x;
}

template <typename T>
void spin_average(std::vector<T> &vec, int FLAVORS){
  for(int orb=0; orb<FLAVORS/2; ++orb){
    int fup=2*orb, fdn=2*orb+1, Ncoef=vec.size()/FLAVORS;
    for(int n=0; n<Ncoef; ++n){
      vec[fup*Ncoef+n]=0.5*(vec[fup*Ncoef+n]+vec[fdn*Ncoef+n]);
      vec[fdn*Ncoef+n]=vec[fup*Ncoef+n];
    }
  }
}


int main(int argc, char ** argv){
  std::string basename=std::string(argv[1]);
  bool write_hr = false;
  if(argc==3 && std::string(argv[2]) == "hr"){
    write_hr=true; //write results in human readable format for quick plotting
    std::cout << "recognized option \"hr\"" << std::endl;
  }
  double BETA;
  double N;
  int N_ORDER;
  int FLAVORS;
  int N_MEAS;
  int MEASURE_gw;
  int MEASURE_fw;
  int MEASURE_gl;
  int MEASURE_fl;
  int MEASURE_g2w;
  int MEASURE_hw;
  int MEASURE_nn;
  int MEASURE_nnt;
  int N_w;
  int N_W;
  int N_w2;
  int N_l;
  int N_nn;
  int PARAMAGNETIC;
  std::ofstream str;
  std::vector<double> Greens;
  std::vector<double> density;
  std::vector<std::complex<double> >gw;//also needed for evaluation
  std::vector<std::complex<double> >fw;//of the vertex function
  {
    alps::hdf5::archive ar(basename+".out.h5", alps::hdf5::archive::READ);
    ar>>alps::make_pvp("/parameters/BETA", BETA);
    ar>>alps::make_pvp("/parameters/N", N);
    ar>>alps::make_pvp("/parameters/N_ORDER", N_ORDER);
    ar>>alps::make_pvp("/parameters/FLAVORS", FLAVORS);
    ar>>alps::make_pvp("/parameters/N_MEAS", N_MEAS);
    ar>>alps::make_pvp("/parameters/MEASURE_gw", MEASURE_gw);
    ar>>alps::make_pvp("/parameters/MEASURE_fw", MEASURE_fw);
    ar>>alps::make_pvp("/parameters/MEASURE_gl", MEASURE_gl);
    ar>>alps::make_pvp("/parameters/MEASURE_fl", MEASURE_fl);
    ar>>alps::make_pvp("/parameters/MEASURE_g2w", MEASURE_g2w);
    ar>>alps::make_pvp("/parameters/MEASURE_hw", MEASURE_hw);
    ar>>alps::make_pvp("/parameters/MEASURE_nn", MEASURE_nn);
    ar>>alps::make_pvp("/parameters/MEASURE_nnt", MEASURE_nnt);
    ar>>alps::make_pvp("/parameters/N_w", N_w);
    ar>>alps::make_pvp("/parameters/N_W", N_W);
    ar>>alps::make_pvp("/parameters/N_w2", N_w2);
    ar>>alps::make_pvp("/parameters/N_l", N_l);
    ar>>alps::make_pvp("/parameters/N_nn", N_nn);
    ar>>alps::make_pvp("/parameters/PARAMAGNETIC", PARAMAGNETIC);

    double N_SWEEPS; ar>>alps::make_pvp("/simulation/results/order/count",N_SWEEPS);
  
    std::vector<double> matrix_size; matrix_size.resize(FLAVORS); ar>>alps::make_pvp("/simulation/results/matrix_size/mean/value",matrix_size);

    std::cout << "simulation details:" << std::endl;
    std::cout << "average sign: " << 1 << std::endl;
    std::cout << "total number of sweeps: " << N_SWEEPS << std::endl;
    std::cout << "total number of measurements: " << N_SWEEPS*N_MEAS << std::endl;
    std::cout << "perturbation order:" << std::endl;

    for(int f=0; f<FLAVORS; ++f)
      std::cout << "orbital " << f/2 << " spin " << f%2 << ": " << matrix_size[f] << std:: endl;

    std::cout << std::endl;
    std::cout << "evaluating observables:" << std::endl;

    std::cout << "evaluating density.." << std::flush;
    density.resize(FLAVORS); ar>>alps::make_pvp("/simulation/results/n/mean/value",density);
    if(PARAMAGNETIC) spin_average(density,FLAVORS);
    std::cout << "done." << std::endl;
  
    if(write_hr){
      str.open("observables.dat"); double total_density(0.);
      for(int orb=0; orb<FLAVORS/2; ++orb){
        int fup=2*orb; int fdn=2*orb+1;
        str << "N" << orb+1 << "up=" << density[fup] << ";" << std::endl;
        str << "N" << orb+1 << "dn=" << density[fdn] << ";" << std::endl;
        str << "N" << orb+1 << "=" << density[fup]+density[fdn] << ";" << std::endl;
        total_density += density[fup]+density[fdn];
      }
      str << "N=" << total_density << ";" << std::endl;
      str.close();
    }
    
    std::vector<double> order; order.resize(FLAVORS*N_ORDER); ar>>alps::make_pvp("/simulation/results/order/mean/value",order);
    if(write_hr){
      std::cout << "evaluating histogram..." << std::flush;
      str.open("histogram.dat");
      for(int n=0; n<N_ORDER; ++n){
        str << n;
        for(int f=0; f<FLAVORS; ++f)
          str << "  " << order[f*N_ORDER+n];
        str << std::endl;
      }
      str.close();
      std::cout << "done." << std::endl;
    }
    std::cout << "evaluating gt..." << std::flush;
    Greens.resize(FLAVORS*(N+1)); 
    ar>>alps::make_pvp("/simulation/results/Greens/mean/value",Greens);
  }//closes ar
    
    if(PARAMAGNETIC) spin_average(Greens, FLAVORS);
    for(int i=0; i<Greens.size(); ++i) Greens[i]*=-1; //different sign convention
    for(int f=0; f<FLAVORS; ++f){
      Greens[f*(N+1)]=-(1-density[f]);
      Greens[f*(N+1)+N]=-density[f];
    }
    {
      alps::hdf5::archive oar(basename+".out.h5", alps::hdf5::archive::WRITE);
      oar << alps::make_pvp("/simulation/results/gt",Greens);
    }
    if(write_hr){
      str.open("gt.dat");
      for(int i=0; i<N+1; ++i){
        str << i*BETA/static_cast<double>(N);
        for(int f=0; f<FLAVORS; ++f)
          str << " " << Greens[f*(N+1)+i];
        str  << std::endl;
      }
      str.close();
    }
    std::cout << "done." << std::endl;
  
  if(MEASURE_nn || write_hr){
    std::cout << "evaluating <nn>..." << std::flush;
    std::vector<double> nn; nn.resize(FLAVORS*(FLAVORS+1)/2);  
    {
      alps::hdf5::archive ar(basename+".out.h5", alps::hdf5::archive::READ);
      ar>>alps::make_pvp("/simulation/results/nn/mean/value",nn);
    }
    if(write_hr){
      str.open("observables.dat",std::ios::app); //write in format that is easy to read from python
      double SzSz=0.0;//total spin
      for(int f1=0; f1<FLAVORS; f1++){
        for(int f2=0; f2<=f1; f2++){
          int n1=f1/2; int s1=f1%2;
          int n2=f2/2; int s2=f2%2;
          double docc=nn[(f1*(f1+1))/2+f2];
          str << "N" << n1+1; if(s1==0) str << "up"; else str << "dn";
          str << "N" << n2+1; if(s2==0) str << "up"; else str << "dn";
          str << "=" << docc << ";" << std::endl;
        }
      }
      for(int n1=0; n1<FLAVORS/2; ++n1){
        for(int n2=0; n2<=n1; ++n2){
          int f1up=2*n1; int f1dn=2*n1+1;
          int f2up=2*n2; int f2dn=2*n2+1;
          double szsz=  nn[f1up*(f1up+1)/2+f2up] - nn[f1up*(f1up+1)/2+f2dn]
                      - nn[f1dn*(f1dn+1)/2+f2up] + nn[f1dn*(f1dn+1)/2+f2dn];
                 szsz*=0.25; SzSz+=szsz;
          str << "s" << n1+1 << "zs" << n2+1 << "z=" << szsz << ";" << std::endl;
        }
      }
      str << "SzSz=" << SzSz << ";" << std::endl;
      str.close();
    }
    std::cout << "done." << std::endl;
  }//MEASURE_nn
  
  if(MEASURE_nnt && write_hr){
    std::cout << "evaluating <Sz(tau)Sz(0)>..." << std::flush;
    std::vector<double> nnt; nnt.resize(sqr(FLAVORS)*(N_nn+1));
    {
      alps::hdf5::archive ar(basename+".out.h5", alps::hdf5::archive::READ);
      ar>>alps::make_pvp("/simulation/results/nnt/mean/value",nnt);
    }
    std::vector<double> szsz_static; szsz_static.resize((FLAVORS/2)*(FLAVORS/2+1)/2, 0.);
    double SzSz_static(0.);
    str.open("szszt.dat");
    str << "#tau  ";
    for(int n1=0; n1<FLAVORS/2; n1++)//orbitals
        for(int n2=0; n2<=n1; n2++)//orbitals
          str << "s" << n1+1 << "zs" << n2+1 << "z val, err, rel. err"; 
    str << std::endl;
    for(int i=0; i<N_nn+1; ++i){
      double tau=i*BETA/static_cast<double>(N_nn);
      str << tau;
      int pos=0;
      double SzSzt(0.);
      for(int n1=0; n1<FLAVORS/2; n1++){//orbitals
        for(int n2=0; n2<=n1; n2++){//orbitals
          int f1up = 2*n1; int f1dn = 2*n1+1;
          int f2up = 2*n2; int f2dn = 2*n2+1;
          double szsz     =  0.25* (  nnt[(f1up*(f1up+1)/2+f2up)*(N_nn+1)+i] - nnt[(f1up*(f1up+1)/2+f2dn)*(N_nn+1)+i]
                                    - nnt[(f1dn*(f1dn+1)/2+f2up)*(N_nn+1)+i] + nnt[(f1dn*(f1dn+1)/2+f2dn)*(N_nn+1)+i] );

          if(i==0 || i==N_nn){ SzSz_static+=0.5*szsz; //trapezoidal rule: factor 1/2 for boundary terms
                               szsz_static[n1*(n1+1)/2+n2]+=0.5*szsz;
                             }
          else{
            SzSz_static+=szsz;
            szsz_static[n1*(n1+1)/2+n2]+=szsz;
          }
          SzSzt+=szsz;
        }
      }
      str << "      " << SzSzt;
      str << std::endl;
    }
    str.close();

    str.open("observables.dat",std::ios::app);
    for(int n1=0; n1<FLAVORS/2; ++n1)
      for(int n2=0; n2<=n1; ++n2)
        str << "s" << n1+1 << "zs" << n2+1 << "z_static=" << BETA*szsz_static[n1*(n1+1)/2+n2]/static_cast<double>(N_nn) << ";" << std::endl;
        str << "SzSz_static=" << BETA*SzSz_static/static_cast<double>(N_nn) << ";" << std::endl;
        str.close();

    std::cout << "done." << std::endl;
  }


  if(MEASURE_gw && N_w >0){
    std::cout << "evaluating gw..." << std::flush;
    std::vector<double> gw_re; gw_re.resize(FLAVORS*N_w);
    std::vector<double> gw_im; gw_im.resize(FLAVORS*N_w);
    {
      alps::hdf5::archive ar(basename+".out.h5", alps::hdf5::archive::READ);
      ar>>alps::make_pvp("/simulation/results/gw_re/mean/value",gw_re);
      ar>>alps::make_pvp("/simulation/results/gw_im/mean/value",gw_im);
    }
    for(int i=0; i<gw_re.size(); ++i) gw.push_back(std::complex<double>(gw_re[i],gw_im[i]));
    if(PARAMAGNETIC)spin_average(gw,FLAVORS);
    {
      alps::hdf5::archive oar(basename+".out.h5", alps::hdf5::archive::WRITE);
      oar << alps::make_pvp("/simulation/results/gw",gw);
    }
    if(write_hr){
      str.open("gw.dat");
      for(int wn=0; wn<N_w; ++wn){
        str << (2*wn+1)*M_PI/BETA;
        for(int f=0; f<FLAVORS; ++f)
          str << " " << gw[f*N_w+wn].real() << " " << gw[f*N_w+wn].imag();
        str  << std::endl;
      }
      str.close();
    }
    std::cout << "done." << std::endl;
    if(MEASURE_fw){
      std::cout << "evaluating sigmaw..." << std::flush;
      std::vector<double> fw_re; fw_re.resize(FLAVORS*N_w);
      std::vector<double> fw_im; fw_im.resize(FLAVORS*N_w);
      {
        alps::hdf5::archive ar(basename+".out.h5", alps::hdf5::archive::READ);
        ar>>alps::make_pvp("/simulation/results/fw_re/mean/value",fw_re);
        ar>>alps::make_pvp("/simulation/results/fw_im/mean/value",fw_im);
      }

      std::vector<std::complex<double> > sigmaw;
      if(PARAMAGNETIC){
        spin_average(fw_re,FLAVORS);
        spin_average(fw_im,FLAVORS);
      }
      for(int i=0; i<fw_re.size(); ++i){ 
        fw.push_back(std::complex<double>(fw_re[i],fw_im[i]));
        sigmaw.push_back(fw[i]/gw[i]);
      }
      {
        alps::hdf5::archive oar(basename+".out.h5", alps::hdf5::archive::WRITE);
        oar << alps::make_pvp("/simulation/results/sigmaw",sigmaw);
        oar << alps::make_pvp("/simulation/results/fw",fw);
      }
      if(write_hr){
        str.open("sigmaw.dat");
        for(int wn=0; wn<N_w; ++wn){
          str << (2*wn+1)*M_PI/BETA;
          for(int f=0; f<FLAVORS; ++f){
            str << " " << real(sigmaw[f*N_w+wn]) << " " << imag(sigmaw[f*N_w+wn]);
          }
          str  << std::endl;
        }
        str.close();
      }
      std::cout << "done." << std::endl;
    }
  }

  if(MEASURE_gl && N_l >0){
    std::cout << "evaluating gw from l..." << std::flush;
    std::vector<double> gl; gl.resize(FLAVORS*N_l); 
    {
      alps::hdf5::archive ar(basename+".out.h5", alps::hdf5::archive::READ);
      ar>>alps::make_pvp("/simulation/results/gl/mean/value",gl);
    }
    if(PARAMAGNETIC) spin_average(gl,FLAVORS);
    std::vector<std::complex<double> >gwl;
    for(int f=0; f<FLAVORS; ++f)
      for(int wn=0; wn<N_w; ++wn){
        std::complex<double> gw_(0.);
          for(int l=0; l<N_l; ++l)
            gw_+=t(wn,l)*sqrt(2.*l+1)*gl[f*N_l+l]; //sqrt(2l+1) has been omitted in the measurement
        gwl.push_back(gw_);
      }
    {
      alps::hdf5::archive oar(basename+".out.h5", alps::hdf5::archive::WRITE);
      oar << alps::make_pvp("/simulation/results/gw_from_l",gwl);
    }
    if(write_hr){
      str.open("gw_from_l.dat");
      for(int wn=0; wn<N_w; ++wn){
        str << (2*wn+1)*M_PI/BETA;
        for(int f=0; f<FLAVORS; ++f)
          str << " " << gwl[f*N_w+wn].real() << " " << gwl[f*N_w+wn].imag();
        str  << std::endl;
      }
      str.close();
    }
    std::cout << "done." << std::endl;
    if(MEASURE_fl){
      std::cout << "evaluating sigmaw from l..." << std::flush;
      std::vector<double> fl; fl.resize(FLAVORS*N_l);

      {
        alps::hdf5::archive ar(basename+".out.h5", alps::hdf5::archive::READ);
        ar>>alps::make_pvp("/simulation/results/fl/mean/value",fl);
      }
      if(PARAMAGNETIC) spin_average(fl,FLAVORS);
      std::vector<std::complex<double> >fwl;
      std::vector<std::complex<double> >sigmawl;
      for(int f=0; f<FLAVORS; ++f)
        for(int wn=0; wn<N_w; ++wn){
          std::complex<double> fw_(0.);
            for(int l=0; l<N_l; ++l)
              fw_+=t(wn,l)*sqrt(2.*l+1)*fl[f*N_l+l]; 
          fwl.push_back(fw_);
          sigmawl.push_back(fw_/gwl[f*N_w+wn]);
        }
      {
        alps::hdf5::archive oar(basename+".out.h5", alps::hdf5::archive::WRITE);
        oar << alps::make_pvp("/simulation/results/fw_from_l",fwl);
        oar << alps::make_pvp("/simulation/results/sigmaw_from_l",sigmawl);
      }
      if(write_hr){
        str.open("sigmaw_from_l.dat");
        for(int wn=0; wn<N_w; ++wn){
          str << (2*wn+1)*M_PI/BETA;
          for(int f=0; f<FLAVORS; ++f){
            str << " " << real(sigmawl[f*N_w+wn]) << " " << imag(sigmawl[f*N_w+wn]);
          }
          str  << std::endl;
        }
        str.close();
      }
      std::cout << "done." << std::endl;
    }
  }


  if(MEASURE_g2w && MEASURE_fw && MEASURE_gw && MEASURE_hw){
    std::cout << "evaluating gamma..." << std::flush;
    const int Nwh=N_w2/2;
    std::vector<double> g2w_re; g2w_re.resize(FLAVORS*FLAVORS*N_w2*N_w2*N_W, 0.);
    std::vector<double> g2w_im; g2w_im.resize(FLAVORS*FLAVORS*N_w2*N_w2*N_W, 0.);
    std::vector<double> hw_re; hw_re.resize(FLAVORS*FLAVORS*N_w2*N_w2*N_W, 0.);
    std::vector<double> hw_im; hw_im.resize(FLAVORS*FLAVORS*N_w2*N_w2*N_W, 0.);
    std::vector<std::complex<double> > gammaw; 
    if(PARAMAGNETIC) gammaw.resize((FLAVORS/2)*FLAVORS*N_w2*N_w2*N_W, 0.);
    else gammaw.resize(FLAVORS*FLAVORS*N_w2*N_w2*N_W, 0.);
    {
      alps::hdf5::archive ar(basename+".out.h5", alps::hdf5::archive::READ);
      ar>>alps::make_pvp("/simulation/results/g2w_re/mean/value",g2w_re);
      ar>>alps::make_pvp("/simulation/results/g2w_im/mean/value",g2w_im);
      ar>>alps::make_pvp("/simulation/results/hw_re/mean/value",hw_re);
      ar>>alps::make_pvp("/simulation/results/hw_im/mean/value",hw_im);
    }
    std::ofstream chi_str;
    std::ofstream chi_irr_str;
    std::ofstream gamma_str;

    if(write_hr){
      chi_str.open("g2w.dat");
      chi_irr_str.open("g2_irr.dat");
      gamma_str.open("gamma.dat");
    }
    std::complex<double> chi, chi_irr, gamma, f, g1, g2, g3, g4;
    std::complex<double> chi_irr_new, chi_irr_sf, gamma_new, gamma_sf, h, chi0;

      for(int n1=0; n1<FLAVORS/2; ++n1)//n1,n2 label orbital indices only
        for(int n2=0; n2<FLAVORS/2; ++n2)
         for(int s1=0; s1<(PARAMAGNETIC ? 1 : 2); ++s1) //in paramagnetic case take spin components 00 and 01 only
          for(int s2=0; s2<2; ++s2){
            int f1=2*n1+s1;
            int f2=2*n2+s2;
            for(int w2n=-Nwh; w2n<Nwh; ++w2n){
              g2=(w2n<0 ? conj(gw[f1*N_w+(-w2n-1)]) : gw[f1*N_w+w2n]);//g2:w2n,n1
              for(int w3n=-Nwh; w3n<Nwh; ++w3n){
                g3=(w3n<0 ? conj(gw[f2*N_w+(-w3n-1)]) : gw[f2*N_w+w3n]);//g3:w4n,n2
                for(int Wn=0; Wn<N_W; ++Wn){
                  int w1n=w2n+Wn; int w4n=w3n+Wn;
                  f = (w1n<0 ? conj(fw[f1*N_w+(-w1n-1)]) : fw[f1*N_w+w1n]);
                  g1= (w1n<0 ? conj(gw[f1*N_w+(-w1n-1)]) : gw[f1*N_w+w1n]);// g1:w1n,f1
                  g4= (w4n<0 ? conj(gw[f2*N_w+(-w4n-1)]) : gw[f2*N_w+w4n]);// g4:w4n,f2
                  chi0=0.0;
                  if(w1n==w2n)chi0+=BETA*g1*g3;
                  if(w2n==w3n && f1==f2) chi0-=BETA*g1*g3; //this gives chi0
                  int f1r=2*n1+(1-s1);
                  int f2r=2*n2+(1-s2);
//                int index  =(FLAVORS*f1 +f2 )*N_w2*N_w2*N_W + (w2n+Nwh)*N_w2*N_W + (w3n+Nwh)*N_W+Wn;
                  int index=Wn*FLAVORS*FLAVORS*N_w2*N_w2 + (FLAVORS*f1+f2) *N_w2*N_w2 + (w2n+Nwh)*N_w2 + w3n+Nwh;
                  if(PARAMAGNETIC){
//                    int index_r=(FLAVORS*f1r+f2r)*N_w2*N_w2*N_W + (w2n+Nwh)*N_w2*N_W + (w3n+Nwh)*N_W+Wn;
                  int index_r=Wn*FLAVORS*FLAVORS*N_w2*N_w2 + (FLAVORS*f1r+f2r) *N_w2*N_w2 + (w2n+Nwh)*N_w2 + w3n+Nwh;
                    chi=0.5*std::complex<double>( g2w_re[index] + g2w_re[index_r],//spin average
                                                  g2w_im[index] + g2w_im[index_r] );
                      h=0.5*std::complex<double>(  hw_re[index] +  hw_re[index_r],//spin average
                                                   hw_im[index] +  hw_im[index_r] );
                  }else{
                    chi=std::complex<double>( g2w_re[index], g2w_im[index] );
                      h=std::complex<double>(  hw_re[index],  hw_im[index] );
                  }
                  chi_irr=g1*h - f*chi; //most accurate results
//                  chi_irr=(g1*h-f*chi0)/(1.0+f); //somewhat less accurate; chi not needed
//                  chi_irr=chi-chi0; //straightforward evaluation; least accurate

                  gamma=chi_irr/(g1*g2*g3*g4);
                  
                  gammaw[index] = gamma;
                  
                  if(write_hr){
                    chi_str << "w: " << (2*w2n+1) << " wp: " << (2*w3n+1) << " W: " << Wn << " " << " z1: " << s1 << " z2: " << s2 << "  "
                            << "   n1: " << n1 << " n2: " << n2 << "   "
                            << chi.real() << " " << chi.imag() << "   " << chi0.real() << " " << chi0.imag() << std::endl;

                    chi_irr_str << "w: " << (2*w2n+1) << " wp: " << (2*w3n+1) << " W: " << Wn << " " << " z1: " << s1 << " z2: " << s2 << "  "
                                << "   n1: " << n1 << " n2: " << n2 << "   "
                                << chi_irr.real() << " " << chi_irr.imag() << std::endl;

                    gamma_str << "w: " << (2*w2n+1) << " wp: " << (2*w3n+1) << " W: " << Wn << " " << " z1: " << s1 << " z2: " << s2 << "  "
                              << "   n1: " << n1 << " n2: " << n2 << "   "
                              << gamma.real() << " " << gamma.imag()  << std::endl;
                  }
                }//Wn
              }//w3n
            }//w2n
          }//s2
    {
      alps::hdf5::archive oar(basename+".out.h5", alps::hdf5::archive::WRITE);
      oar << alps::make_pvp("/simulation/results/gammaw",gammaw);
    }

  std::cout << "done." << std::endl;
  }//MEASURE_g2w

  return 0;
}
