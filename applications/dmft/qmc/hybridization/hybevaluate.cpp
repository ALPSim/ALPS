/****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2012 by Hartmut Hafermann <hafermann@cpht.polytechnique.fr>,
 *                       Emanuel Gull <gull@pks.mpg.de>
 *
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

#include "hyb.hpp"
#include "hybevaluate.hpp"
#include <alps/config.h>
#include <boost/cstdint.hpp>

void evaluate_basics(const alps::results_type<hybridization>::type &results,
                     const alps::parameters_type<hybridization>::type &parms,
                     alps::hdf5::archive &solver_output){

  std::size_t n_orbitals=parms["N_ORBITALS"];
  std::size_t N_meas=parms["N_MEAS"];
  double beta=parms["BETA"];

  boost::uint64_t sweeps = results["Sign"].count()+(boost::uint64_t)parms["THERMALIZATION"];

  if(parms["TEXT_OUTPUT"]|false){
    std::ofstream sim_file("simulation.dat");
    sim_file << "simulation details:" << std::endl;
    sim_file << "average sign: " << results["Sign"].mean<double>() << std::endl;
    sim_file << "total number of Monte Carlo updates: " << sweeps*(int)parms["N_MEAS"] << std::endl;
    sim_file << "total number of Monte Carlo measurements: " << results["Sign"].count() << std::endl;
    if(parms["MEASURE_time"]|true)
      sim_file << "total number of imaginary time measurements: " << results["Sign"].count()*N_meas << std::endl;
    sim_file << "number thermalization sweeps: " << parms["THERMALIZATION"] << std::endl;
    sim_file << "inverse temperature: " << beta << std::endl;
    sim_file << "perturbation order:" << std::endl;
    for(std::size_t i=0;i<n_orbitals;++i){//replace Green function endpoints by corresponding densities
      std::stringstream order_name; order_name<<"order_"<<i;
      double order=results[order_name.str()].mean<double>();
      sim_file << "orbital " << i << ": " << order << std::endl;
    }
    {
      int tot_acc=0,cur_prec = sim_file.precision();
      for (int i=0;i<nacc.size();i++) tot_acc += nacc[i];
      sim_file << std::endl << "|------ Simulation details (master only) after " << sweeps << " sweeps ------|" << std::endl;
      sim_file << "  Total acceptance rate = " << std::setprecision(2) << std::fixed;
      sim_file << (((double)tot_acc)/(sweep_count*N_meas))*100 << "%" << std::endl;
      sim_file << "  Individual acceptance rate for update " << std::endl;
      for (int i=0;i<nacc.size();i++) {
          sim_file << "     " << update_type[i] << " = ";
          sim_file << std::setprecision(2) << std::fixed << (((double)nacc[i])/(sweep_count*N_meas))*100 << "%";
          sim_file << " (proposal rate = ";
          sim_file << std::setprecision(2) << std::fixed << (((double)nprop[i])/(sweep_count*N_meas))*100 << "%)" << std::endl;
      }
      sim_file << "|-----------------------------------------------------------------|" << std::endl;
    }
    sim_file.close();

    std::ofstream obs_file("observables.dat");//equal-time correlators
    for(std::size_t i=0;i<n_orbitals;++i){//replace Green function endpoints by corresponding densities
      std::stringstream density_name; density_name<<"density_"<<i;
      double density=results[density_name.str()].mean<double>();
      obs_file << "n" << i << "=" << density << ";" << std::endl;
      if(parms["MEASURE_nn"]|false){
        for(std::size_t j=0;j<i;++j){
          std::stringstream nn_name; nn_name<<"nn_"<<i<<"_"<<j;
          double nn=results[nn_name.str()].mean<double>();
          obs_file << "n" << i << "n" << j << "=" << nn << ";" << std::endl;
        }
      }
    }
    obs_file.close();

    std::ofstream order_file("orders.dat");
    std::vector<std::vector<double> > order_histogram(n_orbitals);
    std::vector<std::vector<double> > order_histogram_err(n_orbitals);
    for(std::size_t j=0;j<n_orbitals;++j){
      std::stringstream order_name; order_name<<"order_histogram_"<<j;
      order_histogram[j]=results[order_name.str()].mean<std::vector<double> >();
      order_histogram_err[j]=results[order_name.str()].error<std::vector<double> >();
    }
    for(std::size_t j=0;j<n_orbitals;++j){
      for(std::size_t k=0;k<order_histogram[0].size();++k){
        order_file<<k;
        order_file<<" "<<order_histogram[j][k]<<" "<<order_histogram_err[j][k]<<std::endl;
      }
      order_file<<std::endl;
    }
    order_file.close();
  }//text_output

}


void evaluate_time(const alps::results_type<hybridization>::type &results,
                   const alps::parameters_type<hybridization>::type &parms,
                   alps::hdf5::archive &solver_output){

  if(!(parms["MEASURE_time"]|true)) return;
  std::size_t N_t = parms["N_TAU"];
  double beta = parms["BETA"];
  std::size_t n_orbitals = parms["N_ORBITALS"];
  std::size_t n_sites    = 1;
  bool accurate = parms["ACCURATE_COVARIANCE"]|true;

  //Imaginary time Green function
  itime_green_function_t G_tau(N_t+1, n_sites, n_orbitals);
  itime_green_function_t F_tau(N_t+1, n_sites, n_orbitals);
  for(std::size_t i=0;i<n_orbitals;++i){
    std::stringstream g_name; g_name<<"g_"<<i;
    std::stringstream f_name; f_name<<"f_"<<i;
    std::vector<double> G=results[g_name.str()].mean<std::vector<double> >();
    std::vector<double> F=results[f_name.str()].mean<std::vector<double> >();
    for(std::size_t t=0;t<N_t+1;++t){
      G_tau(t,0,0,i)=G[t];
      F_tau(t,0,0,i)=F[t];
    }
    G_tau(0,0,0,i)*=2.; //first and last bin
    G_tau(N_t,0,0,i)*=2.; //have half the size
    F_tau(0,0,0,i)*=2.; //first and last bin
    F_tau(N_t,0,0,i)*=2.; //have half the size
  }

  for(std::size_t i=0;i<n_orbitals;++i){//replace Green function endpoints by corresponding densities
    std::stringstream density_name; density_name<<"density_"<<i;
    double density=results[density_name.str()].mean<double>();
    G_tau(0,0,0,i)  =-1.*(1-density);
    G_tau(N_t,0,0,i)=-1.*(density);
  }

  //store in hdf5
  G_tau.write_hdf5(solver_output, boost::lexical_cast<std::string>(parms["BASEPATH"]|"")+"/G_tau");
  F_tau.write_hdf5(solver_output, boost::lexical_cast<std::string>(parms["BASEPATH"]|"")+"/F_tau");

  // ERROR
  for(std::size_t i=0; i<n_orbitals; i++){
     std::stringstream density_name; density_name<<"density_"<<i;
     double density_err=results[density_name.str()].error<double>();
     std::vector<double> err(N_t+1);
     std::stringstream g_name; g_name<<"g_"<<i;
     err = results[g_name.str()].error<std::vector<double> >();
     err[0] = err[N_t] = density_err;
     std::stringstream data_path;
     data_path << boost::lexical_cast<std::string>(parms["BASEPATH"]|"")+"/G_tau/"<<i<< "/mean/error";
     solver_output<<alps::make_pvp(data_path.str(),err);
     g_name.str(""); g_name<<"f_"<<i;
     err = results[g_name.str()].error<std::vector<double> >();
     data_path.str("");
     data_path << boost::lexical_cast<std::string>(parms["BASEPATH"]|"")+"/F_tau/"<<i<< "/mean/error";
     solver_output<<alps::make_pvp(data_path.str(),err);
  }
    
  //COVARIANCE
  for(std::size_t i=0; i<n_orbitals; i++){
    boost::numeric::ublas::matrix<double> cov(N_t+1, N_t+1);
    std::stringstream g_name; g_name<<"g_"<<i;

#ifdef ALPS_NGS_USE_NEW_ALEA
    {
      typedef alps::accumulator::RealVectorObservable::result_type result_type;
      result_type const & arg = results[g_name.str()].extract<result_type>();
      if (accurate)
        cov = arg.accurate_covariance(arg);
      else
        cov = arg.covariance(arg);
    }
#else
    if (accurate)
      cov=results[g_name.str()].accurate_covariance<std::vector<double> >(results[g_name.str()]);
    else
      cov=results[g_name.str()].covariance<std::vector<double> >(results[g_name.str()]);
#endif

    std::vector<double> data((N_t+1)*(N_t+1));
    for(std::size_t t1=0; t1<=N_t; t1++)
      for(std::size_t t2=0; t2<=N_t; t2++)
        data[t1*(N_t+1)+t2]=cov(t1,t2);

    std::stringstream data_path;
    data_path << boost::lexical_cast<std::string>(parms["BASEPATH"]|"")+"/G_tau/"<<i<< "/mean/covariance";

    solver_output<<alps::make_pvp(data_path.str(), data);
    g_name.str(""); g_name<<"f_"<<i;

#ifdef ALPS_NGS_USE_NEW_ALEA
    {
      typedef alps::accumulator::RealVectorObservable::result_type result_type;
      result_type const & arg = results[g_name.str()].extract<result_type>();
      if (accurate)
        cov = arg.accurate_covariance(arg);
      else
        cov = arg.covariance(arg);
    }
#else
    if (accurate)
      cov=results[g_name.str()].accurate_covariance<std::vector<double> >(results[g_name.str()]);
    else
      cov=results[g_name.str()].covariance<std::vector<double> >(results[g_name.str()]);      
#endif

    for(std::size_t t1=0; t1<=N_t; t1++)
      for(std::size_t t2=0; t2<=N_t; t2++)
         data[t1*(N_t+1)+t2]=cov(t1,t2);
      
    data_path.str("");
    data_path << boost::lexical_cast<std::string>(parms["BASEPATH"]|"")+"/F_tau/"<<i<< "/mean/covariance";
      
    solver_output<<alps::make_pvp(data_path.str(), data);
  }

  if(parms["TEXT_OUTPUT"]|0){
    std::ofstream G_file("Gt.dat");
    for(std::size_t t=0;t<=N_t;++t){
      G_file<<beta*t/N_t;
      for(std::size_t j=0;j<n_orbitals;++j){
        G_file<<" "<<G_tau(t,0,0,j);
      }
      G_file<<std::endl;
    }
    G_file.close();
     std::ofstream F_file("Ft.dat");
    for(std::size_t t=0;t<=N_t;++t){
      F_file<<beta*t/N_t;
      for(std::size_t j=0;j<n_orbitals;++j){
        F_file<<" "<<F_tau(t,0,0,j);
      }
      F_file<<std::endl;
    }
    F_file.close();
  }
}


void evaluate_freq(const alps::results_type<hybridization>::type &results,
                   const alps::parameters_type<hybridization>::type &parms,
                   alps::hdf5::archive &solver_output){

  if(!(parms["MEASURE_freq"]|false)) return;
  //evaluate Matsubara Green's function and self-energy
  double beta = parms["BETA"];
  std::size_t N_w=parms["N_MATSUBARA"];
  std::size_t n_orbitals=parms["N_ORBITALS"];
  std::size_t n_sites = 1;

  matsubara_green_function_t G_omega(N_w, n_sites, n_orbitals);
  matsubara_green_function_t F_omega(N_w, n_sites, n_orbitals);
  matsubara_green_function_t S_omega(N_w, n_sites, n_orbitals);
  for(std::size_t i=0;i<n_orbitals;++i){
    std::stringstream gw_re_name; gw_re_name<<"gw_re_"<<i;
    std::stringstream gw_im_name; gw_im_name<<"gw_im_"<<i;
    std::stringstream fw_re_name; fw_re_name<<"fw_re_"<<i;
    std::stringstream fw_im_name; fw_im_name<<"fw_im_"<<i;
    std::vector<double> Gw_re=results[gw_re_name.str()].mean<std::vector<double> >();
    std::vector<double> Gw_im=results[gw_im_name.str()].mean<std::vector<double> >();
    std::vector<double> Fw_re=results[fw_re_name.str()].mean<std::vector<double> >();
    std::vector<double> Fw_im=results[fw_im_name.str()].mean<std::vector<double> >();
    for(std::size_t w=0;w<N_w;++w){
      std::complex<double> G(Gw_re[w],Gw_im[w]);
      std::complex<double> F(Fw_re[w],Fw_im[w]);
      G_omega(w,0,0,i)=G;
      F_omega(w,0,0,i)=F;
      S_omega(w,0,0,i)=F/G;
    }
  }

  //store in hdf5
  G_omega.write_hdf5(solver_output, boost::lexical_cast<std::string>(parms["BASEPATH"]|"")+"/G_omega");
  F_omega.write_hdf5(solver_output, boost::lexical_cast<std::string>(parms["BASEPATH"]|"")+"/F_omega");
  S_omega.write_hdf5(solver_output, boost::lexical_cast<std::string>(parms["BASEPATH"]|"")+"/S_omega");

  // ERROR
  for(std::size_t i=0; i<n_orbitals; i++){
    std::vector<double> err_g_re(N_w),err_g_im(N_w),err_f_re(N_w),err_f_im(N_w);
    std::vector<std::complex<double> > err(N_w);
    std::stringstream g_name; g_name<<"gw_re_"<<i;
    err_g_re = results[g_name.str()].error<std::vector<double> >();
    g_name.str("");g_name<<"gw_im_"<<i;
    err_g_im = results[g_name.str()].error<std::vector<double> >();
    for (int k=0;k<err.size();k++) err[k] = std::complex<double>(err_g_re[k],err_g_im[k]);
    std::stringstream data_path;
    data_path << boost::lexical_cast<std::string>(parms["BASEPATH"]|"")+"/G_omega/"<<i<< "/mean/error";
    solver_output<<alps::make_pvp(data_path.str(),err);
    g_name.str(""); g_name<<"fw_re_"<<i;
    err_f_re = results[g_name.str()].error<std::vector<double> >();
    g_name.str("");g_name<<"fw_im_"<<i;
    err_f_im = results[g_name.str()].error<std::vector<double> >();
    for (int k=0;k<err.size();k++) err[k] = std::complex<double>(err_f_re[k],err_f_im[k]);
    data_path.str("");
    data_path << boost::lexical_cast<std::string>(parms["BASEPATH"]|"")+"/F_omega/"<<i<< "/mean/error";
    solver_output<<alps::make_pvp(data_path.str(),err);
  }

    
  std::ofstream Gw_file("Gw.dat");
  for(std::size_t n=0;n<N_w;++n){
    Gw_file<<(2.*n+1)*M_PI/beta;
    for(std::size_t j=0;j<n_orbitals;++j){
      Gw_file<<" "<<G_omega(n,0,0,j).real()<<" "<<G_omega(n,0,0,j).imag();
    }
    Gw_file<<std::endl;
  }
  Gw_file.close();

  std::ofstream Fw_file("Fw.dat");
  for(std::size_t n=0;n<N_w;++n){
    Fw_file<<(2.*n+1)*M_PI/beta;
    for(std::size_t j=0;j<n_orbitals;++j){
      Fw_file<<" "<<F_omega(n,0,0,j).real()<<" "<<F_omega(n,0,0,j).imag();
    }
    Fw_file<<std::endl;
  }
  Fw_file.close();

  std::ofstream Sw_file("Sw.dat");
  for(std::size_t n=0;n<N_w;++n){
    Sw_file<<(2.*n+1)*M_PI/beta;
    for(std::size_t j=0;j<n_orbitals;++j){
      Sw_file<<" "<<S_omega(n,0,0,j).real()<<" "<<S_omega(n,0,0,j).imag();
    }
    Sw_file<<std::endl;
  }
  Sw_file.close();
}


void evaluate_legendre(const alps::results_type<hybridization>::type &results,
                       const alps::parameters_type<hybridization>::type &parms,
                       alps::hdf5::archive &solver_output){

  if(!(parms["MEASURE_legendre"]|0)) return;
  double beta = parms["BETA"];
  std::size_t N_l=parms["N_LEGENDRE"];
  std::size_t N_w=parms["N_MATSUBARA"];
  std::size_t N_t = parms["N_TAU"];
  std::size_t n_orbitals = parms["N_ORBITALS"];
  std::size_t n_sites = 1;

  //Legendre Green function (evaluated in Matsubara)
  matsubara_green_function_t G_l_omega(N_w, n_sites, n_orbitals);
  matsubara_green_function_t F_l_omega(N_w, n_sites, n_orbitals);
  matsubara_green_function_t S_l_omega(N_w, n_sites, n_orbitals);
  //Legendre Green function (evaluated in imaginary time)
  itime_green_function_t G_l_tau(N_t+1, n_sites, n_orbitals);
  itime_green_function_t F_l_tau(N_t+1, n_sites, n_orbitals);
  std::vector<std::vector<std::vector<double> > >gc_conv(2, std::vector<std::vector<double> >(n_orbitals, std::vector<double>(N_l,0.)));
  std::vector<std::vector<std::vector<double> > >fc_conv(2, std::vector<std::vector<double> >(n_orbitals, std::vector<double>(N_l,0.)));
  for(std::size_t i=0;i<n_orbitals;++i){
    std::stringstream gl_name; gl_name<<"gl_"<<i;
    std::stringstream fl_name; fl_name<<"fl_"<<i;
    std::vector<double> Gl=results[gl_name.str()].mean<std::vector<double> >();
    std::vector<double> Fl=results[fl_name.str()].mean<std::vector<double> >();
    for(std::size_t wn=0; wn<N_w; ++wn){
      G_l_omega(wn,0,0,i)=0.;
      F_l_omega(wn,0,0,i)=0.;
      for(std::size_t l=0; l<N_l; ++l){
        G_l_omega(wn,0,0,i)+=t(wn,l)*sqrt(2.*l+1)*Gl[l]; //sqrt(2l+1) has been omitted in the measurement
        F_l_omega(wn,0,0,i)+=t(wn,l)*sqrt(2.*l+1)*Fl[l];
      }
      S_l_omega(wn,0,0,i)=F_l_omega(wn,0,0,i)/G_l_omega(wn,0,0,i);
    }
    //Imaginary time Green function from Legendre
    for(std::size_t t=0;t<N_t+1;++t){
      double tau=t*beta/N_t;
      G_l_tau(t,0,0,i)=0.;
      F_l_tau(t,0,0,i)=0.;
      double x=2.0*tau/beta-1.0;
      double pl_2=1; double pl_1=x; double legendre_p;
      for(std::size_t l=0;l<N_l;++l){
        if(l==0) legendre_p=1;
        else if(l==1) legendre_p=x;
        else{
          legendre_p=((2*l-1)*x*pl_1-(l-1)*pl_2)/static_cast<double>(l);//l
          pl_2=pl_1; //l-2
          pl_1=legendre_p; //l-1
        }
        G_l_tau(t,0,0,i)+=(2.*l+1)*Gl[l]*legendre_p/beta;
        F_l_tau(t,0,0,i)+=(2.*l+1)*Fl[l]*legendre_p/beta;
        if(t==N_t/2){
          gc_conv[0][i][l]=G_l_tau(t,0,0,i);
          fc_conv[0][i][l]=F_l_tau(t,0,0,i);
        }
        if(t==N_t){
          gc_conv[1][i][l]=G_l_tau(t,0,0,i);
          fc_conv[1][i][l]=F_l_tau(t,0,0,i);
        }
      }
    }
  }//i

  //store in hdf5
  G_l_omega.write_hdf5(solver_output, boost::lexical_cast<std::string>(parms["BASEPATH"]|"")+"/G_l_omega");
  F_l_omega.write_hdf5(solver_output, boost::lexical_cast<std::string>(parms["BASEPATH"]|"")+"/F_l_omega");
  S_l_omega.write_hdf5(solver_output, boost::lexical_cast<std::string>(parms["BASEPATH"]|"")+"/S_l_omega");
  G_l_tau.write_hdf5(solver_output, boost::lexical_cast<std::string>(parms["BASEPATH"]|"")+"/G_l_tau");
  F_l_tau.write_hdf5(solver_output, boost::lexical_cast<std::string>(parms["BASEPATH"]|"")+"/F_l_tau");

  if(parms["TEXT_OUTPUT"]|false){
    std::ofstream gc_str("Gl_conv.dat");
    std::ofstream fc_str("Fl_conv.dat");
    gc_str << "#lc";
    fc_str << "#lc";
    for(std::size_t i=0;i<n_orbitals;++i){
      gc_str << " Gt_" << i << "(tau=beta/2,lc) Gt_" << i << "(tau=beta)";
      fc_str << " Ft_" << i << "(tau=beta/2,lc) Ft_" << i << "(tau=beta)";
    }
    gc_str << std::endl;
    fc_str << std::endl;
    for(std::size_t l=0;l<N_l;++l){
      gc_str << l;
      fc_str << l;
      for(std::size_t i=0;i<n_orbitals;++i){
        gc_str << " " << gc_conv[0][i][l] << " " << gc_conv[1][i][l];
        fc_str << " " << fc_conv[0][i][l] << " " << fc_conv[1][i][l];
      }
      gc_str << std::endl;
      fc_str << std::endl;
    }
    gc_str.close();
    fc_str.close();
    std::ofstream Gtl_file("Gtl.dat");
    for(std::size_t t=0;t<=N_t;++t){
      Gtl_file<<beta*t/N_t;
      for(std::size_t j=0;j<n_orbitals;++j){
        Gtl_file<<" " <<G_l_tau(t,0,0,j);
      }
      Gtl_file<<std::endl;
    }
    Gtl_file.close();
    std::ofstream Ftl_file("Ftl.dat");
    for(std::size_t t=0;t<=N_t;++t){
      Ftl_file<<beta*t/N_t;
      for(std::size_t j=0;j<n_orbitals;++j){
        Ftl_file<<" " <<F_l_tau(t,0,0,j);
      }
      Ftl_file<<std::endl;
    }
    Ftl_file.close();
    std::ofstream Gw_file("Gwl.dat");
    for(std::size_t t=0;t<N_w;++t){
      Gw_file<<(2.*t+1)*M_PI/beta;
      for(std::size_t j=0;j<n_orbitals;++j){
        Gw_file<<" "<<G_l_omega(t,0,0,j).real()<<" "<<G_l_omega(t,0,0,j).imag();
      }
      Gw_file<<std::endl;
    }
    Gw_file.close();
    std::ofstream Fw_file("Fwl.dat");
    for(std::size_t t=0;t<N_w;++t){
      Fw_file<<(2.*t+1)*M_PI/beta;
      for(std::size_t j=0;j<n_orbitals;++j){
        Fw_file<<" "<<F_l_omega(t,0,0,j).real()<<" "<<F_l_omega(t,0,0,j).imag();
      }
      Fw_file<<std::endl;
    }
    Fw_file.close();
    std::ofstream Sw_file("Swl.dat");
    for(std::size_t t=0;t<N_w;++t){
      Sw_file<<(2.*t+1)*M_PI/beta;
      for(std::size_t j=0;j<n_orbitals;++j){
        Sw_file<<" "<<S_l_omega(t,0,0,j).real()<<" "<<S_l_omega(t,0,0,j).imag();
      }
      Sw_file<<std::endl;
    }
    Sw_file.close();
  }

}


void evaluate_nnt(const alps::results_type<hybridization>::type &results,
                  const alps::parameters_type<hybridization>::type &parms,
                  alps::hdf5::archive &solver_output){

  if(!(parms["MEASURE_nnt"]|0)) return;

  std::size_t N_nn=parms["N_nn"];
  std::size_t n_orbitals=parms["N_ORBITALS"];
  double beta=parms["BETA"];

  std::vector<std::vector<double> > nnt(n_orbitals*(n_orbitals+1)/2);
  int pos=0;
  for(std::size_t i=0;i<n_orbitals;++i){
    for(std::size_t j=0;j<=i;++j){
      std::stringstream nnt_name; nnt_name<<"nnt_"<<i<<"_"<<j;
      nnt[pos]=results[nnt_name.str()].mean<std::vector<double> >();
      solver_output<<alps::make_pvp(nnt_name.str(), nnt[pos++]);
    }
  }

  if(parms["TEXT_OUTPUT"]|0){
    std::ofstream nnt_file("nnt.dat");
    nnt_file << "#tau";
    for(std::size_t i=0;i<n_orbitals;++i)
      for(std::size_t j=0;j<=i;++j)
        nnt_file<<" nnt_"<<i<<j;
    nnt_file << std::endl;
    for(std::size_t n=0;n<=N_nn;++n){
      pos=0;
      double tau=n*beta/N_nn;
      nnt_file << tau;
      for(std::size_t i=0;i<n_orbitals;++i)
        for(std::size_t j=0;j<=i;++j)
          nnt_file << " " << nnt[pos++][n];
      nnt_file << std::endl;
    }
    nnt_file.close();
  }
}


void evaluate_nnw(const alps::results_type<hybridization>::type &results,
                  const alps::parameters_type<hybridization>::type &parms,
                  alps::hdf5::archive &solver_output){

  if(!(parms["MEASURE_nnw"]|0)) return;

  std::size_t N_W=parms["N_W"];
  std::size_t n_orbitals=parms["N_ORBITALS"];
  double beta=parms["BETA"];

  std::vector<std::vector<double> > nnw_re(n_orbitals*(n_orbitals+1)/2);
  int pos=0;
  for(std::size_t i=0;i<n_orbitals;++i){
    for(std::size_t j=0;j<=i;++j){
      std::stringstream nnw_re_name; nnw_re_name<<"nnw_re_"<<i<<"_"<<j;
      nnw_re[pos]=results[nnw_re_name.str()].mean<std::vector<double> >();
      solver_output<<alps::make_pvp(nnw_re_name.str(), nnw_re[pos++]);
    }
  }

  if(parms["TEXT_OUTPUT"]|0){
    std::ofstream nnw_file("nnw.dat");
    nnw_file << "#w";
    for(std::size_t i=0;i<n_orbitals;++i)
      for(std::size_t j=0;j<=i;++j)
        nnw_file<<" nnw_"<<i<<j;
    nnw_file << std::endl;
    for(std::size_t m=0;m<N_W;++m){
      pos=0;
      double wm=2*m*M_PI/beta;
      nnw_file << wm;
      for(std::size_t i=0;i<n_orbitals;++i)
        for(std::size_t j=0;j<=i;++j)
          nnw_file << " " << nnw_re[pos++][m];
      nnw_file << std::endl;
    }
    nnw_file.close();
  }
}


void evaluate_sector_statistics(const alps::results_type<hybridization>::type &results,
                                const alps::parameters_type<hybridization>::type &parms,
                                alps::hdf5::archive &solver_output){

  if(!(parms["MEASURE_sector_statistics"]|0)) return;

  std::size_t n_orbitals=parms["N_ORBITALS"];
  std::ofstream stat_file("sector_statistics.dat");
  stat_file << "#state |n_1={0,1} n_2={0,1} ...> n_i={0,1}: orbital i {empty,occupied}" << std::endl;
  stat_file << "#rel weight (in %)" << std::endl;
  int n_states=1<<n_orbitals;
  std::vector<double> sector_statistics=results["sector_statistics"].mean<std::vector<double> >();
  for(int n=0;n<n_states;++n){
    std::stringstream state; state << "|";
    std::size_t i=0; int m=n;
    while(m>=0 && i<n_orbitals){
      int rem=m%2; ++i; m/=2;
      state << rem;
    }
    state << ">";
    stat_file << n << "\t" << sector_statistics[n]*100. << "\t" << state.str() << std::endl;
  }
  stat_file.close();
}


void evaluate_2p(const alps::results_type<hybridization>::type &results,
                 const alps::parameters_type<hybridization>::type &parms,
                 alps::hdf5::archive &solver_output){
  //write two-particle functions to text file if desired
  //compute the vertex function if needed;
  bool MEASURE_g2w=parms["MEASURE_g2w"]|false;
  bool MEASURE_h2w=parms["MEASURE_h2w"]|false;

  if(!(MEASURE_g2w || MEASURE_h2w)) return;

  //int N_w = parms["N_MATSUBARA"]|0;
  std::size_t N_W = parms["N_W"]|0;
  std::size_t N_w2 = parms["N_w2"]|0;
  std::size_t n_orbitals = parms["N_ORBITALS"]; //number of orbitals
  bool text_output = parms["TEXT_OUTPUT"]|false;
  bool COMPUTE_VERTEX = parms["COMPUTE_VERTEX"]|false;
  double beta=parms["BETA"];

  std::ofstream g2w_str;
  std::ofstream h2w_str;
  std::ofstream gam_str;
  if(text_output){
    if(MEASURE_g2w) g2w_str.open("g2w.dat");
    if(MEASURE_h2w) h2w_str.open("h2w.dat");
    if(COMPUTE_VERTEX) gam_str.open("gammaw.dat");
  }
  std::vector<double> g2w_re;
  std::vector<double> g2w_im;
  std::vector<double> h2w_re;
  std::vector<double> h2w_im;
  std::vector<std::vector<double> >gw_re(n_orbitals);
  std::vector<std::vector<double> >gw_im(n_orbitals);
  std::vector<std::vector<double> >fw_re(n_orbitals);
  std::vector<std::vector<double> >fw_im(n_orbitals);
  std::vector<std::complex<double> > vertex;

  if(COMPUTE_VERTEX){
    vertex.resize(N_w2*N_w2*N_W);
    for(std::size_t i=0;i<n_orbitals;++i){
      std::stringstream gw_re_name; gw_re_name<<"gw_re_"<<i;
      std::stringstream gw_im_name; gw_im_name<<"gw_im_"<<i;
      std::stringstream fw_re_name; fw_re_name<<"fw_re_"<<i;
      std::stringstream fw_im_name; fw_im_name<<"fw_im_"<<i;
      gw_re[i]=results[gw_re_name.str()].mean<std::vector<double> >();
      gw_im[i]=results[gw_im_name.str()].mean<std::vector<double> >();

      if(MEASURE_h2w){
        fw_re[i]=results[fw_re_name.str()].mean<std::vector<double> >();
        fw_im[i]=results[fw_im_name.str()].mean<std::vector<double> >();
      }
    }
  }

  std::complex<double> g1,g2,g3,g4,f1;
  std::complex<double> gg,g2w,h2w,g2w_con,gamma;

  for(std::size_t i=0;i<n_orbitals;++i){
    for(std::size_t j=0;j<=i;++j){
      if(MEASURE_g2w){
        std::stringstream g2w_re_name; g2w_re_name<<"g2w_re_"<<i<<"_"<<j;
        std::stringstream g2w_im_name; g2w_im_name<<"g2w_im_"<<i<<"_"<<j;
        g2w_re=results[g2w_re_name.str()].mean<std::vector<double> >();
        g2w_im=results[g2w_im_name.str()].mean<std::vector<double> >();
      }
      if(MEASURE_h2w){
        std::stringstream h2w_re_name; h2w_re_name<<"h2w_re_"<<i<<"_"<<j;
        std::stringstream h2w_im_name; h2w_im_name<<"h2w_im_"<<i<<"_"<<j;
        h2w_re=results[h2w_re_name.str()].mean<std::vector<double> >();
        h2w_im=results[h2w_im_name.str()].mean<std::vector<double> >();
      }
      int N_wh=N_w2/2;
      for(int w2n=-N_wh;w2n<N_wh;++w2n){
        if(COMPUTE_VERTEX) g2=( w2n<0 ? conj(std::complex<double>(gw_re[i][-w2n-1],gw_im[i][-w2n-1])) : std::complex<double>(gw_re[i][w2n],gw_im[i][w2n]) );//i;w2
        for(int w3n=-N_wh;w3n<N_wh;++w3n){
          if(COMPUTE_VERTEX) g3=( w3n<0 ? conj(std::complex<double>(gw_re[j][-w3n-1],gw_im[j][-w3n-1])) : std::complex<double>(gw_re[j][w3n],gw_im[j][w3n]) );//j;w3
          for(std::size_t Wm=0;Wm<N_W;++Wm){
            int index=Wm*N_w2*N_w2 + (w2n+N_wh)*N_w2 + (w3n+N_wh);
            int w1n=w2n+Wm; int w4n=w3n+Wm;

            if(MEASURE_g2w) g2w = std::complex<double>(g2w_re[index],g2w_im[index]);
            if(MEASURE_h2w) h2w = std::complex<double>(h2w_re[index],h2w_im[index]);

            if(COMPUTE_VERTEX){
              g1=( w1n<0 ? conj(std::complex<double>(gw_re[i][-w1n-1],gw_im[i][-w1n-1])) : std::complex<double>(gw_re[i][w1n],gw_im[i][w1n]) );//i;w1
              g4=( w4n<0 ? conj(std::complex<double>(gw_re[j][-w4n-1],gw_im[j][-w4n-1])) : std::complex<double>(gw_re[j][w4n],gw_im[j][w4n]) );//j;w4
              if(MEASURE_h2w)
                f1=( w1n<0 ? conj(std::complex<double>(fw_re[i][-w1n-1],fw_im[i][-w1n-1])) : std::complex<double>(fw_re[i][w1n],fw_im[i][w1n]) );//i;w1

              gg=0.;
              if(w1n==w2n) gg+=beta*g1*g3;
              if(w2n==w3n && i==j) gg-=beta*g1*g3;

              //evaluate the vertex depending on what has been measured
              if(MEASURE_g2w && MEASURE_h2w) g2w_con=g1*h2w - f1*g2w; //usually most accurate
              else if(MEASURE_h2w)           g2w_con=(g1*h2w-f1*gg)/(1.0+f1); //somewhat less accurate; g2w not needed
              else                           g2w_con=g2w-gg; //straightforward evaluation, least accurate; h2w not needed

              gamma=g2w_con/(g1*g2*g3*g4);
              vertex[index]=gamma;
            }
            if(text_output){

              if(MEASURE_g2w)    g2w_str << "w: " << 2*w2n+1 << " wp: " << 2*w3n+1 << " W: " << Wm << " " << " i: " << i << " j: " << j << "  "
                << g2w.real() << " " << g2w.imag() << std::endl;
              if(MEASURE_h2w)    h2w_str << "w: " << 2*w2n+1 << " wp: " << 2*w3n+1 << " W: " << Wm << " " << " i: " << i << " j: " << j << "  "
                << h2w.real() << " " << h2w.imag() << std::endl;
              if(COMPUTE_VERTEX) gam_str << "w: " << 2*w2n+1 << " wp: " << 2*w3n+1 << " W: " << Wm << " " << " i: " << i << " j: " << j << "  "
                << gamma.real() << " " << gamma.imag() << std::endl;
            }
          }//Wm
        }//w3n
      }//w2n
      //write to hdf5
      if(COMPUTE_VERTEX){
        std::stringstream data_path; data_path<<"vertex_"<<i<<"_"<<j;
        solver_output<<alps::make_pvp(data_path.str(), vertex);
      }
    }//j
  }//i
}

