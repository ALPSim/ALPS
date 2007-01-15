/*****************************************************************************
*
* ALPS DMFT Project
*
* Copyright (C) 2005-2006 by Philipp Werner <werner@comp-phys.org>,
*                            Emanuel Gull <gullc@comp-phys.org>
*                            Matthias Troyer <troyer@comp-phys.org>
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

/* $Id$ */

/// @file hirschfyesim.h
/// @brief the actual Hirsch-Fye simulation

#include "hirschfyesim.h"
#include "xml.h"
#include <alps/alea.h>
#include <cmath>
#include "../blas_classes/vector.h"
using namespace std;
using namespace alps;

/// compute largest element of a matrix
double norm_max(dense_matrix & m) {

  double max(0);
  for (dense_matrix::const_iterator1 it = m.begin1(); it != m.end1(); ++it)
    if (std::fabs((*it)) > max)
	  max = std::fabs((*it));
	  
  return max;

}

/// compute Green's function for given spin configuration from Green0
void update_from_zero(dense_matrix & Green, dense_matrix & Green0, std::vector<int> & spins, double l) {

  int N(spins.size());
  std::vector<double> e_ls(N);
  for(int i=0;i<N;++i){
  	e_ls[i]=exp(l*spins[i]);
  }
  // calculate a = 1 + (1-g)*(exp(v')-1) 
  dense_matrix a(N, N);
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
	  a(i,j) = -Green0(i,j)*(e_ls[j]-1);
	}
	a(i,i) += e_ls[i]; 
  }
  // solve a*Green=Green0
  // gesv solves A*X=B, input for B is Green0, output (=solution X) is Green   
  Green = Green0; 
  //std::cout<<"green0: "<<Green0<<std::endl;
  boost::numeric::bindings::lapack::gesv(a,Green);
  //std::cout<<"updatefz: a"<<a<<"Green0: "<<Green0<<" Green: "<<Green<<std::endl;
}


HirschFyeRun::HirschFyeRun(const alps::ProcessList& where,const alps::Parameters& p,int node)
  : alps::scheduler::MCRun(where,p,node),  
    sweeps(0),
    thermalization_sweeps(static_cast<int>(parms["THERMALIZATION"])),
    total_sweeps(static_cast<int>(parms["SWEEPS"])),
    beta(static_cast<double>(parms["BETA"])),
    N(static_cast<int>(parms["N"])),	
    u(static_cast<double>(parms["U"])),
	N_check(1000), // value is automatically adjusted during the simulation
	check_counter(0),
	tolerance(static_cast<double>(parms["TOLERANCE"])),	
	G_up(static_cast<int>(parms["N"])),
	Green0_up(static_cast<int>(parms["N"]), static_cast<int>(parms["N"])),
	Green_up(static_cast<int>(parms["N"]), static_cast<int>(parms["N"])),
	G_down(static_cast<int>(parms["N"])),
	Green0_down(static_cast<int>(parms["N"]), static_cast<int>(parms["N"])),
	Green_down(static_cast<int>(parms["N"]), static_cast<int>(parms["N"])),			
    spins(static_cast<int>(parms["N"])),
    start_time(time(NULL)),
    max_time(static_cast<int>(parms.value_or_default("MAX_TIME", 86400))),
    measure_fourpoint(static_cast<bool>(parms.value_or_default("MEASURE_FOURPOINT_FUNCTION", false))),
    fp_interval(static_cast<int>(parms.value_or_default("FOURPOINT_INTERVAL", 1000)))

{
  // initialize G0
  itime_green_function_t multiple_G0(N+1, 1, 2);
  multiple_vector_type G0;
  //outside HF solver we use sizeof(G0)=N+1, inside we only compute with N imag
  //time positions.
  G0.first.resize(N);
  G0.second.resize(N);
  std::istringstream in(parms["G0"]);
  read_itime(in,multiple_G0);
  G0=multiple_G0.to_multiple_vector();
  
  std::valarray<double> G0_up(G0.first.size());
  std::valarray<double> G0_down(G0.first.size());
  for (uint i=0; i<G0.first.size(); i++) {
    G0_up[i] = -G0.first[G0.first.size()-1-i];
    G0_down[i] = -G0.second[G0.second.size()-1-i];
    std::cout<<i<<"\t"<<G0_up[i]<<" "<<G0_down[i]<<std::endl;
  }
  // initialize lambda
  double tmp = exp(beta/N*u/2);
  lambda = log(tmp + sqrt(tmp*tmp-1));

  // initialize bath Green's function matrices
  green_matrix_from_vector(G0_up, Green0_up);
  green_matrix_from_vector(G0_down, Green0_down);
  //std::cout<<"Green0 down: "<<Green0_down<<std::endl;
  //std::cout<<"Green0 up: "<<Green0_up<<std::endl;
  // choose random spin configuration
  for (uint i=0; i<spins.size(); i++)
    spins[i] = (random_01()<0.5 ? 1 : -1);

 
  // initialize Green's function matrices
  update_from_zero(Green_up, Green0_up, spins, lambda);
  update_from_zero(Green_down, Green0_down, spins, -lambda);

 
  // create measurement objects
  measurements << RealVectorObservable("G_meas_up");
  measurements << RealVectorObservable("G_meas_down");  
  measurements << SimpleRealVectorObservable("G_fourpoint_uu");  
  measurements << SimpleRealVectorObservable("G_fourpoint_dd");  
  measurements << SimpleRealVectorObservable("G_fourpoint_ud");  
  measurements << SimpleRealVectorObservable("G_fourpoint_du");  
  
  //reset the measurment objects

  measurements.get<RealVectorObservable>("G_meas_up").reset(true);
  measurements.get<RealVectorObservable>("G_meas_down").reset(true);
  measurements.get<SimpleRealVectorObservable>("G_fourpoint_uu").reset(true);
  measurements.get<SimpleRealVectorObservable>("G_fourpoint_dd").reset(true);
  measurements.get<SimpleRealVectorObservable>("G_fourpoint_ud").reset(true);
  measurements.get<SimpleRealVectorObservable>("G_fourpoint_du").reset(true);
  
}


// in principle, need only save spins (and do update_from_zero after loading)
void HirschFyeRun::load(alps::IDump& dump) {
  dump >> sweeps;
  if(!where.empty()) {// skip reading the spins if we are just evaluating
    dump >> Green_up >> Green_down >> spins; 	 
  }
  else
    measurements.compact(); // and throw away time series
}


void HirschFyeRun::save(alps::ODump& dump) const
{
  dump << Green_up << Green_down << spins;
}


bool HirschFyeRun::is_thermalized() const
{
  return (sweeps >= thermalization_sweeps);
}

double HirschFyeRun::work_done() const
{
  if(time(NULL)-start_time> max_time){
    std::cout<<"we ran out of time!"<<std::endl<<"sweeps: "<<sweeps-thermalization_sweeps<<std::endl;
    return 1;
  }

  return (is_thermalized() ? (sweeps-thermalization_sweeps)/double(total_sweeps) : 0.);
}


void HirschFyeRun::dostep()
{

  // increment sweep count
  sweeps++;
  check_counter++;
  // do local update
  for (int i=0; i<N; i++){
    update_single_spin(random_01, Green_up, Green_down, spins, lambda);
	  }
  // calculate Greens function	
  //std::valarray<double> G_tmp_up(N), G_tmp_down(N);		
  green_vector_from_matrix(G_up, Green_up);
  green_vector_from_matrix(G_down,Green_down);	

  //std::cout<<"in dostep"<<std::endl;
  // measure Greens function
  if(is_thermalized()){
  	measurements.get<RealVectorObservable>("G_meas_up") << G_up;
  	measurements.get<RealVectorObservable>("G_meas_down") << G_down;
  }
#ifdef FOURPOINT
  if(is_thermalized() && measure_fourpoint && (sweeps %fp_interval==0)){
    int s=Green_up.size1();
    int s3=s*s*s;
    int s2=s*s;
    static blas::vector G_ijkl_uu(s*s*s, 0);
    static blas::vector G_ijkl_dd(s*s*s, 0);
    static blas::vector G_ijkl_du(s*s*s, 0);
    static blas::vector G_ijkl_ud(s*s*s, 0);
    static int fpsteps=0;
    fpsteps++;
    /*for(int i=0;i<N;++i){
      for(int j=0;j<N;++j){
        G_ij_up[N*i+j]=Green_up  (i,j); //G_up(tau_i - tau_j)
        G_ij_dn[N*i+j]=Green_down(i,j); //G_down(tau_i - tau_j)
       //std::cout<<i<<" "<<j<<"\t"<<G_ij_up[N*i+j]<<" "<<G_ij_dn[N*i+j]<<std::endl;
      }
    }*/
    //std::cout<<"size s is: "<<s<<std::endl;
      for(int j=0;j<s;++j){
        double guji=Green_up(j, 0);
        double gdji=Green_down(j, 0);
        for(int k=0;k<s;++k){
          double gujk=Green_up(j, k);
          double gdjk=Green_down(j, k);
          for(int l=0;l<s;++l){
          double gulk=Green_up(l, k);
          double gdlk=Green_down(l, k);
          double guli=Green_up(l, 0);
          double gdli=Green_down(l, 0);
            //std::cout<<i<<" "<<j<<" "<<k<<" "<<l<<std::endl;
            G_ijkl_uu(j*s2+k*s+l)+=guji*gulk-gujk*guli; //check review formula 154
            G_ijkl_dd(j*s2+k*s+l)+=gdji*gdlk-gdjk*gdli;
            G_ijkl_ud(j*s2+k*s+l)+=guji*gdlk;;
            G_ijkl_du(j*s2+k*s+l)+=gdji*gulk;
          }
        }
      }
    int meas_int=100;
    if(fpsteps% meas_int == 0){
      //std::cout<<G_ijkl_uu[0]<<" "<<G_ijkl_ud[0]<<" "<<G_ijkl_du[0]<<" "<<G_ijkl_dd[0]<<std::endl;
      //std::cout<<Green_up(0,0)<<" "<<Green_up(1,1)<<" "<<Green_up(1,0)<<" "<<Green_up(0,1)<<G_ijkl_uu[16+1]<<std::endl;
      std::valarray<double> G_ijkl_uu_va(s*s*s);
      std::valarray<double> G_ijkl_dd_va(s*s*s);
      std::valarray<double> G_ijkl_ud_va(s*s*s);
      std::valarray<double> G_ijkl_du_va(s*s*s);
      for(int i=0;i<s*s*s;++i){
        G_ijkl_uu_va[i]=G_ijkl_uu(i)/meas_int;
        G_ijkl_ud_va[i]=G_ijkl_ud(i)/meas_int;
        G_ijkl_du_va[i]=G_ijkl_du(i)/meas_int;
        G_ijkl_dd_va[i]=G_ijkl_dd(i)/meas_int;
      }

      measurements.get<SimpleRealVectorObservable>("G_fourpoint_uu")<<G_ijkl_uu_va;
      measurements.get<SimpleRealVectorObservable>("G_fourpoint_dd")<<G_ijkl_dd_va;
      measurements.get<SimpleRealVectorObservable>("G_fourpoint_ud")<<G_ijkl_ud_va;
      measurements.get<SimpleRealVectorObservable>("G_fourpoint_du")<<G_ijkl_du_va;
     
      G_ijkl_uu=0;
      G_ijkl_du=0;
      G_ijkl_ud=0;
      G_ijkl_dd=0;
    }
    //std::cout<<Green_up(0,0)<<" "<<Green_up(1,1)<<" "<<Green_up(1,0)<<" "<<Green_up(0,1)<<std::endl;
    //std::cout<<"mean of 0: "<<measurements.get<SimpleRealVectorObservable>("G_fourpoint_uu").mean()[0]<<" "<<measurements.get<SimpleRealVectorObservable>("G_fourpoint_uu").mean()[s3]<<std::endl;
  }
#endif
  //only for equilibration measurement
  if(is_thermalized() && sweeps % 10000 ==0){
    static std::ofstream equ_file("equilibration.dat");
    equ_file<<sweeps<<"\t"<<measurements.get<RealVectorObservable>("G_meas_up").mean()[0]<<"\t"<<measurements.get<RealVectorObservable>("G_meas_up").error()[0]<<std::endl;
    //std::cout<<G_ijkl_uu[0]<<" "<<G_ijkl_ud[0]<<" "<<G_ijkl_du[0]<<" "<<G_ijkl_dd[0]<<std::endl;
  }
	
  // check if precision has deteriorated
  if (check_counter%N_check == 0) {
  
    dense_matrix Green_diff_up(Green_up), Green_diff_down(Green_down);

    update_from_zero(Green_up, Green0_up, spins, lambda);
    update_from_zero(Green_down, Green0_down, spins, -lambda);

    Green_diff_up -= Green_up;
	Green_diff_down -= Green_down;

    if (norm_max(Green_diff_up) > tolerance || norm_max(Green_diff_down) > tolerance)
	  N_check /= 2;
	else 
	  N_check *= 2;
	  
	check_counter = 0;

  }

}


itime_green_function_t HirschFyeRun::get_result() const {

  multiple_vector_type g;
  alps::RealVectorObsevaluator G_meas_up=get_measurements()["G_meas_up"];
  alps::RealVectorObsevaluator G_meas_down=get_measurements()["G_meas_down"];  
  
  std::valarray<double> g_up_mean_vector(G_meas_up.mean());
  std::valarray<double> g_down_mean_vector(G_meas_down.mean());
  
  for(uint i=0;i<g_up_mean_vector.size();++i){
  	g.first.push_back(g_up_mean_vector[i]);
  	g.second.push_back(g_down_mean_vector[i]);
  }
  //also take into account G(beta)
  g.first.push_back(-1.-g_up_mean_vector[0]);
  g.second.push_back(-1.-g_down_mean_vector[0]);
  itime_green_function_t multiple_G(N+1,1,2);
  multiple_G.from_multiple_vector(g);
   return multiple_G;

}


///take a Green's function in vector form G(\tau-\tau') and create a matrix G(\tau, \tau') out of it
void HirschFyeRun::green_matrix_from_vector(const std::valarray<double> &green_vector, dense_matrix & green_matrix)const{
  green_matrix.clear();
	for (int i=0; i<N; i++) {
    for (int j=i; j<N; j++) {	  
      green_matrix(j,i) = green_vector[j-i];
    }
  }
  for (int i=0; i<N; i++) {
    for (int j=i+1; j<N; j++) {//consider antisymmetry (Fermions)
      green_matrix(i,j) -= green_vector[N-(j-i)];
 	}
  }
}
///take a Green's function in matrix form and transform it back into vector form.
void HirschFyeRun::green_vector_from_matrix(std::valarray<double> &green_vector, const dense_matrix & green_matrix)const{
  green_vector = 0.;
 
  for (int i=0; i<N; i++) {
    for (int j=i; j<N; j++) {
      green_vector[j-i] += green_matrix(j,i);
	  
	  if (j>i) {//consider antisymmetry (Fermions)
        green_vector[N-(j-i)] -= green_matrix(i,j);
	  }
	  
    }
  }
  green_vector /= N;
}


itime_green_function_t HirschFyeSim::get_result() const {

  multiple_vector_type g;
  int N(static_cast<int>(parms["N"]));
  alps::ObservableSet meas = get_measurements();
  alps::RealVectorObsevaluator G_meas_up=meas["G_meas_up"];
  alps::RealVectorObsevaluator G_meas_down=meas["G_meas_down"];  
#ifdef FOURPOINT
  if(parms.defined("MEASURE_FOURPOINT_FUNCTION")){
    int s=N;
    alps::RealVectorObsevaluator G_fourp_uu=meas["G_fourpoint_uu"];  
    alps::RealVectorObsevaluator G_fourp_dd=meas["G_fourpoint_dd"];  
    alps::RealVectorObsevaluator G_fourp_ud=meas["G_fourpoint_ud"];  
    alps::RealVectorObsevaluator G_fourp_du=meas["G_fourpoint_du"];  
    std::ofstream fourp_uu("fourpoint_uu.dat");
    std::ofstream fourp_dd("fourpoint_dd.dat");
    std::ofstream fourp_ud("fourpoint_ud.dat");
    std::ofstream fourp_du("fourpoint_du.dat");
    std::valarray<double> mean_array_uu(G_fourp_uu.mean());
    std::valarray<double> mean_array_dd(G_fourp_dd.mean());
    std::valarray<double> mean_array_ud(G_fourp_ud.mean());
    std::valarray<double> mean_array_du(G_fourp_du.mean());
    std::valarray<double> error_array_uu(G_fourp_uu.error());
    std::valarray<double> error_array_dd(G_fourp_dd.error());
    std::valarray<double> error_array_ud(G_fourp_ud.error());
    std::valarray<double> error_array_du(G_fourp_du.error());
    int s3=s*s*s;
    int s2=s*s;
    std::cout<<"fourpoint size: "<<mean_array_uu.size()<<" "<<s*s*s*s<<std::endl;
    for(int i=0;i<s;++i){
      for(int j=0;j<s;++j){
        for(int k=0;k<s;++k){
          for(int l=0;l<s;++l){
            fourp_uu<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<mean_array_uu[i*s3+j*s2+k*s+l]<<" "<<error_array_uu[i*s3+j*s2+k*s+l]<<std::endl;
            fourp_dd<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<mean_array_dd[i*s3+j*s2+k*s+l]<<" "<<error_array_dd[i*s3+j*s2+k*s+l]<<std::endl;
            fourp_ud<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<mean_array_ud[i*s3+j*s2+k*s+l]<<" "<<error_array_ud[i*s3+j*s2+k*s+l]<<std::endl;
            fourp_du<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<mean_array_du[i*s3+j*s2+k*s+l]<<" "<<error_array_du[i*s3+j*s2+k*s+l]<<std::endl;
          }
        }
      }
    }
  }
#endif //FOURPOINT
  std::valarray<double> g_up_mean_vector(G_meas_up.mean());
  std::valarray<double> g_down_mean_vector(G_meas_down.mean());
  //std::cout<<"measured: "<<std::endl;
  //also take into account G(beta)
  g.first.push_back(-1.+g_up_mean_vector[0]);
  g.second.push_back(-1.+g_down_mean_vector[0]);
  for(uint i=0;i<g_up_mean_vector.size();++i){
  	g.first.push_back(-g_up_mean_vector[g_up_mean_vector.size()-1-i]);
  	g.second.push_back(-g_down_mean_vector[g_up_mean_vector.size()-1-i]);
        //std::cout<<i<<"\t"<<-g_up_mean_vector[g_up_mean_vector.size()-1-i]<<"\t"<<G_meas_up.error()[g_up_mean_vector.size()-1-i]<<std::endl;
  }
      
  itime_green_function_t multiple_G(N+1, 1, 2);
  multiple_G.from_multiple_vector(g);
  return multiple_G;
}

