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

#include"hybretintfun.hpp"
//construct a retarded interaction function. ntime: number of time slices.
//noffidag_orbitals: number of offdiagonal orbitals. ndiag_orbitals: number
//of diagonal orbitals. 
ret_int_fun::ret_int_fun(const alps::params &p):
green_function<double>(p["N_TAU"]+1, 1, 2)//2 "flavors": K and K'
{
  bool use_retarded_interaction=p.defined("RET_INT_K");
  if(!use_retarded_interaction) return;

  if(!p.defined("N_TAU") || (int)(p["N_TAU"])==0) throw std::invalid_argument("define parameter N_TAU, the number of retarded interaction time slices!");
  beta_=p["BETA"];
  
  //read in Green's function from a file
  read_interaction_K_function(p);
  interaction_K_function_sanity_check();

  extern int global_mpi_rank;
  if(global_mpi_rank==0){
    std::cout<<*this<<std::endl;
  }
}

void ret_int_fun::interaction_K_function_sanity_check(void){
for(std::size_t i=0; i<ntime();++i)
  for(std::size_t j=0; j<nflavor();++j)
    if(operator()(i,j) < 0. && j==0) throw std::invalid_argument("Problem with retarded interaction function: RET_INT_K(\\tau) < 0. K should always be positive!");
if(operator()(0,0)!=0.) throw std::invalid_argument("Problem with retarded interaction function: RET_INT_K(\\tau=0) must be zero.");
}

//this routine reads in the retarded interaction function, either from a text file or from an hdf5 file (for easy passing of binary data).
//In  case of text files the file format is index - hyb_1 - hyb2 - hyb3 - ... in columns that go from time=0 to time=beta. Note that
//the retarded interaction function is in imaginary time and always positive both for negative and positive times. It is also symmetric.
void ret_int_fun::read_interaction_K_function(const alps::params &p){
  if(!p.defined("RET_INT_K")) throw(std::invalid_argument(std::string("Parameter RET_INT_K missing, filename for retarded interaction function not specified.")));
  std::string fname=p["RET_INT_K"].cast<std::string>();
  if(p.defined("K_IN_HDF5") && p["K_IN_HDF5"].cast<bool>()){//attempt to read from h5 archive
    alps::hdf5::archive ar(fname, alps::hdf5::archive::READ);
    std::vector<double> tmp(ntime());
    ar>>alps::make_pvp("/Ret_int_K",tmp);
      for(std::size_t i=0; i<ntime(); i++)
        operator()(i,0)=tmp[i];
    tmp.clear();
    ar>>alps::make_pvp("/Ret_int_Kp",tmp);
      for(std::size_t i=0; i<ntime(); i++)
        operator()(i,1)=tmp[i];
    tmp.clear();
  }
  else{//read from text file
    std::ifstream infile(fname.c_str());
    if(!infile.good()){
      throw(std::invalid_argument(std::string("could not open retarded interaction file (text format) ") + p["RET_INT_K"]));
    }
    for (std::size_t i=0; i<ntime(); i++) {
      if(!infile.good()){
        throw(std::invalid_argument(std::string("could not read retarded interaction file (text format) ") + p["RET_INT_K"] + ". probably wrong number of lines"));
      }
      double dummy;
      infile >> dummy;
      //std::cout<<i<<" "<<ntime()<<" "<<"dummy: "<<dummy<<std::endl;
      for (std::size_t j=0; j<nflavor(); j++){
        if(!infile.good()){
          throw(std::invalid_argument(std::string("could not read retarded interaction file (text format) ") + p["RET_INT_K"] + ". probably wrong number of columns"));
        }
        double delta;
        infile >> delta;
        operator()(i,j)=delta;
      }
    }
  }
}

std::ostream &operator<<(std::ostream &os, const ret_int_fun &K){
  os<<"the retarded interaction function and derivative are: "<<std::endl;
  for(int i=0;i<10;++i){ std::cout<<i<<" "; for(std::size_t j=0;j<K.nflavor();++j){ std::cout<<K(i,j)<<" ";} std::cout<<std::endl; }
  os<<"... *** etc *** ...\n";
  os<<K.ntime()-1<<" "; for(std::size_t j=0;j<K.nflavor();++j){ std::cout<<K(K.ntime()-1,j)<<" ";} std::cout<<std::endl;
  return os;
}
//linear interpolation of the retarded interaction function 
double ret_int_fun::interpolate(double time) const{
  time=std::abs(time);

  double n = time/beta_*(ntime()-1);
  int n_lower = (int)n; // interpolate linearly between n_lower and n_lower+1

  return operator()(n_lower,0) + (n-n_lower)*(operator()(n_lower+1,0)-operator()(n_lower,0));
}

//linear interpolation of the retarded interaction function.
double ret_int_fun::interpolate_deriv(double time) const{
  //if(time<-beta_ || time > beta_) std::cout<< "!!" << std::endl;
  if(time<0.) return -1.*interpolate_deriv(-time);

  double n = time/beta_*(ntime()-1);
  int n_lower = (int)n; // interpolate linearly between n_lower and n_lower+1

  return operator()(n_lower,1) + (n-n_lower)*(operator()(n_lower+1,1)-operator()(n_lower,1));
}
