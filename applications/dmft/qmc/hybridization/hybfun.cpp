/****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2012 by Emanuel Gull <gull@pks.mpg.de>,
 *                   
 *  based on an earlier version by Philipp Werner and Emanuel Gull
 *
 *
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#include"hybfun.hpp"
//construct a hybridization function. ntime: number of time slices.
//noffidag_orbitals: number of offdiagonal orbitals. ndiag_orbitals: number
//of diagonal orbitals. 
hybfun::hybfun(const alps::params &p):
green_function<double>(p["N_TAU"]+1, 1, p["N_ORBITALS"])
{
  if(!p.defined("N_TAU") || (int)(p["N_TAU"])==0) throw std::invalid_argument("define parameter N_TAU, the number of hybridization time slices!");
  beta_=p["BETA"];

  //read in Green's function from a file
  read_hybridization_function(p);
  hybridization_function_sanity_check();

  extern int global_mpi_rank;
  if(global_mpi_rank==0){
    std::cout<<*this<<std::endl;
  }
}

void hybfun::hybridization_function_sanity_check(void){
for(std::size_t i=0; i<ntime();++i)
  for(std::size_t j=0; j<nflavor();++j)
    if(operator()(i,j)>0.) {
      std::cerr << "ERROR: Delta(t="<<i<<"; f="<<j<<") = " << operator()(i,j) << "  is positive." << std::endl;
      std::cerr << "Note: small positive values might be due to noise, in that case try to enhance the MAX_TIME." << std::endl << std::flush;
      throw std::invalid_argument("Problem with hybridization function: Delta(\\tau) > 0. Delta should always be negative!");
    }
}

//this routine reads in the hybridization function, either from a text file or from an hdf5 file (for easy passing of binary data).
//In case of text files the file format is index - hyb_1 - hyb2 - hyb3 - ... in columns that go from time=0 to time=beta. Note that
//the hybridization function is in imaginary time and always positive between zero and \beta.
void hybfun::read_hybridization_function(const alps::params &p){
  if(!p.defined("DELTA")) throw(std::invalid_argument(std::string("Parameter DELTA missing, filename for hybridization function not specified.")));
  std::string fname=p["DELTA"].cast<std::string>();
  if(p.defined("DELTA_IN_HDF5") && p["DELTA_IN_HDF5"].cast<bool>()){//attempt to read from h5 archive
    alps::hdf5::archive ar(fname, alps::hdf5::archive::READ);
    if(p.defined("DMFT_FRAMEWORK") && p["DMFT_FRAMEWORK"].cast<bool>()){//read in as green_function
      read_hdf5(ar,"/Delta");
    }
    else{//plain hdf5
      std::vector<double> tmp(ntime());
      for(std::size_t j=0; j<nflavor(); j++){
        std::stringstream path; path<<"/Delta_"<<j;
        ar>>alps::make_pvp(path.str(),tmp);
        for(std::size_t i=0; i<ntime(); i++)
          operator()(i,j)=tmp[i];
      }
      tmp.clear();
    }
  }
  else{//read from text file
    std::ifstream infile(fname.c_str());
    if(!infile.good()){
      throw(std::invalid_argument(std::string("could not open hybridization file (text format) ") + p["DELTA"] + " for Delta function"));
    }
    for (std::size_t i=0; i<ntime(); i++) {
      if(!infile.good()){
        throw(std::invalid_argument(std::string("could not read hybridization file (text format) ") + p["DELTA"] + ". probably wrong number of lines"));
      }
      double dummy;
      infile >> dummy;
      //std::cout<<i<<" "<<ntime()<<" "<<"dummy: "<<dummy<<std::endl;
      for (std::size_t j=0; j<nflavor(); j++){
        if(!infile.good()){
          throw(std::invalid_argument(std::string("could not read hybridization file (text format) ") + p["DELTA"] + ". probably wrong number of columns"));
        }
        double delta;
        infile >> delta;
        operator()(i,j)=delta;
      }
    }
  }
}

std::ostream &operator<<(std::ostream &os, const hybfun &hyb){
  os<<"the hybridization function is: "<<std::endl;
  for(int i=0;i<10;++i){ std::cout<<i<<" "; for(std::size_t j=0;j<hyb.nflavor();++j){ std::cout<<hyb(i,j)<<" ";} std::cout<<std::endl; }
  os<<"... *** etc *** ...\n";
  os<<hyb.ntime()-1<<" "; for(std::size_t j=0;j<hyb.nflavor();++j){ std::cout<<hyb(hyb.ntime()-1,j)<<" ";} std::cout<<std::endl;
  return os;
}

//linear interpolation of the hybridization function. 
double hybfun::interpolate(double time, int orbital) const{

  ///TODO: do a spline here, write a better/faster interpolation routine.
  double sign=1;
  if (time<0) {
    time += beta_;
    sign=-1;
  }

  //this is the overall flip of Delta (in comparison to F)
  sign*=-1;
  //the code takes Delta as input, but internally works with F(tau)=-Delta(beta-tau)
  time=beta_-time;

  double n = time/beta_*(ntime()-1);
  int n_lower = (int)n; // interpolate linearly between n_lower and n_lower+1

  return sign*(operator()(n_lower,orbital) + (n-n_lower)*(operator()(n_lower+1,orbital)-operator()(n_lower,orbital)));
}
