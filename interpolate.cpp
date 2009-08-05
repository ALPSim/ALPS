/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2003-2007 by Sebastian Fuchs <fuchs@comp-phys.org>
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

#include "cluster.h"
#include <alps/alea.h>



std::vector<std::map<double, std::vector<double> > > read_spectra(const alps::Parameters& parms) 
{
  int n_zone = parms["FLAVORS"];
  int n_sites = parms["SITES"];
  std::vector<std::map<double, std::vector<double> > > spectra(n_zone);
  std::ifstream data(parms["SPECTRA"].c_str());
  while (data) {
    int k;
    int z;
    double omega;
    double A;
    data >> k >> z >> omega >> A; 
    if (spectra[z].find(omega)==spectra[z].end()) {
      std::vector<double> kd(n_sites);
      spectra[z][omega] = kd;
    }
    spectra[z][omega][k] = A;
  }
  return spectra;
}



std::map<double, double> interpolate(const std::map<double, std::vector<double> >& spectra, 
				     const ClusterHelper::vector_type k,
				     const ClusterHelper& helper)
{
  std::map<double, double> spectrum;
  for (std::map<double, std::vector<double> >::const_iterator it=spectra.begin(); it!=spectra.end(); ++it) {
    spectrum[it->first] = helper.interpolate(k, it->second);
  }
  return spectrum;
}



void interpolate_spectra(const alps::Parameters& parms)
{
  std::vector<std::map<double, std::vector<double> > > spectra = read_spectra(parms);
  ClusterHelper helper(parms);
  if (helper.dimension() != 3) 
    boost::throw_exception(std::invalid_argument("implement interpolation for d!=3!"));	
  for (int z=0; z<spectra.size(); ++z) {
    std::stringstream name1;
    name1 << "Gamma_R" << "_" << z;
    std::ofstream gamma_r_str(name1.str().c_str());
    gamma_r_str << "# kx\tky\tkz\tomega\tspectrum\n";
    for (int i=0; i<=100; ++i) {
      ClusterHelper::vector_type k(3);
      k[0] = double(i)/100.*M_PI;
      k[1] = k[0];
      k[2] = k[0];
      std::map<double, double> spec = interpolate(spectra[z], k, helper);
      for (std::map<double, double>::const_iterator it=spec.begin(); it!=spec.end(); ++it) {
	gamma_r_str << k[0] << "\t" << k[1] << "\t"  << k[2] << "\t" 
		    << it->first << "\t" << it->second << "\n"; 
      }
      gamma_r_str << "\n";
    }
    gamma_r_str.close();
    std::stringstream name2;
    name2 << "R_X" << "_" << z;
    std::ofstream r_x_str(name2.str().c_str());
    r_x_str << "# kx\tky\nkz\tomega\tspectrum\n";
    for (int i=100; i>=0; --i) {
      ClusterHelper::vector_type k(3);
      k[0] = M_PI;
      k[1] = double(i)/100*M_PI;
      k[2] = k[1];
      std::map<double, double> spec = interpolate(spectra[z], k, helper);
      for (std::map<double, double>::const_iterator it=spec.begin(); it!=spec.end(); ++it) {
	r_x_str << k[0] << "\t" << k[1] << "\t"  << k[2] << "\t" 
		<< it->first << "\t" << it->second << "\n"; 
      }
      r_x_str << "\n";
    }
    r_x_str.close();
    std::stringstream name3;
    name3 << "X_Gamma" << "_" << z;
    std::ofstream x_gamma_str(name3.str().c_str());
    x_gamma_str << "# kx\tky\nkz\tomega\tspectrum\n";
    for (int i=100; i>=0; --i) {
      ClusterHelper::vector_type k(3);
      k[0] = double(i)/100*M_PI;
      k[1] = 0;
      k[2] = 0.;
      std::map<double, double> spec = interpolate(spectra[z], k, helper);
      for (std::map<double, double>::const_iterator it=spec.begin(); it!=spec.end(); ++it) {
	x_gamma_str << k[0] << "\t" << k[1] << "\t"  << k[2] << "\t" 
		    << it->first << "\t" << it->second << "\n"; 
      }
      x_gamma_str << "\n";
    }
    x_gamma_str.close();
  }
}

      

  
int main(int argc, char** argv)
{
#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif
    int i=1;
    alps::Parameters new_parms;
    while (i<argc-1 && argv[i][0]=='-') {
      new_parms[argv[i]+1]=argv[i+1];
      i+=2;
    }
    if (argc-i!=1) {
      std::cerr << "Usage: \n" << argv[0] 
		<<  "inputfile \n";
      exit(1);
    }
    alps::Parameters parms;
    std::ifstream is(argv[i]);
    if(!is.is_open()) {
      std::cerr << "parameter file argv[" << i<< " ] " <<argv[i] << "is not open! exiting!" << std::endl; 
      abort(); 
    }
    is>>parms;
    is.close();
    interpolate_spectra(parms);
#ifndef BOOST_NO_EXCEPTIONS
  }
  catch (std::exception& e)
    {
      std::cerr << "Caught exception: " << e.what() << "\n";
      std::exit(-5);
    }
#endif
}



