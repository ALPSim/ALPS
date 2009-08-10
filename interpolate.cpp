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



void interpolate(std::ofstream& out_str, const ClusterHelper& helper,
		 const std::map<double, std::vector<double> >& spectra, 
		 const std::vector<ClusterHelper::vector_type>& kvec,
		 const double offset, const double sign, 
		 const bool antiferromagnet)
{
  for (std::vector<ClusterHelper::vector_type>::const_iterator itk=kvec.begin(); itk!=kvec.end(); ++itk) {
    ClusterHelper::vector_type k = *itk;
    std::map<double, double> spectrum;
    for (std::map<double, std::vector<double> >::const_iterator it=spectra.begin(); it!=spectra.end(); ++it) {
      spectrum[it->first] = helper.interpolate(k, it->second, antiferromagnet);
    }
    double l = sqrt(k[0]*k[0] + k[1]*k[1] + k[2]*k[2]);
    for (std::map<double, double>::const_iterator it=spectrum.begin(); it!=spectrum.end(); ++it) {
      out_str << k[0] << "\t" << k[1] << "\t"  << k[2] << "\t" << offset+sign*l << "\t" 
	      << it->first << "\t" << it->second << "\n"; 
    }
    out_str << "\n";
  }
}



void interpolate_spectra(const alps::Parameters& parms)
{
  std::vector<std::map<double, std::vector<double> > > spectra = read_spectra(parms);
  ClusterHelper helper(parms);
  const bool antiferromagnet = parms["ANTIFERROMAGNET"];
  if (helper.dimension() != 3) 
    boost::throw_exception(std::invalid_argument("implement interpolation for d!=3!"));	
  for (int z=0; z<spectra.size(); ++z) {
    std::stringstream name1;
    name1 << "Dispersion" << "_" << z;
    std::ofstream out_str(name1.str().c_str());
    std::vector<ClusterHelper::vector_type> kvec(101);
    if (antiferromagnet) {
      for (int i=0; i<=100; ++i) {
	kvec[i].resize(3);
	kvec[i][0] = double(i)/100.*M_PI/2.;
	kvec[i][1] = kvec[i][0];
	kvec[i][2] = kvec[i][0];
      }
      interpolate(out_str, helper, spectra[z], kvec, 0, 1, antiferromagnet);
      for (int i=0; i<=100; ++i) {
	kvec[i].resize(3);
	kvec[i][0] = double(i)/100.*M_PI/2.;
	kvec[i][1] = 0;
	kvec[i][2] = 0;
      }
      interpolate(out_str, helper, spectra[z], kvec, 0, -1, antiferromagnet);
    }
    else {
      for (int i=0; i<=100; ++i) {
	kvec[i].resize(3);
	kvec[i][0] = double(i)/100.*M_PI;
	kvec[i][1] = kvec[i][0];
	kvec[i][2] = kvec[i][0];
      }
      interpolate(out_str, helper, spectra[z], kvec, 0, 1, antiferromagnet);
      for (int i=0; i<=100; ++i) {
	kvec[i].resize(3);
	kvec[i][0] = double(i)/100.*M_PI;
	kvec[i][1] = 0;
	kvec[i][2] = 0; 
      }
      interpolate(out_str, helper, spectra[z], kvec, 0, -1, antiferromagnet);
      for (int i=0; i<=100; ++i) {
	kvec[i].resize(3);
	kvec[i][0] = M_PI;
	kvec[i][1] = double(i)/100.*M_PI;
	kvec[i][2] = kvec[i][1];
      }
      interpolate(out_str, helper, spectra[z], kvec, M_PI, -1, antiferromagnet);
    }
    out_str.close();
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



