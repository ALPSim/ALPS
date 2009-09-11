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



std::vector<std::map<double, std::vector<std::complex<double> > > > read_selfenergies(const alps::Parameters& parms) 
{
  int n_zone = parms["FLAVORS"];
  int n_sites = parms["SITES"];
  std::vector<std::map<double, std::vector<std::complex<double> > > > selfenergies(n_zone);
  std::ifstream data(parms["SELFENERGIES"].c_str());
  while (data) {
    int k;
    int z;
    double omega;
    double Sigma_real;
    double Sigma_imag;
    data >> k >> z >> omega >> Sigma_real >> Sigma_imag; 
    if (selfenergies[z].find(omega)==selfenergies[z].end()) {
      std::vector<std::complex<double> > kd(n_sites);
      selfenergies[z][omega] = kd;
    }
    selfenergies[z][omega][k] = std::complex<double>(Sigma_real, Sigma_imag);
  }
  return selfenergies;
}




void interpolate(std::string name, const boost::shared_ptr<DCATransformer> transform_ptr,
		 const std::map<double, std::vector<std::complex<double> > >& selfenergies, 
		 const std::vector<ClusterHelper::vector_type>& kvec,
		 const double offset, const int sign, const double mu, 
		 const bool antiferromagnet)
{
  std::ofstream out_str(name.c_str());
  double l0 = 0;
  for (int d=0; d<transform_ptr->dimension(); ++d) 
    l0 += kvec[0][d]*kvec[0][d];
  l0 =sqrt(l0);
  for (std::vector<ClusterHelper::vector_type>::const_iterator itk=kvec.begin(); itk!=kvec.end(); ++itk) {
    ClusterHelper::vector_type k = *itk;
    std::map<double, std::complex<double> > spectrum;
    for (std::map<double, std::vector<std::complex<double> > >::const_iterator it=selfenergies.begin(); it!=selfenergies.end(); ++it) {
      spectrum[it->first] = transform_ptr->interpolate(k, it->second, antiferromagnet);
    }
    double l = 0;
    for (int d=0; d<transform_ptr->dimension(); ++d) 
      l += k[d]*k[d];
    l =sqrt(l);
    for (std::map<double, std::complex<double> >::iterator it=spectrum.begin(); it!=spectrum.end(); ++it) {
      std::complex<double> G = 1./(it->first + mu - transform_ptr->epsilon(k) - it->second);
      for (int d=0; d<transform_ptr->dimension(); ++d) 
	out_str << k[d] << "\t" ;
      out_str << (offset+sign*(l-l0))/M_PI << "\t" << it->first << "\t" 
	      << G.real() << "\t" << G.imag() << "\n"; 
    }
    out_str << "\n";
  }
  out_str.close();
}



void interpolate_spectra(const alps::Parameters& parms)
{
  std::vector<std::map<double, std::vector<std::complex<double> > > > selfenergies = read_selfenergies(parms);
  boost::shared_ptr<DCATransformer> transform_ptr;
  DCATransformer::generate_transformer(parms, transform_ptr);
  const double mu = parms["MU"];
  const bool antiferromagnet = parms.value_or_default("ANTIFERROMAGNET", false);
  for (int z=0; z<selfenergies.size(); ++z) {
    std::stringstream name;
    name << "Dispersion" << "_" << z;
    std::vector<ClusterHelper::vector_type> kvec(101);
    switch (transform_ptr->dimension()) {
    case 1 :
      for (int i=0; i<=100; ++i) {
	kvec[i].resize(1);
	kvec[i][0] = double(i)/100.*M_PI;
	if (antiferromagnet) 
	  kvec[i][0] /= 2.;
      }
      interpolate(name.str(), transform_ptr, selfenergies[z], kvec, 0., 1, mu, antiferromagnet);
      break;
    case 2 : {
      for (int i=0; i<=100; ++i) {
	kvec[i].resize(2);
	kvec[i][0] = double(i)/100.*M_PI;
	kvec[i][1] = kvec[i][0];
	if (antiferromagnet)
	  kvec[i][0] /= 2.;
      }
      std::stringstream name1;
      name1 << name.str() << "_1";
      interpolate(name1.str(), transform_ptr, selfenergies[z], kvec, 0., 1, mu, antiferromagnet);
      for (int i=0; i<=100; ++i) {
	kvec[i].resize(2);
	kvec[i][0] = double(i)/100.*M_PI;
	kvec[i][1] = 0;
      }
      std::stringstream name2;
      name2 << name.str() << "_2";
      interpolate(name2.str(), transform_ptr, selfenergies[z], kvec, 0., -1, mu, antiferromagnet);
      if (antiferromagnet) {
	for (int i=0; i<=100; ++i) {
	  kvec[i].resize(2);
	  double del_k = double(i)/100.*M_PI/2.;
	  kvec[i][0] = M_PI/2 + del_k;
	  kvec[i][1] = M_PI/2 - del_k;
	}
	std::stringstream name3;
	name3 << name.str() << "_3";
	interpolate(name3.str(), transform_ptr, selfenergies[z], kvec, M_PI/sqrt(2), 1, mu, antiferromagnet);
      }
      else {
	for (int i=0; i<=100; ++i) {
	  kvec[i].resize(2);
	  kvec[i][0] = M_PI;
	  kvec[i][1] = double(i)/100.*M_PI;
	}
	std::stringstream name3;
	name3 << name.str() << "_3";
	interpolate(name3.str(), transform_ptr, selfenergies[z], kvec, -M_PI, -1, mu, antiferromagnet);
      }
      break;
    }
    case 3 :
      if (antiferromagnet) {
	/*for (int i=0; i<=100; ++i) {
	  kvec[i].resize(3);
	  kvec[i][0] = double(i)/100.*M_PI/2.;
	  kvec[i][1] = kvec[i][0];
	  kvec[i][2] = kvec[i][0];
	}
	interpolate(name1.str(), transform_ptr, selfenergies[z], kvec, 0, 1, mu, antiferromagnet);
	for (int i=0; i<=100; ++i) {
	  kvec[i].resize(3);
	  kvec[i][0] = double(i)/100.*M_PI/2.;
	  kvec[i][1] = 0;
	  kvec[i][2] = 0;
	}
	interpolate(name1.str(), transform_ptr, selfenergies[z], kvec, 0, -1, mu, antiferromagnet);*/
      }
      else {
	for (int i=0; i<=100; ++i) {
	  kvec[i].resize(3);
	  kvec[i][0] = double(i)/100.*M_PI;
	  kvec[i][1] = kvec[i][0];
	  kvec[i][2] = kvec[i][0];
	}
	std::stringstream name1;
	name1 << name.str() << "_1";
	interpolate(name1.str(), transform_ptr, selfenergies[z], kvec, 0, 1, mu, antiferromagnet);
	for (int i=0; i<=100; ++i) {
	  kvec[i].resize(3);
	  kvec[i][0] = double(i)/100.*M_PI;
	  kvec[i][1] = 0;
	  kvec[i][2] = 0; 
	}
	std::stringstream name2;
	name2 << name.str() << "_2";
	interpolate(name2.str(), transform_ptr, selfenergies[z], kvec, 0, -1, mu, antiferromagnet);
	for (int i=0; i<=100; ++i) {
	  kvec[i].resize(3);
	  kvec[i][0] = M_PI;
	  kvec[i][1] = double(i)/100.*M_PI;
	  kvec[i][2] = 0;
	}
	std::stringstream name3;
	name3 << name.str() << "_3";
	interpolate(name3.str(), transform_ptr, selfenergies[z], kvec, M_PI, -1, mu, antiferromagnet);
	for (int i=0; i<=100; ++i) {
	  kvec[i].resize(3);
	  kvec[i][0] = M_PI;
	  kvec[i][1] = M_PI;
	  kvec[i][2] = double(i)/100.*M_PI;
	}
	std::stringstream name4;
	name4 << name.str() << "_4";
	interpolate(name4.str(), transform_ptr, selfenergies[z], kvec, 2*M_PI, -1, mu, antiferromagnet);
      }
      break;
    default :
      boost::throw_exception(std::invalid_argument("QCD?"));	
    }
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



