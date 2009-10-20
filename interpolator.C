/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2009 by Sebastian Fuchs <fuchs@comp-phys.org>
*                       Matthias Troyer <troyer@comp-phys.org>
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

#include "interpolator.h"
#include "cluster.h"
#include <alps/alea.h>




void Interpolator::interpolate() 
{
  for (int z=0; z<n_zone; ++z) {
    switch (transform_ptr->dimension()) {
    case 1 :
      interpolate_1D(z);
      break;
    case 2 : {
      interpolate_2D(z);
      break;
    }
    case 3 : {
      interpolate_3D(z);
      break;
    }
    default :
      boost::throw_exception(std::invalid_argument("Dimension>3: QCD?"));	
    }
  }
}




ParaInterpolator::ParaInterpolator(const alps::Parameters& parms)  : 
  Interpolator(parms), selfenergies(n_zone)
{
  const int n_sites = parms["SITES"];
  std::ifstream data(parms["SELFENERGIES"].c_str());
  if (!data.is_open()) 
    boost::throw_exception(std::invalid_argument("cannot read self-energies"));	
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
}




AntiferroInterpolator::AntiferroInterpolator(const alps::Parameters& parms) :
  Interpolator(parms), selfenergies(n_zone)  
{
  const int n_sites = parms["SITES"];
  std::ifstream data(parms["SELFENERGIES"].c_str());
  if (!data.is_open())
    boost::throw_exception(std::invalid_argument("cannot read self-energies"));	
  while (data) {
    int k1;
    int k2;
    int z;
    double omega;
    double Sigma_real;
    double Sigma_imag;
    data >> k1 >> k2 >> z >> omega >> Sigma_real >> Sigma_imag; 
    if (selfenergies[z].find(omega)==selfenergies[z].end()) {
      std::vector<matrix_t> kd(n_sites);
      for (std::vector<matrix_t>::iterator itm=kd.begin(); itm!=kd.end(); ++itm)
	itm->resize(2,2);
      selfenergies[z][omega] = kd;
    }
    const int k = (k1/2)*2;
    const int i = k1 % 2;
    const int j = k2 % 2;
    //if (i==j) {
      selfenergies[z][omega][k]  (i,j) = std::complex<double>(Sigma_real, Sigma_imag);
      selfenergies[z][omega][k+1](i,j) = std::complex<double>(Sigma_real, Sigma_imag);
      /*}
    else {
      selfenergies[z][omega][k]  (i,j) = std::complex<double>(0., 0.);
      selfenergies[z][omega][k+1](i,j) = std::complex<double>(0., 0.);
      }*/
  }
}




void ParaInterpolator::interpolate_spectrum(const int counter, const int zone, 
					    const std::vector<vector_t>& kvec,
					    const double offset, const int sign)
{
  std::stringstream name_str;
  name_str << name << zone << "_" << counter;
  std::ofstream out_str(name_str.str().c_str());
  for (std::vector<vector_t>::const_iterator itk=kvec.begin(); itk!=kvec.end(); ++itk) {
    vector_t k = *itk;
    std::map<double, std::complex<double> > interpolated_selfenergy;
    for (selfenergy_t::const_iterator it=selfenergies[zone].begin(); it!=selfenergies[zone].end(); ++it) {
      interpolated_selfenergy[it->first] = transform_ptr->interpolate(k, it->second, false);
    }
    double l = 0;
    for (int d=0; d<transform_ptr->dimension(); ++d) 
      l += (k[d]-kvec[0][d])*(k[d]-kvec[0][d]);
    l = sqrt(l);
    for (std::map<double, std::complex<double> >::const_iterator it=interpolated_selfenergy.begin(); 
	 it!=interpolated_selfenergy.end(); ++it) {
      std::complex<double> idelta(0., delta);
      std::complex<double> G = it->first + mu - transform_ptr->epsilon(k) - it->second + idelta;
      G = 1./G;
      for (int d=0; d<transform_ptr->dimension(); ++d) 
	out_str << k[d] << "\t" ;
      out_str << (offset+sign*l)/M_PI << "\t" 
	      << it->first << "\t" 
	      << G.real() << "\t" 
	      << G.imag() << "\n";
    }
    out_str << "\n";
  }
  out_str.close();
}


  

void AntiferroInterpolator::interpolate_spectrum(const int counter, const int zone, 
						 const std::vector<vector_t>& kvec,
						 const double offset, const int sign)
{
  std::stringstream name_str;
  name_str << name << zone << "_" << counter;
  std::ofstream out_str(name_str.str().c_str());
  for (std::vector<vector_t>::const_iterator itk=kvec.begin(); itk!=kvec.end(); ++itk) {
    vector_t k1 = *itk;
    vector_t k2 = k1;
    for (int d=0; d<transform_ptr->dimension(); ++d) 
      k2[d] += M_PI;
    std::map<double, matrix_t> interpolated_selfenergy;
    for (selfenergy_t::const_iterator it=selfenergies[zone].begin(); it!=selfenergies[zone].end(); ++it) {
      interpolated_selfenergy[it->first] = transform_ptr->interpolate(k1, it->second, true);
    }
    double l = 0;
    for (int d=0; d<transform_ptr->dimension(); ++d) 
      l += (k1[d]-kvec[0][d])*(k1[d]-kvec[0][d]);
    l = sqrt(l);
    for (std::map<double, matrix_t>::const_iterator it=interpolated_selfenergy.begin(); 
	 it!=interpolated_selfenergy.end(); ++it) {
      matrix_t G = -it->second;
      std::complex<double> idelta(0., delta);
      G(0,0) += it->first + mu - transform_ptr->epsilon(k1) + idelta;
      G(1,1) += it->first + mu - transform_ptr->epsilon(k2) + idelta; 
      G = ClusterHelper::Matrix<std::complex<double> >::inverse(G);
      std::complex<double> Gs = G(0,0) + G(1,1);
      for (int d=0; d<transform_ptr->dimension(); ++d) 
	out_str << k1[d] << "\t" ;
      out_str << (offset+sign*l)/M_PI << "\t" << it->first << "\t" 
	      << Gs.real() << "\t" << Gs.imag() << "\n";
    }
    out_str << "\n";
  }
  out_str.close();
}


  

void ParaInterpolator::interpolate_1D(const int zone) 
{
  std::vector<ClusterHelper::vector_type> kvec(101);
  for (int i=0; i<=100; ++i) {
    kvec[i].resize(1);
    kvec[i][0] = double(i)/100.*M_PI;
  }
  interpolate_spectrum(0, zone, kvec, 0., 1);
}




void AntiferroInterpolator::interpolate_1D(const int zone) 
{
  std::vector<ClusterHelper::vector_type> kvec(101);
  for (int i=0; i<=100; ++i) {
    kvec[i].resize(1);
    kvec[i][0] = double(i)/200.*M_PI;
  }
  interpolate_spectrum(0, zone, kvec, 0., 1);
}




void ParaInterpolator::interpolate_2D(const int zone) 
{
  std::vector<ClusterHelper::vector_type> kvec(101);
  for (int i=0; i<=100; ++i) {
    kvec[i].resize(2);
    kvec[i][0] = double(i)/100.*M_PI;
    kvec[i][1] = kvec[i][0]; 
  }
  interpolate_spectrum(0, zone, kvec, 0., 1);
  for (int i=0; i<=100; ++i) {
    kvec[i][0] = double(i)/100.*M_PI;
    kvec[i][1] = 0.;
  }
  interpolate_spectrum(1, zone, kvec, 0., -1);
  for (int i=0; i<=100; ++i) {
    kvec[i][0] = M_PI;
    kvec[i][1] = double(i)/100.*M_PI;
  }
  interpolate_spectrum(2, zone, kvec, -M_PI, -1);
}




void AntiferroInterpolator::interpolate_2D(const int zone) 
{
  std::vector<ClusterHelper::vector_type> kvec(101);
  for (int i=0; i<=100; ++i) {
    kvec[i].resize(2);
    kvec[i][0] = double(i)/200.*M_PI;
    kvec[i][1] = kvec[i][0]; 
  }
  interpolate_spectrum(0, zone, kvec, 0., 1);
  for (int i=0; i<=100; ++i) {
    kvec[i][0] = double(i)/100.*M_PI;
    kvec[i][1] = 0.;
  }
  interpolate_spectrum(1, zone, kvec, 0., -1);
  for (int i=0; i<=100; ++i) {
    double del_k = double(i)/100.*M_PI/2.;
    kvec[i][0] = M_PI/2 + del_k;
    kvec[i][1] = M_PI/2 - del_k;
  }
  interpolate_spectrum(2, zone, kvec, M_PI/sqrt(2), 1);
}




void ParaInterpolator::interpolate_3D(const int zone) 
{
  std::vector<ClusterHelper::vector_type> kvec(101);
  for (int i=0; i<=100; ++i) {
    kvec[i].resize(3);
    kvec[i][0] = double(i)/100.*M_PI;
    kvec[i][1] = kvec[i][0];
    kvec[i][2] = kvec[i][0];
  }
  interpolate_spectrum(0, zone, kvec, 0, 1);
  for (int i=0; i<=100; ++i) {
    kvec[i][0] = double(i)/100.*M_PI;
    kvec[i][1] = 0.;
    kvec[i][2] = 0.; 
  }
  interpolate_spectrum(1, zone, kvec, 0, -1);
  for (int i=0; i<=100; ++i) {
    kvec[i][0] = M_PI;
    kvec[i][1] = double(i)/100.*M_PI;
    kvec[i][2] = 0.;
  }
  interpolate_spectrum(2, zone, kvec, -M_PI, -1);
  for (int i=0; i<=100; ++i) {
    kvec[i][0] = M_PI;
    kvec[i][1] = M_PI;
    kvec[i][2] = double(i)/100.*M_PI;
  }
  interpolate_spectrum(3, zone, kvec, -2*M_PI, -1);
}




void AntiferroInterpolator::interpolate_3D(const int zone) 
{
  std::vector<ClusterHelper::vector_type> kvec(101);
  for (int i=0; i<=100; ++i) {
    kvec[i].resize(3);
    kvec[i][0] = double(i)/200.*M_PI;
    kvec[i][1] = kvec[i][0];
    kvec[i][2] = kvec[i][0];
  }
  interpolate_spectrum(0, zone, kvec, 0, 1);
  for (int i=0; i<=100; ++i) {
    kvec[i][0] = double(i)/100.*M_PI;
    kvec[i][1] = 0.;
    kvec[i][2] = 0.; 
  }
  interpolate_spectrum(1, zone, kvec, 0, -1);
  for (int i=0; i<=100; ++i) {
    kvec[i][0] = M_PI;
    kvec[i][1] = double(i)/200.*M_PI;
    kvec[i][2] = 0.;
  }
  interpolate_spectrum(2, zone, kvec, -M_PI, -1);
  for (int i=0; i<=100; ++i) {
    double del_k = double(i)/100.*M_PI/4.;
    kvec[i][0] = M_PI/2. + del_k;
    kvec[i][1] = M_PI/2. + del_k;
    kvec[i][2] = M_PI/2. - 2*del_k;
  }
  interpolate_spectrum(3, zone, kvec, sqrt(3)/2.*M_PI, 1);
  for (int i=0; i<=100; ++i) {
    double del_k = double(i)/100.*M_PI/4.;
    kvec[i][0] = 3./4.*M_PI + del_k;
    kvec[i][1] = 3./4.*M_PI - del_k;
    kvec[i][2] = 0.;
  }
  interpolate_spectrum(4, zone, kvec, (sqrt(3)/2.+sqrt(3./8.))*M_PI, 1);
}

