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

#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H


#include "cluster.h"
#include <alps/alea.h>



class Interpolator 
{

public:

  typedef ClusterHelper::vector_type vector_t;

  Interpolator(const alps::Parameters& parms) :
    mu(parms["MU"]), name("Dispersion_"), n_zone(parms["FLAVORS"])
  {
    DCATransformer::generate_transformer(parms, transform_ptr);
  }

  virtual ~Interpolator() {}
  

  void interpolate();
   
protected:

  boost::shared_ptr<DCATransformer> transform_ptr;
  const double mu;
  const std::string name;
  const unsigned int n_zone;
  
  virtual void interpolate_spectrum(const int counter, const int zone, 
				    const std::vector<vector_t>& kvec,
				    const double offset, const int sign) = 0;
  virtual void interpolate_1D(const int zone) = 0;
  virtual void interpolate_2D(const int zone) = 0;
  virtual void interpolate_3D(const int zone) = 0;
 
};




class ParaInterpolator : public Interpolator 
{
public:
  
  typedef std::map<double, std::vector<std::complex<double> > > selfenergy_t;

  ParaInterpolator(const alps::Parameters& parms);
  
private:
  
  virtual void interpolate_spectrum(const int counter, const int zone, 
				    const std::vector<vector_t>& kvec,
				    const double offset, const int sign);  
  virtual void interpolate_1D(const int zone);
  virtual void interpolate_2D(const int zone);
  virtual void interpolate_3D(const int zone);
  
  std::vector<selfenergy_t> selfenergies;
};


      

class AntiferroInterpolator : public Interpolator 
{
public:

  typedef ClusterHelper::Matrix<std::complex<double> >::Type matrix_t;
  typedef std::map<double, std::vector<matrix_t> > selfenergy_t;
  
  AntiferroInterpolator(const alps::Parameters& parms);
  
private:
  
  virtual void interpolate_spectrum(const int counter, const int zone, 
				    const std::vector<vector_t>& kvec,
				    const double offset, const int sign);  
  virtual void interpolate_1D(const int zone);
  virtual void interpolate_2D(const int zone);
  virtual void interpolate_3D(const int zone);
  
  std::vector<selfenergy_t> selfenergies;
};


      
#endif
