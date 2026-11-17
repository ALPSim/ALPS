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
#include "alps/parameter.h"



class Interpolator 
{

public:

  typedef ClusterHelper::vector_type vector_t;

  Interpolator(const alps::Parameters& parms);
  virtual ~Interpolator() {}
  virtual void local_green_function() = 0;
  void interpolate();

protected:

  boost::shared_ptr<ClusterTransformer> transform_ptr;
  const double mu;
  const double delta;
  const double delta_green;
  const std::string name;
  const unsigned int n_zone;
  const unsigned int n_k;
  
  virtual void interpolate_spectrum(const unsigned int counter, const unsigned int, 
				    const std::vector<vector_t>& kvec,
				    const double offset, const int sign) = 0;

  void bare_dispersion(const unsigned int counter, const unsigned int, 
		       const std::vector<vector_t>& kvec,
		       const double offset, const int sign);

  virtual void interpolate_1D(const unsigned int) = 0;
  virtual void interpolate_2D(const unsigned int) = 0;
  virtual void interpolate_3D(const unsigned int) = 0;
 
};




class ParaInterpolator : public Interpolator 
{
public:
  
  typedef std::map<double, std::vector<std::complex<double> > > selfenergy_t;

  ParaInterpolator(const alps::Parameters& parms);

  virtual void local_green_function();
  
private:
  
  virtual void interpolate_spectrum(const unsigned int counter, const unsigned int, 
				    const std::vector<vector_t>& kvec,
				    const double offset, const int sign);  

  virtual void interpolate_1D(const unsigned int);
  virtual void interpolate_2D(const unsigned int);
  virtual void interpolate_3D(const unsigned int);
  
  std::vector<selfenergy_t> selfenergies;
};


      

class AntiferroInterpolator : public Interpolator 
{
public:

  typedef ClusterHelper::Matrix<std::complex<double> >::Type matrix_t;
  typedef std::map<double, std::vector<matrix_t> > selfenergy_t;
  
  AntiferroInterpolator(const alps::Parameters& parms);
  
  virtual void local_green_function() { ; }

private:
  
  virtual void interpolate_spectrum(const unsigned int counter, const unsigned int, 
				    const std::vector<vector_t>& kvec,
				    const double offset, const int sign);  

  virtual void interpolate_1D(const unsigned int);
  virtual void interpolate_2D(const unsigned int);
  virtual void interpolate_3D(const unsigned int);
  
  std::vector<selfenergy_t> selfenergies;
};


      
#endif
