/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1999-2006 by Matthias Troyer <troyer@comp-phys.org>,
*                            Fabian Stoeckli <fabstoec@student.ethz.ch>
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

#ifndef ALPS_APPLICATIONS_MC_SPIN_ON_H_2
#define ALPS_APPLICATIONS_MC_SPIN_ON_H_2

#include "tinyvec.h"
#include "matrices.h"
#include <iostream>
#include <boost/random/uniform_real.hpp>

template<unsigned int N, class RealType = double>
class uniform_on_sphere {
  /* result is a new spin orientation  */
  typedef TinyVector<RealType, N> result_type;

private:
  result_type container_;

public:
  uniform_on_sphere() {}
  template<class RNG> const result_type & operator()(RNG& rng) {
    RealType sqsum = 0;
    for (int i=0; i<N; i++) {
      RealType val = boost::normal_distribution<>()(rng);
      container_[i] = val;
      sqsum += val * val;
    }
#ifndef BOOST_NO_STDC_NAMESPACE
    using std::sqrt;
#endif
    // renormalise to length 1
    container_ /= sqrt(sqsum);
    return container_;
  }
};

template<class RealType>
class uniform_on_sphere<1,RealType>
{
public:
  typedef TinyVector<RealType,1> result_type;

  uniform_on_sphere() {};

  template <class RNG> result_type operator()(RNG& rng)
  {
    return result_type(rng()<0.5 ? -1. : 1.);
  }
};

template<class RealType>
class uniform_on_sphere<2,RealType>
{
public:
  typedef TinyVector<RealType,2> result_type;

  uniform_on_sphere() {};

  template <class RNG> result_type operator()(RNG& rng)
  {
    double phi=boost::uniform_real<>(0,2.*M_PI)(rng);
    return result_type(cos(phi),sin(phi));
  }
};

template<class RealType>
class uniform_on_sphere<3,RealType>
{
public:
  typedef TinyVector<RealType,3> result_type;

  uniform_on_sphere() {}

  template <class RNG> result_type operator()(RNG& rng)
  {
    boost::uniform_real<> float_random1(0.,2.*M_PI);
    double phi=float_random1(rng);
    double z=boost::uniform_real<>(-1.,1.)(rng);
    double r = sqrt(1.-z*z);
    return result_type(r*cos(phi),r*sin(phi),z);
  }
};




template <int N> class ONMoment;    

/**
 * Computes the change in the bond energy between m1 and m2, when m1 is changed
 * according to dir. J is the couplings matrix between the bound that connects
 * m1 to m2.
 */
template <class M>
inline double bond_energy_change(const ONMoment<M::dim>& m1,
                        const ONMoment<M::dim>& m2, M J, 
                        typename ONMoment<M::dim>::update_type dir);

/**
 * Computes the change in the site energy, when m is changed to dir.
 * 
 * @param m1 the current spin orientation 
 * @param h  the external magnetic field
 * @param dir the update
 * @return    the difference in the site energy
 */
template <int N>
inline double site_energy_change(const ONMoment<N> m,
               const TinyVector<double,N>& h,
               typename ONMoment<N>::update_type dir);
      
/**
 * Computes the change in the onsite energy when m is changed according to dir
 * 
 * @param m    the current spin direction 
 * @param D    the selfinteraction matrix for the current position
 * @param dir  describes the proposed update
 */
// template <class M>
// inline double onsite_energy_change(const ONMoment<M::dim> m, const M& D,
//                         typename ONMoment<M::dim>::update_type dir);

/**
 * Computes the onsite energy of the current spin. 
 * 
 * @param m the current spin orientation
 * @param D the selfinteraction matrix for the current position.
 */
template <class M>
inline double onsite_energy(const ONMoment<M::dim> m, const M& D);

/**
 * Returns the bond energy between two neighbouring nodes
 * 
 * @param m1 the spin vector at the current position
 * @param m2 the spin vector at a neighbouring site
 * @param D  the coupling matrix between the two sites.
 */ 
template <class M>
inline double bond_energy(const ONMoment<M::dim> m1, const ONMoment<M::dim> m2,
                        const M& J);

/**
 * Returns the site energy of the spin in the magnetic field
 *
 * @param m the current spin orientation
 * @param h the outer magnetic field
 */
template <int N>
inline double site_energy(const ONMoment<N> m, const TinyVector<double,N>& h);

template <int N>
class ONMoment {
  public:
    /* vector in N dimensions to save the current spin orientation */
    TinyVector<double, N> state_;
  
    /* return directly the rotated new vector, not the 
       description of the rotation.   */
    typedef TinyVector<double,N> update_type;
    
    typedef TinyVector<double,N> magnetization_type;
    
    static const int dim=N;
    
     /* initialize all components with the same value */
    ONMoment() : state_(1./sqrt(double(N))) {}
      
    template<class RNG> update_type random_update(RNG& rng) {
      update_type axis = uniform_on_sphere<N>()(rng);
      /* by default: return the axis
         should have the same distribution as the rotated vector */
      return axis;
    }
      
    void update(TinyVector<double, N> dir) {
      state_ -= 2.0*dot(dir,state_)*dir;
    }
      
    double project(typename ONMoment::update_type dir) {
      return double(dot(dir,state_));
    }

    friend alps::ODump& operator << (alps::ODump& dump, 
                        const ONMoment<N>& m){
      return dump << m.state_;
    }

    friend std::ostream& operator << (std::ostream& out, const ONMoment<N>& m)
    { return out << m.state_; }
      
    friend alps::IDump& operator >> (alps::IDump& dump, ONMoment<N>& m) {
      return dump >> m.state_;
    }
 
    double mag_h(TinyVector<double,N> h_) const 
    // it is assumed that h has unit length
    { return dot(h_,state_); }

    const magnetization_type& magnetization() const 
    { return state_; }
};

/* ONMoment<1> is implemented as IsingMoment Model */

template<>
class ONMoment<2> {
  public:
    TinyVector<double, 2> state_;
  
  public:
    typedef TinyVector<double,2> update_type;
    typedef TinyVector<double,2> magnetization_type;

    static const int dim=2;
    
    ONMoment() : state_(1./sqrt(2.0)) { }
      
    template<class RNG> update_type random_update(RNG& rng)
    {  return uniform_on_sphere<2>()(rng); }
      
    void update(TinyVector<double, 2> dir) {
      state_ -= 2.0*dot(dir,state_)*dir;
    }
      
    double project(update_type dir) {
      return double(dot(dir,state_));
    }

    friend alps::ODump& operator << (alps::ODump& dump, 
                        const ONMoment<2>& m) {
      return dump << m.state_;
    }

    friend std::ostream& operator << (std::ostream& out, const ONMoment<2>& m)
    { return out << m.state_; }

    friend alps::IDump& operator >> (alps::IDump& dump, ONMoment<2>& m) {
      return dump >> m.state_;
    }
 
    double mag_h(TinyVector<double,2> h_) const
    { return dot(h_,state_); }

    const magnetization_type& magnetization() const { 
      return state_;
    }
};

template<>
class ONMoment<3> {
  public:
    TinyVector<double, 3> state_;
  
    typedef TinyVector<double,3> update_type;
    typedef TinyVector<double,3> magnetization_type;
    
    static const int dim=3;

    ONMoment() : state_(1./sqrt(3.)) { }

/* here follow different versions of update types 
   uncomment the corresponding one. */
/*    
    // rotations a la Gould & Tobochnik
    // choose a vector with normally distributed components, add
    // it to the current string and normalize to 1.
    template<class RNG> update_type random_update(RNG& rng)
    { 
      double x = boost::normal_distribution<>()(rng);
      double y = boost::normal_distribution<>()(rng);
      double z = boost::normal_distribution<>()(rng);
      update_type tmp(x,y,z);
      tmp+=state_;
      double length = std::sqrt(dot(tmp,tmp));
      tmp/=length;
      return tmp;
    }
*/

    // direct random update -- choose a uniformely distributed 
    // unit vector in 3D 
    template<class RNG> update_type random_update(RNG& rng)
    {  return uniform_on_sphere<3>()(rng); }
/*
    // choose a uniformely distributed rotation axis and a uniformely
    // distributed rotation angle.
    // rotate the current spin around this axis by the chosen angle
    template<class RNG> update_type random_update(RNG& rng) {
      double phi = boost::uniform_real<>(0, 2*M_PI)(rng);
      // choose uniformly distributed rotation axis
      update_type axis = uniform_on_sphere<3>()(rng);
      // rotate current state around axis with angle phi
      double c = cos(phi);
      double s = sin(phi);
      double t = 1-c;
      // rotation matrix according to http://www.makegames.com/3drotation/
      double x = axis[0];
      double y = axis[1];
      double z = axis[2];
      double xNEW, yNEW, zNEW;
      xNEW = (t*x*x+c)*state_[0]+(t*x*y-s*z)*state_[1]+(t*x*z+s*y)*state_[2];
      yNEW = (t*x*y+s*z)*state_[0]+(t*y*y+c)*state_[1]+(t*y*z-s*x)*state_[2];
      zNEW = (t*x*z-s*y)*state_[0]+(t*y*z+s*x)*state_[1]+(t*z*z+c)*state_[2];
      update_type tmp(xNEW, yNEW, zNEW);
      return tmp;
    }
  */  
    void update(TinyVector<double,3> dir) {
      state_ -= 2.*dot(dir,state_)*dir;
    }

    double project(update_type dir) {
      return double(dot(dir,state_));
    }

    friend alps::ODump& operator << (alps::ODump& dump, 
                        const ONMoment<3>& m){
      return dump << m.state_;
    }

    friend std::ostream& operator << (std::ostream& out, const ONMoment<3>& m)
    { return out << m.state_; }

    friend alps::IDump& operator >> (alps::IDump& dump, ONMoment<3>& m) {
      return dump >> m.state_;
    }
 
    double mag_h(TinyVector<double,3> h_) const
    { return dot(h_,state_); }

    const magnetization_type& magnetization() const { 
      return state_;
    }
};

/**
 * Computes the energy of the spin  in the magnetic field.
 * 
 * \param m  The current spin orientation
 * \param h  The magnetic field
 */
template <int N>
inline double site_energy(const ONMoment<N> m, const TinyVector<double,N>& h)
{ return -dot(m.state_,h); };

template <int N>
inline double site_energy_change(const ONMoment<N> m,
               const TinyVector<double,N>& h,
               typename ONMoment<N>::update_type dir)
  // new direction is (m.state_ - 2.0*dot(dir,m.state_)*dir)
{ return (2.*dot(h,dir)*dot(m.state_,dir)); }

template <class M>
inline double bond_energy(const ONMoment<M::dim> m1, const ONMoment<M::dim> m2,
                          const M& J)
  /* Bond energy between the two spins */
{ return -J.vec_mat_vec(m1.state_,m2.state_); }

template <class M>
inline double bond_energy_change(const ONMoment<M::dim>& m1, 
                        const ONMoment<M::dim>& m2, M J,
                        typename ONMoment<M::dim>::update_type dir) {
  TinyVector<double,M::dim> tmp = dir*dot(dir,m1.state_);
  return 2.*J.vec_mat_vec(tmp, m2.state_);
}

template <int N>
inline double onsite_energy_change(const ONMoment<N> m,
               const SquareMatrix<double,N>& D, 
               typename ONMoment<N>::update_type dir)
{
  TinyVector<double,N> tmp = dir*(2.*dot(dir,m.state_));
  TinyVector<double,N> lhs = m.state_;
  lhs *= 2.0;
  lhs -= tmp;
  
  double res = D.vec_mat_vec(tmp,lhs);
  return res;
}

template <int N>
inline double onsite_energy_change(const ONMoment<N> m,
               const DiagMatrix<double,N>& D,
               typename ONMoment<N>::update_type dir)
{
  TinyVector<double,N> tmp = dir*(2.*dot(dir,m.state_));
  TinyVector<double,N> lhs = m.state_;
  lhs *= 2.0;
  lhs -= tmp;
  return D.vec_mat_vec(tmp,lhs);
}

template <int N>
inline double onsite_energy_change(const ONMoment<N>,
               const MIdMatrix<double,N>&,
               typename ONMoment<N>::update_type)
{  return 0.0; }

// If D is a multiple of the identity matrix, the onsite_energy is the 
// determinant of the matrix, as the spins always have length 1.
template <int N>
inline double onsite_energy(const ONMoment<N> m, SquareMatrix<double,N>& D)
{ return -D.vec_mat_vec(m.state_); }

template <int N>
inline double onsite_energy(const ONMoment<N> m, DiagMatrix<double,N>& D) 
{ return -D.vec_mat_vec(m.state_); }

template <int N>
inline double onsite_energy(const ONMoment<N>, MIdMatrix<double,N>&)
{ return -1.0; }

// O(3) models are often referred to by common name Heisenberg
typedef ONMoment<3> HeisenbergMoment;
typedef ONMoment<2> XYMoment;

// optimized implementation for O(1) and O(2) moments
template<> class ONMoment<1> : public IsingMoment {};
//template<> class ONMoment<2> : public XYMoment {};

#endif
