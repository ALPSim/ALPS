/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2007      by Sebastian Fuchs <fuchs@comp-phys.org>
 * Copyright (C) 2007-1013 by Emanuel Gull <egull@umich.edu>
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

#pragma once

#include "cluster_graph.h"
#include "cluster_helper.h"
#include "cluster_transform.h"
#include "DCA_transform.h"
#include "simple_DCA_transform.h"

class AntiFerroDCATransformer : public DCATransformer
{
public:
  AntiFerroDCATransformer(const alps::params& params) : DCATransformer(params) {}

  //self-consistency transform operator
  matsubara_green_function_t operator()(const matsubara_green_function_t &G_omega,
                                                                                 const matsubara_green_function_t &G0_omega,
                                                                                 matsubara_green_function_t &selfenergy,
                                                                                 double &mu, const double h, const double beta, const bool use_selfenergy);
  matsubara_green_function_t operator()(const matsubara_green_function_t &G_omega,
                                                                                 const matsubara_green_function_t &G0_omega,
                                                                                 double &mu, const double h, const double beta);
  nambu_matsubara_green_function_t operator()(const nambu_matsubara_green_function_t &G_omega,
                                                                                             const nambu_matsubara_green_function_t &G0_omega,
                                                                                             nambu_matsubara_green_function_t &selfenergy,
                                                                                             double &mu, const double h, const double eta1, const double eta2, const double beta, const bool use_selfenergy);
  itime_green_function_t symmetrize(const itime_green_function_t& G, const bool paramagnet) const
  { throw std::runtime_error("spin symmetrization in AFM not implemented for itime.");}//return spin_symmetrizer(G, paramagnet); }
  matsubara_green_function_t symmetrize(const matsubara_green_function_t& G, const bool paramagnet) const
  { return spin_symmetrizer(G, paramagnet); }
  matsubara_green_function_t transform_into_real_space(const matsubara_green_function_t& G) const
  {
      if(G.shape()!= blockdiagonal) throw std::logic_error("transform of GF to real space: shape incorrect.");
      return transform_coordinates<std::complex<double>, Backward>(G);
    }
  matsubara_green_function_t transform_into_k_space(const matsubara_green_function_t& G) const
  {
      if(G.shape()!= nondiagonal) throw std::logic_error("transform of GF to k space: shape incorrect.");
      return transform_coordinates<std::complex<double>, Forward>(G);
    }
  itime_green_function_t transform_into_real_space(const itime_green_function_t& G) const
  {
      if(G.shape()!= blockdiagonal) throw std::logic_error("transform of GF to real space: shape incorrect.");
      return transform_coordinates<double, Backward>(G);
    }
  itime_green_function_t transform_into_k_space(const itime_green_function_t& G) const
  {
      if(G.shape()!= nondiagonal) throw std::logic_error("transform of GF to k space: shape incorrect.");
      return transform_coordinates<double, Forward>(G);
    }

  nambu_matsubara_green_function_t transform_into_real_space(const nambu_matsubara_green_function_t& G) const {throw std::logic_error("AFM NAMBU not implemented.");}
  //{ return transform_coordinates<std::complex<double>, Backward>(G); }
  nambu_matsubara_green_function_t transform_into_k_space(const nambu_matsubara_green_function_t& G) const {throw std::logic_error("AFM NAMBU not implemented.");}
  //{ return transform_coordinates<std::complex<double>, Forward>(G); }
  nambu_itime_green_function_t transform_into_real_space(const nambu_itime_green_function_t& G) const {throw std::logic_error("AFM NAMBU not implemented.");}
  //{ return transform_coordinates<double, Backward>(G); }
  nambu_itime_green_function_t transform_into_k_space(const nambu_itime_green_function_t& G) const {throw std::logic_error("AFM NAMBU not implemented.");}
  //{ return transform_coordinates<double, Forward>(G); }


  virtual matsubara_green_function_t initial_G0(const alps::params& parms, shape_t shape=blockdiagonal);
  static matsubara_green_function_t green_function_inverse(const matsubara_green_function_t& G);
  static matsubara_green_function_t green_function_multiply(const matsubara_green_function_t& G1, const matsubara_green_function_t& G2);
  static matsubara_green_function_t green_function_add(const matsubara_green_function_t& G1, const matsubara_green_function_t& G2);
  static matsubara_green_function_t green_function_subtract(const matsubara_green_function_t& G1, const matsubara_green_function_t& G2);
  static matsubara_green_function_t green_function_divide(const matsubara_green_function_t& G1, const std::complex<double> d);

  std::map<double, ClusterHelper::Matrix<std::complex<double> >::Type>
  local_green_function_from_selfenergies(const std::map<double,
                                                                                   std::vector<ClusterHelper::Matrix<std::complex<double> >::Type> >& selfenergies,
                                                                                   const double mu, const double delta);

private:
  matsubara_green_function_t spin_symmetrizer(const matsubara_green_function_t& G, const bool paramagnet) const;
  template<typename T, class D>
  green_function<T> transform_coordinates(const green_function<T>& G) const;

  struct Forward;
  struct Backward;

};