  /****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2012 by Emanuel Gull <gull@pks.mpg.de>,
 *                   
 *  based on an earlier version by Philipp Werner and Emanuel Gull
 *
 *
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
 *
 *****************************************************************************/
#ifndef HYB_CONFIG_HPP
#define HYB_CONFIG_HPP

#include"hybfun.hpp"
#include"hybsegment.hpp"
#include"hybmatrix.hpp"
//this class knows everything to do with hybridization / bath operators.
//determinants of hybridization matrices, rank one updates, and so on.
//this is where \det \Delta is implemented.
class hybridization_configuration{
public:
  hybridization_configuration(const alps::params &p);
  double hyb_weight_change_insert(const segment &new_segment, int orbital);
  void insert_segment(const segment &new_segment, int orbital);
  void insert_antisegment(const segment &new_antisegment, int orbital);
  double hyb_weight_change_remove(const segment &new_segment, int orbital);
  void remove_segment(const segment &new_segment, int orbital);
  void remove_antisegment(const segment &new_antisegment, int orbital);
    
  void dump();

  void measure_G(std::vector<std::vector<double> > &G, double sign) const;
  void measure_Gw(std::vector<std::vector<double> > &Gwr,std::vector<std::vector<double> > &Gwi,std::vector<std::vector<double> > &Fwr,std::vector<std::vector<double> > &Fwi, const std::vector<std::map<double,double> > &F_prefactor, double sign) const;
  void measure_G2w(std::vector<std::vector<std::complex<double> > > &G2w, std::vector<std::vector<std::complex<double> > > &F2w, int N_w2, int N_w_aux, const std::vector<std::map<double,double> > &F_prefactor) const;
  void measure_Gl(std::vector<std::vector<double> > &Gl,std::vector<std::vector<double> > &Fl, const std::vector<std::map<double,double> > &F_prefactor, double sign) const;

  friend std::ostream &operator<<(std::ostream &os, const hybridization_configuration &hyb_config);
private:
  //the hybridization function Delta
  hybfun Delta;
  std::vector<hybmatrix> hybmat_;
};
std::ostream &operator<<(std::ostream &os, const hybridization_configuration &hyb_config);

#endif

