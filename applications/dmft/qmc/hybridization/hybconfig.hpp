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
#ifndef HYB_CONFIG_HPP
#define HYB_CONFIG_HPP

#include<alps/ngs.hpp>
#include"hybfun.hpp"
#include"hybsegment.hpp"
#include"hybmatrix.hpp"
//this class knows everything to do with hybridization / bath operators.
//determinants of hybridization matrices, rank one updates, and so on.
//this is where \det \Delta is implemented.
class hybridization_configuration{
public:
  hybridization_configuration(const alps::params &p);
  hybridization_configuration(const hybridization_configuration &rhs) : Delta(rhs.Delta),hybmat_(rhs.hybmat_) {
    std::cerr << hybmat_.size() << std::endl;
    for (int i=0;i<hybmat_.size();i++)
      hybmat_[i].rebuild_hyb_matrix(i,Delta);
  }
  const hybridization_configuration &operator=(const hybridization_configuration &rhs) {
    hybmat_ = rhs.hybmat_;
    Delta = rhs.Delta;
    for (int i=0;i<hybmat_.size();i++)
      hybmat_[i].rebuild_hyb_matrix(i,Delta);
    return *this;
  }
  
  ~hybridization_configuration() {
//    std::cerr << "Deleting hybconfig\n";
  }
  
  double hyb_weight_change_insert(const segment &new_segment, int orbital);
  void insert_segment(const segment &new_segment, int orbital);
  void insert_antisegment(const segment &new_antisegment, int orbital);
  double hyb_weight_change_remove(const segment &new_segment, int orbital);
  void remove_segment(const segment &new_segment, int orbital);
  void remove_antisegment(const segment &new_antisegment, int orbital);
    
  void dump();
  void rebuild();
  void rebuild(int orbital);
  void rebuild(std::vector<int> orbital);

  void measure_G(std::vector<std::vector<double> > &G, std::vector<std::vector<double> > &F, const std::vector<std::map<double,double> > &F_prefactor, double sign) const;
  void measure_Gw(std::vector<std::vector<double> > &Gwr,std::vector<std::vector<double> > &Gwi,std::vector<std::vector<double> > &Fwr,std::vector<std::vector<double> > &Fwi, const std::vector<std::map<double,double> > &F_prefactor, double sign) const;
  void measure_G2w(std::vector<std::vector<std::complex<double> > > &G2w, std::vector<std::vector<std::complex<double> > > &F2w, int N_w2, int N_w_aux, const std::vector<std::map<double,double> > &F_prefactor) const;
  void measure_Gl(std::vector<std::vector<double> > &Gl,std::vector<std::vector<double> > &Fl, const std::vector<std::map<double,double> > &F_prefactor, double sign) const;
  double full_weight() const;

  friend std::ostream &operator<<(std::ostream &os, const hybridization_configuration &hyb_config);
private:
  //the hybridization function Delta
  hybfun Delta;
  std::vector<hybmatrix> hybmat_;
};
std::ostream &operator<<(std::ostream &os, const hybridization_configuration &hyb_config);

#endif

