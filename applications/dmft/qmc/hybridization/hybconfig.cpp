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

#include"hybconfig.hpp"

hybridization_configuration::hybridization_configuration(const alps::params &p):
  Delta(p),
  hybmat_((int)(p["N_ORBITALS"]), p)
{
}

void hybridization_configuration::dump() {
    for (int i=0;i<hybmat_.size();i++)
        std::cerr << "Weight for orbital " << i << " : " << hybmat_[i].full_weight() << std::endl;
}

void hybridization_configuration::rebuild() {
  for (int i=0;i<hybmat_.size();i++)
    hybmat_[i].rebuild_hyb_matrix(i,Delta);
}

void hybridization_configuration::rebuild(int orbital) {
  hybmat_[orbital].rebuild_hyb_matrix(orbital,Delta);
}

void hybridization_configuration::rebuild(std::vector<int> orbital) {
  for (int i=0;i<orbital.size();i++)
    hybmat_[orbital[i]].rebuild_hyb_matrix(orbital[i],Delta);
}

double hybridization_configuration::hyb_weight_change_insert(const segment &new_segment, int orbital){
  return hybmat_[orbital].hyb_weight_change_insert(new_segment, orbital, Delta); //hand this off to the determinant matrix
}
void hybridization_configuration::insert_segment(const segment &new_segment, int orbital){
  //std::cout<<clmagenta<<"before insert recompute "<<cblack<<std::endl;
  //if(hybmat_[orbital].size()>0) hybmat_[orbital].rebuild_hyb_matrix(orbital, Delta);
  //std::cout<<clmagenta<<"done before insert recompute "<<cblack<<std::endl;
  hybmat_[orbital].insert_segment(new_segment, orbital); //hand this off to the determinant matrix
  //std::cout<<clmagenta<<"after insert recompute "<<cblack<<std::endl;
  //hybmat_[orbital].rebuild_hyb_matrix(orbital, Delta);
  //std::cout<<clmagenta<<"done after insert recompute "<<cblack<<std::endl;
}
double hybridization_configuration::hyb_weight_change_remove(const segment &new_segment, int orbital){
  return hybmat_[orbital].hyb_weight_change_remove(new_segment, orbital, Delta); //hand this off to the determinant matrix
}
void hybridization_configuration::remove_segment(const segment &new_segment, int orbital){
  //std::cout<<clmagenta<<"before remove recompute "<<cblack<<std::endl;
  //hybmat_[orbital].rebuild_hyb_matrix(orbital, Delta);
  //std::cout<<clmagenta<<"done before remove recompute "<<cblack<<std::endl;
  hybmat_[orbital].remove_segment(new_segment, orbital); //hand this off to the determinant matrix
  //std::cout<<clmagenta<<"after remove recompute "<<cblack<<std::endl;
  //if(hybmat_[orbital].size()>0) hybmat_[orbital].rebuild_hyb_matrix(orbital, Delta);
  //std::cout<<clmagenta<<"done after remove recompute "<<cblack<<std::endl;
}
void hybridization_configuration::remove_antisegment(const segment &new_antisegment, int orbital){
  hybmat_[orbital].remove_segment(new_antisegment, orbital); //hand this off to the determinant matrix
  //std::cout<<clmagenta<<"after as remove recompute "<<cblack<<std::endl;
  if(hybmat_[orbital].size()>0) hybmat_[orbital].rebuild_hyb_matrix(orbital, Delta);
  //std::cout<<clmagenta<<"done after as remove recompute "<<cblack<<std::endl;
}
void hybridization_configuration::insert_antisegment(const segment &new_antisegment, int orbital){
  hybmat_[orbital].insert_segment(new_antisegment, orbital); //hand this off to the determinant matrix
  //std::cout<<clmagenta<<"after as insert recompute "<<cblack<<std::endl;
  //hybmat_[orbital].rebuild_hyb_matrix(orbital, Delta);
  //std::cout<<clmagenta<<"done after as insert recompute "<<cblack<<std::endl;
}
void hybridization_configuration::measure_G(std::vector<std::vector<double> > &G, std::vector<std::vector<double> > &F, const std::vector<std::map<double,double> > &F_prefactor, double sign) const{
  for(std::size_t orbital=0;orbital<hybmat_.size();++orbital){
    hybmat_[orbital].measure_G(G[orbital], F[orbital], F_prefactor[orbital], sign);
  }
}
void hybridization_configuration::measure_Gw(std::vector<std::vector<double> > &Gwr, std::vector<std::vector<double> > &Gwi, std::vector<std::vector<double> > &Fwr, std::vector<std::vector<double> > &Fwi, const std::vector<std::map<double,double> > &F_prefactor, double sign) const{
  for(std::size_t orbital=0;orbital<hybmat_.size();++orbital){
    hybmat_[orbital].measure_Gw(Gwr[orbital], Gwi[orbital], Fwr[orbital], Fwi[orbital], F_prefactor[orbital], sign);
  }
}
void hybridization_configuration::measure_G2w(std::vector<std::vector<std::complex<double> > >&G2w, std::vector<std::vector<std::complex<double> > > &F2w, int N_w2, int N_w_aux, const std::vector<std::map<double,double> > &F_prefactor) const{
  for(std::size_t orbital=0;orbital<hybmat_.size();++orbital){
    hybmat_[orbital].measure_G2w(G2w[orbital], F2w[orbital], N_w2, N_w_aux, F_prefactor[orbital]);
  }
}
void hybridization_configuration::measure_Gl(std::vector<std::vector<double> > &Gl, std::vector<std::vector<double> > &Fl, const std::vector<std::map<double,double> > &F_prefactor, double sign) const{
  for(std::size_t orbital=0;orbital<hybmat_.size();++orbital){
    hybmat_[orbital].measure_Gl(Gl[orbital], Fl[orbital], F_prefactor[orbital], sign);
  }
}
double hybridization_configuration::full_weight() const{
  double weight=1.;
  for(std::size_t orbital=0;orbital<hybmat_.size();++orbital){
    weight*=hybmat_[orbital].full_weight();
  }
  return weight;
}
std::ostream &operator<<(std::ostream &os, const hybridization_configuration &hyb_config){
  for(std::size_t i=0;i<hyb_config.hybmat_.size();++i){
    os<<cblue<<"------- "<<"orbital: "<<i<<" ------"<<cblack<<std::endl;
    os<<hyb_config.hybmat_[i]<<std::endl;
  }
  return os;
}
