/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Sebastian Fuchs <fuchs@theorie.physik.uni-goettingen.de>
 *                              Matthias Troyer <troyer@comp-phys.org>
 *
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

#include "interaction_expansion.hpp"

///The basic updates for the InteractionExpansion algorithm: adding and removing vertices.
///This is the heart of InteractionExpansion's code.
void InteractionExpansion::interaction_expansion_step(void)
{
  int pert_order=vertices.size();   //current order of perturbation series
  double metropolis_weight=0.;
  //double oldsign=sign;
  static unsigned int i=0; ++i;    
  if(random()<0.5){  //trying to ADD vertex
    if(vertices.size()>=max_order) 
      return; //we have already reached the highest perturbation order
    metropolis_weight=try_add();
    if(fabs(metropolis_weight)> random()){
      measurements["VertexInsertion"]<<1.;
      perform_add();
      sign*=metropolis_weight<0?-1:1;
    }else{
      measurements["VertexInsertion"]<<0.;
      reject_add();
    }
  }else{ // try to REMOVE a vertex
    pert_order=vertices.size(); //choose a vertex
    if(pert_order < 1)
      return;     //we have an empty list or one with just one vertex
    //this might be the ideal place to do some cleanup, e.g. get rid of the roundoff errors and such.
    int vertex_nr=(int)(random() * pert_order);
    metropolis_weight=try_remove(vertex_nr); //get the determinant ratio. don't perform fastupdate yet
    if(fabs(metropolis_weight)> random()){ //do the actual update
      measurements["VertexRemoval"]<<1.;
      perform_remove(vertex_nr);
      sign*=metropolis_weight<0?-1:1;
    }else{
      measurements["VertexRemoval"]<<0.;
      reject_remove();
    }
  }//end REMOVE
  weight=metropolis_weight;
}

///Every now and then we have to recreate M from scratch to avoid roundoff
///error. This is done by iserting the vertices starting from zero.
void InteractionExpansion::reset_perturbation_series()
{
  std::vector<inverse_m_matrix> M2(M); //make a copy of M
  vertex_array vertices_backup;
  for(unsigned int i=0;i<vertices.size();++i){
    vertices_backup.push_back(vertices[i]);
  }
  vertices.clear();
  sign=1;
  for(spin_t z=0;z<n_flavors;++z){
    resize(M[z].matrix(),0,0);
  }
  green_matsubara = bare_green_matsubara;
  green_itime     = bare_green_itime;
  //recompute M from scratch
  for(unsigned int i=0;i<vertices_backup.size();++i){
    vertices.push_back(vertices_backup[i]);
    perform_add();
  }
  for(unsigned int z=0;z<M2.size();++z){
    double max_diff=0;
    for(unsigned int j=0;j<num_cols(M2[z].matrix());++j){
      for(unsigned int i=0;i<num_rows(M2[z].matrix());++i){
        double diff=M[z].matrix()(i,j)-M2[z].matrix()(i,j);
        if(std::abs(diff)>max_diff) max_diff=std::abs(diff);
      }
    }
    if(max_diff > 1.e-8)
      std::cout<<"WARNING: roundoff errors in flavor: "<<z<<" max diff "<<max_diff<<std::endl;
  }
}

