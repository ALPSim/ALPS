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

#include "interaction_expansion.hpp"


//the following functions work for a single site hubbard model, one site, one flavor.
double HalfFillingHubbardInteractionExpansionRun::try_add()
{
  assert(n_flavors==1);
  spin_t flavor=0;
  //construct a creator and an annihilator in flavor 0:
  double t = beta*random_01(); 
  double abs_w = beta*onsite_U*n_site; // onsite_U is const

  unsigned int site = random_int(n_site);
  M[flavor].creators().push_back(creator(up, site, t, n_matsubara));
  M[flavor].annihilators().push_back(annihilator(up, site, t, n_matsubara));
  M[flavor].alpha().push_back(random_01()<0.5?alpha:1-alpha); //symmetrized version

  //keep track of vertex list
  vertices.push_back(vertex(0, site, M[0].creators().size()-1, M[0].annihilators().size()-1, 0, 0, 0, 0, abs_w)); 
  //second part of list is ignored (symmetry)
  //perform fastupdate up for weight
  double lambda=fastupdate_up(0, true); // true means compute_only_weight
  double metropolis_weight=abs_w/(vertices.size())*lambda*lambda;
  //return weight
  return metropolis_weight;
}


void HalfFillingHubbardInteractionExpansionRun::perform_add()
{
  //perform the fastupdate up move
  fastupdate_up(0,false);
}


void HalfFillingHubbardInteractionExpansionRun::reject_add()
{
  //get rid of the operators
  M[0].creators().pop_back();
  M[0].annihilators().pop_back();
  M[0].alpha().pop_back();
  //get rid of the vertex from vertex list
  vertices.pop_back();
}


double HalfFillingHubbardInteractionExpansionRun::try_remove(unsigned int vertex_nr)
{
  assert(num_rows(M[0].matrix()) == num_cols(M[0].matrix()));
  //get weight
  double lambda_1 = fastupdate_down(vertex_nr, 0, true);  // true means compute_only_weight
  double pert_order=num_rows(M[0].matrix());
  //return weight
  return  pert_order/(beta*onsite_U*n_site)*lambda_1*lambda_1;
}


void HalfFillingHubbardInteractionExpansionRun::perform_remove(unsigned int vertex_nr)
{
  //perform fastupdate down
  fastupdate_down(vertex_nr, 0, false);  // false means really perform, not only compute weight
  //get rid of operators
  M[0].creators().pop_back();
  M[0].annihilators().pop_back();
  M[0].alpha().pop_back();
  //get rid of vertex list entries
  vertices.pop_back();
}


void HalfFillingHubbardInteractionExpansionRun::reject_remove()
{
  //do nothing
  return;
}


double HubbardInteractionExpansionRun::try_add()
{
  assert(n_flavors==2);
  double t = beta*random_01(); 
  double abs_w = beta*onsite_U*n_site; // onsite_U is const
  unsigned int site = random_int(n_site);
  double alpha0 = random_01()<0.5?alpha:1-alpha;
  double alpha1 = 1 - alpha0;
  spin_t flavor0=0;
  spin_t flavor1=1;
  M[0].creators().push_back(creator(flavor0, site, t, n_matsubara));
  M[0].annihilators().push_back(annihilator(flavor0, site, t, n_matsubara));
  M[0].alpha().push_back(alpha0); //symmetrized version
  M[1].creators().push_back(creator(flavor1, site, t, n_matsubara));
  M[1].annihilators().push_back(annihilator(flavor1, site, t,n_matsubara));
  M[1].alpha().push_back(alpha1); //symmetrized version
  //keep track of vertex list
  vertices.push_back(vertex(flavor0, site, M[0].creators().size()-1, M[0].annihilators().size()-1, 
                        flavor1, site, M[1].creators().size()-1, M[1].annihilators().size()-1, abs_w)); 
  //perform fastupdate up for weight
  double lambda0=fastupdate_up(flavor0, true); // true means compute_only_weight
  double lambda1=fastupdate_up(flavor1, true);
  //std::cout<<"lambda: "<<lambda<<std::endl;
  double metropolis_weight=-abs_w/(vertices.size())*lambda0*lambda1;
  //return weight
  return metropolis_weight;
}


void HubbardInteractionExpansionRun::perform_add()
{
  //perform the fastupdate up move
  fastupdate_up(0,false);
  fastupdate_up(1,false);
}


void HubbardInteractionExpansionRun::reject_add()
{
  //get rid of the operators
  M[0].creators().pop_back();
  M[0].annihilators().pop_back();
  M[0].alpha().pop_back();
  M[1].creators().pop_back();
  M[1].annihilators().pop_back();
  M[1].alpha().pop_back();
  //get rid of the vertex from vertex list
  vertices.pop_back();
}


double HubbardInteractionExpansionRun::try_remove(unsigned int vertex_nr)
{
  assert(num_rows(M[0].matrix()) == num_cols(M[0].matrix()));
  //get weight
  double lambda0 = fastupdate_down(vertex_nr, 0, true);  // true means compute_only_weight
  double lambda1 = fastupdate_down(vertex_nr, 1, true);  
  double pert_order=num_rows(M[0].matrix());
  //return weight
  return  -pert_order/(beta*onsite_U*n_site)*lambda0*lambda1;
}


void HubbardInteractionExpansionRun::perform_remove(unsigned int vertex_nr)
{
  //perform fastupdate down
  fastupdate_down(vertex_nr, 0, false);  // false means really perform, not only compute weight
  fastupdate_down(vertex_nr, 1, false);  // false means really perform, not only compute weight
  //get rid of operators
  M[0].creators().pop_back();
  M[0].annihilators().pop_back();
  M[0].alpha().pop_back();
  M[1].creators().pop_back();
  M[1].annihilators().pop_back();
  M[1].alpha().pop_back();
  //get rid of vertex list entries
  vertices.pop_back();
}


void HubbardInteractionExpansionRun::reject_remove()
{
  //do nothing
  return;
}


//A term U n_i n_j of the series expansion. i != j.  E.g. for
//orbitals with G0(c_i, c_j)=0.
double MultiBandDensityHubbardInteractionExpansionRun::try_add()
{
  assert(n_site==1 && n_flavors >1);
  spin_t flavor1=(spin_t)(random_01()*n_flavors);
  spin_t flavor2;
  
  do{ flavor2=(spin_t)(random_01()*n_flavors);
  } while (U(flavor1, flavor2)==0);
  
  double t=beta*random_01();

  double abs_w=beta*U(flavor1, flavor2)/2;
  double alpha1 = random_01()<0.5 ? alpha : 1-alpha;
  double alpha2 = 1-alpha1;
  site_t site1=0;
  site_t site2=0;
  M[flavor1].creators().    push_back(creator(flavor1,site1,t, n_matsubara));
  M[flavor1].annihilators().push_back(annihilator(flavor1,site1,t,n_matsubara));
  M[flavor1].alpha().       push_back(alpha1); //symmetrized version
  M[flavor2].creators().    push_back(creator(flavor2,site2,t, n_matsubara));
  M[flavor2].annihilators().push_back(annihilator(flavor2,site2,t,n_matsubara));
  M[flavor2].alpha().       push_back(alpha2); //symmetrized version
  
  vertices.push_back(vertex(flavor1, site1, M[flavor1].creators().size()-1, M[flavor1].annihilators().size()-1, 
                flavor2, site2, M[flavor2].creators().size()-1, M[flavor2].annihilators().size()-1, abs_w));
  
  double lambda1=fastupdate_up(flavor1, true);
  double lambda2=fastupdate_up(flavor2, true);
  double sym_factor = U.n_nonzero();
  //minus sign because we're not working with down holes but down electrons -> everything picks up a minus sign.
  double metropolis_weight=-abs_w*sym_factor/(vertices.size())*lambda1*lambda2;
  return metropolis_weight;
}


void MultiBandDensityHubbardInteractionExpansionRun::perform_add()
{    
  //find flavors
  spin_t flavor1=vertices.back().flavor1();
  spin_t flavor2=vertices.back().flavor2();
  //perform the fastupdate up move
  fastupdate_up(flavor1,false);
  fastupdate_up(flavor2,false);
}


void MultiBandDensityHubbardInteractionExpansionRun::reject_add()
{
  //find flavors
  spin_t flavor1=vertices.back().flavor1();
  spin_t flavor2=vertices.back().flavor2();
  //get rid of the operators
  M[flavor1].creators().pop_back();
  M[flavor1].annihilators().pop_back();
  M[flavor1].alpha().pop_back();
  M[flavor2].creators().pop_back();
  M[flavor2].annihilators().pop_back();
  M[flavor2].alpha().pop_back();
  //get rid of the vertex from vertex list
  vertices.pop_back();
}


double MultiBandDensityHubbardInteractionExpansionRun::try_remove(unsigned int vertex_nr)
{
  spin_t flavor1=vertices[vertex_nr].flavor1();
  spin_t flavor2=vertices[vertex_nr].flavor2();
  if(flavor1==flavor2) {
    std::cerr<<"bug: flavor1 and flavor2 are equal, we'd require a two vertex removal move for that!"<<std::endl; 
    abort();
  }
  //find operator positions in that flavor
  unsigned int operator_nr_1=vertices[vertex_nr].c_dagger_1();
  unsigned int operator_nr_2=vertices[vertex_nr].c_dagger_2();
  //get weight
  double abs_w=vertices[vertex_nr].abs_w();
  double lambda_1 = fastupdate_down(operator_nr_1, flavor1, true);
  double lambda_2 = fastupdate_down(operator_nr_2, flavor2, true);
  double pert_order=vertices.size();
  
  double sym_factor = U.n_nonzero();
  double metropolis_weight = -pert_order/abs_w/sym_factor*lambda_1*lambda_2;

  return  metropolis_weight; 
}


void MultiBandDensityHubbardInteractionExpansionRun::perform_remove(unsigned int vertex_nr)
{
  //find flavors
  spin_t flavor1=vertices[vertex_nr].flavor1();
  spin_t flavor2=vertices[vertex_nr].flavor2();
  //find operator positions in that flavor
  unsigned int operator_nr_1=vertices[vertex_nr].c_dagger_1();
  unsigned int operator_nr_2=vertices[vertex_nr].c_dagger_2();
  //perform fastupdate down
  fastupdate_down(operator_nr_1, flavor1, false);
  fastupdate_down(operator_nr_2, flavor2, false);
  assert(num_rows(M[flavor1].matrix()) == num_cols(M[flavor1].matrix()));
  //take care of vertex list
  for(int i=vertices.size()-1;i>=0;--i){
    //this operator pointed to the last row/column of M[flavor1], which has just been moved to vertex_nr.
    if(vertices[i].flavor1()==flavor1 && vertices[i].c_dagger_1()==num_rows(M[flavor1].matrix())){ 
      vertices[i].c_dagger_1()=operator_nr_1;
      vertices[i].c_1()=operator_nr_1;
      break;
    }
    if(vertices[i].flavor2()==flavor1 && vertices[i].c_dagger_2()==num_rows(M[flavor1].matrix())){
      vertices[i].c_dagger_2()=operator_nr_1;
      vertices[i].c_2()=operator_nr_1;
      break;
    }
  }
  for(int i=vertices.size()-1;i>=0;--i){
    if(vertices[i].flavor1()==flavor2 && vertices[i].c_dagger_1()==num_rows(M[flavor2].matrix())){
      vertices[i].c_dagger_1()=operator_nr_2;
      vertices[i].c_1()=operator_nr_2;
      break;
    }
    if(vertices[i].flavor2()==flavor2 && vertices[i].c_dagger_2()==num_rows(M[flavor2].matrix())){
      vertices[i].c_dagger_2()=operator_nr_2;
      vertices[i].c_2()=operator_nr_2;
      break;
    }
  }
  vertices[vertex_nr]=vertices.back();

  //get rid of operators
  M[flavor1].creators().pop_back();
  M[flavor1].annihilators().pop_back();
  M[flavor1].alpha().pop_back();
  M[flavor2].creators().pop_back();
  M[flavor2].annihilators().pop_back();
  M[flavor2].alpha().pop_back();
  //get rid of vertex list entries
  vertices.pop_back();

}


void MultiBandDensityHubbardInteractionExpansionRun::reject_remove()
{
  //do nothing
  return;
}


