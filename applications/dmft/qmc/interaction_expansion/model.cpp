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
 * THIS SOFTWARE NEEDS AN APPROPRIATE LICENSE BLOCK HERE
 *****************************************************************************/

#include "interaction_expansion.hpp"


//the following functions work for a single site hubbard model, one site, one zone.
double HalfFillingHubbardInteractionExpansionRun::try_add()
{
  assert(n_zone==1);
  spin_t zone=0;
  //construct a creator and an annihilator in zone 0:
  double t = beta*random_01(); 
  double abs_w = beta*onsite_U*n_site; // onsite_U is const
  uint site = random_int(n_site);
  M[zone].creators().push_back(creator(up, site, t, n_matsubara));
  M[zone].annihilators().push_back(annihilator(up, site, t, n_matsubara));
  M[zone].alpha().push_back(random_01()<0.5?alpha:1-alpha); //symmetrized version
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
  spin_t zone=0;
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


double HalfFillingHubbardInteractionExpansionRun::try_remove(uint vertex_nr)
{
  //get weight
  double lambda_1 = fastupdate_down(vertex_nr, 0, true);  // true means compute_only_weight
  double pert_order=M[0].size();
  //return weight
  return  pert_order/(beta*onsite_U*n_site)*lambda_1*lambda_1;
}


void HalfFillingHubbardInteractionExpansionRun::perform_remove(uint vertex_nr)
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
  assert(n_zone==2);
  double t = beta*random_01(); 
  double abs_w = beta*onsite_U*n_site; // onsite_U is const
  uint site = random_int(n_site);
  double alpha0 = random_01()<0.5?alpha:1-alpha;
  double alpha1 = 1 - alpha0;
  spin_t zone0=0;
  spin_t zone1=1;
  M[0].creators().push_back(creator(zone0, site, t, n_matsubara));
  M[0].annihilators().push_back(annihilator(zone0, site, t, n_matsubara));
  M[0].alpha().push_back(alpha0); //symmetrized version
  M[1].creators().push_back(creator(zone1, site, t, n_matsubara));
  M[1].annihilators().push_back(annihilator(zone1, site, t,n_matsubara));
  M[1].alpha().push_back(alpha1); //symmetrized version
  //keep track of vertex list
  vertices.push_back(vertex(zone0, site, M[0].creators().size()-1, M[0].annihilators().size()-1, 
			    zone1, site, M[1].creators().size()-1, M[1].annihilators().size()-1, abs_w)); 
  //perform fastupdate up for weight
  double lambda0=fastupdate_up(zone0, true); // true means compute_only_weight
  double lambda1=fastupdate_up(zone1, true);
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


double HubbardInteractionExpansionRun::try_remove(uint vertex_nr)
{
  //get weight
  double lambda0 = fastupdate_down(vertex_nr, 0, true);  // true means compute_only_weight
  double lambda1 = fastupdate_down(vertex_nr, 1, true);  
  double pert_order=M[0].size();
  //return weight
  return  -pert_order/(beta*onsite_U*n_site)*lambda0*lambda1;
}


void HubbardInteractionExpansionRun::perform_remove(uint vertex_nr)
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
  assert(n_site==1 && n_zone >1);
  spin_t zone1=(spin_t)(random_01()*n_zone);
  spin_t zone2;
  
  do{ zone2=(spin_t)(random_01()*n_zone);
  } while (U(zone1, zone2)==0);
  
  double t=beta*random_01();

  double abs_w=beta*U(zone1, zone2)/2;
  double alpha1 = random_01()<0.5 ? alpha : 1-alpha;
  double alpha2 = 1-alpha1;
  site_t site1=0;
  site_t site2=0;
  M[zone1].creators().    push_back(creator(zone1,site1,t, n_matsubara));
  M[zone1].annihilators().push_back(annihilator(zone1,site1,t,n_matsubara));
  M[zone1].alpha().       push_back(alpha1); //symmetrized version
  M[zone2].creators().    push_back(creator(zone2,site2,t, n_matsubara));
  M[zone2].annihilators().push_back(annihilator(zone2,site2,t,n_matsubara));
  M[zone2].alpha().       push_back(alpha2); //symmetrized version
  
  vertices.push_back(vertex(zone1, site1, M[zone1].creators().size()-1, M[zone1].annihilators().size()-1, 
			    zone2, site2, M[zone2].creators().size()-1, M[zone2].annihilators().size()-1, abs_w));
  
  double lambda1=fastupdate_up(zone1, true);
  double lambda2=fastupdate_up(zone2, true);
  double sym_factor = U.n_nonzero();
  //minus sign because we're not working with down holes but down electrons -> everything picks up a minus sign.
  double metropolis_weight=-abs_w*sym_factor/(vertices.size())*lambda1*lambda2;
  return metropolis_weight;
}


void MultiBandDensityHubbardInteractionExpansionRun::perform_add()
{	
  //find zones
  spin_t zone1=vertices.back().zone1();
  spin_t zone2=vertices.back().zone2();
  //perform the fastupdate up move
  fastupdate_up(zone1,false);
  fastupdate_up(zone2,false);
}


void MultiBandDensityHubbardInteractionExpansionRun::reject_add()
{
  //find zones
  spin_t zone1=vertices.back().zone1();
  spin_t zone2=vertices.back().zone2();
  //get rid of the operators
  M[zone1].creators().pop_back();
  M[zone1].annihilators().pop_back();
  M[zone1].alpha().pop_back();
  M[zone2].creators().pop_back();
  M[zone2].annihilators().pop_back();
  M[zone2].alpha().pop_back();
  //get rid of the vertex from vertex list
  vertices.pop_back();
}


double MultiBandDensityHubbardInteractionExpansionRun::try_remove(uint vertex_nr)
{
  spin_t zone1=vertices[vertex_nr].zone1();
  spin_t zone2=vertices[vertex_nr].zone2();
  if(zone1==zone2) {
    std::cerr<<"bug: zone1 and zone2 are equal, we'd require a two vertex removal move for that!"<<std::endl; 
    abort();
  }
  //find operator positions in that zone
  uint operator_nr_1=vertices[vertex_nr].c_dagger_1();
  uint operator_nr_2=vertices[vertex_nr].c_dagger_2();
  //get weight
  double abs_w=vertices[vertex_nr].abs_w();
  double lambda_1 = fastupdate_down(operator_nr_1, zone1, true);
  double lambda_2 = fastupdate_down(operator_nr_2, zone2, true);
  double pert_order=vertices.size();
  
  double sym_factor = U.n_nonzero();
  double metropolis_weight = -pert_order/abs_w/sym_factor*lambda_1*lambda_2;

  return  metropolis_weight; 
}


void MultiBandDensityHubbardInteractionExpansionRun::perform_remove(uint vertex_nr)
{
  //find zones
  spin_t zone1=vertices[vertex_nr].zone1();
  spin_t zone2=vertices[vertex_nr].zone2();
  //find operator positions in that zone
  uint operator_nr_1=vertices[vertex_nr].c_dagger_1();
  uint operator_nr_2=vertices[vertex_nr].c_dagger_2();
  //perform fastupdate down
  fastupdate_down(operator_nr_1, zone1, false);
  fastupdate_down(operator_nr_2, zone2, false);
  //take care of vertex list
  for(int i=vertices.size()-1;i>=0;--i){
    //this operator pointed to the last row/column of M[zone1], which has just been moved to vertex_nr.
    if(vertices[i].zone1()==zone1 && vertices[i].c_dagger_1()==M[zone1].size()){ 
      vertices[i].c_dagger_1()=operator_nr_1;
      vertices[i].c_1()=operator_nr_1;
      break;
    }
    if(vertices[i].zone2()==zone1 && vertices[i].c_dagger_2()==M[zone1].size()){
      vertices[i].c_dagger_2()=operator_nr_1;
      vertices[i].c_2()=operator_nr_1;
      break;
    }
  }
  for(int i=vertices.size()-1;i>=0;--i){
    if(vertices[i].zone1()==zone2 && vertices[i].c_dagger_1()==M[zone2].size()){
      vertices[i].c_dagger_1()=operator_nr_2;
      vertices[i].c_1()=operator_nr_2;
      break;
    }
    if(vertices[i].zone2()==zone2 && vertices[i].c_dagger_2()==M[zone2].size()){
      vertices[i].c_dagger_2()=operator_nr_2;
      vertices[i].c_2()=operator_nr_2;
      break;
    }
  }
  vertices[vertex_nr]=vertices.back();

  //get rid of operators
  M[zone1].creators().pop_back();
  M[zone1].annihilators().pop_back();
  M[zone1].alpha().pop_back();
  M[zone2].creators().pop_back();
  M[zone2].annihilators().pop_back();
  M[zone2].alpha().pop_back();
  //get rid of vertex list entries
  vertices.pop_back();

}


void MultiBandDensityHubbardInteractionExpansionRun::reject_remove()
{
  //do nothing
  return;
}


