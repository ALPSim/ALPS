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

extern "C" void zaxpy_(const unsigned int *size, const void *alpha, const void *x, const int *incx, void *y, const int *incy);

/// @brief This function performs the InteractionExpansion update of adding one vertex to the set of vertices. If the 
///        Green's function for the measurement has to be tracked it calls the function that does that.
/// @param vertex_nr specify which of the U n_up n_down you want to look at.
/// @param compute_only_weight Do not perform the fastupdate, but only return the Mnn entry (1/lambda), which is used for the acceptance weight.
/// @param track_green_matsubara Track Green function in Matsubara frequencies or imaginary time if required.

double InteractionExpansionRun::fastupdate_up(const int flavor, bool compute_only_weight)
{
  unsigned int noperators=M[flavor].size(); 
  //current size of M: number of vertices -1. We need to add the last vertex. 
  //A pointer to the creator and annihilator is already stored in creators_ and annihilators_ at position M[flavor].size();
  double Green0_n_n=green0_spline(M[flavor].creators()[noperators],M[flavor].annihilators()[noperators]); 
  double Green0_n_j[noperators];
  double Green0_j_n[noperators];
  double lastcolumn[noperators];
  double lastrow[noperators];
  double lambda;
  //compute the Green's functions with interpolation
  for(unsigned int i=0;i<noperators;++i){
    Green0_n_j[i]=green0_spline(M[flavor].creators()[noperators],M[flavor].annihilators()[i]);
    Green0_j_n[i]=green0_spline(M[flavor].creators()[i],M[flavor].annihilators()[noperators]);
  }
  //compute the last row
  M[flavor].right_multiply(Green0_j_n, lastcolumn);
  //compute lambda
  double ip = (noperators==0?0:inner_prod(Green0_n_j, lastcolumn, noperators));
  lambda = Green0_n_n - ip + M[flavor].alpha()[noperators];
  //return weight if we have nothing else to do
  if(compute_only_weight){
    return lambda;
  }
  //compute last column
  M[flavor].left_multiply(Green0_n_j, lastrow);
  //compute norm of vector and mean for roundoff error check
  if(noperators > 0){
    scale(1./lambda, lastrow, noperators);
    M[flavor].add_outer_prod(lastcolumn, lastrow);
    scale(1./lambda, lastcolumn, noperators);
    //std::cout<<lambda<<" "<<M[flavor]<<std::endl;
  }
  //add row and column to M
  M[flavor].resize(noperators+1);
  for(unsigned int i=0;i<noperators;++i){
    M[flavor](i, noperators)=-lastcolumn[i];
  }
  for(unsigned int i=0;i<noperators;++i){
    M[flavor](noperators, i)=-lastrow[i];
  }
  M[flavor](noperators, noperators)=1./lambda;
  return lambda;
}



///Fastupdate formulas, remove order by one (remove a vertex). If necessary
///also take track of the Green's function.
double InteractionExpansionRun::fastupdate_down(const int operator_nr, const int flavor, bool compute_only_weight)
{
  //perform updates according to formula 21.1, 21.2
  unsigned int noperators=M[flavor].size();  //how many operators do we have in total?
  if(compute_only_weight){
    return M[flavor](operator_nr,operator_nr);
  }
  //swap rows and colums of M <-> move selected vertex to the end.
  for(unsigned int i=0;i<noperators;++i){
    double tmp=M[flavor](i,noperators-1);
    M[flavor](i,noperators-1)=M[flavor](i,operator_nr);
    M[flavor](i,operator_nr)=tmp;
  }
  for(unsigned int i=0;i<noperators;++i){
    double tmp=M[flavor](noperators-1,i);
    M[flavor](noperators-1,i)=M[flavor](operator_nr,i);
    M[flavor](operator_nr, i)=tmp;
  }
  //swap creator and annihilator
  std::swap(M[flavor].creators()[operator_nr],     M[flavor].creators()[noperators-1]);
  std::swap(M[flavor].annihilators()[operator_nr], M[flavor].annihilators()[noperators-1]);
  std::swap(M[flavor].alpha()[operator_nr],        M[flavor].alpha()[noperators-1]);
  double Mnn=M[flavor](noperators-1,noperators-1);
  //now perform fastupdate of M
  double lastrow[noperators-1];
  double lastcolumn[noperators-1];
  for(unsigned int j=0;j<noperators-1;++j){
    lastrow[j]=M[flavor](noperators-1, j);
    lastcolumn[j]=M[flavor](j,noperators-1);
  }
  scale(-1./Mnn, lastcolumn, noperators-1);
  M[flavor].resize(noperators-1);  //lose the last row and last column, reduce size by one, but keep contents.
  M[flavor].add_outer_prod(lastcolumn, lastrow);
  return Mnn;  //the determinant ratio det D_{k-1}/det D_{k}
}




