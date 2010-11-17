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
  std::vector<double> Green0_n_j(noperators);
  std::vector<double> Green0_j_n(noperators);
  std::vector<double> lastcolumn(noperators);
  std::vector<double> lastrow(noperators);
  double lambda;
  //compute the Green's functions with interpolation
  for(unsigned int i=0;i<noperators;++i){
    Green0_n_j[i]=green0_spline(M[flavor].creators()[noperators],M[flavor].annihilators()[i]);
    Green0_j_n[i]=green0_spline(M[flavor].creators()[i],M[flavor].annihilators()[noperators]);
  }
  //compute the last row
  if(noperators!=0)
    M[flavor].right_multiply(&(Green0_j_n[0]), &(lastcolumn[0]));
  //compute lambda
  double ip = (noperators==0?0:inner_prod(&(Green0_n_j[0]), &(lastcolumn[0]), noperators));
  lambda = Green0_n_n - ip + M[flavor].alpha()[noperators];
  //return weight if we have nothing else to do
  if(compute_only_weight){
    return lambda;
  }
  //compute last column
  if(noperators!=0)
    M[flavor].left_multiply(&(Green0_n_j[0]), &(lastrow[0]));
  //compute norm of vector and mean for roundoff error check
  if(noperators > 0){
    scale(1./lambda, &(lastrow[0]), noperators);
    M[flavor].add_outer_prod(&(lastcolumn[0]), &(lastrow[0]));
    scale(1./lambda, &(lastcolumn[0]), noperators);
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
  std::vector<double> lastrow(noperators-1);
  std::vector<double> lastcolumn(noperators-1);
  for(unsigned int j=0;j<noperators-1;++j){
    lastrow[j]=M[flavor](noperators-1, j);
    lastcolumn[j]=M[flavor](j,noperators-1);
  }
  if(noperators>1)
    scale(-1./Mnn, &(lastcolumn[0]), noperators-1);
  M[flavor].resize(noperators-1);  //lose the last row and last column, reduce size by one, but keep contents.
  if(noperators>1)
    M[flavor].add_outer_prod(&(lastcolumn[0]), &(lastrow[0]));
  return Mnn;  //the determinant ratio det D_{k-1}/det D_{k}
}




