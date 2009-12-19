/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2001-2003 by Fabien Alet <alet@comp-phys.org>,
*                            Matthias Troyer <troyer@comp-phys.org>
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

/* $Id: SSE.Directed.cpp,v 1.8 2005/05/11 21:14:59 troyer Exp $ */

#include "SSE.hpp"
#include <lpkit.h>

// When checking for the solution given by the linear solver is OK, 
// we allow a small deviation of myeps. 10^{-7} is typically sufficient.
double myeps=1e-7;

// For debugging - not documented
bool print=0;

bool SSE::check_out(double weight[4]) 
  // Check out if the solution found fillfuls all requirements
  // Return true if there is a bounce
{
  bool bounce=0;
  for (int i=0;i<4;++i) 
    // Check for non zero bounces
    if (the_sol[i][i]!=0.) { 
      bounce=1; 
      if (print)
      std::cout << "Bounce is possible : " << the_sol[i][i] << "\n";
    }

  if (print)
    for (int j=0;j<4;++j) { 
      double Sum=0.;
      for (int i=0;i<4;++i) { 
        // Check if only one path is allowed - can induce subtle ergodicity probkems  
        if (the_sol[i][j]==1.)
          std::cerr << "*** Matrix element is 1 - check ergodicity\n";

        // Check if detailed balance is respected
        // Here this is not true anymore with the worm weight
        // This check is not valid anymore
        // TODO : change this
        if (the_sol[j][i]*weight[i]-the_sol[i][j]*weight[j]>myeps)
          std::cerr << "Detailed Balance Error " << i << " " << j << "\n";

        // Check is solution is negative
        if (the_sol[i][j]<0.)
          std::cerr << "Negative Error " << i<< " " << j << "\n";
        Sum+=the_sol[i][j];
      }
      // Check if probabilities sum up to 1
      if (weight[j]!=0. && (fabs(1.-Sum)>myeps) )
        std::cerr << "Sum Error " << j << "\n";
    }
  return bounce;
}

bool SSE::find_sol(double weight_before[4],double weight_after[4],double worm_weight_before[4], double worm_weight_after[4])
  
  // Find the solution scattering matrix, given in input the vertex and worm
  // weights for all possible paths - solution is stored in thesol array

  // weight_after[i] is the weight of the vertex {\it after} leg i is flipped
  // weight_before[i] is the weight of the vertex {\it before} leg i was flipped  // worm_weight_after[i] is the weight of the worms {\it after} leg i is flipped
  // worm_weight_before[i] is the weight of the worms {\it before} leg i was flipped
  // ************************ CAUTION ************************

  // Here flip means state (k) -> state (k-1) for legs 0 and 1
  //                 state (k) -> state (k+1) for legs 2 and 3
  //      (Convention 1 : +1 worms oriented north)
  //  *********************************************************
{
  /*
  cout << "WeightBefore : ";
  for (uint8_t rr=0;rr<4;++rr) { cout << weightBefore[rr] << " ";} 
  cout << "\nWeightAfter : ";
  for (uint8_t rr=0;rr<4;++rr) { cout << weightAfter[rr] << " ";} 
  cout << "\nWormWeightBefore : ";
  for (uint8_t rr=0;rr<4;++rr) { cout << WormweightBefore[rr] << " ";} 
  cout << "\nWormWeightAfter : ";
  for (uint8_t rr=0;rr<4;++rr) { cout << WormweightAfter[rr] << " ";} 
  cout <<"\n";
  */

  // An useful definition
  int TM[4][4]; 
  for (int i=0;i<4;++i) TM[i][i]=i+1;
  TM[0][1]=5; TM[0][2]=7; TM[0][3]=9; TM[1][0]=6; TM[1][2]=11; 
  TM[1][3]=13; TM[2][0]=8; TM[2][1]=12; TM[2][3]=15; TM[3][0]=10; 
  TM[3][1]=14; TM[3][2]=16;
  // End of definition

  lprec *lp1;
  
  lp1=make_lp(0,16); // 16 variables, 0 constraints 
  
  // Probabilities constraints (16)
  for (int i=1;i<=16;i++) {
    set_lowbo(lp1,i,0.); 
    set_upbo(lp1,i,1.);
  }

  
  // Normalization of probabilities constraints (4) \sum_j P(i->j) = 1
  // except if of course weight[i]=0
  for (int i=0;i<4;++i) { 
    if (weight_after[i]!=0.) {   
      double constraint[17]={0.};
      for  (int j=0;j<4;++j)
        constraint[TM[j][i]]=1.;
      add_constraint(lp1,constraint, EQ, 1.);
    }
  }


// Detailed Balance constraints (6) P(i->j)/P(j->i)=Wj/Wi
// coming from i

  for (int i=0;i<4;++i)
    for (int j=i+1;j<4;++j) { 
      double constraint[17]={0.};
      constraint[TM[i][j]]=1.;
      constraint[TM[j][i]]=-(worm_weight_before[i]*weight_before[i])/(weight_after[j]*worm_weight_after[j]);
      add_constraint(lp1,constraint, EQ, 0.);
    }
  
// Weight 0 constraints
// If weight[i]=0, then for all j, we have P(j->i) = P(i->j) = 0
  for (int i=0;i<4;++i) 
    if  (weight_after[i]==0.) 
      for (int j=0;j<4;++j) { 
        set_upbo(lp1,TM[i][j],0.);
        set_upbo(lp1,TM[j][i],0.);
      }

  if ((parms.defined("FORCE_SYMMETRY_CONSTRAINT")))
    { 
  // Symmetry constraints  
  /* If 2 incoming legs have the same weights, then for each outgoing leg
     different from these 2, the probabilities should be the same. The reverse
     situation (2 outgoing legs have same weight) should be included in 
     detailed balance. */

   for (int inc1=0;inc1<4;++inc1) 
     for (int inc2=0;inc2<4;++inc2)
        if (inc1!=inc2 && worm_weight_before[inc1]*weight_before[inc1]==worm_weight_before[inc2]*weight_before[inc2])
          for (int out=0;out<4;++out) 
            if (out!=inc1 && out!=inc2) { 
              double constraint[17]={0.};
              constraint[TM[out][inc1]]=1.; 
              constraint[TM[out][inc2]]=-1.;
              add_constraint(lp1,constraint, EQ, 0.);
            }

   // reverse situation (just in case)
   
   for (int out1=0;out1<4;++out1) 
     for (int out2=0;out2<4;++out2)
        if (out1!=out2 && worm_weight_after[out1]*weight_after[out1]==worm_weight_after[out2]*weight_after[out2])
          for (int inc=0;inc<4;++inc)
            if (inc!=out1 && inc!=out2) { 
              double constraint[17]={0.};
              constraint[TM[out1][inc]]=1.; 
              constraint[TM[out2][inc]]=-1.;
              add_constraint(lp1,constraint, EQ, 0.);
            }
    }
// Function to minimize (here Trace)
   double fm[17]={0.}; 
   fm[1]=fm[2]=fm[3]=fm[4]=1.;
   set_obj_fn(lp1,fm);  

// Now Solve
  int solu=solve(lp1);
  if (solu!=0) 
    std::cerr << "*** CAUTION !!! No solution found !!! " << solu << "\n";

  // Get Solution
  //  double the_sol[5][5];
  for (int i=0;i<4;++i)
    for (int j=0;j<4;++j) 
      the_sol[i][j]=lp1->best_solution[lp1->rows+TM[i][j]]; 


  // Print the solution
 /*   
     for (uint32_t r1=0;r1<4;++r1) {
        for (uint32_t r2=0;r2<4;++r2) 
          cout << the_sol[r1][r2] << " ";
        cout << endl;
      }
 */

  // Check it out
  bool bb=check_out(weight_after);
  
  // 
  delete_lp(lp1); 
  return bb;
}


double SSE::return_weight_after_flip(uint32_t bond_type, uint8_t VS[4],uint8_t leg_to_flip,bool convention)
// Return the matrix element of the vertex if the leg "leg_to_flip"
// is modified
// VS means Vertex State : the status of the legs of the current vertex
  {
  double ans=0.; convention=1;
  uint8_t state_number_0=number_states_for_site_type_[site_type_for_bond_type_[bond_type].first];
  uint8_t state_number_1=number_states_for_site_type_[site_type_for_bond_type_[bond_type].second];

  if (convention) 
    switch(leg_to_flip) {
    case 0: 
     if (VS[0]!=0) 
        ans=matrix_element[bond_type][VS[0]-1][VS[1]][VS[2]][VS[3]];
      break;
    case 1: 
      if (VS[1]!=0)
        ans=matrix_element[bond_type][VS[0]][VS[1]-1][VS[2]][VS[3]];
      break;
    case 2: 
      if (VS[2]!=state_number_0-1)
        ans=matrix_element[bond_type][VS[0]][VS[1]][VS[2]+1][VS[3]];
      break;
    case 3:
      if (VS[3]!=state_number_1-1)
        ans=matrix_element[bond_type][VS[0]][VS[1]][VS[2]][VS[3]+1];
      break;
    default:
      boost::throw_exception(std::logic_error("default encountered in SSE::return_weight_after_flip"));
  }

  /**************** The following is for other convention *************/
  /**************** Not used in the current implementation ************/

  else 
    switch(leg_to_flip) {
    case 0: 
      if (VS[0]!=state_number_0-1)
        ans=matrix_element[bond_type][VS[0]+1][VS[1]][VS[2]][VS[3]];
      break;
    case 1: 
      if (VS[1]!=state_number_1-1) 
        ans=matrix_element[bond_type][VS[0]][VS[1]+1][VS[2]][VS[3]];
      break;
    case 2: 
      if (VS[2]!=0)
        ans=matrix_element[bond_type][VS[0]][VS[1]][VS[2]-1][VS[3]];
      break;
    case 3:
      if (VS[3]!=0) 
        ans=matrix_element[bond_type][VS[0]][VS[1]][VS[2]][VS[3]-1];
      break;
    default:
      boost::throw_exception(std::logic_error("default encountered in SSE::return_weight_after_flip"));
    }

  /********************************************************************/
  /********************************************************************/
  return ans;
}


double SSE::return_weight_before_flip(uint32_t bond_type,uint8_t VS[4],uint8_t leg_that_was_flipped,bool convention)
// Return the matrix element of the vertex before the leg "leg_that_was_flipped" was modified
// VS means Vertex State : the status of the legs of the current vertex
{
double ans=0.; 
convention=1;
uint8_t state_number_0=number_states_for_site_type_[site_type_for_bond_type_[bond_type].first];
uint8_t state_number_1=number_states_for_site_type_[site_type_for_bond_type_[bond_type].second];

if (convention)
  switch(leg_that_was_flipped) {
  case 0: 
    if (VS[0]!=0)
      ans=matrix_element[bond_type][VS[0]-1][VS[1]][VS[2]][VS[3]];
    break;
  case 1:
    if (VS[1]!=0) 
      ans=matrix_element[bond_type][VS[0]][VS[1]-1][VS[2]][VS[3]];
    break;
  case 2: 
    if (VS[2]!=state_number_0-1) 
      ans=matrix_element[bond_type][VS[0]][VS[1]][VS[2]+1][VS[3]];
    break; 
  case 3: 
    if (VS[3]!=state_number_1-1) 
      ans=matrix_element[bond_type][VS[0]][VS[1]][VS[2]][VS[3]+1];
    break; 
  default:
    boost::throw_exception(std::logic_error("default encountered in SSE::return_weight_before_flip"));
  }

  /**************** The following is for other convention *************/
  /**************** Not used in the current implementation ************/

else
  switch(leg_that_was_flipped) {
  case 0: 
    if (VS[0]!=state_number_0-1)
      ans=matrix_element[bond_type][VS[0]+1][VS[1]][VS[2]][VS[3]];
    break;
  case 1: 
    if (VS[1]!=state_number_1-1)
      ans=matrix_element[bond_type][VS[0]][VS[1]+1][VS[2]][VS[3]];
    break;
  case 2:
    if (VS[2]!=0)
      ans=matrix_element[bond_type][VS[0]][VS[1]][VS[2]-1][VS[3]];
    break;
  case 3:
    if (VS[3]!=0) 
      ans=matrix_element[bond_type][VS[0]][VS[1]][VS[2]][VS[3]-1];
    break; 
  default:
    boost::throw_exception(std::logic_error("default encountered in SSE::return_weight_before_flip"));
  }
  /********************************************************************/
  /********************************************************************/
   
return ans;
}

double SSE::return_worm_weight_before(uint32_t bond_type,uint8_t VS[4],uint8_t leg_that_was_flipped, bool convention)
// Return the matrix element of the worm before the leg "leg_that_was_flipped" was modified
// VS means Vertex State : the status of the legs of the current vertex
{
  double ans=0.; convention=1;

if (convention) 
  switch(leg_that_was_flipped) {
  case 0: 
      ans=worm_plus(VS[0]-1,site_type_for_bond_type_[bond_type].first);
    break; 
  case 1: 
      ans=worm_plus(VS[1]-1,site_type_for_bond_type_[bond_type].second);
    break;
  case 2: 
      ans=worm_minus(VS[2]+1,site_type_for_bond_type_[bond_type].first);
    break; 
  case 3:
      ans=worm_minus(VS[3]+1,site_type_for_bond_type_[bond_type].second);
    break;
  default:
    boost::throw_exception(std::logic_error("default encountered in SSE::return_weight_before_flip"));
  }
  /**************** The following is for other convention *************/
  /**************** Not used in the current implementation ************/
else
  switch(leg_that_was_flipped) {
  case 0: 
      ans=worm_minus(VS[0]+1,site_type_for_bond_type_[bond_type].first);
    break;
  case 1: 
      ans=worm_minus(VS[1]+1,site_type_for_bond_type_[bond_type].second);
    break;
  case 2: 
      ans=worm_plus(VS[2]-1,site_type_for_bond_type_[bond_type].first);
    break;
  case 3: 
      ans=worm_plus(VS[3]-1,site_type_for_bond_type_[bond_type].second);
    break;
  default:
    boost::throw_exception(std::logic_error("default encountered in SSE::return_weight_before_flip"));
  }
  /********************************************************************/
  /********************************************************************/

// If we force the standard directed loop equations, we force the worm factor
// to be constant
if (ans!=0. && FORCE_STANDARD_DIRECTED_LOOPS) 
  ans=1.;

return ans;    
}

double SSE::return_worm_weight_after(uint32_t bond_type,uint8_t VS[4],uint8_t leg_to_flip, bool convention)
// Return the matrix element of the worm if the leg "leg_to_flip"
// is modified
// VS means Vertex State : the status of the legs of the current vertex

{
double ans=0.; 
convention=1;

if (convention) 
  switch(leg_to_flip) {
  case 0:
    ans=worm_minus(VS[0],site_type_for_bond_type_[bond_type].first);
    break; 
  case 1:
    ans=worm_minus(VS[1],site_type_for_bond_type_[bond_type].second);
    break; 
  case 2: 
    ans=worm_plus(VS[2],site_type_for_bond_type_[bond_type].first);
    break;
  case 3: 
    ans=worm_plus(VS[3],site_type_for_bond_type_[bond_type].second);
    break; 
  default:
    boost::throw_exception(std::logic_error("default encountered in SSE::return_weight_after_flip"));
  }
else
  switch(leg_to_flip) {
    case 0: 
      ans=worm_plus(VS[0],site_type_for_bond_type_[bond_type].first);
      break;
    case 1:
      ans=worm_plus(VS[1],site_type_for_bond_type_[bond_type].second);
      break; 
    case 2:
      ans=worm_minus(VS[2],site_type_for_bond_type_[bond_type].first);
      break;
    case 3: 
      ans=worm_minus(VS[3],site_type_for_bond_type_[bond_type].second);
      break;
    default:
      boost::throw_exception(std::logic_error("default encountered in SSE::return_weight_after_flip"));
  }

// If we force the standard directed loop equations, we force the worm factor
// to be constant
if (ans!=0. && FORCE_STANDARD_DIRECTED_LOOPS) 
  ans=1.;
  
return ans;
}

bool SSE::calculate_scattering_matrix(uint32_t bond_type,uint8_t VS[4],bool convention)
  // Given a current vertex (VS is the Vertex state), calculates the different worm and vertex weights, feed them in the solver, and return true if there is a bounce 
{
  double weigA[4]; // matrix weight after
  double weigB[4]; // matrix weight before
  double WWA[4]; // worm weight after
  double WWB[4]; // worm weight before
  for (uint8_t i=0;i<4;++i) { 
    weigA[i]=return_weight_after_flip(bond_type,VS,i,convention);
    weigB[i]=return_weight_before_flip(bond_type,VS,i,convention);
    WWA[i]=return_worm_weight_after(bond_type,VS,i,convention); 
    WWB[i]=return_worm_weight_before(bond_type,VS,i,convention);
  }
  

  // Check if at least one path gives rise to non-zero vertex
  bool possible=false; 
  for (unsigned short int i=0;i<4 && !possible;++i) 
    if (weigA[i]!=0.) 
      possible=true;
      
  return (possible ? find_sol(weigB,weigA,WWB,WWA) : false);
}


void SSE::calculate_heat_bath_matrix(uint32_t bond_type,uint8_t VS[4],bool convention)
  // Calculate the scattering matrix with the conventional heat bath choice
  // Scattering matrix is put in the_sol array
{
  double weigA[4]; // matrix weight after
  for (uint8_t i=0;i<4;++i)
    weigA[i]=return_weight_after_flip(bond_type,VS,i,convention);

  // Check if at least one path gives rise to non-zero vertex
  bool possible=false; 
  for (unsigned short int i=0;i<4 && !possible;++i) 
    if (weigA[i]!=0.) 
      possible=true;

  double Sum=0.;
  for (int i=0;i<4;++i)
    Sum+=weigA[i];
  for (int j=0;j<4;++j)
    for (int i=0;i<4;++i)
      the_sol[i][j]=possible ? weigA[i]/Sum : 0.;
}


bool SSE::calculate_all_scattering_matrices(uint32_t bond_type,bool convention)
  // Loops over all the vertices, calculates their scattering matrix 
  // Checks if it has not been calculated before by symmetry 
  // Currently only the symmetry right-left is implemented
  // Time reversal symetries is not needed if one pick one single convention
{
  bool totb=0; 
  uint8_t state_number_0=number_states_for_site_type_[site_type_for_bond_type_[bond_type].first];
  uint8_t state_number_1=number_states_for_site_type_[site_type_for_bond_type_[bond_type].second];

  // Big array collecting the vertex that have already been done
  boost::multi_array<bool,5> already_done(boost::extents[2][state_number_0][state_number_1][state_number_0][state_number_1]); 
  for (unsigned int l1=0;l1<state_number_0;++l1) 
    for (unsigned int l2=0;l2<state_number_1;++l2)
      for (unsigned int l3=0;l3<state_number_0;++l3)
        for (unsigned int l4=0;l4<state_number_1;++l4) { 
          already_done[0][l1][l2][l3][l4]=false;
          already_done[1][l1][l2][l3][l4]=false;
        }

  for (unsigned int l1=0;l1<state_number_0;++l1) 
    for (unsigned int l2=0;l2<state_number_1;++l2)
      for (unsigned int l3=0;l3<state_number_0;++l3)
        for (unsigned int l4=0;l4<state_number_1;++l4) { 
          uint8_t VS[4]={l1,l2,l3,l4};
          for (int i=0;i<4;++i) 
            for (int j=0;j<4;++j) 
              the_sol[i][j]=0.;

          /************ Symetries *************/
          // symmetry right-left (symmetry 1)
          // if already done, looks for the corresponding solution
          if (state_number_0==state_number_1) {
            // do this only if right and left legs are identical
            // (here same number of states on each leg)
            // NOTE : this could be different for complicated models
            // where the number of states are the same but the legs
            // somehow different
            if (already_done[convention][l2][l1][l4][l3])
            put_in(bond_type,VS,1);
          }
          /*****************************************************************
           ********* Do these symetries for different conventions **********
           *****************************************************************
       
           // Time-reversal symmetry (symmetry 2)
           // if already done, looks for the corresponding solution
           if (AlreadyDone[!convention][l3][l4][l1][l2])
             { }
           // right-left + Time reversal (symmetry 3)
           // if already done, looks for the corresponding solution
            if (AlreadyDone[!convention][l4][l3][l2][l1])
             { }
      
           ****************************************************************
           ****************************************************************
           ****************************************************************/

    // If not already done :
    // Checks if we want heat bath or directed loop solutions
    if (!(FORCE_HEATBATH)) { 
      // Calculates the directed loop scattering matrix, check for bounce
      if (calculate_scattering_matrix(bond_type,VS,convention))
         totb=1;
    }
    else { 
      // Calculates the heat bath scattering matrix, bounce is always present
      totb=1;
      calculate_heat_bath_matrix(bond_type,VS,convention);
    }

    // Put the solution in the correct place
    put_in(bond_type,VS,0);
    // Mark this vertex as already done
    already_done[convention][l1][l2][l3][l4]=1;  
      }
// Returns true if at least one bounce appeared 
 return totb;
}


void SSE::put_in(uint32_t bond_type,uint8_t VS[4],uint8_t symmetry)
  // Puts the solution at the correct place
  // If the calculation has already been done by symmetry, looks for the corresponding solution, unrolls it and put it at the correct place
{
  if (symmetry) { 
    // already done, can use symmetries to obtain elements
    
    // array that gives the symetric legs 
    // according to the corresponding symmetry (1,2,3)

    int Sym[3][4]={{1,0,3,2},{2,3,0,1},{3,2,1,0}};

    // The symetric Vertex State we are looking for
    uint8_t SymVS[4]; 
    for (uint8_t r=0;r<4;++r) 
      SymVS[r]=VS[Sym[symmetry-1][r]];

    
    for (int inc=0;inc<4;++inc)  { 
      int ent=Sym[symmetry-1][inc];
      // Unrolls the cumulated probabilities of the symetric vertex
      double befa=proba_worm[bond_type][SymVS[0]][SymVS[1]][SymVS[2]][SymVS[3]][0][ent];
      double befb=proba_worm[bond_type][SymVS[0]][SymVS[1]][SymVS[2]][SymVS[3]][1][ent]-befa;
      double befc=proba_worm[bond_type][SymVS[0]][SymVS[1]][SymVS[2]][SymVS[3]][2][ent]-befa-befb;
      double befd=1.-befa-befb-befc;
       
      if (symmetry==1) { // right-left
        proba_worm[bond_type][VS[0]][VS[1]][VS[2]][VS[3]][0][inc]=befb;
        proba_worm[bond_type][VS[0]][VS[1]][VS[2]][VS[3]][1][inc]=befa+befb;
        proba_worm[bond_type][VS[0]][VS[1]][VS[2]][VS[3]][2][inc]=befa+befb+befd;
        proba_worm[bond_type][VS[0]][VS[1]][VS[2]][VS[3]][3][inc]=1.;
      }
/*****************************************************************************/
/**** These two symetries are not used in the current convention choice ******/
/*****************************************************************************/
      if (symmetry==2) { // time
        proba_worm[bond_type][VS[0]][VS[1]][VS[2]][VS[3]][0][inc]=befc;
        proba_worm[bond_type][VS[0]][VS[1]][VS[2]][VS[3]][1][inc]=befc+befd;
        proba_worm[bond_type][VS[0]][VS[1]][VS[2]][VS[3]][2][inc]=befa+befc+befd;
        proba_worm[bond_type][VS[0]][VS[1]][VS[2]][VS[3]][3][inc]=1.;
      }
      if (symmetry==3) { // time + right-left
        proba_worm[bond_type][VS[0]][VS[1]][VS[2]][VS[3]][0][inc]=befd;
        proba_worm[bond_type][VS[0]][VS[1]][VS[2]][VS[3]][1][inc]=befc+befd;
        proba_worm[bond_type][VS[0]][VS[1]][VS[2]][VS[3]][2][inc]=befb+befc+befd;
        proba_worm[bond_type][VS[0]][VS[1]][VS[2]][VS[3]][3][inc]=1.;
      }
/*****************************************************************************/
/*****************************************************************************/
    }
  } // end of already done
  else  // not already done
    // Put the cumulated probabilities at the correct place
    for (int inc=0;inc<4;++inc) {
      proba_worm[bond_type][VS[0]][VS[1]][VS[2]][VS[3]][0][inc]=the_sol[0][inc];
      proba_worm[bond_type][VS[0]][VS[1]][VS[2]][VS[3]][1][inc]=the_sol[0][inc]+the_sol[1][inc];
      proba_worm[bond_type][VS[0]][VS[1]][VS[2]][VS[3]][2][inc]=the_sol[0][inc]+the_sol[1][inc]+the_sol[2][inc];
      proba_worm[bond_type][VS[0]][VS[1]][VS[2]][VS[3]][3][inc]=1.;
    }
}


void SSE::calculate_proba_worm()
// Loops over all bond types and calculate all their scattering matrices
{
  // resize the worm probability array
   proba_worm.resize( boost::extents[number_of_bond_types][maximum_number_of_states][maximum_number_of_states][maximum_number_of_states][maximum_number_of_states][4][4]);

   bool convention=true;
   for (uint32_t i=0;i<number_of_bond_types;++i)
     calculate_all_scattering_matrices(i,convention); 
}
