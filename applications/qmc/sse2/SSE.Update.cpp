/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2001-2005 by Fabien Alet <alet@comp-phys.org>,
*                            Matthias Troyer <troyer@comp-phys.org>
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

/* $Id: SSE.Update.cpp,v 1.25 2005/07/20 04:41:30 troyer Exp $ */

#include "SSE.hpp"

/*************************************************************************
 **************************** Diagonal Update ****************************
 ************************************************************************/


// Do a diagonal update
void SSE::diagonal_update()
{
        vector<uint8_t> states=site_state;
        std::vector<unsigned int> last_level;
        
        std::valarray<double> initm(num_sites());
        if (measure_site_compressibility_) {
                last_level.resize(num_sites(),std::numeric_limits<unsigned int>::max());
                for (unsigned int  i=0;i<num_sites();++i) {
                        initm[i] = matrix_element_n_[site_type(i)][states[i]];
                        localint[i]=operator_string.size()*initm[i];
                        localint2[i]=operator_string.size()*initm[i]*initm[i];
                }
        }
        double proba=0.;


        if(variable_length_representation) {
                vector<vertex_type>::iterator oit=operator_string.begin();
                vector<vertex_type>::iterator bit=buffer_string.begin();
                bool deletion_occured=false;               
                 
                for (int i=0; i<=current_number_of_non_identity*2;++i) {
                        proba=0.;
                        
                        if (i%2==0)  {//"Gap" between succeeding operators    
              if (current_number_of_non_identity<norder_max) {  
                uint32_t random_bond=random_int(0,num_bonds()-1);
                bond_descriptor this_bond=bond(random_bond);
                proba=ensemble_weight_fraction(current_number_of_non_identity);
                proba*=proba_diagonal[bond_type[this_bond]][states[source(this_bond)]][states[target(this_bond)]]
                  /(current_number_of_non_identity+1);
                if (proba!=0. && (proba>=1. || random_real()<proba)) {
                    bit->vertex_type=1;
                bit->bond_number=random_bond;
                bit->leg[0]=bit->leg[2]=states[source(this_bond)];                              
                bit->leg[1]=bit->leg[3]=states[target(this_bond)];
                ++bit;
                ++current_number_of_non_identity;
                ++i;                    
                    } 
              }
                        } else { //Operator
              deletion_occured=false;
              if (oit->diagonal()) {
                if (current_number_of_non_identity>norder_min) {
                  proba=1./ensemble_weight_fraction(current_number_of_non_identity-1);
                  proba*=current_number_of_non_identity
                /proba_diagonal[bond_type[bond(oit->bond_number)]][oit->leg[0]][oit->leg[1]];
                  
                  if ( proba>=1 || random_real()<proba ) {
                --current_number_of_non_identity;
                      --i;
                ++oit;  //Proceed to next operator in original string
                    deletion_occured=true;
                  }     
                }
              } else { // Non-Diagonal

                int s=source(bond(oit->bond_number));
                int t=target(bond(oit->bond_number));
                
                /*    if (measure_site_compressibility_) {
                  if (last_level[s]!=std::numeric_limits<unsigned int>::max()) {
                  localint[s]+=(level-last_level[s])*(matrix_element_n_[site_type(s)][states[s]]-initm[s]);
                  localint2[s]+=(level-last_level[s])*(matrix_element_n_[site_type(s)][states[s]]*matrix_element_n_[site_type(s)][states[s]]-initm[s]*initm[s]);
                  }
                  if (last_level[t]!=std::numeric_limits<unsigned int>::max()) {
                  localint[t]+=(level-last_level[t])*(matrix_element_n_[site_type(t)][states[t]]-initm[t]);
                  localint2[t]+=(level-last_level[t])*(matrix_element_n_[site_type(t)][states[t]]*matrix_element_n_[site_type(t)][states[t]]-initm[t]*initm[t]);
                  }
                  last_level[s]=last_level[t]=level;
                  }*/
                states[s]=oit->leg[2];
                states[t]=oit->leg[3];
              }
              if (!deletion_occured) {   //If nothing was deleted (offdiagonal or by probability, we copy the op.
                *bit=*oit;
                ++oit;
                ++bit;
              }
            }
                        diagonal_update_iteration_end(i);   //Used to update histograms, weights (WL) etc.
                }
                operator_string=buffer_string;
                buffer_string=identity_string;             
        diagonal_update_sweep_end();
                
        } else {
                
                // Loops over all vertices
                int level=0;
                
                for ( vector<vertex_type>::iterator it=operator_string.begin();it!=operator_string.end();++it,++level) {
                        if (!(it->non_diagonal())) { // Identity or Diagonal
                                if (it->identity()) { // Identity operator
                                        if (current_number_of_non_identity<norder_max) {
                                                
                                                // Get a random bond
                                                uint32_t random_bond=random_int(0,num_bonds()-1);
                                                bond_descriptor this_bond=bond(random_bond);
                                                
                                                // Get the probability part from the use of WL/optimized ensembles
                                                proba=ensemble_weight_fraction(current_number_of_non_identity);
                                                
                                                
                                                // Now, multiply with the usual matrix element.
                                                proba*=proba_diagonal[bond_type[this_bond]][states[source(this_bond)]][states[target(this_bond)]]/(cutoff_L-current_number_of_non_identity);
                                                
                                                if (proba!=0. && (proba>=1. || random_real()<proba)) {
                                                        // Insert a diagonal element
                                                        ++current_number_of_non_identity;
                                                        it->vertex_type=1;
                                                        it->bond_number=random_bond;                                                        
                                                        it->leg[0]=it->leg[2]=states[source(this_bond)]; 
                                                        it->leg[1]=it->leg[3]=states[target(this_bond)]; 
                                                }
                                        }
                                }
                                else { // Diagonal
                                        if (current_number_of_non_identity>norder_min) {
                                                // Get the removal probability part from the use of WL/OE
                                                proba=1./ensemble_weight_fraction(current_number_of_non_identity-1);
                                                
                                                // Now, multiply with the usual matrix element
                                                proba*=(cutoff_L-current_number_of_non_identity+1) / proba_diagonal[bond_type[bond(it->bond_number)]][it->leg[0]][it->leg[1]];
                                                if ( proba>=1 || random_real()<proba ) {
                                                        // Remove it
                                                        --current_number_of_non_identity;
                                                        it->vertex_type=0;
                                                        it->bond_number=0;
                                                }
                                        }
                                }
                        }
                        else { // Non-Diagonal
                                int s=source(bond(it->bond_number));
                                int t=target(bond(it->bond_number));
                                if (measure_site_compressibility_) {
                                        if (last_level[s]!=std::numeric_limits<unsigned int>::max()) {
                                                localint[s]+=(level-last_level[s])*(matrix_element_n_[site_type(s)][states[s]]-initm[s]);
                                                localint2[s]+=(level-last_level[s])*(matrix_element_n_[site_type(s)][states[s]]*matrix_element_n_[site_type(s)][states[s]]-initm[s]*initm[s]);
                                        }
                                        if (last_level[t]!=std::numeric_limits<unsigned int>::max()) {
                                                localint[t]+=(level-last_level[t])*(matrix_element_n_[site_type(t)][states[t]]-initm[t]);
                                                localint2[t]+=(level-last_level[t])*(matrix_element_n_[site_type(t)][states[t]]*matrix_element_n_[site_type(t)][states[t]]-initm[t]*initm[t]);
                                        }
                                        last_level[s]=last_level[t]=level;
                                }
                                states[s]=it->leg[2];
                                states[t]=it->leg[3];
                        }    
                        diagonal_update_iteration_end(0);
                }
        diagonal_update_sweep_end();
    }
}

/****************************************************************************
 *********************************  Link legs ******************************
****************************************************************************/

// Just after the diagonal update, link (non-Identity) vertices to their 
// following neighbour at each site. Necessary for worm update.
 
void SSE::link_legs()
{

//std::cout<<"in Link_legs\n";

  // for each site, the current vertex sitting there
  // A value of cutoff_L indicates no vertex at this site (for the moment)
  vector<uint32_t> current_vertex(num_sites(),cutoff_L);
  // for each site, the first vertex sitting there
  // A value of cutoff_L indicates no vertex at this site (for the moment)
  vector<uint32_t> first_vertex(num_sites(),cutoff_L);
//  cerr << "in link_legs\n";
  std::vector<vertex_type>::iterator it; 
  uint32_t i=0;
  // Loop over all vertices
  for (it=operator_string.begin();it!=operator_string.end();++it,++i) {
    if (!(it->identity())) { // Non-Identity vertices
          //          cerr << "bond "<< it->bond_number<<": not id\n";
      for (uint8_t rr=0;rr<2;++rr) { // For both legs 
        // Get the site on which this vertex sits
              
        sites_size_type site = rr ? target(bond(it->bond_number)) : source (bond(it->bond_number));
        //cerr << "rr " << (int)rr <<", site "<<site<<"\n";
        
        if (current_vertex[site]==cutoff_L)
          // First vertex on this site
        { first_vertex[site]=i;
          //    cerr << "first vertex on site\n";
        
        }
        else {// There was already a vertex on this site
                // Get the previous vertex on this site
            //          cerr << "already a vertex on site \n";
                
          uint32_t previous_vertex=current_vertex[site];         
          // This boolean indicates if it was previous vertex' left (0) or 
          // right (leg) at this site
          bool previous_vertex_first_leg=(source(bond(operator_string[previous_vertex].bond_number))!=site);
          // Link the previous vertex' top legs to the current one
          operator_string[previous_vertex].linked_vertices[2+previous_vertex_first_leg]=make_pair(i,bool(rr)); 
          // Link the current vertex' bottom leg to the previous one
          it->linked_vertices[rr]=make_pair(previous_vertex,previous_vertex_first_leg);
          
        }
        // Put this vertex as the current vertex on this site
        current_vertex[site]=i;
      }
    }
  }
  //cerr << "after first try\n";
  //  cerr << "cutoffL:"<<cutoff_L<<"\n";

  // Boundary conditions in imaginary time
  // Link last and first vertices on all sites
  for (uint32_t ii=0;ii<num_sites();++ii) { 
    if (current_vertex[ii]!=cutoff_L) { 
      uint32_t last = current_vertex[ii];
      uint32_t first = first_vertex[ii];
      bool last_legs  = (source(bond(operator_string[last].bond_number))!=ii);
      bool first_legs = (source(bond(operator_string[first].bond_number))!=ii);
      operator_string[last].linked_vertices[2+last_legs]=make_pair(first,first_legs);
      operator_string[first].linked_vertices[first_legs]=make_pair(last,last_legs);
    }
    else // There was no vertex on this site - Pick a random state

      site_state[ii]=random_int(0,site_number_of_states[ii]-1);
  }
//cerr << "after link legs loop\n";

  //print_operator_string();
  
}

// Helper function to print the Operator string
void SSE::print_operator_string()
{
  std::vector<vertex_type>::iterator it; uint32_t i=0;
  for (it=operator_string.begin();it!=operator_string.end();++it,++i) { 
    if (!(it->identity())) {
      cout << "Vertex " << i << " between [" << source(bond(it->bond_number)) 
           << ", " << target(bond(it->bond_number)) << "] : ";
      for (uint8_t rr=0;rr<4;++rr)
           cout << "Links " << rr << "={ "<< it->linked_vertices[rr].first << "," <<  it->linked_vertices[rr].second << "} ";
      cout << "\nState legs ";  
      for (uint8_t rr=0;rr<4;++rr)
        cout << (int) it->leg[rr] << " ";
      cout << "\n";
    } 
  }
}

/*************************************************************************
 **************************** Worm Start ********************************
 ************************************************************************/

std::pair<uint32_t,uint8_t> SSE::initial_vertex() 
{
        //  cerr << "start initial vertex\n";
  uint32_t current_vertex;
  uint8_t current_leg_number;
  if (!measure_green_function_) {
    current_vertex=random_int(0,cutoff_L-1); 
    //cerr << "ugly loop 1\n";
    // go to the first non trivial vertex
    while(operator_string[current_vertex].identity()) { 
      if ((++current_vertex)==cutoff_L)
        current_vertex=0;
    }
    //cerr << "current_vertex " << current_vertex << ", bond_number: "<<operator_string[current_vertex].bond_number  <<"\n";
    
    current_leg_number=random_int(0,3);
  }
  else {
    initial_site=random_int(0,num_sites()-1);
    current_vertex=initial_level=random_int(0,cutoff_L-1);
    bool searchup = random()<0.5;
  
    if (searchup) {
      while(operator_string[current_vertex].identity() || 
          (source(bond(operator_string[current_vertex].bond_number))!=initial_site &&
          target(bond(operator_string[current_vertex].bond_number))!=initial_site)) { 
        if (++current_vertex==cutoff_L)
          current_vertex=0;
        if (current_vertex==initial_level)
          return std::make_pair(current_vertex,4); // no vertex found
      }
    }
    else {
      while(operator_string[current_vertex].identity() || 
          (source(bond(operator_string[current_vertex].bond_number))!=initial_site &&
          target(bond(operator_string[current_vertex].bond_number))!=initial_site)) { 
        if (current_vertex==0)
          current_vertex=cutoff_L-1;
        else
          --current_vertex;
        if (current_vertex==initial_level)
          return std::make_pair(current_vertex,5); // no vertex found
      }
      ++initial_level;
      if (initial_level==cutoff_L)
        initial_level=0;
    }
    if (source(bond(operator_string[current_vertex].bond_number))==initial_site)
      current_leg_number=0 + (searchup ? 0 : 2);
    else if (target(bond(operator_string[current_vertex].bond_number))==initial_site)
      current_leg_number=1 + (searchup ? 0 : 2);
    else
      boost::throw_exception(std::logic_error("Invalid starting vertex"));
//    std::cerr << "Initial level: " << initial_level << "\n";
//    std::cerr << "Initial site: " << initial_site << "\n";
//    std::cerr << "Initial vertex: " << current_vertex << "\n";
//    std::cerr << "Initial leg: " <<int( current_leg_number) << "\n";
  }
  return std::make_pair(current_vertex,current_leg_number);
}

/*************************************************************************
 **************************** Worm Update *******************************
 ************************************************************************/

// Worm update
void SSE::worm_update()
{ 
  // pick a random slot and leg
  uint32_t current_vertex; 
  uint8_t current_leg_number;
  
//cerr << "doing a worm update\n";

  boost::tie(current_vertex,current_leg_number) = initial_vertex();
  
  if (current_leg_number>3) {
    // no vertex: just measure 
    if (measure_green_function_) {
      double matel = 0.5*(
      matrix_element_raise_[site_type(initial_site)][site_state[initial_site]]*
      matrix_element_raise_[site_type(initial_site)][site_state[initial_site]]+
      matrix_element_lower_[site_type(initial_site)][site_state[initial_site]]*
      matrix_element_lower_[site_type(initial_site)][site_state[initial_site]]);
      if (!measurement_origin_)
        green[distance(initial_site,initial_site)] += matel;
      else if (initial_site==measurement_origin_.get())
        green[measurement_origin_.get()] += matel;
    }
    return;
  }
  
  // pick a convention
  // Convention 1 = +1 worm go up, -1 worm down
  //            0 = +1 worm go down, -1 worm up
  // bool Convention=(random_real()<0.5);
  // For the moment, only convention 1 is implemented
  // bool convention=true;

  // boolean current_leg : 0 for left legs, 1 for right legs
  bool current_leg=((current_leg_number==1) || (current_leg_number==3));
  
  // boolean is_upper_leg : 0 for bottom legs, 1 for upper legs
  bool is_upper_leg=(current_leg_number>1);
  
  //************ The following might be useful for other conventions *****
  //************ Please keep it like this (excuse my french) *************
  //******* NON !! is_upper_leg dit vraiment si on est au nord ou non ****
  // Convention 0  
  // bool is_upper_leg=(current_leg_number<2);
  // Generic *** Check This ****
  // bool is_upper_leg=(convention==(current_leg_number>1));
  //**********************************************************************
 
  // Initial definitions - useful to check for the end of the worm
  uint32_t initial_vertex=current_vertex;
  bool initial_leg=current_leg;
  uint32_t next_vertex=current_vertex;
  bool next_leg=current_leg;
  bool initial_is_upper_leg=is_upper_leg;
    
    /************** For different conventions only ****************
     *************  Please keep it like this **********************
     **************************************************************

    bool plus_worm=(random_real()<0.5);

    if (plus_worm)
      { if (convention==1)
        { if (is_upper_leg==1)
          { go to upper-vertex;}
        }
        else 
          { if (is_upper_leg==0)
            { go to lower vertex;}
          }
      }
    else {
      if (convention==1)
        { if (is_upper_leg==0)
          { go to lower vertex;}
        }
      else 
        { if (is_upper_leg==1)
          { go to upper vertex;}
          }
      }
    
    ***************************************************************
    ***************************************************************/

  // Will the worm be accepted ?
  bool do_it=true;

    // Don't propagate (+/-) 1 worms if they are to act on the max(min)imum
    //  number of states at this site
    // for different conventions, change this to if (!(Plus_Worm))

    // TODO : For 2 states models : always construct a worm


  if (is_upper_leg) { // Propagation of a -1 worm
    do_it=(operator_string[initial_vertex].leg[current_leg_number]!=0);
    if (measure_green_function_)
      worm_weight = matrix_element_lower_[site_type(initial_site)][operator_string[initial_vertex].leg[current_leg_number]]; 
  }
  else { // Propagation of a +1 worm 
    bond_descriptor b = bond(operator_string[initial_vertex].bond_number);
    do_it=(operator_string[initial_vertex].leg[current_leg_number]!=(site_number_of_states[current_leg ? target(b) : source (b)]-1));
    if (measure_green_function_)
      worm_weight = matrix_element_raise_[site_type(initial_site)][operator_string[initial_vertex].leg[current_leg_number]]; 
  }
  if (measure_green_function_)
    worm_weight *= worm_weight;

    // for model: is matrixelement of op applied to state nonzero? TOCHANGE

  if (!do_it)
    return;

  if (measure_green_function_) {
    if (!measurement_origin_)
      green[distance(initial_site,initial_site)] += worm_weight;
    else if (initial_site==measurement_origin_.get())
      green[measurement_origin_.get()] += worm_weight;
  }    
  //cerr << "vor worm main loop\n";

  // main loop
  while (true) {
    // Add 1 to worm size
    worm_size+=1;
    // Modify entrance leg
    // This has to be modified for different convention
    //std::cerr << "Modifying entrance leg: " << next_vertex << " " << int(current_leg_number) << "\n";
    if (is_upper_leg) 
      --operator_string[next_vertex].leg[current_leg_number];
    else
      ++operator_string[next_vertex].leg[current_leg_number];
        
    // Modify vertex
    // the modify_vertex function modify the exit leg (given by 
    // return_exit_leg and the vertex type
    operator_string[next_vertex].modify_vertex(return_exit_leg(next_vertex,current_leg_number),next_leg,is_upper_leg);
    //std::cerr << "Passed through: " << next_vertex << " " << next_leg << "\n";

    // Eventual End
    if (next_vertex==initial_vertex && next_leg==initial_leg && is_upper_leg==initial_is_upper_leg)
      break;
    
    
    // Otherwise get the vertex to come
    current_vertex=next_vertex; 
    current_leg=next_leg;

    // Get the next vertex and the new values for next_leg and
    // current_leg_number
    boost::tie(next_vertex,next_leg)=operator_string[current_vertex].linked_vertices[current_leg+is_upper_leg*2];
    current_leg_number=next_leg+(1-is_upper_leg)*2;
    //std::cerr << "Following to: " << next_vertex << " " << next_leg << "\n";
          
    // Change states of sites if the worm go accross imaginary time
    // boundary conditions **** Improve this ??? ****
            
    int s = current_leg ? target(bond(operator_string[current_vertex].bond_number))
                          : source(bond(operator_string[current_vertex].bond_number));
    
        
    bool crossed=false;
    if (is_upper_leg) {
      if (next_vertex<=current_vertex) {
        site_state[s]++;
        //std::cerr << "Boundary up\n";
        crossed = (next_vertex >= initial_level || current_vertex < initial_level);
      } 
      else 
        crossed = (next_vertex >= initial_level && current_vertex < initial_level);
    }
    else {
      if (next_vertex>=current_vertex) {
        site_state[s]--;
        //std::cerr << "Boundary down\n";
        crossed=(next_vertex < initial_level || current_vertex >= initial_level);
      }
      else 
        crossed=(next_vertex < initial_level && current_vertex >= initial_level);
    }

    // Change position of is_upper_leg boolean for the next vertex
    is_upper_leg=!is_upper_leg;

    // eventual end of the worm here
    if ( next_vertex==initial_vertex && next_leg==initial_leg  && is_upper_leg==initial_is_upper_leg)
      break;
      
    if(crossed && measure_green_function_) {
      if (!measurement_origin_)
        green[distance(initial_site,s)] += worm_weight;
      else if (initial_site==measurement_origin_.get())
        green[s] += worm_weight;
    }    
  }
}


// Given the current vertex, and the incoming leg, looks into the 
// worm scattering probability array to return an exit leg
uint8_t SSE::return_exit_leg(uint32_t vertex,uint8_t incomingleg)
{
  uint8_t* MP= operator_string[vertex].leg;
  boost::multi_array<double, 2>  prob = proba_worm[bond_type[bond(operator_string[vertex].bond_number)]][MP[0]][MP[1]][MP[2]][MP[3]];

  double aa=random_real(); 
  if (aa<prob[0][incomingleg])
    return 0;
  else if (aa<prob[1][incomingleg])
    return 1;
  else if (aa<prob[2][incomingleg])
    return 2;
  else
    return 3;
}



// A Monte Carlo step
void SSE::do_update()
{
  worm_size=0.;

  // Eventually increase the cutoff
  /*
  if (current_number_of_non_identity>0.8*cutoff_L) {
     increase_cutoff_L(uint32_t(1+0.1*cutoff_L));
     cerr << "increase_cutoff_L not implemented.\n";
  }*/
  
  diagonal_update();
//std::cout<<"Nach diag update\n";
  link_legs();
//  std::cout<<"Nach link legs\n";

  if (current_number_of_non_identity || measure_green_function_) { // do worm update
 
    if (parms.defined("NUMBER_OF_WORMS_PER_SWEEP")) { 

      for (uint32_t kk=0;kk<number_of_worms_per_sweep;++kk)
        worm_update();
      worm_size/=number_of_worms_per_sweep;
    }
    else {
      if (!(is_thermalized())) {
             //  std::cout<<"Worms not_thermalized\n";

      // If in the thermalization part, do enough worms such that
      // the total worm size is larger than 2n 
      // (n=number of non-Identity vertices)
      
        int yy=0;
        while (worm_size<2*current_number_of_non_identity) { 
                  
          worm_update(); 
          yy++;
        }

        // If in the last 20% part of thermalization, make an average
        // on how many worms are needed to realize this condition
        if ((simulation_phase==1 && logf_step==logf_steps_total-1) || (simulation_phase==2 && block_sweeps>0.8*OE_nb_thermalization)) { 
          nb_worms_thermalization+=yy; 
          count_worms_thermalization++;
        }
      }
      else { // in the measurement part
      
        // For the very first step, calculate the average number of worms 
        // necessary and modify the value of number_of_worms_per_sweep
        // Then don't touch it
        if (count_worms_thermalization) {
          number_of_worms_per_sweep=int(nb_worms_thermalization/count_worms_thermalization);
          count_worms_thermalization=0; 
          green=0.;
         //           cout << "Number of worms per sweep : " << number_of_worms_per_sweep << endl;
        }
        
        // Do number_of_worms_per_sweep worm updates
        for (uint32_t kk=0;kk<number_of_worms_per_sweep;++kk)
          worm_update();
        worm_size/=number_of_worms_per_sweep;
      }
        }
  }
}
