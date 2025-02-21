/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2001-2005 by Fabien Alet <alet@comp-phys.org>,
*                            Matthias Troyer <troyer@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
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

/* $Id: SSE.cpp,v 1.41 2005/09/16 16:20:51 honecker Exp $ */

#include <alps/osiris/comm.h>
#include "SSE.hpp"
#include <fstream>

void SSE::print_copyright(std::ostream& out)
{
  out << "Quantum Monte simulations using the generalized directed loop algorithm v. 1.1\n"
      << "  available from http://alps.comp-phys.org/\n"
      << "  copyright (c) 2001-2005 by Fabien Alet <alet@comp-phys.org>,\n"
      << "                             Synge Todo <wistaria@comp-phys.org>,\n"
      << "                             and Matthias Troyer <troyer@comp-phys.org>\n"
      << "  see F. Alet, S. Wessel and M. Troyer, Phys. Rev. E 71, 036706 (2005) for details.\n\n";
}

SSE::SSE(const alps::ProcessList& w, const alps::Parameters& myparms,int n)
        : QMCRun<>(w,myparms,n),
          steps_done_total(0),
          measurements_done(0),
          nb_worms_thermalization(0),
          count_worms_thermalization(0),
          site_state(num_sites(),0),
          operator_string(0),
          current_number_of_non_identity(0),
          boundary_crossing(alps::get_or_default(alps::boundary_crossing_t(),graph(),alps::boundary_crossing())),
      all_done(0),
          block_sweeps(0),
          thermalized(false),
      logf_step(1),
          thistime(0.),realtime(0.),
          upwalker(1)

{
    initialize_simulation();
        create_observables();
        initialize_current_simulation();  
    std::cout<<"Simulation initialized\n";
}

// Check if simulation is finished
double SSE::work_done() const
{
        if (!is_thermalized())
                return 0.;
        return ((all_done) ? 1. :0.);
}

// It is desirable to be able to increase the number of SWEEPS
bool SSE::change_parameter(const std::string& name, const StringValue& value) {
        uint64_t new_sweeps = 0;
        if(name=="SWEEPS")
                new_sweeps = uint64_t(value);
        if(name=="MCS")
                new_sweeps = uint64_t(value);
        if(name=="Steps")
                new_sweeps = uint64_t(value);
        // Is it sensoble to do it ?
        if(new_sweeps > 0) {
                nb_steps = new_sweeps;
                return true;
        }
        // Otherwise we cannot do it
        return false;
}

// Do a MC step - and do measurements (if thermalized)
void SSE::dostep()
{

        do_update();

        steps_done_total++;
        block_sweeps++;
        
    if (simulation_phase==1)
        do_WL_step();
    else     
        do_OE_step();
}


void SSE::do_WL_step() {
    if (is_thermalized()) {
        do_measurements();  // Physical measurements. We can only measure after g has been fixed.
                measurements_done++;            
                if (measurements_done>=nb_measurement_steps && all_done!=1.) {
                    finish_measurements();
            store_histo_vars();
                        cerr << "[" << norder_min << "-" << norder_max << "]  all done." << endl;
                        all_done=1.;
                }
    } else {
        if (block_sweeps >= block_sweeps_total) {
            cerr << "Successfully done one sweep with "<<block_sweeps_total<<" sweeps\n";         
            block_sweeps=0;
            g.subtract(); 
              double histomin=histo.min();
              double histoflatness=histo.flatness();
              cerr << "operator_string length: " << operator_string.size() << "non-ids: " << current_number_of_non_identity << "\n";
              cerr << "[" << norder_min << "-" << norder_max << "]  step "  << logf_step << ", ln[f]="<< logf 
                       << " : flatness="<< histoflatness << " ratio=" << histomin/minimum_histogram <<endl;
              if (histoflatness<flatness_treshold && histomin>=minimum_histogram) {
                    histo.fill(0.0); 
                cerr<<"current_logf_step: "<<logf_step<<"\n"; 
                    if (logf_step==logf_steps_total) {  
                          thermalized=true;
                    measurements_done=0;
                    // Update the final WL weights:
                    for (unsigned int i=norder_min;i<=norder_max;++i) {
                        final_weight_fract[i]=exp(g[i]-g[i+1]);  // = w(n+1)/w(n)
                      }
                          cerr << "WL: [" << norder_min << "-" << norder_max << "]  continuing using final weights..." << endl;   
                    } else {
                    ++logf_step; 
                    logf/=2.0;
                    if (use_zhou_bhatt) {
                        block_sweeps_total=static_cast<unsigned int>(block_sweeps_total*1.41);
                        minimum_histogram*=2;
                          }
                    }
              }       
                }
    }
}


void SSE::do_OE_step() {
    if (all_done==1.) 
        return;
    if (is_thermalized()) {
        do_measurements();  // Physical measurements - we can always measure (as soon as it is thermalized), since det. bal. is fulfilled.
        measurements_done++;
    } else if (block_sweeps>=OE_nb_thermalization)  // Wait some time before starting measurements    
        thermalized=true;    
    
    if (block_sweeps >= block_sweeps_total) { // In this implementation, only 1 optimization step is performed per OE run. I.e., we stop now:
        finish_measurements();
        store_histo_vars();
        all_done=1.;            
    }
}




void SSE::store_histo_vars() {
  if (simulation_phase==1) {  // Only in this case we store g. For OE sims, we already started with an \
                  //  external g from the LOGG_FILE input parameter
      valarray<double> gval=g.getvalarray(norder_min,norder_max);
      measurements["logg"] << gval;
  }
  measurements["Histogram"] << histo.getvalarray(norder_min,norder_max);
  measurements["HistoUp"] << histoup.getvalarray(norder_min,norder_max);
  for (int i=0; i<timeup.size();++i) {
           measurements["Time Up"]<<timeup[i];
    measurements["RealTime Up"]<<realtimeup[i];
  }     
  for (int i=0; i<timedown.size();++i) {
        measurements["Time Down"]<<timedown[i];
        measurements["RealTime Down"]<<realtimedown[i];
  }
}




void SSE::initialize_current_simulation() {

  final_weight_fract.resize(norder_max-norder_min+1,norder_min);  
  histo.resize(norder_max-norder_min+1,norder_min);
  histoup.resize(norder_max-norder_min+1,norder_min);
  timeup.resize(0);
  timedown.resize(0);
  realtimeup.resize(0);
  realtimedown.resize(0);
  if (parms.defined("NUMBER_OF_WORMS_PER_SWEEP")) 
      number_of_worms_per_sweep=parms.value_or_default("NUMBER_OF_WORMS_PER_SWEEP",1);
  

  simulation_phase=parms.value_or_default("SIMULATION_PHASE",1);
 
  if (simulation_phase==1) { // WL sim
    g.resize(norder_max-norder_min+1,norder_min);
      logf_steps_total=parms.value_or_default("NUMBER_OF_WANG_LANDAU_STEPS",12);     
      use_zhou_bhatt=parms.value_or_default("USE_ZHOU_BHATT_METHOD",1.);
    if (use_zhou_bhatt) {
        block_sweeps_total=g.size();
            logf=log(double(parms.value_or_default("INITIAL_MODIFICATION_FACTOR",exp(1.))));
            minimum_histogram=1./logf;
            flatness_treshold=parms.value_or_default("FLATNESS_TRESHOLD",1E10);
      } else {
            block_sweeps_total=parms.value_or_default("WL_BLOCK_SWEEPS",10000);
            logf=log(double(parms.value_or_default("INITIAL_INCREASE_FACTOR",exp(cutoff_L*log((double)num_sites())/block_sweeps_total))));
            minimum_histogram=0.;
            flatness_treshold=parms.value_or_default("FLATNESS_TRESHOLD",0.2);
      }
      
  } else { // OE sim
    std::string logg_filename((std::string)parms.value_or_default("LOGG_FILENAME","logg.dat"));
        block_sweeps_total = parms.value_or_default("OE_BLOCK_SWEEPS", 20000);
    OE_nb_thermalization=parms.value_or_default("OE_NB_THERMALIZATION", block_sweeps_total/10);
      g.resize(norder_max-norder_min+1,norder_min);
      cerr << "Loading logg from file "<<logg_filename<<"\n";
      ifstream myfile (logg_filename.c_str());
      if (myfile.is_open()) {
            int i=norder_min;
            double tmplogg=0.;
            while ( myfile>>tmplogg) {
                  if (i<=norder_max)
                    g[i] = tmplogg;
                  else {
                cerr << "WARNING: Logg file larger than norder_max! Disregarding further input\n";
                break;
                  }
                  ++i;
            }
            myfile.close();
      } else {
        cout << "Unable to open logg file!\n";
        exit(0);
      }
    //Convert logg to weight fractions:
    //final_weight_fract[i] gives the *ratio* of w(n+1)/w(n) !!!!!
      for (int i = norder_min; i<norder_max;++i) {
              final_weight_fract[i]= exp(-g[i+1]+g[i]); //g(n) is in Log!!          
      }
  }
}

// This is not optimal, since this function is called for each operator and does lots of checks...
double SSE::ensemble_weight_fraction(int norder) {
  if (simulation_phase==2 || is_thermalized())
    return final_weight_fract[norder];
  else 
    return ((norder<norder_min)?(1.):(exp(g[norder]-g[norder+1])));

}

// This is executed after a complete diagonal update sweep was performed.
void SSE::diagonal_update_sweep_end() { 
  // IMPORTANT: Independently of the op string representation, increment the time at the end of complete traversals,
  // since this determines how long the rest of the (off-diagonal) updates will take. And this performance is independent of the representation.
  if (is_thermalized()) { // WL or OE sim are thermalized
    thistime++;
    realtime+=(double)current_number_of_non_identity/cutoff_L;  // Time rescaled by the current op string length
     
    // IMPORTANT: It turns out that even the histos of a WL sim in variable length should continuously be updated, and not after sweep end.
    // Uncomment all the following lines if you still want to try the non-optimal case of updating after complete sweeps only.
    // If so, remember to *comment* in diagonal_update_iteration_end, as indicated.
    
    // histo[current_number_of_non_identity]++;
    // if (upwalker) histoup[current_number_of_non_identity]++;    
   }
   // else if (simulation_phase==1) {    
    //if (variable_length_representation) {
    //    g[current_number_of_non_identity]+=logf;
    //    histo[current_number_of_non_identity]++;
    //}
// }

}

//This is called at the end of a diagonal update step *for each* operator
void SSE::diagonal_update_iteration_end(int i) { 
  // Test for upper/lower boundary after continuously, independent of the representation. If this is done differently, variable length gives a fraction 
  // with a large, discontinuous "jump" from f(n)>0 to 0 at n=norder_max
  if (upwalker && current_number_of_non_identity==norder_max) {
    upwalker=0;
        if (is_thermalized()) { 
          timeup.push_back((double)thistime);
      realtimeup.push_back(realtime);
      thistime=0;
          realtime=0.;
        }
  } else if (!upwalker && current_number_of_non_identity==norder_min) {
    upwalker=1; 
        if (is_thermalized()) {
      timedown.push_back((double)thistime);
          realtimedown.push_back(realtime);
          thistime=0;
      realtime=0.;
        }
  }
  // IMPORTANT: If you have uncommented as indicated in diagonal_update_sweep_end, you alsa HAVE to comment the following lines until the end of this 
  // function. It's a simple case: For fixed length, it does not matter at all (except more data). For variable length, you either optimize the local
  // current (if you leave the code) or you optimize the "global" current (ie. measured after complete string traversals) (if you change the code).
  // Norbert Stoop, Dec. 2006
  if (is_thermalized()) {  // WL or OE sim are thermalized, so measure those thingies...
        histo[current_number_of_non_identity]++;
    if (upwalker) histoup[current_number_of_non_identity]++;
  } else if (simulation_phase==1) { // In the WL phase, we always "measure" the histo and g
    g[current_number_of_non_identity]+=logf;
        histo[current_number_of_non_identity]++;
  } 
}



// save simulation
void SSE::save(alps::ODump& dump) const
{
  dump <<  steps_done_total << site_state << operator_string << cutoff_L <<  current_number_of_non_identity << number_of_worms_per_sweep << nb_worms_thermalization << count_worms_thermalization<<block_sweeps << measurements_done << block_sweeps_total << thermalized << upwalker << timeup << timedown << realtimeup << realtimedown << simulation_phase;
  g.save(dump);
  histo.save(dump);
  histoup.save(dump);
 
  if (simulation_phase==1) { // WL simulation
      dump << minimum_histogram << flatness_treshold<<logf;
  }
  
}

// load simulation
void SSE::load(alps::IDump& dump)
{
  dump >> steps_done_total>> site_state>> operator_string>> cutoff_L>> current_number_of_non_identity>> number_of_worms_per_sweep>> nb_worms_thermalization>> count_worms_thermalization>> block_sweeps>> measurements_done>> block_sweeps_total>> thermalized>> upwalker>> timeup>> timedown>> realtimeup>> realtimedown>>simulation_phase;
  g.load(dump);
  histo.load(dump);
  histoup.load(dump);
 
  if (simulation_phase==1) { // WL simulation
      dump>> minimum_histogram >> flatness_treshold >> logf;
  }
}



int main(int argc, char** argv)
{
        
#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif
   return alps::scheduler::start(argc,argv,alps::scheduler::SimpleMCFactory<SSE>());
#ifndef BOOST_NO_EXCEPTIONS
  }
  catch (std::exception& exc) {
    std::cerr << exc.what() << "\n";
      alps::comm_exit(true);
      return -1;
    }
  catch (...) {
    std::cerr << "Fatal Error: Unknown Exception!\n";
    return -2;
  }
#endif
        

}
