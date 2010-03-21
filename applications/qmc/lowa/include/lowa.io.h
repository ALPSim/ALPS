/****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Applications
*
* Copyright (C) 2006-2010 by Lode Pollet <lpollet@physics.harvard.edu>,
*                            Ping Nang Ma <pingnang@itp.phys.ethz.ch>,
*                            Matthias Troyer <troyer@itp.phys.ethz.ch>
*
* This software is part of the ALPS libraries, published under the ALPS
* Library License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
*
* You should have received a copy of the ALPS Library License along with
* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
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

/* $Id: lowa.io.h 3520 2010-03-03 16:55:00Z tamama $ */


/*
 *
 * 1) Code modification      -- partly done
 * 2) Replacing raw pointers -- not done yet
 *
 */



#ifndef ALPS_APPLICATIONS_LOWA_IO_H
#define ALPS_APPLICATIONS_LOWA_IO_H


namespace alps {
namespace applications {


void Lowa::save_lowa()
{
  outFile.open(filename1_trial.c_str(), std::ios::binary);
  outFile << std::setprecision(20);
  export_lowa_simulation(outFile);
  outFile.close();
  int renameinfo1 = std::rename(filename1_trial.c_str(),filename1.c_str());
}



void Lowa::load_lowa()
{
  inFile.open(filename1.c_str(), std::ios::binary);
  if (inFile.good()) { // continuing a previous simulation
    std::cout << "\n# Retrieving old configuration....";
    import_lowa_simulation(inFile);
    std::cout << "\n# Retrieving old configuration done. Testing...";
    if (!test_conf()) {std::cin.get();};
    std::cout << "\n# Testing OK...\n";
    update_N_and_E();
    std::cout << "# Number of particles : " << _nrpart << std::endl;
    std::cout << "# Worm head           : " << worm_head << std::endl;
    std::cout << "# Worm tail           : " << worm_tail << std::endl;
    std::cout << "# wormdiag - wormdir  : " << is_worm_diagonal << "\t" << is_worm_moving_forward << std::endl;
    std::cout << "# wormrising          : " << is_worm_rising << std::endl;
  }
  else { // starting a new simulation
#ifdef HARDCORE
    for (site_type s=0; s < _N; ++s)  {  _state[s] = ((random_real()<0.5) ? 0 : 1);  }
#else
    for (site_type s=0; s < _N; ++s)  {  _state[s] = 0;  }
#endif
#ifdef TRAPPEDSYSTEM
    for (site_type i = 0; i < _N; ++i) {
      if (is_border_site(i) == 1) {
        _state[i] = 0;
        _mu_eff[i] = std::abs(_mu_eff[i])*10000;
      }
    }
#endif
    // insert dummy kinks on all sites, keeping the code ler and making a line of where to measure diag properties
    for (site_type i = 0; i < _N; ++i) {
      kink_type new_elem(zcoord,_state[i],_state[i], i, i, beta);
      site_it[i]       = _kinks_on_lattice[i].insert(site_it[i], new_elem);
      dummy_kink_it[i] = site_it[i];
#ifdef DEBUGMODE
      std::cout << "\nInserting dummy : " << i << "\t " << site_it[i]->time() << "\t" << site_it[i]->before();
#endif
    }
    for (site_type i = 0; i < _N; ++i) {
      for (dim_type z = 0; z < zcoord; ++z)  {  site_it[i]->set_assoc(z, site_it[nb(i,z)]);  }
    }
#ifdef DEBUGMODE
    //print_configuration(std::cout);
#endif
  }
  inFile.close();
}


void Lowa::export_lowa_simulation(std::ostream& os)
{
  os << counter_MEASURE << "\t" << counter_TEST << "\n";
  os << MCstep << "\t" << MCstep_total << "\n";
  os << worm_head << std::endl;
  os << worm_tail << std::endl;
  os << _state_enclosed_within_wormpair_at_creation << std::endl;
  os << *worm_it << std::endl;
  for (site_type s = 0; s < _N; s++) {
    int n = 0;
    os << _kinks_on_lattice[s].size() << std::endl;
    kinks_iterator_type it = dummy_kink_it[s];
    os << *it << std::endl;
    it->set_name(n); n++;
    ++it; if (it == _kinks_on_lattice[s].end()) it = _kinks_on_lattice[s].begin();
    while (it != dummy_kink_it[s]) {
      os << *it << std::endl;
      it->set_name(n); n++;
      ++it; if (it == _kinks_on_lattice[s].end()) it = _kinks_on_lattice[s].begin();
    }
  }
  // save associations
  for (site_type s = 0; s < _N; s++) {
    kinks_iterator_type it = dummy_kink_it[s];
    for (int j = 0; j < zcoord; j++) os << (it->assoc(j))->name() << "\t";
    os << std::endl;
    ++it; if (it == _kinks_on_lattice[s].end()) it = _kinks_on_lattice[s].begin();
    while (it != dummy_kink_it[s]) {
      for (int j = 0; j < zcoord; j++) os << (it->assoc(j))->name() << "\t";
      os << std::endl;
      ++it; if (it == _kinks_on_lattice[s].end()) it = _kinks_on_lattice[s].begin();
    }
  }

  for (site_type s = 0; s < _N; s++) {
    os << av_dnsmat[s] << "\t" << av_dnsmat_inf[s] << std::endl;
  }
  os << _Z_dnsmat << std::endl;

  os << is_worm_diagonal << "\t" << is_worm_moving_forward << "\t" << is_worm_rising << std::endl;
  os << _nrvertex << std::endl;
#ifdef BOUND_WORM
  os << M__ << std::endl;
  os << correct_lattice_structure << std::endl;
#endif
  os << _Ekin << "\t" << _Epot << "\t" << new_measurement << std::endl;
  os << no_of_accepted_worm_insertions << "\t" << no_of_proposed_worm_insertions << std::endl;
  os << sweeps << std::endl;
  return;
}


void Lowa::import_lowa_simulation(std::istream& is)
{
  kinks_iterator_type** v;
  v = new kinks_iterator_type* [_N];
  kink_type a_worm_elem(zcoord);
  is >> counter_MEASURE >> counter_TEST;
  is >> MCstep >> MCstep_total;

  MCstep = 0.;  MCstep_total = 0. ;

  is >> worm_head;
  is >> worm_tail;
  is >> _state_enclosed_within_wormpair_at_creation;
  is >> a_worm_elem;
  for (site_type s = 0; s < _N; s++) {
    int64_t opsize_n;
    is >> opsize_n;
    v[s] = new kinks_iterator_type [opsize_n];
    for (int32_t i = 0; i < opsize_n; i++) {
      kink_type elem(zcoord);
      is >> elem;
      _kinks_on_lattice[s].push_back(elem);
      if (elem == a_worm_elem) { worm_it = _kinks_on_lattice[s].end(); --worm_it;}
    }
    int32_t i = 0;
    for (kinks_iterator_type it=_kinks_on_lattice[s].begin(); it != _kinks_on_lattice[s].end(); ++it, i++) {
      v[s][i]=it;
    }
    dummy_kink_it[s] = _kinks_on_lattice[s].begin();
    site_it[s] = dummy_kink_it[s];
  }
  // load associations
  for (site_type i = 0; i < _N; i++) {
    int ass;
    kinks_iterator_type it = dummy_kink_it[i];
    for (int j = 0; j < zcoord; j++) {
      site_type s = nb(i,j);
      is >> ass;
      it->set_assoc(j, v[s][ass] );
    }
    ++it; if (it == _kinks_on_lattice[i].end()) it = _kinks_on_lattice[i].begin();
    while (it != dummy_kink_it[i]) {
      for (int j = 0; j < zcoord; j++) {
        site_type s = nb(i,j);
        is >> ass;
        it->set_assoc(j, v[s][ass] );
      }
      ++it; if (it == _kinks_on_lattice[i].end()) it = _kinks_on_lattice[i].begin();
    }
  } // ... associations

  for (site_type s = 0; s < _N; s++) {
    is >> av_dnsmat[s] >> av_dnsmat_inf[s];
  }
  is >> _Z_dnsmat;

  is >> is_worm_diagonal >>  is_worm_moving_forward >> is_worm_rising;
  is >> _nrvertex;
#ifdef BOUND_WORM
  is >> M__;
  is >> correct_lattice_structure;
#endif
  is >> _Ekin >> _Epot >> new_measurement;
  is >> no_of_accepted_worm_insertions >> no_of_proposed_worm_insertions;
  is >> sweeps;

  calc_N_and_E();
  for (int i = 0; i < _N; i++) {
    delete [] v[i];
  }
  delete [] v;
  site_it[worm_head.to()] = worm_it;
#ifdef DEBUGMODE
  std::cout << "\n dummy :\n ";
  for (site_type i = 0; i < _N; i++) {
    std::cout << dummy_kink_it[i]->time() << "\t";
  }
  std::cout << "\nWorm iterator : " << worm_it->time() << "\t" << &*worm_it << "\n";
  print_configuration(std::cout);
#endif
}


std::ostream& Lowa::print_params(std::ostream& out) const
{
  out << "\n Parameters : ";
  out << "\n -------------";
  out << "\n label           :  " << label;
  out << "\n dim             :  " << dim;
  out << "\n v0       [Er]   :  "; for(int k =0; k < dim; k++) out << v0[k] << "\t";
  out << "\n _mu      [nK]   :  " << _mu;
  out << "\n vc       [nK]   :  "; for(int k =0; k < dim; k++) out << vc[k] << "\t";
  out << "\n waist    [um]   :  "; for(int k =0; k < dim; k++) out << waist[k] << "\t";
  out << "\n lambda   [nm]   :  "; for(int k =0; k < dim; k++) out << lambda[k] << "\t";
  out << "\n beta     [nK-1] :  " << beta;
  out << "\n Ls              :  "; for(int k =0; k < dim; k++) out << Ls[k] << "\t";
  out << "\n _N              :  " << _N;
  out << "\n _statemax       :  " << _statemax;
  out << "\n _Eoffset        :  " << _Eoffset;
  //out << "\n Ntherm        :  " << Ntherm;
  //out << "\n Nloop         :  " << Nloop;
  if (Ncan >= 0) {
    out << "\n Ncan            :  " << Ncan;
  }
  else {
    out << "\n Grand can. ";
  }
  out << "\n Nmeasure : " << Nmeasure;
  out << "\n Ntest    : " << Ntest;
  out << "\n Nsave    : " << Nsave;
  out << std::endl << std::endl << std::endl;
}


std::ostream& Lowa::print_configuration(std::ostream& out)
{
  out << "\n\n Printing operator string";
  for (site_type i = 0; i < _N; i++) {
    out << "\nSite : " << i;
    int n = _kinks_on_lattice[i].size();
    int ii = 0;
    for (kinks_iterator_type it = _kinks_on_lattice[i].begin(); it != _kinks_on_lattice[i].end(); ++it, ++ii) {
      it->print();
      if (ii > n) {
        std::cout << "\n Error with std::list!\n";
        char ch; std::cin >> ch;
      }
    }
    out << "\n--------------------\n\n\n";
  }
#ifdef DEBUGMODE
  std::cout << "\n Worm site : " << worm_head.to() << "\t time " << worm_head.time() << "\t dir "<< is_worm_moving_forward << "\t rising " << is_worm_rising << "\n";
#endif
  return out;
}


}
}


#endif
