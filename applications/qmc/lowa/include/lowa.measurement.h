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

/* $Id: lowa.measurement.h 3520 2010-03-03 16:55:00Z tamama $ */


/*
 *
 * 1) Code modification      -- partly done
 * 2) Replacing raw pointers -- not done yet
 *
 */


#ifndef ALPS_APPLICATIONS_LOWA_MEASUREMENT_H
#define ALPS_APPLICATIONS_LOWA_MEASUREMENT_H


namespace alps {
namespace applications {


void Lowa::calc_N_and_E()
{
  _Ekin = -tx_plus(_site0)*(_nrvertex/beta);   // obviously not correct for non-isotropic and site dependent t and U... Tama: Please change it!
  _Epot = 0.;

  _nrpart=0;  
  _meas_nrpart=0;
  for (site_type j = 0; j < _N; ++j) 
  {
    fock_basis_type n = dummy_kink_it[j]->before();
    _state[j] = n;

    _nrpart += n;
    if (!is_doublon_treated_as_hole)  {  _meas_nrpart += n; }
    else                              {  _meas_nrpart += (n%2); }

    _Epot += (0.5*U(j)*n*(n-1) + _mu_eff[j]*n);
  }
}


void Lowa::statebinning()
{
  for (site_type index=0; index < _nrbin; ++index)  {  _proj_binstate[index] = 0;  _cs_binstate[index] = 0;  }

  for (site_type p=0, j=0; j < Ls[1]; ++j) {
    for (site_type i=0; i < Ls[0]; ++i, ++p) {
      site_type cur_proj_bin = _proj_binnr[p];

      if (!is_doublon_treated_as_hole)
      {
        for (site_type k=0; k < Ls[2]; ++k) {  _proj_binstate[cur_proj_bin] += _state[p + k*_Nxy];  }
         _cs_binstate[cur_proj_bin] += _state[p + (Ls[2]/2)*_Nxy];
      }
      else
      {
        for (site_type k=0; k < Ls[2]; ++k) {  _proj_binstate[cur_proj_bin] += (_state[p + k*_Nxy] % 2);  }
        _cs_binstate[cur_proj_bin] += (_state[p + (Ls[2]/2)*_Nxy] % 2);
      }

    }
  }

  for (site_type index=0; index < _nrbin; ++index) {
    _proj_binstate[index] /= _proj_binfreq[index];
    _cs_binstate[index]   /= _proj_binfreq[index];
  }
}


void Lowa::take_diagonal_measurements()
{
  outFile.open(filename_N.c_str(),std::ios::app);       outFile << get_Npart() << "\n";               outFile.close();
  outFile.open(filename_mN.c_str(),std::ios::app);      outFile << get_Nmeaspart() << "\n";           outFile.close();
  outFile.open(filename_dns0.c_str(),std::ios::app);    outFile << static_cast<int>(get_density(_site0))     << "\n";   outFile.close();
  outFile.open(filename_mdns0.c_str(),std::ios::app);   outFile << static_cast<int>(get_measdensity(_site0)) << "\n";   outFile.close();

  outFile.open(filename_proj_cymdns.c_str(),std::ios::app);
  for (site_type index=0; index < _nrbin; ++index)  {  outFile << _proj_binstate[index] << "\t";  }
  outFile << "\n";
  outFile.close();

  outFile.open(filename_cs_cymdns.c_str(),std::ios::app);
  for (site_type index=0; index < _nrbin; ++index)  {  outFile << _cs_binstate[index] << "\t";  }
  outFile << "\n";
  outFile.close();

#ifdef MEASURE_TIME_SERIES_DENSITY
  alps::hdf5::oarchive oa_mdns(filename_mdns.c_str());

  std::string cur_desc_mdns = "TimeSeries/Density/Set" + boost::lexical_cast<std::string>(sweeps);

  oa_mdns << alps::make_pvp("No_of_datasets",sweeps);
  oa_mdns << alps::make_pvp(cur_desc_mdns,_state,_N);

/*
  std::string cur_filename_mdns_ASCII = "./timeseries_density_measurements/" + filename_mdns + "_" + ss.str() + ".o";
  outFile.open(cur_filename_mdns_ASCII.c_str(),std::ios::app);
  for (site_type index=0; index < _N; ++index)  {  outFile << static_cast<int>(get_measdensity(index)) << "\t";  }
  outFile << "\n";
  outFile.close();
*/
#endif
            
  //measurements["Total Particle Number (Actual)"]     << get_Npart();
  //measurements["Total Particle Number (Measured)"]   << get_Nmeaspart();
  //measurements["Density at center (Actual)"]         << _state[_site0];
  //measurements["Density at center (Measured)"]       << (_state[_site0]%2);
  //measurements["Kinetic Energy"]                     << get_kinetic_energy();
  //measurements["Potential Energy"]                   << get_potential_energy();
  //measurements["Energy"]                             << get_energy();

  //measurements["Columnn integrated Density (cylindrically binned)"]    << _proj_binstate;
  //measurements["Cross sectional Density (cylindrically binned)"]       << _cs_binstate;
}


void Lowa::take_offdiagonal_measurements()
{
#ifdef MEASURE_TIME_SERIES_DENSITY_MATRIX
  // *** replace this by hdf5 measurements here
  outFile.open(filename_mdnsmat.c_str(),std::ios::app);
  for (site_type index=0; index < _N; ++index)  {  outFile << (static_cast<obs_type>(av_dnsmat[index])/_Z_dnsmat) << "\t";  }
  outFile << "\n";
  outFile.close();

  outFile.open(filename_mdnsmatinf.c_str(),std::ios::app);
  for (site_type index=0; index < _N; ++index)  {  outFile << (static_cast<obs_type>(av_dnsmat_inf[index])/_Z_dnsmat) << "\t";  }
  outFile << "\n";
  outFile.close();
#endif  
}


}
}

#endif
