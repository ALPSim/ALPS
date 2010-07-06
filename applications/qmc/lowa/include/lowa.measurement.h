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
 * 1) Code modification                                                  -- partly done
 * 2) Replacing raw pointers                                             -- not done yet
 * 3) Replacing ugly loops to std::for_each, std::copy, std::transform, with the help of boost::lambda         -- partly done
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
  std::fill(&_proj_binstate[0], &_proj_binstate[_nrbin], 0);
  std::fill(&_cs_binstate[0], &_cs_binstate[_nrbin], 0);

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

  std::transform(&_proj_binstate[0], &_proj_binstate[_nrbin], &_proj_binfreq[0], &_proj_binstate[0], boost::lambda::_1 / boost::lambda::_2);
  std::transform(&_cs_binstate[0], &_cs_binstate[_nrbin], &_proj_binfreq[0], &_cs_binstate[0], boost::lambda::_1 / boost::lambda::_2);
}


void Lowa::update_off_diag_obs()
{
  if (measure_time_series_density_matrix)
  {
    std::for_each(&av_dnsmat[0], &av_dnsmat[_N], boost::lambda::_1 /= _Z_dnsmat);
    std::for_each(&av_dnsmat_inf[0], &av_dnsmat_inf[_N], boost::lambda::_1 /= _Z_dnsmat);
  }
}


void Lowa::take_diagonal_measurements()
{
  // These are useful for checking at intermediate times
  outFile.open(filename_N.c_str(),std::ios::app);       outFile << get_Npart() << "\n";               outFile.close();
  outFile.open(filename_dns0.c_str(),std::ios::app);    outFile << static_cast<int>(get_density(_site0))     << "\n";   outFile.close();

  if (is_doublon_treated_as_hole)
  {
    outFile.open(filename_mN.c_str(),std::ios::app);      outFile << get_Nmeaspart() << "\n";           outFile.close();
    outFile.open(filename_mdns0.c_str(),std::ios::app);   outFile << static_cast<int>(get_measdensity(_site0)) << "\n";   outFile.close();
  }

  outFile.open(filename_proj_cymdns.c_str(),std::ios::app);
  std::copy(&_proj_binstate[0], &_proj_binstate[_nrbin], std::ostream_iterator<time_type>(outFile, "\t"));
  outFile << "\n";
  outFile.close();

  outFile.open(filename_cs_cymdns.c_str(),std::ios::app);
  std::copy(&_cs_binstate[0], &_cs_binstate[_nrbin], std::ostream_iterator<time_type>(outFile, "\t"));
  outFile << "\n";
  outFile.close();

  if (measure_time_series_density)
  {
    alps::hdf5::oarchive oa_mdns(filename_mdns.c_str());

    std::string cur_desc_mdns = "TimeSeries/Density/Set" + boost::lexical_cast<std::string>(sweeps);

    oa_mdns << alps::make_pvp("No_of_datasets",sweeps);
    oa_mdns << alps::make_pvp(cur_desc_mdns,_state,_N);
  }
  else
  {
    _Z_dns += 1.;
    std::transform(&av_dns[0], &av_dns[_N], &_state[0], &av_dns[0], boost::lambda::_1 + boost::lambda::_2);

    outFile.open(filename_dns_trial.c_str(),std::ios::out);
    std::transform(&av_dns[0], &av_dns[_N], std::ostream_iterator<obs_type>(outFile, "\n"), boost::lambda::_1 / _Z_dns);
    outFile.close();
    int renameinfo1 = std::rename(filename_dns_trial.c_str(),filename_dns.c_str());
  }
}


void Lowa::take_offdiagonal_measurements()
{
  if (measure_time_series_density_matrix)
  {
    // finite TOF...
    alps::hdf5::oarchive oa_mdnsmat(filename_mdnsmat.c_str());

    std::string cur_desc_mdnsmat     = "TimeSeries/DensityMatrix/Set" + boost::lexical_cast<std::string>(sweeps_green);
   
    oa_mdnsmat << alps::make_pvp("No_of_datasets",sweeps_green);
    oa_mdnsmat << alps::make_pvp(cur_desc_mdnsmat,av_dnsmat,_N);


    // infinite TOF...
    alps::hdf5::oarchive oa_mdnsmat_inf(filename_mdnsmat_inf.c_str());

    std::string cur_desc_mdnsmat_inf = "TimeSeries/DensityMatrixInfinity/Set" + boost::lexical_cast<std::string>(sweeps_green);

    oa_mdnsmat_inf << alps::make_pvp("No_of_datasets",sweeps_green);
    oa_mdnsmat_inf << alps::make_pvp(cur_desc_mdnsmat_inf,av_dnsmat_inf,_N);
  }
  else
  {
    // finite TOF...
    outFile.open(filename_dnsmat_trial.c_str(),std::ios::out);
    std::transform(&av_dnsmat[0], &av_dnsmat[_N], std::ostream_iterator<obs_type>(outFile, "\n"), boost::lambda::_1 / _Z_dnsmat);
    outFile.close();
    int renameinfo1 = std::rename(filename_dnsmat_trial.c_str(),filename_dnsmat.c_str());

    // infinite TOF...
    outFile.open(filename_dnsmat_inf_trial.c_str(),std::ios::out);
    std::transform(&av_dnsmat_inf[0], &av_dnsmat_inf[_N], std::ostream_iterator<obs_type>(outFile, "\n"), boost::lambda::_1 / _Z_dnsmat);
    outFile.close();
    int renameinfo2 = std::rename(filename_dnsmat_inf_trial.c_str(),filename_dnsmat_inf.c_str());
  }
}


}
}

#endif
