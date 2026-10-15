/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1998-2018 ALPS Collaboration
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

// Legacy ALPS model construction and interpolation from src/apps/maxent/src/default_model.hpp
// at d4533fd30c65656cf5ab7fa448dd6b921b39aa6f; original notice:
/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2010 by Sebastian Fuchs <fuchs@comp-phys.org>
*                       Thomas Pruschke <pruschke@comp-phys.org>
*                       Matthias Troyer <troyer@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

// Adapted from CQMP/Maxent src/default_model.cpp at
// e727650a1f0a427e95efebc409d907b4c84011ec; see ../README.md.

#include "default_model.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include "maxent_string.hpp"



///This class deals with tabulated model functions
TabFunction::TabFunction(const alps::params& p, std::string const& name){
  std::string p_name = p[name].as<std::string>();
  std::ifstream defstream(p_name.c_str());
  if (!defstream)
    throw std::invalid_argument("could not open default model file: "+p[name].as<std::string>());
  std::string line;
  while (getline(defstream, line)) {
    std::stringstream linestream(line);
    linestream >> std::ws;
    if (linestream.eof() || linestream.peek() == '#')
      continue;
    double om, D;
    if (!(linestream >> om >> D) || !std::isfinite(om) || !std::isfinite(D))
      throw std::invalid_argument("expected finite omega and model value in: " + p_name);
    if (!Omega_.empty() && om <= Omega_.back())
      throw std::invalid_argument("default-model frequencies must be strictly increasing: " + p_name);
    // Retain ALPS's allowance for trailing text after the two numeric columns.
    Omega_.push_back(om);
    Def_.push_back(D);
  }
  if (Omega_.size() < 2)
    throw std::invalid_argument("default-model table needs at least two rows: " + p_name);
  double omega_max = p["OMEGA_MAX"];
  double omega_min = p.exists("OMEGA_MIN") ? p["OMEGA_MIN"] : -omega_max; //we had a 0 here in the bosonic case. That's not a good idea if you're continuing symmetric functions like chi(omega)/omega. Change omega_min to zero manually if you need it.
  if (Omega_[0]!=omega_min || Omega_.back()!=omega_max){
    std::cout<<"Omega[ 0] "<<Omega_[0]<<" omega min: "<<omega_min<<std::endl;
    std::cout<<"Omega[-1] "<<Omega_.back()<<" omega max: "<<omega_max<<std::endl;
  }
}

// Preserve ALPS interpolation and extend the end segments outside the table.
double TabFunction::operator()(const double omega) {
  std::vector<double>::const_iterator ub = std::upper_bound(Omega_.begin(), Omega_.end(), omega);
  const auto index = std::clamp<std::size_t>(ub - Omega_.begin(), 1, Omega_.size() - 1);

  double om1 = Omega_[index-1];
  double om2 = Omega_[index];
  double D1 = Def_[index-1];
  double D2 = Def_[index];
  return -(D2-D1)/(om2-om1)*(om2-omega)+D2;
}


GeneralDefaultModel::GeneralDefaultModel(const alps::params& p, std::shared_ptr<Model> mod)
: DefaultModel(p)
, Mod(mod)
, ntab(5001)
, xtab(ntab) {
  double sum = norm();
  for (int o=0; o<ntab; ++o) {
    xtab[o] *= 1./sum;
  }
}

///given a number x between 0 and 1, find the frequency omega belonging to x.
double GeneralDefaultModel::omega(const double x) const {
  //range check for x
  if(x>1. || x<0.)
    throw std::logic_error("parameter x is out of bounds!");
  std::vector<double>::const_iterator ub = std::upper_bound(xtab.begin(), xtab.end(), x);
  int omega_index = ub - xtab.begin();
  if (ub==xtab.end())
    omega_index = xtab.end() - xtab.begin() - 1;
  double om1 = omega_min + (omega_index-1)*(omega_max-omega_min)/(ntab-1);
  double om2 = omega_min + omega_index*(omega_max-omega_min)/(ntab-1);
  double x1 = xtab[omega_index-1];
  double x2 = xtab[omega_index];
  return -(om2-om1)/(x2-x1)*(x2-x)+om2;
}

/// returns the value of the model function at frequency omega
double GeneralDefaultModel::D(const double omega) const {
  return (*Mod)(omega);
}

//I have no idea what this does.
double GeneralDefaultModel::x(const double t) const {
  if(t>1. || t<0.) throw std::logic_error("parameter t is out of bounds!");
  int od = (int)(t*(ntab-1));
  if (od==(ntab-1))
    return 1.;
  double x1 = xtab[od];
  double x2 = xtab[od+1];
  return -(x2-x1)*(od+1-t*ntab)+x2;
}

double GeneralDefaultModel::norm() {
  double sum = 0;
  xtab[0] = 0.;
  //this is an evaluation on an equidistant grid; sum integrated by trapezoidal rule
  double delta_omega = (omega_max - omega_min) / (ntab - 1);
  for (int o = 1; o < ntab; ++o) {
    // This weighted-endpoint form gives an exact zero on symmetric grids.
    // That matters for models with a deliberately discontinuous zero branch
    // and avoids compiler-dependent FMA rounding at the branch point.
    double omega1 = ((ntab - o) * omega_min + (o - 1) * omega_max) / (ntab - 1);
    double omega2 = ((ntab - 1 - o) * omega_min + o * omega_max) / (ntab - 1);
    sum += ((*Mod)(omega1) + (*Mod)(omega2)) / 2. * delta_omega;
    xtab[o] = sum;
  }
  return sum;
}

std::shared_ptr<DefaultModel> make_default_model(const alps::params& parms, std::string const& name){
  std::string p_name = parms[name].as<std::string>();
  to_lower(p_name);
  if (p_name == "flat") {
    std::cout << "Using flat default model" << std::endl;
    return std::shared_ptr<DefaultModel>(new FlatDefaultModel(parms));
  }
  else if (p_name == "gaussian") {
    std::cout << "Using Gaussian default model" << std::endl;
    std::shared_ptr<Model> Mod(new Gaussian(parms));
    return std::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
  else if (p_name == "twogaussians" || p_name == "two gaussians") {
    std::cout << "Using sum of two Gaussians default model" << std::endl;
    std::shared_ptr<Model> Mod(new TwoGaussians(parms));
    return std::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
  else if (p_name == "shifted gaussian" || p_name == "shiftedgaussian") {
    std::cout << "Using shifted Gaussian default model" << std::endl;
    std::shared_ptr<Model> Mod(new ShiftedGaussian(parms));
    return std::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
  else if (p_name == "double gaussian" || p_name == "doublegaussian") {
    std::cout << "Using double Gaussian default model" << std::endl;
    std::shared_ptr<Model> Mod(new DoubleGaussian(parms));
    return std::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
  else if (p_name == "general double gaussian" || p_name == "generaldoublegaussian") {
    std::cout << "Using general double Gaussian default model" << std::endl;
    std::shared_ptr<Model> Mod(new GeneralDoubleGaussian(parms));
    return std::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
  else if (p_name == "linear rise exp decay" || p_name == "linearriseexpdecay") {
    std::cout << "Using linear rise exponential decay default model" << std::endl;
    std::shared_ptr<Model> Mod(new LinearRiseExpDecay(parms));
    return std::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
  else if (p_name == "quadratic rise exp decay" || p_name == "quadraticriseexpdecay") {
    std::cout << "Using quadratic rise exponential decay default model" << std::endl;
    std::shared_ptr<Model> Mod(new QuadraticRiseExpDecay(parms));
    return std::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
  else if (p_name == "lorentzian") {
    std::cout << "Using Lorentzian default model" << std::endl;
    std::shared_ptr<Model> Mod(new Lorentzian(parms));
    return std::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
  else if (p_name == "twolorentzians" || p_name == "two lorentzians") {
    std::cout << "Using sum of two Lorentzians default model" << std::endl;
    std::shared_ptr<Model> Mod(new TwoLorentzians(parms));
    return std::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
  else if (p_name == "shifted lorentzian" || p_name == "shiftedlorentzian") {
    std::cout << "Using shifted Lorentzian default model" << std::endl;
    std::shared_ptr<Model> Mod(new ShiftedLorentzian(parms));
    return std::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
  else if (p_name == "double lorentzian" || p_name == "doublelorentzian") {
    std::cout << "Using double Lorentzian default model" << std::endl;
    std::shared_ptr<Model> Mod(new DoubleLorentzian(parms));
    return std::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
  else { 
    std::cout << "Using tabulated default model" << std::endl;
    std::shared_ptr<Model> Mod(new TabFunction(parms, name));
    return std::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
}
