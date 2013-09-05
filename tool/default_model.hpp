/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2010 by Sebastian Fuchs <fuchs@comp-phys.org>
*                       Thomas Pruschke <pruschke@comp-phys.org>
*                       Matthias Troyer <troyer@comp-phys.org>
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


#ifndef ALPS_TOOL_DEFAULT_MODEL_HPP
#define ALPS_TOOL_DEFAULT_MODEL_HPP

#include <math.h>
#include <alps/parameter.h>
#include <alps/ngs.hpp>
#include <boost/shared_ptr.hpp>


//Note the slightly crooked structure here:

//class DefaultModel
// -> FlatDefaultModel : public DefaultModel
// -> GeneralDefaultModel : public DefaultModel
//  and GeneralDefaultModel contains a 'Model' object.

//class Model
// -> Gaussian : public Model 
//    -> ShiftedGaussian : public Gaussian
//      -> DoubleGaussian : public ShiftedGaussian
//        -> GeneralDoubleGaussian : public ShiftedGaussian
// ->  TabFunction : public Model
// ->  LinearRiseExpDecay : public Model
// ->  QuadraticRiseExpDecay : public Model

class DefaultModel 
{
public:
//  DefaultModel(const alps::Parameters& p) :
  DefaultModel(const alps::params& p) :
    omega_max(p["OMEGA_MAX"]),
    omega_min(static_cast<double>(p["OMEGA_MIN"]|(-omega_max))), //we had a 0 here in the bosonic case. That's not a good idea if you're continuing symmetric functions like chi(omega)/omega. Change omega_min to zero manually if you need it.
    blow_up_(p["BLOW_UP"]|1.)
  { //std::cout<<"found omega_min: "<<omega_min<<std::endl;
    //std::cout<<"found omega_max: "<<omega_max<<std::endl;
    //std::cout<<"found blowup:    "<<blow_up_<<std::endl;
  }

  virtual ~DefaultModel(){}
  
  //omega returns a frequency point, given x between 0 and 1.
  virtual double omega(const double x) const = 0;
  
  //D returns the derivative of the integrated default model
  virtual double D(const double omega) const = 0;
  
  //returns the integrated default model
  virtual double x(const double t=0) const = 0;
  
  //equidistant mapping from [0,1] to [omega_min, omega_max]
  double omega_of_t(const double t) const { return omega_min + (omega_max-omega_min)*t; }
  //equidistant mapping from [omega_min, omega_max] to [0,1]
  double t_of_omega(const double omega) const { return (omega-omega_min)/(omega_max-omega_min); }
  double blow_up() const { return blow_up_; }

protected:
  const double omega_max;
  const double omega_min;
  const double blow_up_;
};



class FlatDefaultModel : public DefaultModel
{
public:

//    FlatDefaultModel(const alps::Parameters& p) : DefaultModel(p) {}
  FlatDefaultModel(const alps::params& p) : DefaultModel(p) {}

  double omega(const double x) const {
    return x/blow_up()*(omega_max-omega_min) + omega_min;
  }

  double D(const double) const {
    return 1./(omega_max-omega_min);
  }

  double x(const double t) const {
    return t*blow_up();
  }

};



class Model
{
public:
  virtual double operator()(const double omega)=0;
  virtual ~Model(){}
};



class Gaussian : public Model 
{
public:
//  Gaussian(const alps::Parameters& p) : sigma(static_cast<double>(p["SIGMA"])) {}
  Gaussian(const alps::params& p) : sigma(static_cast<double>(p["SIGMA"])) {}
  
  virtual double operator()(const double omega) {
    return std::exp(-omega*omega/2./sigma/sigma)/sqrt(2*M_PI)/sigma;
  }

private:
  const double sigma;
};

class TwoGaussians : public Model
{
public:
//    TwoGaussians(const alps::Parameters& p) : sigma1(static_cast<double>(p["SIGMA1"])),
    TwoGaussians(const alps::params& p) : sigma1(static_cast<double>(p["SIGMA1"])),
    sigma2(static_cast<double>(p["SIGMA2"])),
    shift1(static_cast<double>(p["SHIFT1"]|0.0)),
    shift2(static_cast<double>(p["SHIFT2"])),
    norm1(static_cast<double>(p["NORM1"]|0.5)) {}
    
    virtual double operator()(const double omega) {
        return norm1*std::exp(-(omega-shift1)*(omega-shift1)/2./sigma1/sigma1)/sqrt(2*M_PI)/sigma1+(1.0-norm1)*std::exp(-(omega-shift2)*(omega-shift2)/2./sigma2/sigma2)/sqrt(2*M_PI)/sigma2;
    }
    
private:
    const double sigma1,sigma2,shift1,shift2,norm1;
};


class ShiftedGaussian : public Gaussian
{
public:
//  ShiftedGaussian(const alps::Parameters& p) :
  ShiftedGaussian(const alps::params& p) :
    Gaussian(p), shift(static_cast<double>(p["SHIFT"])){}

  double operator()(const double omega) {
    return Gaussian::operator()(omega-shift);
  }
  
protected:
  const double shift;
};


class DoubleGaussian : public ShiftedGaussian
{
public:
//  DoubleGaussian(const alps::Parameters& p) :
  DoubleGaussian(const alps::params& p) :
    ShiftedGaussian(p){}

  double operator()(const double omega) {
    return 0.5*(Gaussian::operator()(omega-shift) + Gaussian::operator()(omega+shift));
  }
};

class LinearRiseExpDecay : public Model{
public:
//  LinearRiseExpDecay(const alps::Parameters &p): lambda_(p["LAMBDA"]){}
  LinearRiseExpDecay(const alps::params &p): lambda_(p["LAMBDA"]){}
  double operator()(const double omega) {
    return lambda_*lambda_*omega*std::exp(-lambda_*omega);
  }

private:
  const double lambda_;
};

class QuadraticRiseExpDecay : public Model{
public:
//  QuadraticRiseExpDecay(const alps::Parameters &p): lambda_(p["LAMBDA"]){}
  QuadraticRiseExpDecay(const alps::params &p): lambda_(p["LAMBDA"]){}
  double operator()(const double omega) {
    return (lambda_*lambda_*lambda_)/2.*(omega*omega)*std::exp(-lambda_*omega);
  }
  
private:
  const double lambda_;
};

class GeneralDoubleGaussian : public ShiftedGaussian
{
public:
//  GeneralDoubleGaussian(const alps::Parameters& p) :
  GeneralDoubleGaussian(const alps::params& p) :
    ShiftedGaussian(p), bnorm(static_cast<double>(p["BOSE_NORM"])) {}
  
  double operator()(const double omega) {
    if (omega > 0)
      return Gaussian::operator()(omega);
    else
      return bnorm*Gaussian::operator()(omega+shift);
  }
  
private:
  const double bnorm; 
};


class TabFunction : public Model
{
public:
//  TabFunction(const alps::Parameters& p, std::string const& name) //: index(0)
  TabFunction(const alps::params& p, std::string const& name) //: index(0)
  {
    std::string p_name = p[name].cast<std::string>();
    std::ifstream defstream(p_name.c_str());
    if (!defstream)
      boost::throw_exception(std::invalid_argument("could not open default model file: "+p[name]));
    double om, D;
    while (defstream >> om >> D) {
      Omega.push_back(om);
      Def.push_back(D);
      defstream.ignore(1000,'\n'); // Anything beyond is considered as junk
    }
    double omega_max = p["OMEGA_MAX"]; 
    double omega_min(static_cast<double>(p["OMEGA_MIN"]|-omega_max)); //we had a 0 here in the bosonic case. That's not a good idea if you're continuing symmetric functions like chi(omega)/omega. Change omega_min to zero manually if you need it.
    //double omega_min = (p["KERNEL"] == "bosonic") ? 0. : 
    //     static_cast<double>(p.value_or_default("OMEGA_MIN", -omega_max));
    if (Omega[0]!=omega_min || Omega[Omega.size()-1]!=omega_max){
      std::cout<<"Omega[ 0] "<<Omega[0]<<" omega min: "<<omega_min<<std::endl;
      std::cout<<"Omega[-1] "<<Omega[Omega.size()-1]<<" omega max: "<<omega_max<<std::endl;
//      boost::throw_exception(std::invalid_argument("invalid omega range for default model"));
    }
  }
  
  double operator()(const double omega) {
    std::vector<double>::const_iterator ub = std::upper_bound(Omega.begin(), Omega.end(), omega);
    int index = ub - Omega.begin();
    if (ub==Omega.end())
      index = Omega.end()-Omega.begin()-1;
    double om1 = Omega[index-1];
    double om2 = Omega[index];
    double D1 = Def[index-1];
    double D2 = Def[index];
    return -(D2-D1)/(om2-om1)*(om2-omega)+D2;      
  }
   
private:
  std::vector<double> Omega;
  std::vector<double> Def;
};



class GeneralDefaultModel : public DefaultModel
{
public:
  
//  GeneralDefaultModel(const alps::Parameters& p, boost::shared_ptr<Model> mod)
  GeneralDefaultModel(const alps::params& p, boost::shared_ptr<Model> mod)
   : DefaultModel(p)
   , Mod(mod)
   , ntab(5001)
   , xtab(ntab) 
  {
    double sum = 0;
    xtab[0] = 0.;
    //this is an evaluation on an equidistant grid; sum integrated by trapezoidal rule
    double delta_omega = (omega_max-omega_min)/(ntab-1);
    for (int o=1; o<ntab; ++o) {
      double omega1 = omega_min + (o-1)*delta_omega;
      double omega2 = omega_min + o*delta_omega;
      sum += ((*Mod)(omega1)+(*Mod)(omega2))/2.*delta_omega;
      xtab[o] = sum;
    }
    for (int o=0; o<ntab; ++o) {
      xtab[o] *= blow_up()/sum;
    }
    /*std::cout<<std::setprecision(14)<<"tabulated a model, these are the values: "<<std::endl;
    for(int i=0;i<ntab;++i){
      std::cout<<i<<" "<<xtab[i]<<std::endl;
    }
    std::cout<<"total sum is: "<<sum<<std::endl;*/
    
    /*for(int N=0;N<=1000;++N){
      double d=1./1000.*N;
      std::cout<<omega(d)<<" "<<D(omega(d))<<" "<<x(d)<<" "<<d<<std::endl;
    }*/
  }
  
  double omega(const double x) const {
    if(!(x<=blow_up() && x>=0.)) throw std::logic_error("parameter x is out of bounds!"); //DNDEBUG switches off debug assertions
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
  
  //this returns the value of the model function at frequency omega
  double D(const double omega) const {
    return (*Mod)(omega);
  }
  
  //I have no idea what this does.
  double x(const double t) const {
    if(t>1. || t<0.) throw std::logic_error("parameter t is out of bounds!");
    int od = (int)(t*(ntab-1));
    if (od==(ntab-1)) 
      return blow_up();
    double x1 = xtab[od];
    double x2 = xtab[od+1];
    return -(x2-x1)*(od+1-t*ntab)+x2;      
  }
  
private:
  boost::shared_ptr<Model> Mod;
  const int ntab;
  std::vector<double> xtab; //xtab has an equidistantly tabulated discretized model function
};



//inline boost::shared_ptr<DefaultModel> make_default_model(const alps::Parameters& parms, std::string const& name)
inline boost::shared_ptr<DefaultModel> make_default_model(const alps::params& parms, std::string const& name)
{
    std::string p_name = parms[name]|"flat";
  if (p_name == "flat") {
    if (alps::is_master())
      std::cerr << "Using flat default model" << std::endl;
    return boost::shared_ptr<DefaultModel>(new FlatDefaultModel(parms));
  }
  else if (p_name == "gaussian") {
    if (alps::is_master())
      std::cerr << "Using Gaussian default model" << std::endl;
    boost::shared_ptr<Model> Mod(new Gaussian(parms));
    return boost::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
  else if (p_name == "twogaussians") {
      if (alps::is_master())
          std::cerr << "Using sum of two Gaussians default model" << std::endl;
      boost::shared_ptr<Model> Mod(new TwoGaussians(parms));
      return boost::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
  else if (p_name == "shifted gaussian") {
    if (alps::is_master())
      std::cerr << "Using shifted Gaussian default model" << std::endl;
    boost::shared_ptr<Model> Mod(new ShiftedGaussian(parms));
    return boost::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
  else if (p_name == "double gaussian") {
    if (alps::is_master())
      std::cerr << "Using double Gaussian default model" << std::endl;
    boost::shared_ptr<Model> Mod(new DoubleGaussian(parms));
    return boost::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
  else if (p_name == "general double gaussian") {
    if (alps::is_master())
      std::cerr << "Using general double Gaussian default model" << std::endl;
    boost::shared_ptr<Model> Mod(new GeneralDoubleGaussian(parms));
    return boost::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
  else if (p_name == "linear rise exp decay") {
    if (alps::is_master())
      std::cerr << "Using linear rise exponential decay default model" << std::endl;
    boost::shared_ptr<Model> Mod(new LinearRiseExpDecay(parms));
    return boost::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
  else if (p_name == "quadratic rise exp decay") {
    if (alps::is_master())
      std::cerr << "Using quadratic rise exponential decay default model" << std::endl;
    boost::shared_ptr<Model> Mod(new QuadraticRiseExpDecay(parms));
    return boost::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
  else { 
    if (alps::is_master())
      std::cerr << "Using tabulated default model" << std::endl;
    boost::shared_ptr<Model> Mod(new TabFunction(parms, name));
    return boost::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
}



#endif
