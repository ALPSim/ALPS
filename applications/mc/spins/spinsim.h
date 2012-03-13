/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1999-2003 by Matthias Troyer <troyer@comp-phys.org>,
*                            Fabian Stoeckli <fabstoec@student.ethz.ch>
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

/* $Id$ */

#ifndef ALPS_APPLICATIONS_MC_SPIN_SPINS_H_
#define ALPS_APPLICATIONS_MC_SPIN_SPINS_H_

#include "clusterupdate.h"
#include "localupdate.h"
#include "abstractspinsim.h"
#include "ising.h"
#include "on.h"

#include "potts.h"
#include "tinyvec.h"
#include "matrices.h"
#include "helper.h"

#include <boost/property_map/vector_property_map.hpp>

namespace boost {

template <class T> struct property_traits<std::vector<T> > {
  typedef T value_type;
};

}

template <class M>
struct measurement_traits {
   static double component_number() { return 1.;}
   static const bool use_improved_estimator = false;
};

template <class M>
const bool measurement_traits<M>::use_improved_estimator;

template <>
struct measurement_traits<IsingMoment> {
   static double component_number() { return 1.;}
   static const bool use_improved_estimator = true;
};

const bool measurement_traits<IsingMoment>::use_improved_estimator;

template <int N>
struct measurement_traits<ONMoment<N> > {
   static double component_number() { return N;}
   static const bool use_improved_estimator = true;
};

template <int N>
const bool measurement_traits<ONMoment<N> >::use_improved_estimator;

template <class M, class MAT>
class SpinSim : public AbstractSpinSim<MAT> {
public:
  typedef M moment_type;
  SpinSim(const alps::ProcessList& w,const alps::Parameters& p,int n);

  void save(alps::ODump& dump) const;
  void load(alps::IDump& dump);
  update_info_type do_update();
  void do_measurements(update_info_type);

private:
  std::vector<moment_type> spins_;
};

/* special case for MIdMatrix (allows cluster updates) */
template <class M>
class SpinSim<M,MIdMatrix<double,M::dim> > 
        : public AbstractSpinSim<MIdMatrix<double,M::dim> > {
public:
  typedef M moment_type;
  SpinSim(const alps::ProcessList& w,const alps::Parameters& p,int n);

  void save(alps::ODump& dump) const;
  void load(alps::IDump& dump);
  update_info_type do_update();
  void do_measurements(update_info_type);
  

private:
  std::vector<moment_type> spins_;
};


template <class M, class MAT>
SpinSim<M, MAT>::SpinSim(const alps::ProcessList& w,const alps::Parameters& p,int n)
 : AbstractSpinSim<MAT>::AbstractSpinSim(w,p,n),
   spins_(this->num_sites())
{
}

template <class M>
SpinSim<M, MIdMatrix<double,M::dim> >::SpinSim(const alps::ProcessList& w,
                const alps::Parameters& p,int n)
 : AbstractSpinSim<MIdMatrix<double,M::dim> >::AbstractSpinSim(w,p,n),
   spins_(this->num_sites())
{
}


template <class M, class MAT>
void SpinSim<M, MAT>::save(alps::ODump& dump) const
{ 
  AbstractSpinSim<MAT>::save(dump);
  dump << spins_;
}

template <class M>
void SpinSim<M, MIdMatrix<double,M::dim> >::save(alps::ODump& dump) const
{ 
  AbstractSpinSim<MIdMatrix<double,M::dim> >::save(dump);
  dump << spins_;
}


template <class M, class MAT>
void SpinSim<M, MAT>::load(alps::IDump& dump) 
{ 
  AbstractSpinSim<MAT>::load(dump);
  dump >> spins_;
  if (AbstractSpinSim<MAT>::where.empty())
    spins_ = std::vector<moment_type>();
}

template <class M>
void SpinSim<M, MIdMatrix<double,M::dim> >::load(alps::IDump& dump) 
{ 
  AbstractSpinSim<MIdMatrix<double,M::dim> >::load(dump);
  dump >> spins_;
  if (AbstractSpinSim<MIdMatrix<double,M::dim> >::where.empty())
    spins_ = std::vector<moment_type>();
}


template <class M, class MAT> 
update_info_type SpinSim<M,MAT>::do_update()
{
    typedef AbstractSpinSim<MAT> parent;
    if (this->general_case_)
      local_update_self(this->graph(),spins_, this->beta_,this->random_01,
          this->couplings_,this->selfinteraction_,this->spinfactor_,this->g_,this->h_);
    else
      local_update(this->graph(),spins_, this->beta_,this->random_01,
          this->couplings_,this->spinfactor_,this->g_,this->h_);
    update_info_type up_info;
    up_info.clustersize = this->num_sites();
    up_info.m2_ = 0;
    return up_info;
}

template<class M>
update_info_type SpinSim<M,MIdMatrix<double,M::dim> >::do_update()
{
  typedef AbstractSpinSim<MIdMatrix<double,M::dim> > parent;
  if (this->cluster_updates_)
    if (this->general_case_)
      return cluster_update(this->graph(),spins_, this->beta_,
          this->random_01, this->couplings_, this->spinfactor_);
    else
      return cluster_update(this->graph(),spins_, this->beta_, 
          this->random_01,this->couplings_,this->spinfactor_);
 else {
    if (this->general_case_)
      local_update_self(this->graph(),spins_,this->beta_,this->random_01,
          this->couplings_,this->selfinteraction_,this->spinfactor_,this->g_,this->h_);
    else
      local_update(this->graph(),spins_,this->beta_,this->random_01,
          this->couplings_,this->spinfactor_,this->g_,this->h_);
    update_info_type up_info;
    up_info.clustersize = this->num_sites();
    up_info.m2_ = 0;
    return up_info;
  }
}
    
    
    
    
template <class M, class MAT> 
void SpinSim<M,MAT>::do_measurements(update_info_type update_info)
{
  double en=0.0;
  typedef AbstractSpinSim<MAT> parent;
  typename parent::bond_iterator bi;
  bi = (this->bonds()).first;
  std::vector<double> bt_energy;
  std::vector<int> bt_count;
  for (; bi!=(this->bonds()).second; ++bi) {
    int bt = this->bond_type(*bi);
    if (bt>=bt_energy.size()) {
      bt_energy.resize(bt+1);
      bt_count.resize(bt+1);
    }
    double be = this->spinfactor_[this->source(*bi)]*this->spinfactor_[this->target(*bi)]*bond_energy(spins_[this->source(*bi)],
            spins_[this->target(*bi)],this->couplings_[*bi]);
    en -= be;
    bt_energy[bt] -= be;
    bt_count[bt]++;
  }
  
  if (this->cluster_updates_)
    this->measurements.template get<alps::RealObservable>("Cluster size") 
        << double(update_info.clustersize)/double(this->num_sites());

  typename parent::site_iterator si;
  si=this->sites().first;
  if (si==this->sites().second)
    return;
  typename M::magnetization_type m= spins_[*si].magnetization();
  m *= this->spinfactor_[*si];
  typename M::magnetization_type mst= spins_[*si].magnetization();
  mst *= this->spinfactor_[*si];
  double m_h = 0.;
  if (this->general_case_)
      en+=onsite_energy(spins_[*si],this->selfinteraction_[*si]);
  if (this->has_magnetic_field_) 
    // use normed h_ - vector to allow faster computation
    m_h = this->spinfactor_[*si]*spins_[*si].mag_h(this->h_normalized);

  if (this->has_magnetic_field_) 
    en+=this->spinfactor_[*si]*site_energy(spins_[*si],this->h_)*this->g_;

  double par=this->parity(*si);
  ++si;
  
  for (; si !=this->sites().second ; ++si) {
    if (this->has_magnetic_field_) 
      en+=this->spinfactor_[*si]*site_energy(spins_[*si],this->h_)*this->g_;
      
    if (this->general_case_)
      en+=onsite_energy(spins_[*si],this->selfinteraction_[*si]);
    
    m+=this->spinfactor_[*si]*spins_[*si].magnetization();
    if (this->has_magnetic_field_) 
      m_h += this->spinfactor_[*si]*spins_[*si].mag_h(this->h_normalized);

    if (this->is_bipartite()) 
      mst+=par*this->parity(*si)*this->spinfactor_[*si]*spins_[*si].magnetization();
  }
  
  std::valarray<double> bt_en(bt_energy.size());
  std::valarray<double> bt_en_density(bt_energy.size());
  for (int i=0;i<bt_energy.size();++i) {
    bt_en[i] = bt_energy[i];
    bt_en_density[i] = bt_energy[i]/this->num_sites();
  }
  
  this->measurements.template get<alps::RealObservable>("Energy") << en;
  this->measurements.template get<alps::RealObservable>("Energy Density") << en/this->num_sites();
  this->measurements.template get<alps::RealVectorObservable>("Bond-type Energy") << bt_en;
  this->measurements.template get<alps::RealVectorObservable>("Bond-type Energy Density") << bt_en_density;
  this->measurements.template get<alps::RealObservable>("Energy^2") << en*en;
  this->measurements.template get<alps::RealObservable>("beta * Energy / sqrt(N)") << this->beta_ *en/sqrt(static_cast<double>(this->num_sites()));
  this->measurements.template get<alps::RealObservable>("(beta * Energy)^2 / N") << this->beta_*this->beta_*en*en/this->num_sites();;

  
  using std::abs;
  double am=abs(m)/this->num_sites();
  double mh=m_h/this->num_sites();
  double chi;
  double m2;
  double m4; 
  double enm2;
  if(this->cluster_updates_ && measurement_traits<M>::use_improved_estimator && this->ferromagnetic_) {
    chi = (update_info.m2_)*(update_info.m2_)*this->beta_/double(update_info.clustersize);
    m2 = (update_info.m2_)*(update_info.m2_)*measurement_traits<M>::component_number()/(update_info.clustersize*double(this->num_sites()));
  }
  else {
    m2=am*am;
    chi=m2*this->beta_*this->num_sites()/measurement_traits<M>::component_number();
  }
   m4=am*am*am*am; 
   enm2=en*am*am;
   
  this->measurements.template get<alps::RealObservable>("|Magnetization|") << am;
  this->measurements.template get<alps::RealObservable>("Magnetization along Field") << mh;
  this->measurements.template get<alps::RealObservable>("Magnetization^2") << m2;
  this->measurements.template get<alps::RealObservable>("E.Magnetization^2") << enm2;
  this->measurements.template get<alps::RealObservable>("Magnetization^4") << m4;
  this->measurements.template get<alps::RealObservable>("E.Magnetization^4") << en*m4;
  this->measurements.template get<alps::RealObservable>("Susceptibility") << chi;
  if (this->is_bipartite()) {
    double amst = abs(mst)/this->num_sites();
    double am2;
    if(this->cluster_updates_ && measurement_traits<M>::use_improved_estimator && this->antiferromagnetic_)
      am2 = (update_info.m2_)*(update_info.m2_)*measurement_traits<M>::component_number()/(update_info.clustersize*double(this->num_sites()));
    else
      am2 = amst*amst;
    this->measurements.template get<alps::RealObservable>("|Staggered Magnetization|") << amst;
    this->measurements.template get<alps::RealObservable>("Staggered Magnetization^2") << am2;
  } 
  
  if (this->print_sweeps_) {
    int good_sweeps = this->sweeps_done_-this->thermalization_sweeps_;
    if (good_sweeps > 0 && (good_sweeps % this->print_sweeps_ ==0)) {
      for (int i=0;i<this->num_sites();++i) {
        std::cout << i << " ";
        std::vector<double> coord = this->coordinate(i);
        for (int j=0;j<coord.size();++j)
          std::cout << coord[j] << " ";
        std::cout << spins_[i] << "\n";
      }
    }
  }
}


    
template <class M> 
void SpinSim<M,MIdMatrix<double,M::dim> >::do_measurements(update_info_type update_info)
{
  typedef AbstractSpinSim<MIdMatrix<double,M::dim> > parent;
  double en=0.0;
  typename parent::bond_iterator bi;
  bi = (this->bonds()).first;
  std::vector<double> bt_energy;
  std::vector<int> bt_count;
  for (; bi!=(this->bonds()).second; ++bi) {
   int bt = this->bond_type(*bi);
    if (bt>=bt_energy.size()) {
      bt_energy.resize(bt+1);
      bt_count.resize(bt+1);
    }
    double be = this->spinfactor_[this->source(*bi)]*this->spinfactor_[this->target(*bi)]*bond_energy(spins_[this->source(*bi)],spins_[this->target(*bi)],
        this->couplings_[*bi]);

    en -= be;
    bt_energy[bt] -= be;
    bt_count[bt]++;

  }
  
  if (this->cluster_updates_)
    this->measurements.template get<alps::RealObservable>("Cluster size") 
        << double(update_info.clustersize)/double(this->num_sites());

  typename parent::site_iterator si;
  
  si=this->sites().first;
  if (si==this->sites().second)
    return;
  typename M::magnetization_type m= spins_[*si].magnetization();
  m *= this->spinfactor_[*si];
  typename M::magnetization_type mst= spins_[*si].magnetization();
  mst *= this->spinfactor_[*si];
  double m_h = 0.;
  if (this->has_magnetic_field_) 
      {
        m_h = this->spinfactor_[*si]*spins_[*si].mag_h(this->h_normalized);
        en+=this->spinfactor_[*si]*site_energy(spins_[*si],this->h_)*this->g_;
      }
  double par=this->parity(*si);
  ++si;
  
  for (; si !=this->sites().second ; ++si) {
    if (this->has_magnetic_field_) 
      en+=this->spinfactor_[*si]*site_energy(spins_[*si],this->h_)*this->g_;
    if (this->general_case_)
      en+=onsite_energy(spins_[*si],this->selfinteraction_[*si]);
    
    m+=this->spinfactor_[*si]*spins_[*si].magnetization();
    if (this->has_magnetic_field_) 
      m_h += this->spinfactor_[*si]*spins_[*si].mag_h(this->h_normalized);
    if (this->is_bipartite()) 
      mst+=par*this->parity(*si)*this->spinfactor_[*si]*spins_[*si].magnetization();
  }

  std::valarray<double> bt_en(bt_energy.size());
  std::valarray<double> bt_en_density(bt_energy.size());
  for (int i=0;i<bt_energy.size();++i) {
    bt_en[i] = bt_energy[i];
    bt_en_density[i] = bt_energy[i]/this->num_sites();
  }
  
  this->measurements.template get<alps::RealObservable>("Energy") << en;
  this->measurements.template get<alps::RealObservable>("Energy Density") << en/this->num_sites();
  this->measurements.template get<alps::RealVectorObservable>("Bond-type Energy") << bt_en;
  this->measurements.template get<alps::RealVectorObservable>("Bond-type Energy Density") << bt_en_density;
  this->measurements.template get<alps::RealObservable>("Energy^2") << en*en;
  this->measurements.template get<alps::RealObservable>("beta * Energy / sqrt(N)") << this->beta_ *en/sqrt(static_cast<double>(this->num_sites()));
  this->measurements.template get<alps::RealObservable>("(beta * Energy)^2 / N") << this->beta_*this->beta_*en*en/this->num_sites();;

  
  using std::abs;
  double am=abs(m)/this->num_sites();
  double mh=m_h/this->num_sites();
  double chi;
  double m2;
  double m4; 
  double enm2;
  if(this->cluster_updates_ && 
          measurement_traits<M>::use_improved_estimator && 
        this->ferromagnetic_) {
    chi = (update_info.m2_)*(update_info.m2_)*this->beta_/double(update_info.clustersize);
    m2 = (update_info.m2_)*(update_info.m2_)*measurement_traits<M>::component_number()/(update_info.clustersize*double(this->num_sites()));
  }
  else {
    m2=am*am; 
    chi=m2*this->beta_*this->num_sites()/measurement_traits<M>::component_number();
  }
  m4=am*am*am*am;
  enm2=en*am*am;
  this->measurements.template get<alps::RealObservable>("|Magnetization|") << am;
  this->measurements.template get<alps::RealObservable>("Magnetization along Field") << mh;
  this->measurements.template get<alps::RealObservable>("Magnetization^2") << m2;
  this->measurements.template get<alps::RealObservable>("E.Magnetization^2") << enm2;
  this->measurements.template get<alps::RealObservable>("Magnetization^4") << m4;
  this->measurements.template get<alps::RealObservable>("E.Magnetization^4") << en*m4;
  this->measurements.template get<alps::RealObservable>("Susceptibility") << chi;
  if (this->is_bipartite()) {
    double amst = abs(mst)/this->num_sites();
    double am2;
    if(this->cluster_updates_ && 
            measurement_traits<M>::use_improved_estimator && 
        this->antiferromagnetic_)
      am2 = (update_info.m2_)*(update_info.m2_)*measurement_traits<M>::component_number()/(update_info.clustersize*double(this->num_sites()));
    else
      am2 = amst*amst;
    this->measurements.template get<alps::RealObservable>("|Staggered Magnetization|") << amst;
    this->measurements.template get<alps::RealObservable>("Staggered Magnetization^2") << am2;
  } 
  if (this->print_sweeps_) {
    int good_sweeps = this->sweeps_done_-this->thermalization_sweeps_;
    if (good_sweeps > 0 && (good_sweeps % this->print_sweeps_ ==0)) {
      std::cout << "Dumping configuration: \n";
      for (int i=0;i<this->num_sites();++i) {
        std::cout << i << " ";
        std::vector<double> coord = this->coordinate(i);
        for (int j=0;j<coord.size();++j)
          std::cout << coord[j] << " ";
        std::cout << spins_[i] << "\n";
      }
    }
  }
}

#endif
