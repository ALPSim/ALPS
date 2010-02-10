/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2004 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Beat Ammon <ammon@ginnan.issp.u-tokyo.ac.jp>,
*                            Andreas Laeuchli <laeuchli@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
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

/* $Id$ */

#ifndef ALPS_ALEA_ABSTRACTBINNING_H
#define ALPS_ALEA_ABSTRACTBINNING_H

#include <alps/osiris/dump.h>
#include <alps/alea/obsvalue.h>

namespace alps {

//=======================================================================
// AbstractBinning
//-----------------------------------------------------------------------

template <class T>
class AbstractBinning {
 public: 
  typedef T value_type;
  typedef typename obs_value_traits<T>::time_type time_type;
  typedef typename obs_value_traits<T>::convergence_type convergence_type;
  AbstractBinning(std::size_t=0) {}

  time_type tau()                  const { boost::throw_exception(std::logic_error("Called non-implemented function of AbstractBinning")); return time_type(); }
  uint32_t max_bin_number()        const { return 0; }
  uint32_t bin_number()            const { return 0; }
  uint32_t filled_bin_number()     const { return 0; }
  uint32_t filled_bin_number2()     const { return 0; }
  uint32_t bin_size()              const { return 0; }
  const value_type& bin_value(uint32_t  ) const {
    boost::throw_exception(std::logic_error("Binning is not supported for this observable"));
    return *(new value_type); // dummy return
  } 
  const value_type& bin_value2(uint32_t  ) const {
    boost::throw_exception(std::logic_error("Binning is not supported for this observable"));
    return *(new value_type); // dummy return
  } 
#ifndef ALPS_WITHOUT_OSIRIS
  void extract_timeseries(ODump& dump) const { dump << 0 << 0 << 0;}
#endif
#ifdef ALPS_HAVE_HDF5
  template<typename E> void read_hdf5 (const E &engine);
  template<typename E> void write_hdf5 (const E &engine)const;
#endif

  bool has_variance() const { return true;} // for now
  void write_scalar_xml(oxstream&) const {}
  template <class IT> void write_vector_xml(oxstream&, IT) const {}

#ifndef ALPS_WITHOUT_OSIRIS
  void save(ODump& /* dump */) const {}
  void load(IDump& dump) 
  { 
    bool thermalized_; 
    if (dump.version() < 306 && dump.version() != 0) 
      dump >> thermalized_;
  }
#endif

#ifdef ALPS_HAVE_HDF5
  void serialize(hdf5::oarchive &, bool = false) const {};
  void serialize(hdf5::iarchive &, bool = false) {};
#endif

  std::string evaluation_method() const { return "simple";}
};

} // end namespace alps

#endif // ALPS_ALEA_ABSTRACTBINNING_H
