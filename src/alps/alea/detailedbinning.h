/***************************************************************************
* ALPS++/alea library
*
* alps/alea/detailedbinning.h     Monte Carlo observable class
*
* $Id$
*
* Copyright (C) 1994-2003 by Matthias Troyer <troyer@comp-phys.org>,
*                            Beat Ammon <ammon@ginnan.issp.u-tokyo.ac.jp>,
*                            Andreas Laeuchli <laeuchli@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>,
*
* This program is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
**************************************************************************/

#ifndef ALPS_ALEA_DETAILEDBINNING_H
#define ALPS_ALEA_DETAILEDBINNING_H

#include <alps/config.h>
#include <alps/alea/observable.h>
#include <alps/alea/simplebinning.h>

#ifdef ALPS_HAVE_VALARRAY
# include <valarray>
#endif

//=======================================================================
// DetailedBinning
//
// detailed binning strategy
//-----------------------------------------------------------------------

namespace alps{

template <class T=double>
class BasicDetailedBinning : public SimpleBinning<T> {
public:
  typedef T value_type;  
  typedef typename obs_value_traits<T>::time_type time_type;
  typedef typename obs_value_traits<T>::size_type size_type;
  typedef typename obs_value_traits<T>::count_type count_type;
  typedef typename obs_value_traits<T>::result_type result_type;

  static const bool has_tau=true;
  static const int magic_id=3;

  BasicDetailedBinning(uint32_t binsize=1, uint32_t binnum=std::numeric_limits<uint32_t>::max());

  void reset(bool=false);
  void operator<<(const T& x);

  
  uint32_t max_bin_number() const { return maxbinnum_;}
  uint32_t bin_number() const;
  uint32_t filled_bin_number() const;
  uint32_t filled_bin_number2() const 
  { return(values2_.size() ? filled_bin_number() : 0);}
  
  void set_bin_number(uint32_t binnum);
  void collect_bins(uint32_t howmany);
  
  uint32_t bin_size() const { return binsize_;}
  void set_bin_size(uint32_t binsize);

  const value_type& bin_value(uint32_t i) const { return values_[i];}
  const value_type& bin_value2(uint32_t i) const { return values2_[i];}
  
  void compact();
  
#ifndef ALPS_WITHOUT_OSIRIS
  virtual void save(ODump& dump) const;
  virtual void load(IDump& dump);
  void extract_timeseries(ODump& dump) const { dump << binsize_ << values_.size() << binentries_ << values_;}
#endif

private:
  uint32_t binsize_;       // number of measurements per bin
  uint32_t minbinsize_;    // minimum number of measurements per bin
  uint32_t maxbinnum_;      // maximum number of bins 
  uint32_t  binentries_; // number of measurements in last bin
  std::vector<value_type> values_; // bin values
  std::vector<value_type> values2_; // bin values of squares
};

template<class T> class DetailedBinning : public BasicDetailedBinning<T>
{
public:
  typedef T value_type;
  static const int magic_id=4;
  DetailedBinning(uint32_t binnum=128) 
  : BasicDetailedBinning<T>(1,binnum==0 ? 128 : binnum) {}
};

template<class T> class FixedBinning : public BasicDetailedBinning<T>
{
public:
  typedef T value_type;
  static const int magic_id=5;
  FixedBinning(uint32_t binsize=1) 
  : BasicDetailedBinning<T>(binsize,std::numeric_limits<uint32_t>::max()) {}
};

typedef BasicSimpleObservable<int32_t,DetailedBinning<int32_t> > IntObservable;
typedef BasicSimpleObservable<double,DetailedBinning<double> > RealObservable;
typedef BasicSimpleObservable<std::complex<double>,DetailedBinning<std::complex<double> > > ComplexObservable;
typedef BasicSimpleObservable<double,FixedBinning<double> > RealTimeSeriesObservable;
typedef BasicSimpleObservable<int32_t,FixedBinning<int32_t> > IntTimeSeriesObservable;

#ifdef ALPS_HAVE_VALARRAY
typedef BasicSimpleObservable< std::valarray<int32_t> , 
                         DetailedBinning<std::valarray<int32_t> > > IntVectorObservable;
typedef BasicSimpleObservable< std::valarray<double> , 
                         DetailedBinning<std::valarray<double> > > RealVectorObservable;
typedef BasicSimpleObservable< std::valarray<std::complex<double> > , 
                         DetailedBinning<std::valarray<std::complex<double> > > > ComplexVectorObservable;
typedef BasicSimpleObservable< std::valarray<int32_t> , 
                         FixedBinning<std::valarray<int32_t> > > IntVectorTimeSeriesObservable;
typedef BasicSimpleObservable< std::valarray<double> , 
                         FixedBinning<std::valarray<double> > > RealVectorTimeSeriesObservable;
typedef BasicSimpleObservable< std::valarray<std::complex<double> > , 
                         FixedBinning<std::valarray<std::complex<double> > > > ComplexVectorTimeSeriesObservable;
#endif
typedef BasicSimpleObservable< alps::multi_array<int32_t,2> , DetailedBinning<alps::multi_array<int32_t,2> > > Int2DArrayObservable;
typedef BasicSimpleObservable< alps::multi_array<double,2> , DetailedBinning<alps::multi_array<double,2> > > Real2DArrayObservable;
typedef BasicSimpleObservable< alps::multi_array<std::complex<double>,2> , DetailedBinning<alps::multi_array<std::complex<double>,2> > > Complex2DArrayObservable;
typedef BasicSimpleObservable< alps::multi_array<int32_t,2> , FixedBinning<alps::multi_array<int32_t,2> > > Int2DArrayTimeSeriesObservable;
typedef BasicSimpleObservable< alps::multi_array<double,2> , FixedBinning<alps::multi_array<double,2> > > Real2DArrayTimeSeriesObservable;
typedef BasicSimpleObservable< alps::multi_array<std::complex<double>,2> , FixedBinning<alps::multi_array<std::complex<double>,2> > > Complex2DArrayTimeSeriesObservable;


template <class T>
inline BasicDetailedBinning<T>::BasicDetailedBinning(uint32_t binsize, uint32_t binnum)
 : SimpleBinning<T>(),
   binsize_(0), minbinsize_(binsize), maxbinnum_(binnum), binentries_(0)    
{
  reset();
}

template <class T>
inline void BasicDetailedBinning<T>::compact()
{
  std::vector<value_type> tmp1;
  std::vector<value_type> tmp2;
  values_.swap(tmp1);
  values2_.swap(tmp2);
  binsize_=minbinsize_;
  binentries_=0;
}


template <class T>
inline void BasicDetailedBinning<T>::reset(bool forthermal)
{
  compact();
  SimpleBinning<T>::reset(forthermal);
}


template <class T>
inline void BasicDetailedBinning<T>::operator<<(const T& x)
{
  if (values_.empty())
  { 
    // start first bin
    values_.push_back(x);
    values2_.push_back(x*x);
    binentries_ = 1;
    binsize_=1;
  }
  else if (values_.size()==1 && binentries_ < minbinsize_)
  {
    // fill first bin
    values_[0]+=x;
    values2_[0]+=x*x;
    binentries_++;
    binsize_++;
  }
  else if (binentries_==binsize_) // have a full bin
  {
    if(values_.size()<maxbinnum_)
    {
      // start a new bin
      values_.push_back(x);
      values2_.push_back(x*x);
      binentries_ = 1;
    }
    else
    {
      // halve the bins and add
      collect_bins(2);
      *this << x; // and just call again
      return;
    }
  }
  else
  {
    values_[values_.size()-1] += x;
    values2_[values_.size()-1] += x*x;
    ++binentries_;
  }
  SimpleBinning<T>::operator<<(x);
}

template <class T>
void BasicDetailedBinning<T>::collect_bins(uint32_t howmany)
{
  if (values_.empty() || howmany<=1)
    return;
    
  uint32_t newbins = (values_.size()+howmany-1)/howmany;
  
  // full bins
  for (uint32_t i=0;i<values_.size()/howmany;++i)
  {
    values_[i]=values_[howmany*i];
    values2_[i]=values2_[howmany*i];
    for (uint32_t j = 1 ; j<howmany;++j)
    {
      values_[i] += values_[howmany*i+j];
      values2_[i] += values2_[howmany*i+j];
    }
  }
  
  // last part of partly full bins
  values_[newbins-1]=values_[howmany*(newbins-1)];
  values2_[newbins-1]=values2_[howmany*(newbins-1)];
  for ( uint32_t i=howmany*(newbins-1)+1;i<values_.size();++i){
    values_[newbins-1]+=values_[i];
    values2_[newbins-1]+=values2_[i];
  }
    
  // how many in last bin?
  binentries_+=((values_.size()-1)%howmany)*binsize_;
  binsize_*=howmany;

  values_.resize(newbins);
  values2_.resize(newbins);
}

template <class T>
void BasicDetailedBinning<T>::set_bin_size(uint32_t minbinsize)
{
  minbinsize_=minbinsize;
  if(binsize_ < minbinsize_ && binsize_ > 0)
    collect_bins((minbinsize-1)/binsize_+1);
}

template <class T>
void BasicDetailedBinning<T>::set_bin_number(uint32_t binnum)
{
  maxbinnum_=binnum;
  if(values_.size() > maxbinnum_)
    collect_bins((values_.size()-1)/maxbinnum_+1);
}

template <class T>
inline uint32_t BasicDetailedBinning<T>::bin_number() const 
{ 
  return values_.size();
}

template <class T>
inline uint32_t BasicDetailedBinning<T>::filled_bin_number() const 
{ 
  if(values_.size()==0) return 0;
  else return values_.size() + (binentries_ ==binsize_ ? 0 : -1);
}

#ifndef ALPS_WITHOUT_OSIRIS
template <class T>
inline void BasicDetailedBinning<T>::save(ODump& dump) const
{
  SimpleBinning<T>::save(dump);
  dump << binsize_ << minbinsize_ << maxbinnum_<< binentries_ << values_
       << values2_;
}

template <class T>
inline void BasicDetailedBinning<T>::load(IDump& dump) 
{
  SimpleBinning<T>::load(dump);
  dump >> binsize_ >> minbinsize_ >> maxbinnum_ >> binentries_ >> values_
       >> values2_;
}
#endif

} // end namespace alps

#endif // ALPS_ALEA_DETAILEDBINNING_H
