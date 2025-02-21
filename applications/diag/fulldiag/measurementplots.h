/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1994-2009 by Matthias Troyer <troyer@comp-phys.org>,
*                            Andreas Honecker <ahoneck@uni-goettingen.de>
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

/* $Id: fulldiag.h 2686 2008-01-08 09:32:34Z honecker $ */

#include <alps/plot.h>
#include <alps/scheduler/measurement_operators.h>

template <class T>
class MeasurementsPlot : public alps::MeasurementLabels
{
public:
  MeasurementsPlot(MeasurementLabels const& labels, std::string const& xn)
   : MeasurementLabels(labels), xname(xn)
  {}
  
  void reset(T xval); 
  void normalize(T z);
  
  template <class ValueType>
  void add (T w, alps::EigenvectorMeasurements<ValueType> const& meas, int j);
  
  void write(std::string const& basename) const;
private:
  void write_one(std::string const& basename, std::string const& name,
                  alps::plot::Set<T> const& sets) const;
             
  void write_one(std::string const& basename, std::string const& name,
                 std::vector<alps::plot::Set<T> > const& sets,
                 std::vector<std::string> const& labels) const;
  
  void write_plot(std::string const& name, alps::plot::Plot<T> const& plot) const;

  std::string xname;
  
  T xvalue;

  std::map<std::string,T> average_sum;
  std::map<std::string,std::vector<T > > local_sum;
  std::map<std::string,std::vector<T > > correlation_sum;
  std::map<std::string,std::vector<T > > structurefactor_sum;
    
  std::map<std::string,alps::plot::Set<T> > average_sets;
  std::map<std::string,std::vector<alps::plot::Set<T> > > local_sets;
  std::map<std::string,std::vector<alps::plot::Set<T> > > correlation_sets;
  std::map<std::string,std::vector<alps::plot::Set<T> > > structurefactor_sets;
  
};

template <class T>
void MeasurementsPlot<T>::reset(T xval)
{
  xvalue = xval;
  
  average_sum.clear();
  local_sum.clear();
  correlation_sum.clear();
  structurefactor_sum.clear();
    
}

template <class T>
template <class ValueType>
void MeasurementsPlot<T>::add (T w, alps::EigenvectorMeasurements<ValueType> const& meas, int j)
{
  typedef std::pair<std::string,std::vector<ValueType> > element_type;
  typedef std::pair<std::string,std::vector<std::vector<ValueType> >  > vector_element_type;
  using alps::numeric::real;
  
  BOOST_FOREACH(element_type const& s, meas.average_values)
    average_sum[s.first] += w * real(s.second[j]);

  BOOST_FOREACH(vector_element_type const& s, meas.local_values) 
  {
    if (local_sum[s.first].size() < s.second[j].size())
      local_sum[s.first].resize(s.second[j].size());
    for (int i=0; i < s.second[j].size(); ++i)
      local_sum[s.first][i] += w * real(s.second[j][i]);
  }

  BOOST_FOREACH(vector_element_type const& s, meas.correlation_values)
  {
    if (correlation_sum[s.first].size() < s.second[j].size())
     correlation_sum[s.first].resize(s.second[j].size());
    for (int i=0; i < s.second[j].size(); ++i)
      correlation_sum[s.first][i] += w * real(s.second[j][i]);
  }

  BOOST_FOREACH(vector_element_type const& s, meas.structurefactor_values)
  {
    if (structurefactor_sum[s.first].size() < s.second[j].size())
      structurefactor_sum[s.first].resize(s.second[j].size());
    for (int i=0; i < s.second[j].size(); ++i)
      structurefactor_sum[s.first][i] += w * real(s.second[j][i]);
  }
}


template <class T>
void MeasurementsPlot<T>::normalize(T z)
{
  typedef std::pair<std::string,T> element_type;
  typedef std::pair<std::string,std::vector<T > >  vector_element_type;

  BOOST_FOREACH(element_type const& s, average_sum)
    average_sets[s.first] << boost::make_tuple(xvalue,s.second/z);

  BOOST_FOREACH(vector_element_type const& s, local_sum) 
  {
    if (local_sets[s.first].size() < s.second.size())
      local_sets[s.first].resize(s.second.size());
    for (int i=0; i < s.second.size(); ++i)
      local_sets[s.first][i] << boost::make_tuple(xvalue,s.second[i]/z);
  }

  BOOST_FOREACH(vector_element_type const& s, correlation_sum)
  {
    if (correlation_sets[s.first].size() < s.second.size())
      correlation_sets[s.first].resize(s.second.size());
    for (int i=0; i < s.second.size(); ++i)
      correlation_sets[s.first][i] << boost::make_tuple(xvalue,s.second[i]/z);
  }

  BOOST_FOREACH(vector_element_type const& s, structurefactor_sum)
  {
    if (structurefactor_sets[s.first].size() < s.second.size())
      structurefactor_sets[s.first].resize(s.second.size());
    for (int i=0; i < s.second.size(); ++i)
      structurefactor_sets[s.first][i] << boost::make_tuple(xvalue,s.second[i]/z);
  }
    
}



template <class T>
void MeasurementsPlot<T>::write(std::string const& basename) const
{
  typedef std::pair<std::string,alps::plot::Set<T> > element_type;

  BOOST_FOREACH(element_type const& s, average_sets)
    write_one(basename, s.first, s.second);

  typedef std::pair<std::string, std::vector<alps::plot::Set<T> > > vector_element_type;

  BOOST_FOREACH(vector_element_type const& s, local_sets)
    if (bond_operator_[s.first])
      write_one(basename, s.first, s.second,bondlabel_);
    else
      write_one(basename, s.first, s.second,sitelabel_);

  BOOST_FOREACH(vector_element_type const& s, correlation_sets)
    write_one(basename, s.first, s.second,distlabel_);

  BOOST_FOREACH(vector_element_type const& s, structurefactor_sets)
    write_one(basename, s.first, s.second,momentumlabel_);
}


template <class T>
void MeasurementsPlot<T>::write_one(std::string const& basename, 
           std::string const& name, alps::plot::Set<T> const& s) const
{
  alps::plot::Plot<T> plot;
  plot.set_labels(xname,name); 
  alps::plot::Set<T> theset(s);
  theset << name; // add legend
  plot << theset; // add set to plot
  write_plot(basename+".measurements."+name,plot);
}

template <class T>
void MeasurementsPlot<T>::write_one(std::string const& basename, 
           std::string const& name, 
           std::vector<alps::plot::Set<T> > const& s,
           std::vector<std::string> const& labels) const
{
  alps::plot::Plot<T> plot;
  plot.set_labels(xname,name); 
  for (int i=0;i<s.size();++i) {
    alps::plot::Set<T> theset(s[i]);
    if (i<labels.size())
      theset << labels[i]; // add legend text
    plot << theset; // add set to plot
  }
  write_plot(basename+".measurements."+name,plot);
}


template <class T>
void MeasurementsPlot<T>::write_plot(std::string const& name, alps::plot::Plot<T> const& plot) const
{
  boost::filesystem::path p(name+".plot.xml");
  alps::oxstream out(p);
  out << plot;
}


