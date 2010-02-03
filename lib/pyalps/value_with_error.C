/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2010 by Ping Nang Ma <pingnang@itp.phys.ethz.ch>,
*                            Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Bela Bauer <bauerb@itp.phys.ethz.ch>
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

/* $Id: nobinning.h 3520 2009-12-11 16:49:53Z gamperl $ */


#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <alps/alea/value_with_error.h>


using namespace boost::python;
using namespace alps::alea;


namespace alps { 
  namespace alea {

    // for pickling support
    template<class T>
    struct value_with_error_pickle_suite : boost::python::pickle_suite
    {
      static boost::python::tuple getinitargs(const value_with_error<T>& v)
      {   
        return boost::python::make_tuple(v.mean(),v.error());
      }   
    };


    // for printing support
    inline static boost::python::str print_value_with_error(value_with_error<double> const & self)
    {
      return boost::python::str(boost::python::str(self.mean()) + " +/- " + boost::python::str(self.error()));
    }

    static boost::python::str print_vector_with_error(value_with_error<std::vector<double> > self)
    {
      boost::python::str s;
      for (std::size_t index=0; index < self.size(); ++index)
      {
        s += boost::python::str(self.at(index));
        s += boost::python::str("\n");
      }
      return s;
    }

    template<class T>
    inline static boost::python::str print_vector_of_value_with_error(std::vector<value_with_error<T> > const & self)
    {
      typename std::vector<value_with_error<T> >::const_iterator it;

      boost::python::str s;
      for (it = self.begin(); it != self.end(); ++it)
      {
        s += print_value_with_error(*it);
        if (it != (self.end()-1))  {  s += '\n';  }
      }
      return s;
    }

    template <class T>
    static boost::python::str print_vector_list(std::vector<T> self)
    {
      boost::python::str s;
      for (typename std::vector<T>::iterator it = self.begin(); it != self.end(); ++it)
      {
        s += boost::python::str(*it);
        s += boost::python::str("\n");
      }
      return s;
    }

  }
}



BOOST_PYTHON_MODULE(pyalea)
{
  class_<value_with_error<double> >("value_with_error",init<optional<double,double> >())
    .add_property("mean", &value_with_error<double>::mean)
    .add_property("error",&value_with_error<double>::error)  

    .def("__repr__", &print_value_with_error)

    .def(+self)
    .def(-self)
    .def("__abs__",&abs<double>)
    .def(self == value_with_error<double>())

    .def(self += value_with_error<double>())
    .def(self += double())
    .def(self -= value_with_error<double>())
    .def(self -= double())
    .def(self *= value_with_error<double>())
    .def(self *= double())
    .def(self /= value_with_error<double>())
    .def(self /= double())

    .def(self + value_with_error<double>())
    .def(self + double())
    .def(double() + self)
    .def(self - value_with_error<double>())
    .def(self - double())
    .def(double() - self)
    .def(self * value_with_error<double>())
    .def(self * double())
    .def(double() * self)
    .def(self / value_with_error<double>())
    .def(self / double())
    .def(double() / self)


    .def("__pow__",&pow<double>)
    .def("sq",&sq<double>)
    .def("cb",&cb<double>)
    .def("sqrt",&sqrt<double>)
    .def("cbrt",&cbrt<double>)
    .def("exp",&exp<double>)
    .def("log",&log<double>)

    .def("sin",&sin<double>)
    .def("cos",&cos<double>)
    .def("tan",&tan<double>)
    .def("asin",&asin<double>)
    .def("acos",&acos<double>)
    .def("atan",&atan<double>)
    .def("sinh",&sinh<double>)
    .def("cosh",&cosh<double>)
    .def("tanh",&tanh<double>)
    .def("asinh",&asinh<double>)
    .def("acosh",&acosh<double>)
    .def("atanh",&atanh<double>)

    .def_pickle(value_with_error_pickle_suite<double>())
    ;

  class_<value_with_error<std::vector<double> > >("vector_with_error",init<optional<std::vector<double>,std::vector<double> > >())
    .add_property("mean",&value_with_error<std::vector<double> >::mean)
    .add_property("error",&value_with_error<std::vector<double> >::error)

    .def("__repr__", &print_vector_with_error)

    .def("__len__",&value_with_error<std::vector<double> >::size)         
    .def("append",&value_with_error<std::vector<double> >::push_back)     
    .def("push_back",&value_with_error<std::vector<double> >::push_back)  
    .def("insert",&value_with_error<std::vector<double> >::insert)    
    .def("pop_back",&value_with_error<std::vector<double> >::pop_back)   
    .def("erase",&value_with_error<std::vector<double> >::erase)        
    .def("clear",&value_with_error<std::vector<double> >::clear)         
    .def("at",&value_with_error<std::vector<double> >::at)

    .def("obtained_from",&obtain_vector_with_error_from_vector_of_value_with_error<double>)

    .def(self + value_with_error<std::vector<double> >())
    .def(self + double())
    .def(double() + self)
    .def(self + std::vector<double>())
    .def(std::vector<double>() + self)
    .def(self - value_with_error<std::vector<double> >())
    .def(self - double())
    .def(double() - self)
    .def(self - std::vector<double>())
    .def(std::vector<double>() - self)
    .def(self * value_with_error<std::vector<double> >())
    .def(self * double())
    .def(double() * self)
    .def(self * std::vector<double>())
    .def(std::vector<double>() * self)
    .def(self / value_with_error<std::vector<double> >())
    .def(self / double())
    .def(double() / self)
    .def(self / std::vector<double>())
    .def(std::vector<double>() / self)

    .def(+self)
    .def(-self)
    .def("__abs__",&abs<std::vector<double> >)

    .def("__pow__",&pow<std::vector<double> >)
    .def("sq",&sq<std::vector<double> >)
    .def("cb",&cb<std::vector<double> >)
    .def("sqrt",&sqrt<std::vector<double> >)
    .def("cbrt",&cbrt<std::vector<double> >)
    .def("exp",&exp<std::vector<double> >)
    .def("log",&log<std::vector<double> >)

    .def("sin",&sin<std::vector<double> >)
    .def("cos",&cos<std::vector<double> >)
    .def("tan",&tan<std::vector<double> >)
    .def("asin",&asin<std::vector<double> >)
    .def("acos",&acos<std::vector<double> >)
    .def("atan",&atan<std::vector<double> >)
    .def("sinh",&sinh<std::vector<double> >)
    .def("cosh",&cosh<std::vector<double> >)
    .def("tanh",&tanh<std::vector<double> >)
    .def("asinh",&asinh<std::vector<double> >)
    .def("acosh",&acosh<std::vector<double> >)
    .def("atanh",&atanh<std::vector<double> >)
    ;


  class_<std::vector<value_with_error<double> > >("vector_of_value_with_error")
    .def(vector_indexing_suite<std::vector<value_with_error<double> > >())

    .def("__repr__", &print_vector_of_value_with_error<double>)

    .def("obtained_from",&obtain_vector_of_value_with_error_from_vector_with_error<double>)

    .def(self + std::vector<value_with_error<double> >())
    .def(self + double())
    .def(double() + self)
    .def(self + std::vector<double>())
    .def(std::vector<double>() + self)
    .def(self - std::vector<value_with_error<double> >())
    .def(self - double())
    .def(double() - self)
    .def(self - std::vector<double>())
    .def(std::vector<double>() - self)
    .def(self * std::vector<value_with_error<double> >())
    .def(self * double())
    .def(double() * self)
    .def(self * std::vector<double>())
    .def(std::vector<double>() * self)
    .def(self / std::vector<value_with_error<double> >())
    .def(self / double())
    .def(double() / self)
    .def(self / std::vector<double>())
    .def(std::vector<double>() / self)

    .def(-self)

    .def("__abs__",&vec_abs<double>)

    .def("__pow__",&vec_pow<double>)
    .def("sq",&vec_sq<double>)
    .def("cb",&vec_cb<double>)
    .def("sqrt",&vec_sqrt<double>)
    .def("cbrt",&vec_cbrt<double>)
    .def("exp",&vec_exp<double>)
    .def("log",&vec_log<double>)

    .def("sin",&vec_sin<double>)
    .def("cos",&vec_cos<double>)
    .def("tan",&vec_tan<double>)
    .def("asin",&vec_asin<double>)
    .def("acos",&vec_acos<double>)
    .def("atan",&vec_atan<double>)
    .def("sinh",&vec_sinh<double>)
    .def("cosh",&vec_cosh<double>)
    .def("tanh",&vec_tanh<double>)
    .def("asinh",&vec_asinh<double>)
    .def("acosh",&vec_acosh<double>)
    .def("atanh",&vec_atanh<double>)
    ;


  class_<std::vector<double> >("vector")
    .def(vector_indexing_suite<std::vector<double> >())

    .def("__repr__", &print_vector_list<double>) 
    ;

}
