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
    template<class T>
    inline static boost::python::str print_value_with_error_element(value_with_error<T> const & self)
    {
      return boost::python::str(boost::python::str(self.mean()) + " +/- " + boost::python::str(self.error()));
    }

    template<class T>
    static boost::python::str print_value_with_error_container(value_with_error_container<T> self)
    {
      boost::python::str s;
      for (typename value_with_error_container<T>::index_type index=0; index < self.size(); ++index)
      {
        s += boost::python::str(boost::python::str(self.mean(index)) + " +/- " + boost::python::str(self.error(index)) + "\n");
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
  class_<value_with_error<double> >("value_with_error",init<optional<value_with_error<double>::value_type,value_with_error<double>::value_type> >())
    .add_property("mean", &value_with_error<double>::mean)
    .add_property("error",&value_with_error<double>::error)  

    .def("__repr__", &print_value_with_error_element<double>)

    .def(+self)
    .def(-self)
    .def("__abs__", &value_with_error<double>::abs)
    .def(self == value_with_error<double>())

    .def(self += value_with_error<double>())
    .def(self += value_with_error<double>::value_type())
    .def(self -= value_with_error<double>())
    .def(self -= value_with_error<double>::value_type())
    .def(self *= value_with_error<double>())
    .def(self *= value_with_error<double>::value_type())
    .def(self /= value_with_error<double>())
    .def(self /= value_with_error<double>::value_type())

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

    .def("__pow__",&value_with_error_pow<double>)
    .def("__sq__",&value_with_error_sq<double>)
    .def("cb",&value_with_error_cb<double>)
    .def("sqrt",&value_with_error_sqrt<double>)
    .def("cbrt",&value_with_error_cbrt<double>)
    .def("exp",&value_with_error_exp<double>)
    .def("log",&value_with_error_log<double>)

    .def("sin",&value_with_error_sin<double>)
    .def("cos",&value_with_error_cos<double>)
    .def("tan",&value_with_error_tan<double>)
    .def("asin",&value_with_error_asin<double>)
    .def("acos",&value_with_error_acos<double>)
    .def("atan",&value_with_error_atan<double>)
    .def("sinh",&value_with_error_sinh<double>)
    .def("cosh",&value_with_error_cosh<double>)
    .def("tanh",&value_with_error_tanh<double>)
    .def("asinh",&value_with_error_asinh<double>)
    .def("acosh",&value_with_error_acosh<double>)
    .def("atanh",&value_with_error_atanh<double>)

    .def_pickle(value_with_error_pickle_suite<double>())
    ;

/*
  class_<value_with_error_container<double> >("value_with_error_vector",init<optional<std::vector<double>,std::vector<double> > >())
    .def("__getitem__",&value_with_error_container<double>::getitem)
    .def("__getslice__",&value_with_error_container<double>::getslice)
    .def("__setitem__",&value_with_error_container<double>::setitem)
    .def("__setslice__",&value_with_error_container<double>::setslice)
    .def("__delitem__",&value_with_error_container<double>::delitem)
    .def("__delslice__",&value_with_error_container<double>::delslice)

    .def("__repr__", &print_value_with_error_container<double>)
    .def("mean",&value_with_error_container<double>::mean)
    .def("error",&value_with_error_container<double>::error)

    .def("__len__",&value_with_error_container<double>::size)
    .def("append",&value_with_error_container<double>::append)
    .def("extend",&value_with_error_container<double>::extend)
    .def("push_back",&value_with_error_container<double>::push_back)
    .def("fill", &value_with_error_container<double>::fill)
    .def("insert",&value_with_error_container<double>::insert)
    .def("pop_back",&value_with_error_container<double>::pop_back)
    .def("erase",&value_with_error_container<double>::erase) 
    .def("unfill",&value_with_error_container<double>::unfill)
    .def("clear",&value_with_error_container<double>::clear)

    .def(self + value_with_error_container<double>())
    .def(self + double())
    .def(double() + self)
    .def(self - value_with_error_container<double>())
    .def(self - double())
    .def(double() - self)
    .def(self * value_with_error_container<double>())
    .def(self * double())
    .def(double() * self)
    .def(self / value_with_error_container<double>())
    .def(self / double())
    .def(double() / self)

    .def(+self)
    .def(-self)
    .def("__abs__",&value_with_error_container_abs<double>)

    .def("__pow__",&value_with_error_container_pow<double>)
    .def("sq",&value_with_error_container_sq<double>)
    .def("cb",&value_with_error_container_cb<double>)
    .def("sqrt",&value_with_error_container_sqrt<double>)
    .def("cbrt",&value_with_error_container_cbrt<double>)
    .def("exp",&value_with_error_container_exp<double>)
    .def("log",&value_with_error_container_log<double>)

    .def("sin",&value_with_error_container_sin<double>)
    .def("cos",&value_with_error_container_cos<double>)
    .def("tan",&value_with_error_container_tan<double>)
    .def("asin",&value_with_error_container_asin<double>)
    .def("acos",&value_with_error_container_acos<double>)
    .def("atan",&value_with_error_container_atan<double>)
    .def("sinh",&value_with_error_container_sinh<double>)
    .def("cosh",&value_with_error_container_cosh<double>)
    .def("tanh",&value_with_error_container_tanh<double>)
    .def("asinh",&value_with_error_container_asinh<double>)
    .def("acosh",&value_with_error_container_acosh<double>)
    .def("atanh",&value_with_error_container_atanh<double>)

    ;

  class_<std::vector<double> >("float_vector")
    .def(vector_indexing_suite<std::vector<double> >())

    .def("__repr__", &print_vector_list<double>) 
    ;
    */
}

