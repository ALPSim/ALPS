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

    static boost::python::str print_value_with_error_vector(value_with_error<std::vector<double> > self)
    {
      boost::python::str s;
      for (std::size_t index=0; index < self.size(); ++index)
      {
        s += boost::python::str(self[index]);
        s += boost::python::str("\n");
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

  class_<value_with_error<std::vector<double> > >("value_with_error_vector",init<optional<std::vector<double>,std::vector<double> > >())
    //.def(vector_indexing_suite<value_with_error<std::vector<double> > >())   // what is left is just the proxy thing

    .add_property("mean",&value_with_error<std::vector<double> >::mean)
    .add_property("error",&value_with_error<std::vector<double> >::error)

    .def("__repr__", &print_value_with_error_vector)


    .def("__len__",&value_with_error<std::vector<double> >::size)          // will be depreciated after vector indexing suite is okay
    .def("append",&value_with_error<std::vector<double> >::push_back)      // will be depreciated after vector indexing suite is okay
    //.def("extend",&value_with_error<std::vector<double> >::extend)       // will be depreciated after vector indexing suite is okay
    .def("push_back",&value_with_error<std::vector<double> >::push_back)   // will be depreciated after vector indexing suite is okay
    //.def("fill", &value_with_error<std::vector<double> >::fill)          // will be depreciated after vector indexing suite is okay
    //.def("insert",&value_with_error<std::vector<double> >::insert)       // will be depreciated after vector indexing suite is okay
    //.def("pop_back",&value_with_error<std::vector<double> >::pop_back)   // will be depreciated after vector indexing suite is okay
    //.def("erase",&value_with_error<std::vector<double> >::erase)         // will be depreciated after vector indexing suite is okay
    //.def("unfill",&value_with_error<std::vector<double> >::unfill)       // will be depreciated after vector indexing suite is okay
    //.def("clear",&value_with_error<std::vector<double> >::clear)         // will be depreciated after vector indexing suite is okay


    .def(self + value_with_error<std::vector<double> >())
    .def(self + double())
    .def(double() + self)
    .def(self - value_with_error<std::vector<double> >())
    .def(self - double())
    .def(double() - self)
    .def(self * value_with_error<std::vector<double> >())
    .def(self * double())
    .def(double() * self)
    .def(self / value_with_error<std::vector<double> >())
    .def(self / double())
    .def(double() / self)

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


  class_<std::vector<double> >("float_vector")
    .def(vector_indexing_suite<std::vector<double> >())

    .def("__repr__", &print_vector_list<double>) 
    ;

}
