/*****************************************************************************
*
* ALPS Project Applications: Directed Worm Algorithm 
*
* Copyright (C) 2013 by Lode Pollet      <pollet@phys.ethz.ch>, 
*                       Ping Nang Ma     <pingnang@phys.ethz.ch>,
*                       Matthias Troyer  <troyer@phys.ethz.ch>  
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


#include "bandstructure.hpp"



BOOST_PYTHON_MODULE(bandstructure_c) {

boost::python::class_<bandstructure>("bandstructure", boost::python::init<double, double, double, double, unsigned int>())
  .def(boost::python::init<boost::python::object, boost::python::object, double, double, unsigned int>())
  .def("__repr__", static_cast<std::string (bandstructure::*)()>(&bandstructure::representation))

  .def("t"     , static_cast<std::vector<double> (bandstructure::*)()>(&bandstructure::get_t))
  .def("U"     , static_cast<double              (bandstructure::*)()>(&bandstructure::get_U))
  .def("Ut"    , static_cast<std::vector<double> (bandstructure::*)()>(&bandstructure::get_Ut))

  .def("norm"  , static_cast<std::vector<double> (bandstructure::*)()>(&bandstructure::get_norm))

  .def("q"     , static_cast<std::vector<double> (bandstructure::*)(unsigned int)>(&bandstructure::get_q))
  .def("wk2"   , static_cast<std::vector<double> (bandstructure::*)(unsigned int)>(&bandstructure::get_wk2))

  .def("wk2_c" , static_cast<double              (bandstructure::*)()>(&bandstructure::get_wk2_c))
  .def("wk2_d" , static_cast<double              (bandstructure::*)()>(&bandstructure::get_wk2_d)) 
  ;

}
