/*****************************************************************************
*
* ALPS Project Applications: Directed Worm Algorithm 
*
* Copyright (C) 2013 by Matthias Troyer  <troyer@phys.ethz.ch> ,
*                       Lode Pollet      <pollet@phys.ethz.ch> ,
*                       Ping Nang Ma     <pingnang@phys.ethz.ch> 
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



#include <alps/ngs/boost_python.hpp>
#include <boost/bind.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <alps/python/numpy_array.hpp>

#include "../worldlines.hpp"
#include "../bandstructure.hpp"


BOOST_PYTHON_MODULE(dwa_c) {

boost::python::class_<std::vector<unsigned int> >("std_vector_unsigned_int")
  .def(boost::python::vector_indexing_suite<std::vector<unsigned int> >())
  ;
boost::python::class_<std::vector<double> >("std_vector_double")
  .def(boost::python::vector_indexing_suite<std::vector<double> >())
  ;
boost::python::class_<std::vector<unsigned short> >("std_vector_unsigned_short")
  .def(boost::python::vector_indexing_suite<std::vector<unsigned short> >())
  ;

boost::python::class_<std::vector<std::vector<unsigned int> > >("std_vector_std_vector_double")
  .def(boost::python::vector_indexing_suite<std::vector<std::vector<unsigned int> > >())
  ;
boost::python::class_<std::vector<std::vector<double> > >("std_vector_std_vector_double")
  .def(boost::python::vector_indexing_suite<std::vector<std::vector<double> > >())
  ;
boost::python::class_<std::vector<std::vector<unsigned short> > >("std_vector_std_vector_unsigned_short")
  .def(boost::python::vector_indexing_suite<std::vector<std::vector<unsigned short> > >())
  ;

boost::python::class_<kink>("kink", boost::python::init<unsigned int>())
  .def(boost::python::init<unsigned int, double, unsigned short>())

  .def("__repr__", &kink::representation)
  
  .def("siteindicator", &kink::siteindicator)
  .def("time"         , &kink::time)
  .def("state"        , &kink::state)
  ;

boost::python::class_<worldlines::location_type>("location_type", boost::python::no_init);

boost::python::class_<worldlines>("worldlines", boost::python::init<>())
  .def(boost::python::init<unsigned int>())

  .def("__repr__", &worldlines::representation)

  .def("load", static_cast<void (worldlines::*)(std::string const &)>(&worldlines::load))
  .def("save", static_cast<void (worldlines::*)(std::string const &) const>(&worldlines::save))  

  .def("open_worldlines", &worldlines::open_worldlines)

  .def("worldlines_siteindicator" , &worldlines::worldlines_siteindicator)
  .def("worldlines_time"          , &worldlines::worldlines_time)
  .def("worldlines_state"         , &worldlines::worldlines_state)

  .def("num_sites", &worldlines::num_sites)
  .def("num_kinks", &worldlines::num_kinks)

  .def("states", &worldlines::states)

  .def("location", &worldlines::location)

  .def("state_before", &worldlines::state_before)
  .def("state",        &worldlines::state)

  .def("is_valid", static_cast<bool (worldlines::*)(unsigned short)>(&worldlines::is_valid))
  ;  

boost::python::class_<wormpair>("wormpair", boost::python::init<>())
  .def(boost::python::init<worldlines::location_type, kink, bool, bool>())

  .def("__repr__", &wormpair::representation)

  .def("wormhead", &wormpair::wormhead)
  .def("wormtail", &wormpair::wormtail)

  .def("wormhead_site",    &wormpair::site)
  .def("wormhead_time",    &wormpair::time)
  .def("wormhead_forward", &wormpair::forward)

  .def("wormtail_site", &wormpair::wormtail_site) 
  .def("wormtail_time", &wormpair::wormtail_time)

  .def("next_partnersite", &wormpair::next_partnersite) 
  .def("next_time"       , &wormpair::next_time)

  .def("wormhead_turns_around"                          , &wormpair::wormhead_turns_around)
  .def("wormhead_moves_to_new_time"                     , &wormpair::wormhead_moves_to_new_time)
  .def("wormhead_inserts_vertex_and_jumps_to_new_site"  , &wormpair::wormhead_inserts_vertex_and_jumps_to_new_site)
  .def("wormhead_deletes_vertex_and_jumps_to_new_site"  , &wormpair::wormhead_deletes_vertex_and_jumps_to_new_site)
  .def("wormhead_relinks_vertex_and_jumps_to_new_site"  , &wormpair::wormhead_relinks_vertex_and_jumps_to_new_site)
  .def("wormhead_crosses_vertex"                        , &wormpair::wormhead_crosses_vertex)
  .def("wormhead_annihilates_wormtail"                  , &wormpair::wormhead_annihilates_wormtail)
  ;

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
