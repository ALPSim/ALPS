#include <alps/testing/stream_fixture.hpp>
/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2003-2005 by Matthias Troyer <troyer@itp.phys.ethz.ch>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/parser/xslt_path.h>
#include <alps/model.h>
#include <fstream>
#include <iostream>

TEST(ModelSerialization, Example5) {
  alps::testing::StreamFixture transcript;
  { // Flush serialization objects before checking the captured stream.
    // create the library from an XML file
    std::ifstream in(alps::search_xml_library_path("models.xml"));
    alps::ModelLibrary lib(in);

    // get operators in one bond term

    std::cout << "Operator names:\n";
    std::set<std::string> names = lib.get_hamiltonian("spin",alps::Parameters(),true).bond_term().operator_names();
    std::copy(names.begin(),names.end(),std::ostream_iterator<std::string>(std::cout,"\n"));

    std::cout << "\nSplit terms:\n\n";

    typedef std::vector<boost::tuple<alps::Term,alps::SiteOperator,alps::SiteOperator > > V;
    V  ops = lib.get_hamiltonian("spin",alps::Parameters(),true).bond_term().split();
    for (V::const_iterator it=ops.begin(); it!=ops.end();++it)
      std::cout << "Prefactor: " << it->get<0>() << "\nSite 1: " << it->get<1>().term() << "\nSite 2: " << it->get<2>().term() << "\n\n";
  }
  transcript.expect_output(ALPS_TEST_SOURCE_DIR "/example5.output");
}
