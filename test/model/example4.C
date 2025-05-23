/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2003-2004 by Matthias Troyer <troyer@itp.phys.ethz.ch>
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

/* $Id$ */

#include <alps/model.h>
#include <fstream>
#include <iostream>

alps::multi_array<alps::Expression,2> bondmatrix(const alps::ModelLibrary lib, const std::string& name, const alps::Parameters& p=alps::Parameters())
{
  alps::HamiltonianDescriptor<short> ham=lib.get_hamiltonian(name,p,true);
  // ham.set_parameters(p);
  int dim=ham.basis().site_basis().num_states();
  
  // get site and bond terms
  alps::multi_array<alps::Expression,2> sitematrix = 
    alps::get_matrix(alps::Expression(),ham.site_term(),ham.basis().site_basis());
  alps::multi_array<alps::Expression,4> bondtensor = 
    alps::get_matrix(alps::Expression(),ham.bond_term(),ham.basis().site_basis(), ham.basis().site_basis());
    
  // add site terms to bond terms
  for (int i=0;i<dim;++i)
    for (int j=0;j<dim;++j)
      for (int k=0;k<dim;++k) {
        bondtensor[i][j][i][k]+=sitematrix[j][k];
        alps::simplify(bondtensor[i][j][i][k]);
        bondtensor[j][i][k][i]+=sitematrix[j][k];
        alps::simplify(bondtensor[j][i][k][i]);
      }
      
  //convert tensor into matrix
  alps::multi_array<alps::Expression,2> bondmatrix(boost::extents[dim*dim][dim*dim]);
  for (int i=0;i<dim;++i)
    for (int j=0;j<dim;++j)
      for (int k=0;k<dim;++k)
        for (int l=0;l<dim;++l)
          bondmatrix[i+j*dim][k+l*dim]=bondtensor[i][j][k][l];
  return bondmatrix;
}

int main()
{

#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif
    // create the library from an XML file
    std::ifstream in("../../lib/xml/models.xml");
    alps::ModelLibrary lib(in);

    // calculate bond matrices 
    alps::Parameters parms;
    
    std::cout << "HHardcoreBoson = \n" << bondmatrix(lib,"hardcore boson") << "\n\n";
    parms["Nmax"]=2;
    std::cout << "HBoson = \n" << bondmatrix(lib,"boson Hubbard",parms)  << "\n\n";
    parms["local_S"]="1/2";
    std::cout << "HSpinHalf = \n" << bondmatrix(lib,"spin")  << "\n\n";
    parms["local_S"]=1;
    std::cout << "HSpinOne = \n" << bondmatrix(lib,"spin",parms)  << "\n\n";
    

#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& e)
{
  std::cerr << "Caught exception: " << e.what() << "\n";
  exit(-1);
}
catch (...)
{
  std::cerr << "Caught unknown exception\n";
  exit(-2);
}
#endif
  return 0;
}
