/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Sebastian Fuchs <fuchs@theorie.physik.uni-goettingen.de>
 *                              Matthias Troyer <troyer@comp-phys.org>
 *
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

#include "interaction_expansion.hpp"

using namespace std;

ostream& operator <<(ostream &os, const vertex_array &vertices)
{
  for(std::size_t i=0;i<vertices.size();++i){
    os<<(vertices[i])<<std::endl;
  }
  return os;
}

ostream& operator <<(ostream &os, const vertex &v)
{
  std::cout <<"z1: "<<v.flavor1()<<" c1^dag: "<<v.c_dagger_1()
            <<" c1: "<<v.c_1()<<" z2: "<<v.flavor2()
            <<" c2^dag: "<<v.c_dagger_2()<<" c2: "
            <<v.c_2()<<"\t"<<v.abs_w();
  return os;
}

ostream &operator<<(ostream &os, const creator &c)
{
  std::cout<<c.flavor()<<" "<<c.s()<<" "<<c.t();
  return os;
}

ostream& operator << (ostream& os, const simple_hist &h)
{
  for(unsigned long i=0;i<h.size();++i){
    std::cout<<i<<"\t"<<h[i]<<std::endl;
  }
  return os;
}

void InteractionExpansionRun::load(alps::IDump &/*dump*/)
{
  std::cerr<<"loading: doing nothing."<<std::endl;
}

void InteractionExpansionRun::save(alps::ODump &/*dump*/) const
{
  std::cerr<<"saving: doing nothing"<<std::endl;
}

std::ostream & operator<<(std::ostream &os, const inverse_m_matrix &M)
{
    os << M.matrix() << std::endl;
  
  std::cout<<"creators: ";
  for(unsigned int i=0;i<M.creators().size();++i){
    os<<M.creators()[i]<<"\t";
  }
  os<<std::endl<<"annihils: ";
  for(unsigned int i=0;i<M.creators().size();++i){
    os<<M.annihilators()[i]<<"\t";
  }
  os<<std::endl;
  return os;
}
