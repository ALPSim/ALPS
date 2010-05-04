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
 * THIS SOFTWARE NEEDS AN APPROPRIATE LICENSE BLOCK HERE
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
  std::cout <<"z1: "<<v.zone1()<<" c1^dag: "<<v.c_dagger_1()
            <<" c1: "<<v.c_1()<<" z2: "<<v.zone2()
            <<" c2^dag: "<<v.c_dagger_2()<<" c2: "
            <<v.c_2()<<"\t"<<v.abs_w();
  return os;
}

ostream &operator<<(ostream &os, const creator &c)
{
  std::cout<<c.zone()<<" "<<c.s()<<" "<<c.t();
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

std::ostream & operator<<(std::ostream &os, const resizeable_matrix &M)
{
  os<<"[ \n";
  for(uint i=0;i<M.size();++i){
    os<<"[";
    for(uint j=0;j<M.size();++j){
      os<<M(i,j)<<" ";
    }
    os<<"] \n";
  }
  os<<"] "<<std::endl;
  
  std::cout<<"creators: ";
  for(uint i=0;i<M.creators().size();++i){
    os<<M.creators()[i]<<"\t";
  }
  os<<std::endl<<"annihils: ";
  for(uint i=0;i<M.creators().size();++i){
    os<<M.annihilators()[i]<<"\t";
  }
  os<<std::endl;
  return os;
}
