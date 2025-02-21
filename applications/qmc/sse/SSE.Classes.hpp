/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2001-2003 by Fabien Alet <alet@comp-phys.org>,
*                            Matthias Troyer <troyer@comp-phys.org>
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

#ifndef SSE_CLASSES_HPP
#define SSE_CLASSES_HPP

#include <alps/osiris/dump.h>

// Main object
template <class StateType=uint8_t>
struct Vertex
{
  typedef StateType state_type;
  Vertex() : vertex_type(0) {} 
  
  Vertex(Vertex const& v) 
   : vertex_type(v.vertex_type)
   , bond_number(v.bond_number)
  {
    leg[0]=v.leg[0];
    leg[1]=v.leg[1];
    leg[2]=v.leg[2];
    leg[3]=v.leg[3];
    linked_vertices[0]=v.linked_vertices[0];
    linked_vertices[1]=v.linked_vertices[1];
    linked_vertices[2]=v.linked_vertices[2];
    linked_vertices[3]=v.linked_vertices[3];
  }

  Vertex const& operator=(Vertex const& v) 
  {
    vertex_type=v.vertex_type;
    bond_number=v.bond_number;
    leg[0]=v.leg[0];
    leg[1]=v.leg[1];
    leg[2]=v.leg[2];
    leg[3]=v.leg[3];
    linked_vertices[0]=v.linked_vertices[0];
    linked_vertices[1]=v.linked_vertices[1];
    linked_vertices[2]=v.linked_vertices[2];
    linked_vertices[3]=v.linked_vertices[3];
    return *this;
  }

  bool operator==(Vertex const& v)
  {
    return vertex_type==v.vertex_type
      && bond_number==v.bond_number
      && leg[0]==v.leg[0]
      && leg[1]==v.leg[1]
      && leg[2]==v.leg[2]
      && leg[3]==v.leg[3];
  }
  
  
  // Type of vertex
  state_type vertex_type; // (0=Id, 1=Diagonal, 2=Non Diagonal)
  // The bond where this vertex sits
  alps::uint32_t bond_number; 

  // States at the four legs
  state_type leg[4]; // 0=down left,1=down right,2=up left,3=up right
  
  // The vertices linked to the current one
  // uint32_t is the vertex number in Operator String
  // bool specifies whether it is target's left or right leg
  std::pair < alps::uint32_t, bool > linked_vertices[4];

  bool identity() { return (vertex_type==0); }
  bool diagonal() { return (vertex_type==1); }
  bool non_diagonal() { return (vertex_type==2); }

  // Modify vertex 
  void modify_vertex(state_type,bool&,bool&); 
};

template <class StateType>
inline void Vertex<StateType>::modify_vertex(state_type exit_leg,bool& incoming_leg,bool& is_upper_leg)
  // Modify the exit leg, changes vertex type
{
  
  // Modify exit leg according to the convention used
  // This has to be modified for other convention
  is_upper_leg=(exit_leg>1);
  if (is_upper_leg)
    ++leg[exit_leg]; 
  else 
    --leg[exit_leg];

  // Modify the vertex type
  if ((leg[0]==leg[2]) && (leg[1]==leg[3]))
    { vertex_type=1;} else {vertex_type=2;}

  // Modify the bool incoming_leg (left or right) for next vertex
  incoming_leg=((exit_leg==1) || (exit_leg==3));
}

template <class StateType>
inline alps::ODump& operator<<(alps::ODump& dump, const Vertex<StateType>& v)
  // save vertex
{
  return dump << v.vertex_type << v.bond_number << v.leg[0] << v.leg[1] << v.leg[2] << v.leg[3];
}
 
template <class StateType>
inline alps::IDump& operator>>(alps::IDump& dump, Vertex<StateType>& v)
  // load vertex
{
  return dump >> v.vertex_type >> v.bond_number >> v.leg[0] >> v.leg[1] >> v.leg[2] >> v.leg[3];
}

#endif
