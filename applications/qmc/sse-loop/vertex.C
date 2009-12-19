/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2001-2003 by Matthias Troyer <troyer@comp-phys.org>,
*                            Fabien Alet <alet@comp-phys.org>,
*                            Andreas Lange <alange@phys.ethz.ch>
*
* This software is part of the ALPS Applications, published under the ALPS
* Application License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Application License along with
* the ALPS Applications; see the file LICENSE.txt. If not, the license is also
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

/* $Id$ */

#include "vertex.h"
#include <functional>

vertex_number_type Vertex::n_ = 0; // TODO: remove

std::ostream& operator<<(std::ostream& out, const Vertex& v)
{
  
  if(v.is_diagonal())
      out << "] diagonal";
    else 
      out << "] nondiagonal";
    out << " on bond connecting " << v.site(0) 
        << " and " << v.site(1) << "\n";
  return out;
}


std::ostream& operator<< (std::ostream& ost, const StateVector& stv)
{
  int i;
  for (i=0;i<stv.size();i++)
    {
    if (stv[i])
      ost << 1;
    else ost << 0;
    }
   
  return ost;
}

void VertexString::init() // TODO: remove
{
  Vertex::n_=std::count_if(super_type::begin(),super_type::end(),std::not1(std::mem_fun_ref(&Vertex::is_identity)));
}
