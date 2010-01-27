/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2006 -2010 by Adrian Feiguin <afeiguin@uwyo.edu>
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

#ifndef __RANGE_H__
#define __RANGE_H__

namespace dmtk
{

class Range: public std::slice
{
  public:
    Range(size_t _o, size_t _d, size_t _s): std::slice(_o, (_d-_o+1)/_s, _s) {}
    Range(size_t _o, size_t _d): std::slice(_o, (_d-_o+1), 1) {}
    Range(const Range &r): std::slice(r) {}

    bool operator==(const Range &r)
     { 
       if(begin() == r.begin() && end() == r.end() && stride() == r.stride())
         return true;
       return false;
     }

    size_t begin() const { return std::slice::start(); }
    size_t end() const { return std::slice::start()+std::slice::size()-1; }
};

} // namespace dmtk

#endif // __RANGE_H__
