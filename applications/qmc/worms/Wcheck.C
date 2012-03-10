/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2001-2004 by Matthias Troyer <troyer@comp-phys.org>,
*                            Simon Trebst <trebst@comp-phys.org>
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

#include "WRun.h"
#include <iomanip>

//- Checks, output --------------------------------------------------------
  
void WRun::check_spins()
{
#ifdef Print_spins
  print_spins();
#endif
  cyclic_iterator w0 = worm_head[0].kink();
  cyclic_iterator w1 = worm_head[1].kink();
  if(have_worm && 
      (w0->state() - (w0-1)->state() +
       w1->state() - (w1-1)->state() )) {
      std::cerr <<   *w0 <<  " " << *w1<< "\n";
      std::cerr <<   *(w0-1) <<  " " << *(w1-1) << "\n";
      boost::throw_exception(std::logic_error("Worm end state problem"));
    }
    bool result=false;
    for (int i=0;i<num_sites();i++) {
#ifdef Print_spins
    cyclic_iterator k=first_kink(i);
    std::cerr << "Site: " << i << "\n";
    if (k.valid())
      do {
        std::cerr << "Kinks : " << *k << " Time: " << k->time() << " state: " << int(k->state())
                  << " num_neighbors: " << *(k-1) << ", " << *(k+1) << " / ";
        for (int nb=0;nb<num_neighbors(i);nb++) {
          cyclic_iterator nn=first_kink(neighbor(i,nb));
          if(nn.valid()) {
            cyclic_iterator nnn(kinks[neighbor(i,nb)],k->adjacent(nn));
            std::cerr << *nnn << " ";
            time_struct t=k->time();
            time_struct t1=nnn->time();
            time_struct t2=(nnn+1)->time();
            if( t>=t1 && t2>t1 && t>t2 || t<t1 && (t2>t1 || t2<t))
              std::cerr << " Neighbor problem ";
          }
        }
        if(k->state()!=(k-1)->state() &&
             k->state()!=annihilate((k-1)->state()) && 
             k->state()!=create((k-1)->state()))
               std::cerr << "state_type problem\n";
          std::cerr << "\n";
          for (int jj=0;jj<num_sites();jj++)
            if(jj!=i && k==first_kink(jj))
              std::cerr << "Site problem: " << *k << " " << i << " " << jj << "\n";
    ++k;
    } while(k!=first_kink(i));
#endif
    cyclic_iterator j=first_kink(i);
    if(j.valid())
      do {
              if(j+1!=first_kink(i)&& (j+1)->time()< j->time()) {
                std::cerr << "Bond " << i << " starting at " << *j << " has time problem: " << (j+1)->time() << " " <<  j->time() << "\n";
                result=true;
              }
              ++j;
      } while(j!=first_kink(i));
  }
  if(result)
    boost::throw_exception(std::logic_error("invalid spin configuration"));
}   // WRun::check_spins
   
  
void WRun::print_spins()
{     
  std::cout << parms;
  std::cout << "Spin configuration:\n";
  std::cout << "Wormheads at " << worm_head[0].site() << " " << worm_head[0].time() 
            << " and " << worm_head[1].site() << " " << worm_head[1].time() << std::endl;
  for (int i=0;i<num_sites();i++) {
    if(!kinks[i].empty()) {
      std::cout << "Site: " << i << "\n";
      for (iterator k=kinks[i].begin();k!=kinks[i].end();++k)
        std::cout << "Kink : " << *k << "\n";
    }
  }
}   // WRun::print_spins  

