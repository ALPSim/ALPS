/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2003-2004 by Salvatore R. Manmana <Salva@theo3.physik.uni-stuttgart.de>,
*                            Reinhard M. Noack <Reinhard.Noack@physik.uni-marburg.de>
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

#include "superblock.h"

// Implementation of the 'multiply'-function 
// returning H_super |psi> needed for IETL-Lanczos: 

void Superblock::multiply(const WaveFunction& psi_in, WaveFunction& psi_out) const {  
  using namespace boost::numeric::ublas;
  
  psi_out.clear();
  
  System* psys = System::sysdefine();
  
  std::vector<OpTerm> oppairs(psys -> OpTermList());
  std::vector<OpTerm>::iterator vecopterm_it;
  
  for(vecopterm_it = oppairs.begin(); 
      vecopterm_it != oppairs.end(); 
      ++vecopterm_it){
    
    if((*vecopterm_it).Target() == Hblock){
      
      psi_out += (*vecopterm_it).sign() * (*vecopterm_it).coefficient() *
	prec_prod(leftblock.operators[(*vecopterm_it).Pair().first],
	     trans(prec_prod(rightblock.operators[(*vecopterm_it).Pair().second],
			trans(psi_in))));
      
    } 
  } 
  
}
