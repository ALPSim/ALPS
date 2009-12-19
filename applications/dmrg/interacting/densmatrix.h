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

#ifndef DENSMAT_H
#define DENSMAT_H

#include "headers.h"

// identify left or right block: 
enum LR {Left, Right};

class DensityMatrix {                   
  
 public:
  
  // constructor: returns reduced density matrix of 
  // the left or right system block, respectively with the
  // ground state wave function of the actual DMRG-sweep:
  DensityMatrix(const WaveFunction& psi, LR lr); 
  
  // now diagonalize rho, obtain the new Basis, and calculate 
  // the discarded weight:  
  void newBasis(int cutoff, Matrix& basis, double& discarded_weight);
  
 private:
  Matrix rho;
  
};

#endif 
