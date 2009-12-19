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

#ifndef BLOCK_H
#define BLOCK_H

#include "headers.h"
#include "system.h"

class Block {
  
 public:

  // contains all operators on the block:
  std::map<OpType,Matrix> operators; 
  
  // contains the operators n_i for measuring the local density:
  std::vector<Matrix> measure_n;
  
  // dimension of the matrices of the block-operators:
  int Dim() {return operators[Hblock].size1();}
  
  Block() { } 
  
  Block(vector<OpType> OpList, int dim){
    
    vector<OpType>::iterator vectoroptype_it;
    
    for(vectoroptype_it = OpList.begin(); vectoroptype_it != OpList.end(); ++vectoroptype_it){
      operators.insert(pair<OpType,Matrix>(*vectoroptype_it,Matrix(dim,dim)));
      operators[*vectoroptype_it].clear();
    }
    
    measure_n.resize(0);
    
  }
  
  // create a siteblock -> program in a more elegant way...
  Block(vector<OpType> OpList, int dim, int do_site){
    
    vector<OpType>::iterator vectoroptype_it;
    
    for(vectoroptype_it = OpList.begin(); vectoroptype_it != OpList.end(); ++vectoroptype_it){
      operators.insert(pair<OpType,Matrix>(*vectoroptype_it,Matrix(dim,dim)));
      
      if(do_site == 1){
        System* psys = System::sysdefine();
        operators[*vectoroptype_it] = (psys -> init(*vectoroptype_it)); 
      }
      else
        operators[*vectoroptype_it].clear();
    
    }
    measure_n.resize(0);
    measure_n.push_back(operators[n]);
  }  
  
  Block(const Block& Other){
    *this = Other;
  } 
  
  Block Reflect() const { return *this; }
    
  // function to add two blocks:
  Block addBlocks(Block& Other);
  
  // function to transform a block+site into the new DMRG-basis:
  Block mergeBlocks(Matrix& basis); 
  
};

#endif 

