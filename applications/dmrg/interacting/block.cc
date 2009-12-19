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

#include "block.h"

Block Block::addBlocks(Block& Other){
  
  using namespace boost::numeric::ublas;
  
  Block& This = *this;
  
  // WAY OF COUNTING IN THE BASIS |ij>: |00>, |10>, |20>,...,|01>,|11>,|21>,...
  // |i> : states of the block 'This', |j> states of the block 'Other'
  
  int dim = This.Dim() * Other.Dim();
  
  System* psys = System::sysdefine();
  
  Block newblock(psys -> OpTypeList(), dim);
  
  std::vector<OpTerm> OpPairs(psys -> OpTermList());
  
  std::vector<OpTerm>::iterator optermvec_it;

  // perform all operations needed to 'update' a block, when it shall represent the result
  // of the addition of two blocks:
  
  for(optermvec_it = OpPairs.begin(); 
      optermvec_it != OpPairs.end(); 
      ++optermvec_it){
    

    if((*optermvec_it).Pair().second == delta){
      
      for(int j = 0; j < Other.Dim(); j++)
	
	if((*optermvec_it).Target() != delta)
	  matrix_range<Matrix> 
	    (newblock.operators[(*optermvec_it).Target()], 
	     range(j * This.Dim() , (j+1) * This.Dim()), 
	     range(j * This.Dim() , (j+1) * This.Dim())) 
	    += (*optermvec_it).sign() * (*optermvec_it).coefficient() * 
	    This.operators[(*optermvec_it).Pair().first];
      
    }
    
    else{
      
      for(int j = 0; j < Other.Dim(); j++){	
	for(int jp = 0; jp < Other.Dim(); jp++){
	  
	  if((*optermvec_it).Target() != delta)
	    if(Other.operators[(*optermvec_it).Pair().second](j,jp) != 0){
	      
	      matrix_range<Matrix>
		(newblock.operators[(*optermvec_it).Target()], 
		 range(j * This.Dim() , (j+1) * This.Dim()), 
		 range(jp * This.Dim(), (jp+1) * This.Dim())) 
		+= (*optermvec_it).sign() * (*optermvec_it).coefficient() * 
		Other.operators[(*optermvec_it).Pair().second](j,jp) * 
		This.operators[(*optermvec_it).Pair().first];  
	      
	    }
	}
      }
    }
  } // end optermvec_it-iterator
  
  // the 'delta-operator' is treated extra:
  newblock.operators[delta].clear();
  for(int j = 0; j < dim; j++)
    newblock.operators[delta](j,j) = 1.0;
  
  // update operators for the measurements:
  std::vector<Matrix>::iterator This_n_it;
  
  for(This_n_it = This.measure_n.begin(); This_n_it != This.measure_n.end(); 
      ++This_n_it){
    
    Matrix temp;
    temp.resize(dim,dim);
    temp.clear();
    
    for(int j = 0; j < Other.Dim(); j++)
      matrix_range<Matrix> 
	(temp, 
	 range(j * This.Dim() , (j+1) * This.Dim()), 
	 range(j * This.Dim() , (j+1) * This.Dim())) 
	+= *This_n_it;
    
    newblock.measure_n.push_back(temp);
    
  }

  //(assuming that Other is always a siteblock -> true in the actual version...)
  newblock.measure_n.push_back(newblock.operators[n]);
  
  return newblock;
  
}

Block Block::mergeBlocks(Matrix& basis) {
  
  using namespace boost::numeric::ublas;    
  
  Block& This = *this;

  Block result(This); //copy constructor 
  
  map<OpType,Matrix>::iterator map_it;
  
  for(map_it = result.operators.begin(); 
      map_it != result.operators.end(); 
      ++map_it){
    
    if ((*map_it).first != delta){     
      
      (*map_it).second.resize(basis.size2(),basis.size2());
    
      // the basis change to the new DMRG-basis is done here:
      (*map_it).second = 
	prec_prod(trans(basis),
	 	  Matrix(prec_prod(This.operators[(*map_it).first],basis)));
    } 
  } 
  
  // again, the delta is treated extra:
  result.operators[delta].resize(basis.size2(),basis.size2());
  result.operators[delta].clear();
  for(int i = 0; i < basis.size2(); i++)
    result.operators[delta](i,i) = 1.0;
  
  //measurements:
  std::vector<Matrix>::iterator vector_it;
  std::vector<Matrix>::iterator vector_This_it = This.measure_n.begin();
  
  for(vector_it = result.measure_n.begin(); 
      vector_it != (result.measure_n.end() - 1); 
      ++vector_it){
    
    *vector_it =
      prec_prod(trans(basis),Matrix(
		prec_prod((*vector_This_it),
			  basis)));
    
    ++vector_This_it;
  
  }
  
  // last site has already been transformed, save one step...
  *(result.measure_n.end() - 1) = result.operators[n]; 
  
  return result;
}
