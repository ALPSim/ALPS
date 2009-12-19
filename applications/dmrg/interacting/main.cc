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

#include "headers.h"
#include "system.h"
#include "block.h"
#include "densmatrix.h"
#include "superblock.h"
#include "groundstate.h"


/*
namespace ietl{
  inline void mult(const Superblock& S, 
		   const WaveFunction& psi_in, WaveFunction& psi) { 
    S.multiply(psi_in, psi);
  }
}
*/

using namespace std;

int main() {
  
  // input parameters: length of the chain, number of sweeps
  cerr << "length of the chain, number of iterations? ";
  int length, nsweeps;  
  cin >> length >> nsweeps;
  std::vector<int> cutoff(nsweeps+1);
  
  cerr << "number of states kept for the warmup sweep?";
  cin >> cutoff[0];
  for(int i = 1; i < cutoff.size();i++){
    cerr << "number of states kept for sweep #" << i <<" ?";
    cin >> cutoff[i];
  }
  
  cerr << "use reflection symmetry? (1=yes) ";
  int rsym = 0;
  cin >> rsym;
  
  cerr << "length = " << length << ", nsweeps = " << nsweeps << endl; 
  cout << "#length = " << length << ", nsweeps = " << nsweeps << endl;
  
  //initialize "system" :
  System* psys = System::sysdefine();
  
  psys -> input();
  
  // ground state wave function:
  WaveFunction psi, psi_measure;
  // GS energy:
  double energy = 0.0;                    
  
  //define a siteblock with initial val. for the block operators
  Block siteblock(psys -> OpTypeList(), psys -> site_dim(), 1); 
  Block twosites = siteblock.addBlocks(siteblock);
  Block threesites = twosites.addBlocks(siteblock);
  
  // this vector will contain all block Hamiltonians 
  // computed in the DMRG-steps 1 ... length:
  std::vector<Block> allblocks(0);
  
  std::vector<Matrix> n_left(0);
  std::vector<Matrix> n_right(0);
  std::vector<Matrix>::iterator measure_n_it;
  
  // initial step: start with a siteblock 
  allblocks.push_back(siteblock);

  // start at...
  int Lmax = length - 3;
  double discarded_weight_now = 1.0;
  double discarded_weight_old = 1.0;
  
  // if using reflection symmetry:
  if(rsym == 1)
    Lmax = length/2 - 1;
  
  // Warmup sweep, using the infinite-system algorithm
  for(int i = 0; i < Lmax; i++)
    {                  
      cerr << "Infinite system sweep #" << i << " : \n";
      
      // add site to left block:
      Block leftblock = (*(allblocks.end()-1)).addBlocks(siteblock);
      
      // if no reflection symmetry: 
      Block rightblock = twosites;
      if(i%2 == 1 ) rightblock = threesites;
      
      // use reflection symmetry to define "right block"      
      if(rsym == 1)
	rightblock = leftblock.Reflect();
      
      Superblock S(leftblock, rightblock);
      
      psi.resize(leftblock.Dim(),rightblock.Dim());
      psi.clear();
      energy = GetGroundState(S,psi);
      
      cerr << energy << endl;
      cout << "in infinite system sweep #" << i << " : \n"; 
      
      cout << "ground state energy and wave function: \n energy = " 
	   << energy << endl << "wave function dimension = " 
	   << psi.size1() << " x " << psi.size2() << endl;
      
      cout << endl;
      
      // reduced density matrix of left system block:
      DensityMatrix rho(psi,Left);
      
      // 'DMRG basis':
      Matrix basis;
      discarded_weight_old = discarded_weight_now;
      rho.newBasis(cutoff[0], basis, discarded_weight_now);     
      
      //new block Hamiltonian in 'DMRG basis': block_left+site
      allblocks.push_back(leftblock.mergeBlocks(basis));
      
      if(i == Lmax -1){ //do measurements only when final system size is reached
	cerr << "measure <n_i> after infinite-system sweep \n";
	cout << "measure <n_i> after infinite-system sweep \n";	
	
	int l = 0;
	
	for(measure_n_it = leftblock.measure_n.begin() ; 
	    measure_n_it != leftblock.measure_n.end(); ++measure_n_it){
	  double n_l = 
	    ietl::dot(psi,
		      boost::numeric::ublas::prec_prod((*measure_n_it),
						       psi));
	  cout << l << " " << n_l << endl;
	  l++;
	}
	if(rsym != 1) // for refl.sym. systems, do measurements only on the leftblock
	  for(measure_n_it = rightblock.measure_n.end() - 1; 
	      measure_n_it >= rightblock.measure_n.begin(); --measure_n_it){
	    double n_l = 
	      ietl::dot(psi,
			boost::numeric::ublas::trans(prec_prod((*measure_n_it),
							       trans(psi))));
	    cout << l << " " << n_l << endl;
	    l++;
	  }  
      }
    }
  
  // Finite System Sweeps
  
  cerr << "starting finite system sweeps \n";
  
  for(int swp = 1; swp <= nsweeps; swp++)
    {
      int Lmin;
      std::vector<Block>::iterator allblocks_it_left;
      std::vector<Block>::iterator allblocks_it_right;
      double groundstate_energy;
      // if using reflection symmetry
      if(rsym == 1){
	Lmax = length/2+1;
	Lmin = 2;
	if(swp == 1)
	  allblocks.push_back((*(allblocks.end() - 2)).Reflect());
	else
	  *(allblocks.end() - 1) = (*(allblocks.end() - 3)).Reflect();
      }
      else{ // if not using reflection symmetry
	Lmax = length - 1;
	Lmin = 2;
	if(swp == 1)
	  allblocks.push_back(siteblock);//entry at length-1: siteblock
	else
	  *(allblocks.end() - 1) = siteblock;//entry at length-1: siteblock
      }
      
      allblocks_it_right = allblocks.end() - 1;	
      allblocks_it_left = allblocks_it_right - 2;
      
      // Right to Left:
      cerr << "right to left: \n";
      for(int i = Lmax; i > Lmin; i--)
	{
	  cerr << "Finite system step right -> left #" 
	       << swp <<"." << i << " : \n";
	  Block leftblock = (*allblocks_it_left).addBlocks(siteblock);
	  Block rightblock = (*allblocks_it_right).addBlocks(siteblock);
	  
	  Superblock S(leftblock,rightblock);

	  psi.resize(leftblock.Dim(),rightblock.Dim());
	  psi.clear();
	  energy = GetGroundState(S,psi); 
	  
	  cout << "in finite system step right -> left #" 
	       << swp << "." << i << " : \n"; 
	  cout << "ground state energy and wave function: \n energy = " 
	       << energy << endl << "wave function dimension = " 
	       << psi.size1() << " x " << psi.size2() << endl;	       
	  cout << endl;
	  
	  DensityMatrix rho(psi,Right);
	  
	  Matrix basis;
	  discarded_weight_old = discarded_weight_now;
	  rho.newBasis(cutoff[swp], basis, discarded_weight_now);     
	  
	  --allblocks_it_right;
	  --allblocks_it_left;
	  
	  *allblocks_it_right = rightblock.mergeBlocks(basis);
	  
	  if(discarded_weight_now <= discarded_weight_old){    
	    n_left.clear();
	    n_right.clear();
	    n_left = leftblock.measure_n;
	    n_right = rightblock.measure_n;
	    psi_measure = psi;
	    groundstate_energy = energy;
	  }
	  
	  if(i == Lmin + 1){
	    cerr << "measure <n_i> after right -> left sweep: \n";
	    cout << "measure <n_i> after right -> left sweep: \n";
	    
	    int l = 0;
	    for(measure_n_it = n_left.begin();
		measure_n_it != n_left.end(); ++measure_n_it){
	      
	      double n_l = ietl::dot(psi_measure,boost::numeric::ublas::
				     prec_prod((*measure_n_it),psi_measure));
	      
	      cout << l << " " << n_l << endl;    
	      l++;
	    }
	    for(measure_n_it = n_right.end() - 1;
		measure_n_it >= n_right.begin(); --measure_n_it){
	      
	      double n_l = ietl::dot(psi_measure,
				     trans(boost::numeric::ublas::
					   prec_prod((*measure_n_it),
						     trans(psi_measure))));

	      cout << l << " " << n_l << endl;    
	      l++;
	    }
	  }
	}
      
      cerr << "ground state energy after right -> left sweep #" 
	   << swp << " : " << groundstate_energy << endl;
      cout << "ground state energy after right -> left sweep #" 
	   << swp << " : " << groundstate_energy << endl;
      
      // Left to Right:
            
      cerr << "left to right: \n";
      
      // if not using reflection symmetry: 
      int Lleftmax = length - 3;
      // if using reflection symmetry:
      if(rsym == 1)
	Lleftmax = length/2 - 2;
      
      allblocks_it_left = allblocks.begin();
      allblocks_it_right = allblocks_it_left + 2;
            
      for(int i = 0; i < Lleftmax; i++)
	{
	  cerr << "Finite system step left -> right #" 
	       << swp << "." << i << " : \n"; 
	  Block leftblock = (*allblocks_it_left).addBlocks(siteblock);
	  Block rightblock = (*allblocks_it_right).addBlocks(siteblock);
	  
	  Superblock S(leftblock,rightblock);
	  psi.resize(leftblock.Dim(),rightblock.Dim());
	  psi.clear();
	  energy = GetGroundState(S,psi);	  
	  
	  cout << "in finite system sweep left -> right #" 
	       << swp << "." << i << " : \n"; 
	  cout << "ground state energy and wave function: \n energy = " 
	       << energy << endl << "wave function dimension = " 
	       << psi.size1() << " x " << psi.size2() << endl;
	  cout << endl;
	  
	  DensityMatrix rho(psi,Left);
	  
	  Matrix basis;
	  discarded_weight_old = discarded_weight_now;
	  rho.newBasis(cutoff[swp], basis, discarded_weight_now);     
	  
	  ++allblocks_it_right;
	  ++allblocks_it_left;
	  *allblocks_it_left = leftblock.mergeBlocks(basis);
	  
	  if(discarded_weight_now <= discarded_weight_old){    
	    n_left.clear();
	    n_right.clear();
	    n_left = leftblock.measure_n;
	    n_right = rightblock.measure_n;
	    psi_measure = psi;
	    groundstate_energy = energy;
	  }
	  
	  if(i == Lleftmax - 1){
	    
	    cerr << "measure <n_i> after left -> right sweep: \n";
	    int l = 0;
	    for(measure_n_it = n_left.begin();
		measure_n_it != n_left.end(); ++measure_n_it){
	      
	      double n_l = ietl::dot(psi_measure,boost::numeric::ublas::
				     prec_prod((*measure_n_it),psi_measure));
	      
	      cout << l << " " << n_l << endl;    
	      l++;
	    }
	    for(measure_n_it = n_right.end() - 1;
		measure_n_it >= n_right.begin(); --measure_n_it){
	      
	      double n_l = ietl::dot(psi_measure,
				     boost::numeric::ublas::
				     trans(prec_prod((*measure_n_it),
						     trans(psi_measure))));
	      cout << l << " " << n_l << endl;    
	      l++;
	    }
	    
	  }
	}
      
      cerr << "ground state energy after left -> right sweep #" 
	   << swp << " : " << groundstate_energy << endl;
      cout << "ground state energy after left -> right sweep #" 
	   << swp << " : " << groundstate_energy << endl;
      
    } // end sweeps
  
  cerr << "\n*************** The End ***************\n";
  
  return 0;
  
} // main end
