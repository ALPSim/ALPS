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

#ifndef SYSTEM_H
#define SYSTEM_H

#include "headers.h"

// for now, c and cdagger are treated separatly, although 
// they only are the complex conjugated of each other
// at the moment for the delta, a unit-matrix is used and updated
enum OpType {Hblock, cdagger, c, n, delta};

class OpTerm {
  
 public:
  
  // contains the different pairs of operators which should be multiplied
  // when adding two blocks
  std::pair<OpType,OpType> Pair() {return (*this).pair_;}
  void Pair(std::pair<OpType,OpType> input) {pair_ = input;}
  
  // the matrix, which the product of the above two operators is added on:
  OpType Target() {return (*this).target_;}
  void Target(OpType input) {target_ = input;}

  // a possible way to introduce the fermionic sign
  // (for spinless fermions at half filling this sign 
  // is always positive)
  int sign() {return (*this).sign_;}
  void sign(int input) {sign_ = input;}
  
  // the coefficient (e.g. t, V) of the pair of operators 
  // given by the Hamiltonian:
  double coefficient() {return (*this).coefficient_;}
  void coefficient(double input) {coefficient_ = input;}
  
  OpTerm(){} 
  
 private:
  std::pair<OpType, OpType> pair_;
  OpType target_;
  int sign_;
  double coefficient_;
  
};

class System {
  
 public:
  
  static System* sysdefine();
  
  virtual int site_dim() = 0; 
  virtual Matrix init(OpType& type) = 0; 
  virtual void input() = 0;
  virtual std::vector<OpTerm> OpTermList() = 0;
  virtual std::vector<OpType> OpTypeList() = 0;   
  
 private:
  static System* sysdefine_;
  
 protected:
  System(){}
  virtual ~System(){}
  
};

class tVSystem : public System {
  
 public:  
  
  int site_dim() {
    tvsite_dim = 2;
    return tvsite_dim;
  } 
  
  Matrix init(OpType& type) {
    
    Matrix H(2,2);
    H.clear();
    
    switch(type){

      // Half filling: add the corresponding chemical potential
      // to make the Hamiltonian particle-hole symmetric
    case Hblock :
      H(0,0) = 0.25*V_;
      H(1,1) = -0.5*V_  + 0.25*V_;
      H(0,1) = H(1,0) = 0.0;
      break;
      
    case cdagger :
      H(1,0) = 1.0;
      break;
      
    case c :
      H(0,1) = 1.0;
      break;
      
    case n :
      H(1,1) = 1.0;
      break;
      
    case delta :
      H(0,0) = 1.0;
      H(1,1) = 1.0;
      break;
    default:
      boost::throw_exception(std::logic_error("reached default in switch statement in tVSystem::init"));
    }
    return H;
  }

  void input(){
    cerr << "nearest neighbour hopping t, ";
    cerr << "nearest neighbour interaction V? ";
    double t, V;
    cin >> t >> V;
    
    t_ = t;
    V_ = V;
    
    cerr << "t = " << t << ", V = " << V << endl;
    cout << "#t = " << t << ", V = " << V << endl;
    
  } 
  
  vector<OpTerm> OpTermList() { 
    
    // The different pairs of operators, which have to be multiplied,
    // and the corresponding 'target-operator':

    optermlist.resize(0);
    
    opterm.Pair(std::make_pair(Hblock,delta));
    opterm.coefficient(1.0);
    opterm.sign(1);
    opterm.Target(Hblock);
    optermlist.push_back(opterm);    
    
    opterm.Pair(std::make_pair(delta,Hblock));
    opterm.coefficient(1.0);
    opterm.sign(1);
    opterm.Target(Hblock);
    optermlist.push_back(opterm);
    
    opterm.Pair(std::make_pair(cdagger,c));
    opterm.coefficient(-1.0*t_);
    opterm.sign(1);
    opterm.Target(Hblock);
    optermlist.push_back(opterm);
    
    opterm.Pair(std::make_pair(c,cdagger));
    opterm.coefficient(-1.0*t_);
    opterm.sign(1);
    opterm.Target(Hblock);
    optermlist.push_back(opterm);
    
    opterm.Pair(std::make_pair(n,n));
    opterm.coefficient(V_);
    opterm.sign(1);
    opterm.Target(Hblock);
    optermlist.push_back(opterm);
    
    opterm.Pair(std::make_pair(delta,cdagger));
    opterm.coefficient(1.0);
    opterm.sign(1);
    opterm.Target(cdagger);
    optermlist.push_back(opterm);
    
    opterm.Pair(std::make_pair(delta,c));
    opterm.coefficient(1.0);
    opterm.sign(1);
    opterm.Target(c);
    optermlist.push_back(opterm);
    
    opterm.Pair(std::make_pair(delta,n));
    opterm.coefficient(1.0);
    opterm.sign(1);
    opterm.Target(n);
    optermlist.push_back(opterm);
    
    return optermlist;
  }
  
  vector<OpType> OpTypeList() { 

    // all the operators that are needed to describe a block:

    optypelist.resize(0);
    
    optypelist.push_back(Hblock);
    optypelist.push_back(cdagger);
    optypelist.push_back(c);
    optypelist.push_back(n);
    optypelist.push_back(delta);
    
    return optypelist;
  }
  
 private:
  
  int tvsite_dim;  
  double t_, V_;
  std::vector<OpTerm> optermlist;
  OpTerm opterm;
  std::vector<OpType> optypelist;
};

#endif
