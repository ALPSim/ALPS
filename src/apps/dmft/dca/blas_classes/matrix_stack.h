/*********************************************************
*  BLAS Compatibility headers
*  Matrix Stack class
 * Copyright (C) 2005-2026 by Emanuel Gull
 *
 * ALPS Project: https://alps.comp-phys.org/
 * SPDX-License-Identifier: MIT
 *
*********************************************************/
#ifndef MAT_STACK
#define MAT_STACK
#include<stack>

//new and delete are both TERRIBLY slow. Instead, we recycle our
//matrices by pushing them onto a stack.
template<class M> class matrix_stack{
  public: 
    matrix_stack(int size1, int size2, int max_stack_size){
      size1_=size1;
      size2_=size2;
      max_stack_size_=max_stack_size;
      matrices_=new std::stack<M*>*[size1_];
      for(int i=0;i<size1_;++i){
        matrices_[i]=new std::stack<M*>[size2_];
      }
    }
    ~matrix_stack(){ /*apres nous le deluge - I don't care, only one instance of this class, and it's static.*/ }
    void release(M* mat){
      if(mat->size1()<size1_ && mat->size2()<size2_){
        if(matrices_[mat->size1()][mat->size2()].size()<(uint)max_stack_size_){
          matrices_[mat->size1()][mat->size2()].push(mat);
        }else
          delete mat; //stack is full
      }else //matrix too large for stack
        delete mat;
    }
    M* get(int size1, int size2){
      if(size1>=size1_ || size2>=size2_){
        return new M(size1, size2);
      }else if(matrices_[size1][size2].size()>0){
        M* res=matrices_[size1][size2].top();
        matrices_[size1][size2].pop();
        return res;
      }else{
        return new M(size1, size2);
      }
      abort();
      return new M(size1, size2);
    }
  private:
    int size1_; //max matrix size in x
    int size2_; //max matrix size in y
    int max_stack_size_; //max size of stack
    std::stack<M*>** matrices_; //the actual stack: a stack of M pointers pointing at a matrix M with sizes i,j.
};

#endif
