/****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2012 by Emanuel Gull <gull@pks.mpg.de>,
 *                   
 *  based on an earlier version by Philipp Werner and Emanuel Gull
 *
 *
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/
#ifndef HYB_BLAS_MATRIX
#define HYB_BLAS_MATRIX
#include <boost/numeric/bindings/blas.hpp>
#include <boost/numeric/bindings/lapack.hpp>

//This is the resizable square matrix class. Whoever adapted the version in ALPS f*ed it up badly: this matrix needs to have a memory size DIFFERENT from the current size, or we'll be asking for memory all the time!

//simple matrix that uses BLAS calls for rank one and matrix vector.
class blas_matrix{
public:
  blas_matrix(fortran_int_t size){
    if(size>0)
      values_=new double[size*size];
    else values_=0;
    size_=size;
    memory_size_=size;
  }
  blas_matrix(){
    values_=0;
    size_=0;
    memory_size_=0;
  }
  ~blas_matrix(){
    if(values_!=0)
      delete[] values_;
  }
  blas_matrix(const blas_matrix&M){
    size_=M.size_;
    memory_size_=M.memory_size_;
    if(memory_size_>0){
      values_=new double[M.memory_size_*M.memory_size_];
      memcpy(values_, M.values_, memory_size_*memory_size_*sizeof(double));
    } else values_=0;
  }
  inline double &operator()(const fortran_int_t i, const fortran_int_t j){return *(values_+(i*memory_size_+j));}
  inline const double &operator()(const fortran_int_t i, const fortran_int_t j) const {return *(values_+(i*memory_size_+j));}
  //blas_matrix size
  inline const  fortran_int_t &size()const{return size_;}
  inline fortran_int_t &size(){return size_;}
  inline const fortran_int_t &memory_size()const{return memory_size_;}
  inline fortran_int_t &memory_size(){return memory_size_;}
  //blas_matrix size 
  /*inline void add_outer_product(const blas::vector &v1, const blas::vector &v2, double alpha=1.){
    add_outer_product(&v1(0), &v2(0), alpha);
  }
  inline void add_outer_product(const blas::rsvector &v1, const blas::rsvector &v2, double alpha=1.){
    add_outer_product(&v1(0), &v2(0), alpha);
  }*/
  inline void add_outer_product(fortran_int_t s, const double *v1, const double *v2, double alpha=1.){
    fortran_int_t inc=1;
    if(s>1){
      FORTRAN_ID(dger)(&s, &s, &alpha,v2, &inc, v1, &inc, values_, &memory_size_); 
    }else if(s==1){
      values_[0]+=alpha*v1[0]*v2[0];
    }else
      return;
  }
  /*inline void insert_row_column_last(blas::vector &row, blas::vector &col, double Mkk){
    resize(size_+1);
    fortran_int_t one=1;
    fortran_int_t oldsize=size_-1;
    dcopy_(&oldsize, &(col(0)), &one, &(values_[oldsize     ]), &memory_size_); //copy in row (careful: col. major)
    dcopy_(&oldsize, &(row(0)), &one, &(values_[oldsize*memory_size_]), &one         );   //copy in column
    operator()(oldsize, oldsize)=Mkk;
  }*/
  inline void getrow(fortran_int_t k, double *row) const{
    fortran_int_t one=1;
    dcopy_(&size_, &(values_[k*memory_size_]), &one, row, &one);
  }
  inline void getcol(fortran_int_t k, double *col) const{
    fortran_int_t one=1;
    dcopy_(&size_, &(values_[k]), &memory_size_, col, &one);
  }
  inline void setrow(fortran_int_t k, const double *row){
    fortran_int_t one=1;
    dcopy_(&size_, row, &one, &(values_[k*memory_size_]), &one);
  }
  inline void setcol(fortran_int_t k, const double *col){
    fortran_int_t one=1;
    dcopy_(&size_, col, &one, &(values_[k]), &memory_size_);
  }
  //delete last column
  inline void remove_row_column_last(){
    size_--;
  }
  //swap two columns:
  inline void swap_row_column(fortran_int_t c1, fortran_int_t c2){
    if(c1==c2) return;
    swap_row(c1,c2);
    swap_column(c1,c2);
  }
  inline void swap_column(fortran_int_t c1, fortran_int_t c2){
    if(c1==c2) return;
    FORTRAN_ID(dswap)(&size_, &(values_[c1]), &memory_size_, &(values_[c2]), &memory_size_);
  }
  inline void swap_row(fortran_int_t c1, fortran_int_t c2){
    if(c1==c2) return;
    fortran_int_t one=1;
    FORTRAN_ID(dswap)(&size_, &(values_[c1*memory_size_]), &one, &(values_[c2*memory_size_]), &one);
  }
  inline void right_multiply(const std::vector<double> &v1, std::vector<double> &v2) const{ //perform v2[i]=M[ij]v1[j]
    //call the BLAS routine for blas_matrix vector multiplication:
    char trans='T';
    double alpha=1., beta=0.;    //no need to multiply a constant or add a vector
    fortran_int_t inc=1;
    FORTRAN_ID(dgemv)(&trans, &size_, &size_, &alpha, values_, &memory_size_, &(v1[0]), &inc, &beta, &(v2[0]), &inc);
  }
  /*inline void right_multiply(const vector &v1, vector &v2) const{ //perform v2[i]=M[ij]v1[j]
    //call the BLAS routine for blas_matrix vector multiplication:
    char trans='T';
    double alpha=1., beta=0.;    //no need to multiply a constant or add a vector
    fortran_int_t inc=1;
    FORTRAN_ID(dgemv)(&trans, &size_, &size_, &alpha, values_, &memory_size_, v1.values_, &inc, &beta, v2.values_, &inc);
  }
  inline void right_multiply(const rsvector &v1, rsvector &v2) const{ //perform v2[i]=M[ij]v1[j]
    //call the BLAS routine for blas_matrix vector multiplication:
    char trans='T';
    double alpha=1., beta=0.;    //no need to multiply a constant or add a vector
    fortran_int_t inc=1;
    FORTRAN_ID(dgemv)(&trans, &size_, &size_, &alpha, values_, &memory_size_, v1.values_, &inc, &beta, v2.values_, &inc);
  }
  inline void left_multiply(const vector &v1, vector &v2) const{ //perform v2[i]=v1[j]M[ji]
    //call the BLAS routine for blas_matrix vector multiplication:
    char trans='N';
    double alpha=1., beta=0.;       //no need to multiply a constant or add a vector
    fortran_int_t inc=1;
    FORTRAN_ID(dgemv)(&trans, &size_, &size_, &alpha, values_, &memory_size_, v1.values_, &inc, &beta, v2.values_, &inc);
    
  }
  inline void left_multiply(const rsvector &v1, rsvector &v2) const{ //perform v2[i]=v1[j]M[ji]
    //call the BLAS routine for blas_matrix vector multiplication:
    char trans='N';
    double alpha=1., beta=0.;       //no need to multiply a constant or add a vector
    fortran_int_t inc=1;
    FORTRAN_ID(dgemv)(&trans, &size_, &size_, &alpha, values_, &memory_size_, v1.values_, &inc, &beta, v2.values_, &inc);
  }*/
  void set_to_identity()
  {
    clear();
    for(fortran_int_t i=0;i<size_;++i){
      operator()(i,i)=1.;
    }
  }
  inline double determinant() const{
    //the simple ones...
    if(size_==0) return 1;
    if(size_==1) return values_[0];
    if(size_==2) return operator()(0,0)*operator()(1,1)-operator()(0,1)*operator()(1,0);
    fortran_int_t info=0;
    std::vector<fortran_int_t> ipiv(size_);
    //assert(size_==memory_size_); //otherwise think about plugging ing memory size to lda
    
    std::vector<double> det_blas_matrix(size()*size());
    for(fortran_int_t i=0;i<size();++i){ for (fortran_int_t j=0;j<size();++j){ det_blas_matrix[i*size()+j]=operator()(i,j);}}
    std::vector<double> identity(size()*size(), 0.); for(fortran_int_t i=0;i<size();++i){ identity[i*size()+i]=1.; }
    //LU factorization
    FORTRAN_ID(dgesv)(&size_, &size_, &(det_blas_matrix[0]), &size_, &(ipiv[0]), &(identity[0]), &size_,&info);
    if(info < 0) {
      std::cout << "lapack error in solver" << std::endl;
      std::cout << "INFO:" << info << std::endl;
    } else if(info > 0){
      //check dgesv page: that means that we haven an exactly singular
      //blas_matrix and the det is therefore =0: 
      return 0.;
    }
    //compute the determinant:
    double det=1.;
    //CAREFUL when pivoting: fortran uses array indexing starting
    //from one. But we're in C here -> 'one off' error
    for(fortran_int_t i=0;i<size_;++i){
      if(ipiv[i]-1!=i){
        det*=-det_blas_matrix[i*size()+i];
      }
      else{
        det*=det_blas_matrix[i*size()+i];
      }
    }
    return det;
  }
  void clear(){
    memset(values_, 0, memory_size_*memory_size_*sizeof(double));
  }
  void resize(fortran_int_t size1, fortran_int_t size2){
    if(size1!=size2){std::cerr<<"size1 has to be size2. aborting. "<<std::endl; abort();}
    resize(size1);
  }
  void resize_nocopy(fortran_int_t new_size){
    if(new_size<=memory_size_){
      size_=new_size;
      return;
    }
    if(size_!=0)
      delete[] values_;
    values_=new double[new_size*new_size];
    size_=new_size;
    memory_size_=new_size;
  }
  void resize(fortran_int_t new_size){
    if(new_size==size_) return;
    if(new_size<(fortran_int_t)(size_)){ //down is easy
      if((fortran_int_t)new_size < (fortran_int_t)(memory_size_-30) && (fortran_int_t) new_size > 10){
        double *new_values_=new double[new_size*new_size];
        for(fortran_int_t i=0;i<new_size;++i){
          memcpy(new_values_+i*new_size, values_+i*memory_size_, sizeof(double)*new_size); //for each row: copy the entire row.
        }
        delete[] values_;       //free memory
        values_=new_values_;    //let the blas_matrix point to the new memory location.  
        size_=new_size;
        memory_size_=new_size;
      }
      size_=new_size;
      return;
    } else if(new_size<= memory_size_){ //up is easy as long as we don't have to allocate new memory
      size_=new_size;
    } else{ //get new memory */
      double *new_values_=new double[new_size*new_size];
      for(fortran_int_t i=0;i<size_;++i){
        memcpy(new_values_+i*new_size, values_+i*memory_size_, sizeof(double)*size_); //for each row: copy the entire row.
      }
      delete[] values_;       //free memory
      values_=new_values_;    //let the blas_matrix point to the new memory location.
      size_=new_size;
      memory_size_=new_size;
    }
  }
  blas_matrix operator-(const blas_matrix &M2) const{
    blas_matrix Msum(*this);
    for(int i=0;i<size_;++i){
      for(int j=0;j<size_;++j){
        Msum(i,j)=operator()(i,j)-M2(i,j);
      }
    }
    return Msum;
  }
  double max() const{
    double *rowmax=new double[size_];
    fortran_int_t rowmax_index, max_index, inc=1;
    if(size_<1) return 0;
    for(int i=0;i<size_;++i){
      rowmax_index=idamax_(&size_, values_+i*memory_size_,&inc);
      rowmax[i]=*(values_+i*memory_size_+rowmax_index-1); //fortran convention: start counting from one
    }
    max_index=idamax_(&size_, rowmax ,&inc);
    delete [] rowmax;
    return std::abs(rowmax[max_index-1]);
  }
  void swap(blas_matrix &M2){
    std::swap(size_, M2.size_);
    std::swap(memory_size_, M2.memory_size_);
    std::swap(values_, M2.values_); //just exchange pointers
  }
  void invert(){
    std::vector<double> B(size_*size_, 0.);
    std::vector<fortran_int_t> ipiv(size_,0);
    fortran_int_t info;
    for(fortran_int_t i=0;i<size_;++i) B[i*size_+i]=1.;
    FORTRAN_ID(dgesv)(&size_, &size_, values_, &memory_size_, &(ipiv[0]), &(B[0]), &size_, &info);
    if(info){ throw(std::logic_error("in dgesv: info was not zero.")); }
    
    for(fortran_int_t i=0;i<size_;++i){
      for(fortran_int_t j=0;j<size_;++j){
        operator()(i,j)=B[i*size_+j];
      }
    }
  }
private:
  fortran_int_t size_; //current size of blas_matrix
  fortran_int_t memory_size_; //current size of blas_matrix
  double *values_; //where the actual values are stored
};

/*inline vector operator*(const blas_matrix &M, const vector &v1){
  assert(v1.size()==M.size());
  vector vres(v1.size());
  if(M.size()>0)
    M.right_multiply(v1, vres);
  return vres;
}
inline vector operator*(const vector &v1, const blas_matrix &M){
  assert(v1.size()==M.size());
  vector vres(v1.size());
  if(M.size()>0)
    M.left_multiply(v1, vres);
  return vres;
}*/

inline std::ostream &operator<<(std::ostream &os, const blas_matrix &M){
  os<<"[ ";
  for(int i=0;i<M.size();++i){
    //os<<"[ ";
    for(int j=0;j<M.size();++j){
      os<<M(i,j)<<" ";
    }
    if(i<M.size()-1)
      os<<" ;"<<" ";
  }
  os<<"]"<<" ";
  return os;
}

std::ostream   &operator<<(std::ostream   &os, const blas_matrix &M); //forward declaration

#endif
