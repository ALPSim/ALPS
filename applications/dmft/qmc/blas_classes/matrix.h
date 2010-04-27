/*****************************************************************************
 *
 * ALPS DMFT Project - BLAS Compatibility headers
 *  Square Matrix Class
 *
 * Copyright (C) 2005 - 2009 by 
 *                              Emanuel Gull <gull@phys.columbia.edu>,
 *
 *
 * THIS SOFTWARE NEEDS AN APPROPRIATE LICENSE BLOCK HERE
 *****************************************************************************/

#ifndef BLAS_MATRIX
#define BLAS_MATRIX

#include "./blasheader.h"
#include "./vector.h"
#include "./resizeable_vector.h"
#include "./general_matrix.h"
#include <cmath>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <sys/types.h> //on osx for uint

#ifdef UBLAS
#include <alps/config.h> // needed to set up correct bindings
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
typedef boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> dense_matrix;
typedef boost::numeric::ublas::matrix<std::complex<double>,boost::numeric::ublas::column_major> complex_dense_matrix;
#endif

namespace blas{
  //simple BLAS matrix that uses BLAS calls for rank one and matrix vector.
  std::ostream &operator<<(std::ostream &os, const matrix &M); //forward declaration
  class matrix{
        public:
    matrix(int size){
      if(size>0)
        values_=new double[size*size];
      else values_=0;
      size_=size;
      memory_size_=size;
    }
    matrix(){ 
      values_=0;
      size_=0;
      memory_size_=0;
    }
    ~matrix(){
      if(values_!=0)
        delete[] values_;
    }
    matrix(int size1, int size2){
      assert(size1==size2); //square matrix only
      if(size1>0)
        values_=new double[size1*size2];
      else
        values_=0;
      size_=size1;
      memory_size_=size1;
    }
    const matrix& operator=(const matrix &M){
      if(this==&M) return *this;
      size_=M.size_;
      if(memory_size_!=M.memory_size_){
        memory_size_=M.memory_size_;
        if(values_!=0)
          delete[] values_;
        if(memory_size_>0){
          values_=new double[M.memory_size_*M.memory_size_];
        }else values_=0;
      }
      if(values_!=0)
        memcpy(values_, M.values_, memory_size_*memory_size_*sizeof(double));
      return *this;
    }
    matrix(const matrix&M){
      size_=M.size_;
      memory_size_=M.memory_size_;
      if(memory_size_>0){
        values_=new double[M.memory_size_*M.memory_size_];
        memcpy(values_, M.values_, memory_size_*memory_size_*sizeof(double));
      } else values_=0;
    }
    inline double &operator()(const unsigned int i, const unsigned int j){return *(values_+(i*memory_size_+j));}
    inline const double &operator()(const unsigned int i, const unsigned int j) const {return *(values_+(i*memory_size_+j));}
    //matrix size
    inline const int &size()const{return size_;}
    inline const int &memory_size()const{return memory_size_;}
    //matrix size 
    inline int size1()const{return size_;}
    inline int size2()const{return size_;}
    inline void add_outer_product(const blas::vector &v1, const blas::vector &v2, double alpha=1.){
      add_outer_product(&v1(0), &v2(0), alpha);
    }
    inline void add_outer_product(const blas::rsvector &v1, const blas::rsvector &v2, double alpha=1.){
      add_outer_product(&v1(0), &v2(0), alpha);
    }
    inline void add_outer_product(const double *v1, const double *v2, double alpha=1.){
      int inc=1;
      if(size_>1){
        dger_(&size_, &size_, &alpha,v2, &inc, v1, &inc, values_, &memory_size_); 
      }else if(size_==1){
        values_[0]+=alpha*v1[0]*v2[0];
      }else
        return;
    }
    inline double sum(){ 
      double sum=0;
      for(int i=0;i<size_;++i){
        for(int j=0;j<size_;++j){
          sum+=operator()(i,j);
        }
      }
      return sum;
    }
    inline void insert_row_column_last(blas::vector &row, blas::vector &col, double Mkk){
      resize(size_+1);
      int one=1;
      int oldsize=size_-1;
      dcopy_(&oldsize, &(col(0)), &one, &(values_[oldsize     ]), &memory_size_); //copy in row (careful: col. major)
      dcopy_(&oldsize, &(row(0)), &one, &(values_[oldsize*memory_size_]), &one         );   //copy in column
      operator()(oldsize, oldsize)=Mkk;
    }
    inline void getrow(int k, double *row) const{
      int one=1;
      dcopy_(&size_, &(values_[k*memory_size_]), &one, row, &one);
    }
    inline void getcol(int k, double *col) const{
      int one=1;
      dcopy_(&size_, &(values_[k]), &memory_size_, col, &one);
    }
    inline void setrow(int k, const double *row){
      int one=1;
      dcopy_(&size_, row, &one, &(values_[k*memory_size_]), &one);
    }
    inline void setcol(int k, const double *col){
      int one=1;
      dcopy_(&size_, col, &one, &(values_[k]), &memory_size_);
    }
    inline void insert_row_column(blas::vector &row, blas::vector &col, double Mkk, int k){
      resize(size_+1);
      double *new_values=new double[size_*size_];
      //insert a row and a column at position k. could be done more
      //efficiently, I suppose.
      for(int j=size_-2;j>=0;--j){
        int j_new=j<k?j:j+1;
        for(int i=size_-2;i>=0;--i){
          int i_new=i<k?i:i+1;
          new_values[i_new*size_+j_new]=operator()(i, j);
        }
        new_values[k*size_+j_new]=row(j);
        new_values[j_new*size_+k]=col(j);
      }
      new_values[k*size_+k]=Mkk;
      delete[] values_;
      values_=new_values;
      memory_size_=size_;
    }
    //delete column k.
    inline void remove_row_column(int c){
      int new_size=size_-1;
      blas::matrix new_matrix(new_size);
      for(int i=0, k=0;i<new_size;++i,++k){
        if(k==c) ++k;
        for(int j=0,l=0;j<new_size;++j,++l){
          if(c==l) ++l;
          new_matrix(i,j)=operator()(k, l);
        }
      }
      this->swap(new_matrix);
    }
    //delete last column
    inline void remove_row_column_last(){
      size_--;
    }
    //swap two columns:
    inline void swap_row_column(int c1, int c2){
      if(c1==c2) return;
      int one=1;
      dswap_(&size_, &(values_[c1]), &memory_size_, &(values_[c2]), &memory_size_);
      dswap_(&size_, &(values_[c1*memory_size_]), &one, &(values_[c2*memory_size_]), &one);
    }
    template<class mat> matrix convert_from(const mat &M){ //convert dense or sparse matrices (or any other with  (i,j) and size1() to blas matrix.
      assert(M.size1()==M.size2());
      resize_nocopy(M.size1());
      for(int i=0;i<(int)M.size1();++i){
        for(int j=0;j<(int)M.size1();++j){
          operator()(i,j)=M(i,j);
        }
      }
      return *this;
    }
    inline void right_multiply(const vector &v1, vector &v2) const{ //perform v2[i]=M[ij]v1[j]
      //call the BLAS routine for matrix vector multiplication:
      char trans='T';
      double alpha=1., beta=0.;        //no need to multiply a constant or add a vector
      int inc=1;
      dgemv_(&trans, &size_, &size_, &alpha, values_, &memory_size_, v1.values_, &inc, &beta, v2.values_, &inc);
    }
    inline void right_multiply(const rsvector &v1, rsvector &v2) const{ //perform v2[i]=M[ij]v1[j]
      //call the BLAS routine for matrix vector multiplication:
      char trans='T';
      double alpha=1., beta=0.;        //no need to multiply a constant or add a vector
      int inc=1;
      dgemv_(&trans, &size_, &size_, &alpha, values_, &memory_size_, v1.values_, &inc, &beta, v2.values_, &inc);
    }
    inline void left_multiply(const vector &v1, vector &v2) const{ //perform v2[i]=v1[j]M[ji]
      //call the BLAS routine for matrix vector multiplication:
      char trans='N';
      double alpha=1., beta=0.;       //no need to multiply a constant or add a vector
      int inc=1;
      dgemv_(&trans, &size_, &size_, &alpha, values_, &memory_size_, v1.values_, &inc, &beta, v2.values_, &inc);
      
    }
    inline void left_multiply(const rsvector &v1, rsvector &v2) const{ //perform v2[i]=v1[j]M[ji]
      //call the BLAS routine for matrix vector multiplication:
      char trans='N';
      double alpha=1., beta=0.;       //no need to multiply a constant or add a vector
      int inc=1;
      dgemv_(&trans, &size_, &size_, &alpha, values_, &memory_size_, v1.values_, &inc, &beta, v2.values_, &inc);
      
    }
    //Mres=this*M2
    inline void matrix_right_multiply(const matrix &M2,  matrix &Mres) const{
      //call the BLAS routine for matrix matrix multiplication:
      double one_double=1.;
      double zero_double=0.;
      char notrans='N';
      Mres.clear(); //to satisfy valgrind on hreidar / acml
      dgemm_(&notrans, &notrans, &size_, &size_, &size_, &one_double, M2.values_, &memory_size_, values_, &memory_size_, &zero_double, Mres.values_, &Mres.memory_size_);
    }
    inline void matrix_right_multiply(const blas::general_matrix<double> &M2,  blas::general_matrix<double>&Mres) const{
      //call the BLAS routine for matrix matrix multiplication:
      double one_double=1.;
      double zero_double=0.;
      char notrans='N';
      Mres.clear(); //to satisfy valgrind on hreidar / acml
      //      dgemm_(&notrans, &notrans, &M2.size2_, &size1_, &size2_,
      //      &one_double, M2.values_, &M2.size2_, values_,
      //     &size2_, &beta, Mres.values_, &M2.size2_);
      int M2_size2=M2.size2();
      dgemm_(&notrans, &notrans, &M2_size2, &size_, &size_, 
             &one_double, &(M2(0,0)), &(M2.size2()), values_,
             &memory_size_, &zero_double, &(Mres(0,0)), &(Mres.size2()));
    }
    inline double trace()const{
      double tr=0;
      for(int i=0;i<size_;++i) tr+=operator()(i,i);
      return tr;
    }
    inline double determinant() const{
      //the simple ones...
      if(size_==0) return 1;
      if(size_==1) return values_[0];
      if(size_==2) return values_[0]*values_[3]-values_[1]*values_[2];
      int info=0;
      int *ipiv = new int[size_];
      assert(size_==memory_size_); //otherwise think about plugging ing memory size to lda
      
      blas::matrix det_matrix(*this); //will be changed...
      blas::matrix identity(*this);
      identity.set_to_identity();
      //LU factorization
      lapack::dgesv_(&size_, &size_, det_matrix.values_, &memory_size_, ipiv, identity.values_, &memory_size_,&info);
      if(info < 0) {
        std::cout << "LAPACK ERROR IN APPLY!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        std::cout << "INFO:" << info << std::endl;
      } else if(info > 0){
        //check dgesv page: that means that we haven an exactly singular
        //matrix and the det is therefore =0: 
        return 0.;
      }
      //compute the determinant:
      double det=1.;
      //CAREFUL when pivoting: fortran uses array indexing starting
      //from one. But we're in C here -> 'one off' error
      for(int i=0;i<size_;++i){
        if(ipiv[i]-1!=i){
          det*=-det_matrix(i,i);
        }
        else{
          det*=det_matrix(i,i);
        }
      }
      delete[] ipiv;
      return det;
    }
    inline void transpose(){
      for(int i=0;i<size_;++i){
        for(int j=0;j<i;++j){
          std::swap(operator()(i,j), operator()(j,i));
        }
      }
    }
    void set_to_identity(){
      for(int i=0;i<size_;++i){
        for(int j=0;j<size_;++j){
          operator()(i,j)=0;
        }
        operator()(i,i)=1;
      }
    }
    void clear(){
      memset(values_, 0., memory_size_*memory_size_*sizeof(double));
      /*for(int i=0;i<size_;++i){
        for(int j=0;j<size_;++j){
          operator()(i,j)=0;
        }
      }*/
    }
    bool operator!=(const matrix &M)const {return ! operator==(M);}
    bool operator==(const matrix &M)const {
      if(size_!=M.size_) return false;
      for(int i=0;i<size_;++i){
        for(int j=0;j<size_;++j){
          if(std::fabs(M(i,j) - operator()(i,j))>1.e-5){
            return false;
          }
        }
      }
      return true;
    }
    //multiply matrix by value.
    matrix &operator *=(double lambda){
      //int inc=1;
      //int total_size=memory_size_*memory_size_;
      //dscal_(&total_size, &lambda, values_, &inc);
      for(int i=0;i<size_;++i){
        for(int j=0;j<size_;++j){
          operator()(i,j)*=lambda;
        }
      }
      return *this;
    }
    matrix operator+(const matrix &M2) const{
      matrix Msum(*this);
      Msum+=M2;
      return Msum;
    }
    matrix &operator+=(const matrix &M2) {
      assert(size_==M2.size_);
      for(int i=0;i<size_;++i){
        for(int j=0;j<size_;++j){
          operator()(i,j)+=M2(i,j);
        }
      }
      return *this;
    }
    matrix operator-(const matrix &M2) const{
      matrix Msum(*this);
      for(int i=0;i<size_;++i){
        for(int j=0;j<size_;++j){
          Msum(i,j)=operator()(i,j)-M2(i,j);
        }
      }
      return Msum;
    }
    void resize(int size1, int size2){
      if(size1!=size2){std::cerr<<"size1 has to be size2. aborting. "<<std::endl; abort();}
      resize(size1);
    }
    void resize_nocopy(int new_size){
      if(new_size<=memory_size_){
        size_=new_size;
        //memset(values_, 0, sizeof(double)*new_size*new_size);
        return;
      }
      if(size_!=0)
        delete[] values_;
      values_=new double[new_size*new_size];
      //memset(values_, 0, sizeof(double)*new_size*new_size);
      size_=new_size;
      memory_size_=new_size;
    }
    void resize(int new_size){
      if(new_size==size_) return;
      if(new_size<size_){ //down is easy
        if((int)new_size < (int)(memory_size_-30) && (int) new_size > 10){
          double *new_values_=new double[new_size*new_size];
          for(int i=0;i<new_size;++i){
            memcpy(new_values_+i*new_size, values_+i*memory_size_, sizeof(double)*new_size); //for each row: copy the entire row.
          }
          delete[] values_;       //free memory
          values_=new_values_;    //let the matrix point to the new memory location.  
          size_=new_size;
          memory_size_=new_size;
        }
        size_=new_size;
        return;
      } else if(new_size<= memory_size_){ //up is easy as long as we don't have to allocate new memory
        size_=new_size;
      } else{ //get new memory */
        double *new_values_=new double[new_size*new_size];
        for(int i=0;i<size_;++i){
          memcpy(new_values_+i*new_size, values_+i*memory_size_, sizeof(double)*size_); //for each row: copy the entire row.
        }
        delete[] values_;       //free memory
        values_=new_values_;    //let the matrix point to the new memory location.
        size_=new_size;
        memory_size_=new_size;
      }
    }
    double max() const{
      double *rowmax=new double[size_];
      int rowmax_index, max_index, inc=1;
      for(int i=0;i<size_;++i){
        rowmax_index=idamax_(&size_, values_+i*memory_size_,&inc);
        rowmax[i]=*(values_+i*memory_size_+rowmax_index-1); //fortran convention: start counting from one
      }
      max_index=idamax_(&size_, rowmax ,&inc);
      delete[] rowmax;
      return std::abs(rowmax[max_index-1]);
    }
    void swap(matrix &M2){
      std::swap(size_, M2.size_);
      std::swap(memory_size_, M2.memory_size_);
      std::swap(values_, M2.values_); //just exchange pointers
    }
    void eigenvalues_eigenvectors_symmetric( vector &eigenvalues, matrix &eigenvectors) const{
      //perform dsyev call (LAPACK)
      eigenvectors=*this;
      assert(memory_size_==size_);
      int lwork=-1;
      int info;
      double work_size;
      char jobs='V';
      char uplo='L';
      //get optimal size for work
      lapack::dsyev_(&jobs, &uplo, &size_, eigenvectors.values_, &size_, &(eigenvalues(0)),&work_size, &lwork, &info); 
      lwork=(int)work_size;
      double *work=new double[lwork];
      lapack::dsyev_(&jobs, &uplo, &size_, eigenvectors.values_, &size_, &eigenvalues(0),work, &lwork, &info); 
      delete[] work;
      /*matrix eigenvectors_trans(eigenvectors);
       eigenvectors_trans.transpose();
       std::cout<<*this<<std::endl;
       matrix diag(size_);
       diag.clear();
       for(int i=0;i<size_;++i){
       diag(i,i)=eigenvalues(i);
       }
       std::cout<<eigenvectors<<std::endl;
       std::cout<<diag<<std::endl;
       std::cout<<(eigenvectors_trans*diag*eigenvectors)<<std::endl;
       std::cout<<(eigenvectors*(*this)*eigenvectors_trans)<<std::endl;*/
    }
    //compute M=D*M, where D is a diagonal matrix represented by the
    //vector.
    void multiply_diagonal_matrix(const blas::vector &diagonal_matrix){
      int inc=1;
      for(int i=0;i<size_;++i){
        dscal_(&size_, &diagonal_matrix(i), values_+i*memory_size_, &inc);
      }
    }
    //same, but M=M*D. Slow because of much larger stride.
    void right_multiply_diagonal_matrix(const blas::vector &diagonal_matrix){
      for(int i=0;i<size_;++i){
        dscal_(&size_, &diagonal_matrix(i), values_+i, &memory_size_);
      }
    }
        private:
    int size_; //current size of matrix
    int memory_size_; //current size of matrix
    double *values_; //where the actual values are stored
  };
  
  inline vector operator*(const matrix &M, const vector &v1){
    assert(v1.size()==M.size());
    vector vres(v1.size());
    if(M.size()>0)
      M.right_multiply(v1, vres);
    return vres;
  }
  inline vector operator*(const vector &v1, const matrix &M){
    assert(v1.size()==M.size());
    vector vres(v1.size());
    if(M.size()>0)
      M.left_multiply(v1, vres);
    return vres;
  }
  inline matrix operator*(const matrix &M, const double lambda){
    matrix M2=M;
    M2*=lambda;
    return M2;
  }
  inline matrix operator*(const double lambda, const matrix &M){
    matrix M2=M;
    M2*=lambda;
    return M2;
  }
  inline matrix operator*(const matrix M, const matrix M2){
    matrix Mres(M.size());
    M.matrix_right_multiply(M2, Mres);
    return Mres;
  }
  /*inline std::ostream &operator<<(std::ostream &os, const matrix &M){
   //os<<"[ ";
   for(int i=0;i<M.size();++i){
   //os<<"[ ";
   for(int j=0;j<M.size();++j){
   os<<M(i,j)<<" ";
   }
   os<<std::endl;
   }
   //os<<"]"<<std::endl;
   return os;
   }*/
  
  inline std::ostream &operator<<(std::ostream &os, const matrix &M){
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
} //namespace
#ifdef UBLAS
inline complex_dense_matrix mult(const complex_dense_matrix &A, const blas::matrix &B){ //slow matrix multiplication real & complex.
  if((int)A.size1() !=(int)B.size() || (int)A.size2()!=(int)B.size()){ std::cerr<<"wrong matrix sizes in mult."<<std::endl; exit(1);}
  int size=A.size1();
  //complex_dense_matrix C(B.size1(), B.size());
  complex_dense_matrix D(B.size1(), B.size());
  complex_dense_matrix B_complex(B.size1(), B.size());
  B_complex.clear(); D.clear();
  for(int i=0;i<B.size();++i){ //copy matrix
    for(int j=0;j<B.size();++j){
      B_complex(i,j)=B(i,j);
    }
  }
  //blas call
  std::complex<double> one_complex=1.;
  std::complex<double> zero_complex=0.;
  char notrans='N';
  blas::zgemm_(&notrans, &notrans, &size, &size, &size, &one_complex, (void*)&(A(0,0)), &size, (void*)&(B_complex(0,0)), &size,(void*) &zero_complex,(void*) &(D(0,0)), &size);
  /*C.clear();
  for(int i=0;i<B.size();++i){
    for(int j=0;j<B.size();++j){
      for(int k=0;k<B.size();++k){
        C(i,j)+=A(i,k)*B(k,j);
      }
    }
  }
  std::cout<<D-C<<std::endl;*/
  return D;
}
inline complex_dense_matrix mult(const blas::matrix &A, const complex_dense_matrix &B){ //slow matrix multiplication real & complex.
  if((int)B.size1() !=A.size() || (int)B.size2()!=A.size()){ std::cerr<<"wrong matrix sizes in mult."<<std::endl; exit(1);}
  int size=B.size1();
  //complex_dense_matrix C(A.size(), A.size());
  complex_dense_matrix D(B.size1(), B.size1());
  complex_dense_matrix A_complex(A.size(), A.size());
  D.clear(); A_complex.clear();
  for(int i=0;i<A.size();++i){ //copy matrix
    for(int j=0;j<A.size();++j){
      A_complex(i,j)=A(i,j);
    }
  }
  //blas call
  std::complex<double> one_complex=1.;
  std::complex<double> zero_complex=0.;
  char notrans='N';
  blas::zgemm_(&notrans, &notrans, &size, &size, &size, &one_complex, (void*)&(A_complex(0,0)), &size, (void*)&(B(0,0)), &size,(void*) &zero_complex,(void*) &(D(0,0)), &size);
  /*C.clear();
  for(int i=0;i<A.size();++i){
    for(int j=0;j<A.size();++j){
      for(int k=0;k<A.size();++k){
        C(i,j)+=A(i,k)*B(k,j);
      }
    }
  }*/
  return D;
}
inline complex_dense_matrix mult(const complex_dense_matrix &A, const complex_dense_matrix &B){ //slow matrix multiplication real & complex.
  if(B.size1() !=A.size1() || B.size2()!=A.size1() || A.size1() !=A.size2()){ std::cerr<<"wrong matrix sizes in mult."<<std::endl; exit(1);}
  int size=B.size1();
  complex_dense_matrix C(A.size1(), A.size1());
  complex_dense_matrix D(A.size1(), A.size1());
  C.clear(); D.clear();
  //blas call
  std::complex<double> one_complex=1.;
  std::complex<double> zero_complex=0.;
  char notrans='N';
  blas::zgemm_(&notrans, &notrans, &size, &size, &size, &one_complex, (void*)&(A(0,0)), &size, (void*)&(B(0,0)), &size,(void*) &zero_complex,(void*) &(D(0,0)), &size);
  return D;
/*  C.clear();
  for(int i=0;i<A.size1();++i){
    for(int j=0;j<A.size1();++j){
      for(int k=0;k<A.size1();++k){
        C(i,j)+=A(i,k)*B(k,j);
      }
    }
  }
  return C;*/
}
#endif //UBLAS
  /*class block_diagonal_matrix: public std::vector<matrix>{
public:
        double trace(){
          double tr=0;
          for(int i=0;i<size();++i){
            tr+=operator[](i).trace();
          }
          return tr;
 */ /*class block_diagonal_matrix: public std::vector<matrix>{
public:
        double trace(){
          double tr=0;
          for(int i=0;i<size();++i){
            tr+=operator[](i).trace();
          }
          return tr;
        }
} ;
*/

#endif
