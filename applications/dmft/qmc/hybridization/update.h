/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Philipp Werner <werner@itp.phys.ethz.ch>
 *                              Emanuel Gull <gull@phys.columbia.edu>
 *
 *
 * THIS SOFTWARE NEEDS AN APPROPRIATE LICENSE BLOCK HERE
 *****************************************************************************/
#ifndef ___UPDATE___
#define ___UPDATE___

#include "impurity.h"



template <class G> void construct_matrix(blas_matrix & M, segment_container_t & segments, double BETA,  G& F) {
  int N = segments.size();
  M.resize(N,N);
  int row=-1;
  int column=-1;
  for (segment_container_t::const_iterator it1=segments.begin(); it1!=segments.end(); it1++) {
    row++;
    for (segment_container_t::const_iterator it2=segments.begin(); it2!=segments.end(); it2++) {
      column++;
      
      double argument = it1->t_end()-it2->t_start();
      double sign = 1;
      if (argument<0) {
        argument += BETA;      
        sign = -1;
      }
      M(row,column) = interpolate_F(argument, BETA, F)*sign;
    }
    column = -1;
  }
  
}



// determine F(\tau)
template <class G> inline double interpolate_F(double t, double BETA, G& F) {
  
  double sign=1;
  if (t<0) {
    t += BETA;
    sign=-1;
  }
  
  int N = F.size()-1;
  double n = t/BETA*N;
  int n_lower = (int)n; // interpolate linearly between n_lower and n_lower+1
  
  return sign*(F[n_lower] + (n-n_lower)*(F[n_lower+1]-F[n_lower]));
  
}


// compute distances up/down to the next segment and iterators of these segments
// note: s_down always points to a physical segment, while s_up may point to segments.end() 
template <class S> void compute_intervals(double t, double BETA, double& t_up, double& t_down, S& segments, typename S::iterator & s_up, typename S::iterator& s_down) {  
  
  if (segments.size() == 0) {
    t_up = BETA;
    t_down = BETA;
    s_up = segments.end();
    s_down = segments.end();
  }
  else {
    
    s_up = lower_bound(segments.begin(), segments.end(), t);
    //s_up = segments.lower_bound(times(t,BETA));
    
    if (s_up == segments.begin()) {
      s_down = segments.end(); s_down--;
      if (s_down->t_end() < s_down->t_start())
        t_down = t - s_down->t_end();
      else
        t_down = t + BETA - s_down->t_end();
    }
    else {
      s_down = s_up; s_down--;
      if (s_down->t_end()>s_down->t_start())
        t_down = t - s_down->t_end();
      else 
        t_down = t - (BETA+s_down->t_end());
    }
    
    if(s_up == segments.end()) {
      t_up = BETA - t + segments.begin()->t_start();
    }
    else {
      t_up = s_up->t_start() - t;
    }
    
  }
  
}

// compute overlap between a segment and a list of segments
// requires segment with 0<=t_begin<t_end<=BETA
template <class S> inline double segment_overlap(times segment, const S& other_segments, int other_full_line, double BETA) {
  
  double length = (segment.t_start()<segment.t_end() ? segment.t_end()-segment.t_start() : segment.t_end()-segment.t_start()+BETA);
  double t_final = segment.t_start()+length;
  double t = segment.t_start();
  double t_final_segment;        
  double other_length=0;
  if (other_full_line==1)
    other_length=length;
  else if (other_segments.size()>0){
    typename S::const_iterator it; //this function does NOT change the segments.
    //it =lower_bound(other_segments.begin(), other_segments.end(), t);    
    //find first segment that has starting time AFTER t.
    it=other_segments.begin();
    while(it != other_segments.end() && it->t_start()<t){
        it++;
    }         
    if (it!=other_segments.begin()) {
      it--;
      //find end point of this segment
      t_final_segment = (it->t_start()<it->t_end() ? it->t_end() : it->t_end()+BETA);
      if (t<t_final_segment) {
        //overlap of segment pointed to by this iterator element with 'segment'
        other_length += (t_final_segment<t_final ? t_final_segment-t : t_final-t);
      }
      it++;
      
    }
    while(it!=other_segments.end() && it->t_start()<t_final) { //go through all segments that  have start times before this segment's end time.
      t_final_segment = (it->t_start()<it->t_end() ? it->t_end() : it->t_end()+BETA);
      //additional segment overlap
      other_length += (t_final_segment<t_final ? t_final_segment-it->t_start() : t_final-it->t_start());
      it++;
    }
    // check if last segment overlaps
    it=other_segments.end();
    //get the last segment
    it--;
    //wrap around, check for overlap with last segment
    if (it->t_end()<it->t_start() && t<it->t_end()) {
      other_length += (t_final<it->t_end() ? t_final-t : it->t_end()-t);
    }
  }     
  return other_length;   
}


template <class S> inline double compute_overlap(times segment, S& other_segments, int other_full_line, double BETA) {
  if (segment.t_start()<segment.t_end())
    return segment_overlap(segment, other_segments, other_full_line, BETA);
  else {
    double other_length=0;
    times segment1(0,segment.t_end());
    times segment2(segment.t_start(), BETA);
    other_length += segment_overlap(segment1, other_segments, other_full_line, BETA);
    other_length += segment_overlap(segment2, other_segments, other_full_line, BETA);
    return other_length;
  }
}


// functions required to compute determinant ratios and perform fast matrix updates 

template <class G, class S, class V> double det_rat_up(const times & new_segment, blas_matrix & M, const S& segments_old, G& F, V&Fe, V&Fs, double BETA, double & det_rat_sign, double & overlap) {
  
  typename S::const_iterator it=segments_old.begin();
  for (int i=0; i<(int)segments_old.size(); i++) {
    Fe(i) = interpolate_F(new_segment.t_end()-it->t_start(), BETA, F);
    Fs(i) = interpolate_F(it->t_end()-new_segment.t_start(), BETA, F);
    it++;
  }
  
  double det_rat = interpolate_F(new_segment.t_end()-new_segment.t_start(), BETA, F);
  det_rat -= Fe*M*Fs;
  /*for (int i=0; i<(int)M.size1(); i++) {
   for (int j=0; j<(int)M.size1(); j++) {
   det_rat -= Fe[i]*M(i,j)*Fs[j];
   }
   }*/
  
  // take care of sign changes produced by segments which "wind around"
  if (new_segment.t_end() < new_segment.t_start()) {
    det_rat *= -1;      
    overlap = -1;    
  }
  else {
    overlap = 1;
  }
    
  if (det_rat < 0) {
    det_rat_sign = -1;
    det_rat *= -1;
  }
  else {
    det_rat_sign = 1;
  }
  
  return det_rat;
}

template <class V> void compute_M_up(int k, blas_matrix & M, V& Fs, V& Fe, double det_rat) {
  
  blas_matrix M_new(M.size1()+1,M.size1()+1);
  int i_new, j_new;
  double det_rat_inv=1./det_rat;
  
  // element (k,k)
  M_new(k,k) = det_rat_inv;
  
  // row k and column k
  for (int i=0; i<(int)M.size1(); i++) {
    //set kth row and column  to zero.
    i_new = (i<k ? i : i+1);
    M_new(i_new,k) = 0;
    M_new(k,i_new) = 0;
    
    for (int n=0; n<(int)M.size1(); n++) {
      M_new(i_new,k) -= M(i,n)*Fs(n);
      M_new(k,i_new) -= M(n,i)*Fe(n);    
    } 
    M_new(i_new,k) *= det_rat_inv;
    M_new(k,i_new) *= det_rat_inv;
  }
  
  // remaining elements
  for (int i=0; i<(int)M.size1(); i++) {
    i_new = (i<k ? i : i+1);
    for (int j=0; j<(int)M.size1(); j++) {
      j_new = (j<k ? j : j+1);
      M_new(i_new,j_new) = M(i,j) + det_rat*M_new(i_new,k)*M_new(k,j_new);
    }
  }
    swap(M,M_new);
 // M_new.swap(M);
  /*blas::matrix M2(M);
   blas::vector k_row(M.size());
   blas::vector k_col(M.size());
   blas::vector f_start(M.size());
   blas::vector f_end(M.size());
   
   for(int i=0;i<M.size();++i){
   f_start(i)=Fs(i);
   f_end(i)=Fe(i);
   }
   if(M.size()>0){
   k_row=M*f_start;
   k_col=f_end*M;
   
   k_row*=-1./det_rat;
   k_col*=-1./det_rat;
   M.add_outer_product(k_row, k_col, det_rat);
   //std::cout<<"Up, M before resize: "<<M<<std::endl;
   }
   //this is the bottleneck for the BLAS implementation...
   M.insert_row_column(k_row, k_col, 1./det_rat,  k);
   //M_new.swap(M);
   //if(M_new !=M) std::cout<<"up: "<<k<<" "<<M_new<<" "<<M<<std::endl;*/
  return;
}  


template <class S> double det_rat_down(int k, blas_matrix & M, const S& segments_old, double & det_rat_sign) {
  
  double det_rat = M(k,k);
  
  // take care of sign changes produced by segments which "wind around"
  if (k==(int)segments_old.size()-1) {
    typename S::const_iterator it=segments_old.end(); it--;
    if (it->t_end() < it->t_start())
      det_rat *= -1;      
  }
    
  if (det_rat < 0) {
    det_rat_sign = -1;
    det_rat *= -1;
  }
  else {
    det_rat_sign = 1;
  }
  
  return det_rat;
}


template<class Mat> void compute_M_down(int k, Mat& M) {
  
  blas_matrix M_new(M.size1()-1, M.size1()-1);
  int i_old;
  double Mkk_inv=1./M(k,k);
  for (int i=0; i<(int)M_new.size1(); i++) {
    i_old = (i<k ? i : i+1);
    for (int j=0; j<k; j++) {
      M_new(i,j) = M(i_old, j)-M(i_old,k)*M(k,j)*Mkk_inv;
    }
    for (int j=k; j<(int)M_new.size1(); ){
      j++;
      M_new(i,j-1) = M(i_old, j)-M(i_old,k)*M(k,j)*Mkk_inv; //one year in prison for writing it like this.
    }
  }
  /*double Mkk=M(k,k);
   blas::vector row_k(M.size()-1);
   blas::vector col_k(M.size()-1);
   for(int i=0;i<M.size()-1;++i){
   int i_old=i<k?i:i+1;
   row_k(i)=M(i_old,k);
   col_k(i)=M(k,i_old);
   }
   //std::cout<<"removing row and column"<<std::endl;
   M.remove_row_column(k);
   if(M.size()>0)
   M.add_outer_product(row_k, col_k, -1./Mkk);*/
  //M.swap(M_new);
    swap(M_new,M);
  //if(M != M_new)
  //  std::cout<<"algo is wrong: "<<M<<" "<<M_new<<std::endl;
}

// move segment without changin its length
template <class G, class S> double det_rat_move(times & new_segment, int k, blas_matrix & M, const S& segments_old, G& F, double BETA, double & det_rat_sign, double & overlap) {
  
  double F_i, F_j;
  typename S::const_iterator it1, it2;
  
  double det_rat = M(k,k)*interpolate_F(new_segment.t_end()-new_segment.t_start(), BETA, F);
  
  it1=segments_old.begin();
  for (int i=0; i<M.size1(); i++) {
    if (i != k) {
      F_i = interpolate_F(new_segment.t_end()-it1->t_start(), BETA, F);
      
      it2=segments_old.begin();
      for (int j=0; j<M.size1(); j++) {
        if (j != k) {
          F_j = interpolate_F(it2->t_end()-new_segment.t_start(), BETA, F);
          det_rat -= F_i*(M(k,k)*M(i,j)-M(i,k)*M(k,j))*F_j;
        }
        it2++;
      }
    }
    it1++;
  }
  
  overlap = 1;
  // take care of sign changes produced by segments which "wind around"
  if (k==segments_old.size()-1) {
    it1--;
    // check if last segment has been shifted across beta
    if ((new_segment.t_end()-new_segment.t_start())*(it1->t_end()-it1->t_start())<0) {
      det_rat *= -1;
      overlap = -1;      
    }
  }
    
  if (det_rat < 0) {
    det_rat_sign = -1;
    det_rat *= -1;
  }
  else {
    det_rat_sign = 1;
  }
  
  return det_rat;
}


template <class G, class S> void compute_M_move(times & new_segment, int k, blas_matrix & M, const S& segments_old, G& F, double BETA, double det_rat) {
  
  blas_matrix M_new(M.size1(),M.size1());
  //double argument;
  
  // row k and column k
  double det_rat_inv=1./det_rat;
  for (int i=0; i<M.size1(); i++) {
    if (i!=k) {
      M_new(i,k) = 0;
      M_new(k,i) = 0;
      
      typename S::const_iterator it=segments_old.begin();
      for (int n=0; n<M.size1(); n++) {
        if (n!=k) {
          M_new(i,k) -= det_rat_inv*(M(k,k)*M(i,n)-M(i,k)*M(k,n))*interpolate_F(it->t_end()-new_segment.t_start(), BETA, F);
          M_new(k,i) -= det_rat_inv*(M(k,k)*M(n,i)-M(n,k)*M(k,i))*interpolate_F(new_segment.t_end()-it->t_start(), BETA, F);      
        }
        it++;
      } 
    }
    else {
      M_new(k,k) = M(k,k)*det_rat_inv;
    }
  }
  
  // remaining elements
  for (int i=0; i<(int)M.size1(); i++) {
    if (i!=k) {
      for (int j=0; j<(int)M.size1(); j++) {
        if (j!=k)
          M_new(i,j) = M(i,j) + (-M(i,k)*M(k,j)+det_rat*M_new(i,k)*M_new(k,j))/M(k,k);
      }
    }
  }
    swap(M,M_new);
  //M_new.swap(M);
  return;
}  

// shift end point of segment
template <class G, class S> double det_rat_shift(times & new_segment, int k, blas_matrix & M, S& segments_old, G& F, double BETA, double & det_rat_sign, double & overlap) {
  
  typename S::const_iterator it;
  double det_rat = 0;
  
  it=segments_old.begin();
  for (int i=0; i<(int)M.size1(); i++) {
    det_rat += interpolate_F(new_segment.t_end()-it->t_start(), BETA, F)*M(i,k);
    it++;
  }
  
  overlap = 1;
  // take care of sign changes produced by segments which "wind around"
  if (k==(int)segments_old.size()-1) {
    it--;
    // check if last segment has been shifted across beta
    if ((new_segment.t_end()-new_segment.t_start())*(it->t_end()-it->t_start())<0) {
      det_rat *= -1;
      overlap = -1;      
    }
  }
    
  if (det_rat < 0) {
    det_rat_sign = -1;
    det_rat *= -1;
  }
  else {
    det_rat_sign = 1;
  }
  
  return det_rat;
}


template <class G, class S> void compute_M_shift(times & new_segment, unsigned k, blas_matrix & M, const S& segments_old, G& F, double BETA, double det_rat) {
  std::vector<double> R(M.size1(),0), M_k(M.size1(),0), Fe(M.size1(),0);
  
  typename S::const_iterator it=segments_old.begin();
  for (unsigned i=0; i<M_k.size(); i++) {
    M_k[i] = M(i,k);
    Fe[i] = interpolate_F(new_segment.t_end()-it->t_start(), BETA, F);
    it++;
  }
  
  for (unsigned i=0; i<R.size(); i++) {
    if (i!=k) {
      for (unsigned j=0; j<R.size(); j++)
        R[i] += Fe[j]*M(j,i);  
    } 
  }   
  
  for (unsigned m=0; m<(unsigned)M.size1(); m++) {
    if (m!=k) {
      for (unsigned n=0; n<(unsigned)M.size1(); n++) {
        M(n,m) -= M_k[n]*R[m]/det_rat;
      }
    }
    else {
      for (int n=0; n<M.size1(); n++) {
        M(n,m) = M_k[n]/det_rat;
      }
    }
  }
  
  //std::vector<double> R(M.size1(),0), M_k(M.size1(),0), Fe(M.size1(),0);
  /*blas::vector R(M.size()), M_k(M.size()), Fe(M.size()); 
   typename S::const_iterator it=segments_old.begin();
   for (int i=0; i<(int)M_k.size(); i++) {
   M_k(i) = M(i,k);
   Fe(i) = interpolate_F(new_segment.t_end()-it->t_start(), BETA, F);    
   it++;
   }
   
   M.left_multiply(Fe, R);
   R(k)=1.;
   M.add_outer_product(M_k, R, -1./det_rat);
   */ 
  return;
}  


template <class G, class S, class V> double det_rat_insert_anti(times & anti_segment, blas_matrix & M, S& segments_old, G& F, double BETA, double & det_rat_sign, double & overlap, V& R) {
  
  std::vector<double> F_k(R.size());
  
  typename S::const_iterator it=segments_old.begin();
  for (int i=0; i<(int)F_k.size(); i++) {
    F_k[i]=interpolate_F(anti_segment.t_start()-it->t_start(), BETA, F);
    it++;
  }
  
  double det_rat = -interpolate_F(anti_segment.t_start()-anti_segment.t_end(), BETA, F);
  
  it=segments_old.begin();
  for (int i=0; i<(int)R.size(); i++) {
    R(i)=0;
    for (int l=0; l<(int)R.size(); l++) {  
      R(i) += F_k[l]*M(l,i);
    }
    det_rat += interpolate_F(it->t_end()-anti_segment.t_end(), BETA, F)*R(i);
    it++;
  }
  
  overlap = 1;
  // take care of sign changes produced by segments which "wind around"
  // check if anti-segment winds around
  if (anti_segment.t_end()<anti_segment.t_start()) {
    det_rat *= -1;
    overlap = -1;      
  }
    
  if (det_rat < 0) {
    det_rat_sign = -1;
    det_rat *= -1;
  }
  else {
    det_rat_sign = 1;
  }
  
  return det_rat;
  
}


inline int cycle(int i, int size) {
  return (i>0 ? i-1 : size-1); 
}

template <class G, class S, class V> void compute_M_insert_anti(times & anti_segment, int s, int r, blas_matrix & M, S& segments_old, G& F, double BETA, double det_rat, V& R) {
  
  blas_matrix M_new(M.size1()+1,M.size1()+1);
  //std::vector<double> F_kp1(R.size()), L(R.size());
  blas::vector F_kp1(R.size()), L(R.size());
  
  typename S::const_iterator it=segments_old.begin();
  for (int i=0; i<(int)F_kp1.size(); i++) {
    F_kp1(i)=interpolate_F(it->t_end()-anti_segment.t_end(), BETA, F);
    it++;
  }
  M.right_multiply(F_kp1, L);
  /*for (int i=0; i<(int)L.size(); i++) {
   L[i]=0;
   for (int l=0; l<(int)L.size(); l++) {  
   L[i] += M(i,l)*F_kp1[l];
   }
   }*/
  
  int i_new, j_new;
  int size=(int)M.size1();
  
  // element (k+1,k)
  M_new(r,s) = -1./det_rat;
  
  if (r!=0) { // segments remain in the usual order
    
    // row k+1 and column k
    for (int i=0; i<size; i++) {
      i_new = (i<r ? i : i+1);
      j_new = (i<s ? i : i+1);
      
      M_new(i_new,s) = L(i)/det_rat;
      M_new(r,j_new) = R(i)/det_rat;
    }
    
    // remaining elements
    for (int i=0; i<size; i++) {
      i_new = (i<r ? i : i+1);
      for (int j=0; j<s; j++) {
        M_new(i_new,j) = M(i,j) - L(i)*R(j)/det_rat;
      }
      for (int j=s; j<size; j++) {
        M_new(i_new,j+1) = M(i,j) - L(i)*R(j)/det_rat;
      }
    }
  }
  else { // need to permute indices of R, L, M
    
    // row k+1 and column k
    for (int i=0; i<size; i++) {
      i_new = (i<r ? i : i+1);
      j_new = (i<s ? i : i+1);
      
      M_new(i_new,s) = L(i)/det_rat;
      M_new(r,j_new) = R(cycle(i,size))/det_rat;
    }
    
    // remaining elements
    for (int i=0; i<size; i++) {
      i_new = (i<r ? i : i+1);
      for (int j=0; j<size; j++) {
        j_new = (j<s ? j : j+1);
        M_new(i_new,j_new) = M(i,cycle(j,size)) - L(i)*R(cycle(j,size))/det_rat;
      }
    }  
  }
    swap(M_new,M);
  //M_new.swap(M);
  return;
}

template <class G, class S> double det_rat_remove_anti(times anti_segment, int r, int s, blas_matrix & M,const S& segments_old, G& F, double BETA, double & det_rat_sign) {
  
  // r is the index of the segment which is removed
  // s is the index of the segment which is shifted
  
  typename S::const_iterator it=segments_old.begin();
  typename S::const_iterator its(it), itr(it);
  advance(its, s); 
  advance(itr, r);
  
  double inv_det_rat = -interpolate_F(its->t_end()-itr->t_start(), BETA, F);
  
  for (int i=0; i<(int)segments_old.size(); i++) {
    if (i!=s) {
      inv_det_rat -= interpolate_F(it->t_end()-itr->t_start(), BETA, F)*M(r,i)/M(r,s);
    }
    it++;
  }
  
  // take care of sign changes produced by segments which "wind around"
  if (anti_segment.t_end() < anti_segment.t_start()) {
    inv_det_rat *= -1;
  }
    
  if (inv_det_rat < 0) {
    det_rat_sign = -1;
    inv_det_rat *= -1;
  }
  else {
    det_rat_sign = 1;
  }
  
  return 1/inv_det_rat;
  
}


template<class Mat> void compute_M_remove_anti(Mat & M, int s, int r) {
  
  blas_matrix M_new(M.size1()-1,M.size1()-1);
  int i_old;
  int size=M_new.size1();
  
  if(r!=0) { // order of segments remains unchanged
    for (int i=0; i<size; i++) {
      i_old = (i<r ? i : i+1);
      for (int j=0; j<s; j++) {
        M_new(i,j) = M(i_old,j) - M(i_old, s)*M(r, j)/M(r, s);
      }
      for (int j=s; j<size; j++) {
        M_new(i,j) = M(i_old,j+1) - M(i_old, s)*M(r, j+1)/M(r, s);
      }
    }
  }
  else { // need to permute indices of M
    for (int i=0; i<size; i++) {
      for (int j=0; j<size; j++) {
        M_new(i,cycle(j,size)) = M(i+1,j) - M(i+1, s)*M(r, j)/M(r, s);
      }
    }  
  }
    swap(M,M_new);
  //M_new.swap( M);
  return;
}



#endif
