#include "impurity.h"
#include "update.h"
#include <boost/numeric/bindings/lapack/driver/gesv.hpp>

//these updates should be replaced by the blas updates (dger, dgemv, etc) ASAP

// invert matrix A and calculate its determinant
void invert(alps_matrix & A, double & det) {
  assert(num_rows(A) == num_cols(A));
  using std::fill;
  using alps::numeric::diagonal;
  using boost::numeric::bindings::lapack::gesv;
  if(num_rows(A)==0){ det=0; return;} 
  alps_matrix B(num_rows(A), num_cols(A));
  fill(diagonal(B).first,diagonal(B).second,1.);
  
  alps::numeric::vector<int> ipivot(num_rows(A),0.);
  //LU factorization
  gesv(A,ipivot,B);
  swap(A,B);
  det = 1;
  for (std::size_t i=0; i<num_rows(B); i++) {
    det *= B(i,i);
  }
  det = std::fabs(det);
}

void construct_matrix(alps_matrix & M, segment_container_t & segments, double BETA,  hybridization_t& F) {
  int N = segments.size();
  resize(M,N,N);
  int row=-1;
  int column=-1;
  for (segment_container_t::iterator it1=segments.begin(); it1!=segments.end(); it1++) {
    row++;
    for (segment_container_t::iterator it2=segments.begin(); it2!=segments.end(); it2++) {
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

double construct_inverse(alps_matrix & M, segment_container_t & segments, double BETA,  hybridization_t& F) {
  construct_matrix(M, segments, BETA, F);
  double dummy;
  invert(M, dummy);
  return dummy;
}




// compute distances up/down to the next segment and iterators of these segments
// note: s_down always points to a physical segment, while s_up may point to segments.end() 
void compute_intervals(double t, double BETA, double& t_up, double& t_down, segment_container_t & segments, segment_container_t::iterator& s_up, segment_container_t::iterator& s_down) {  
  
  if (segments.size() == 0) { //no segments there.
    t_up = BETA;
    t_down = BETA;
    s_up = segments.end(); //let all iterators point to the end.
    s_down = segments.end();
  }
  else { //there is at least one segment in this orbital
    
    //s_up = lower_bound(segments.begin(), segments.end(), t); //find the segment that starts just after t. (operator< with double t)
    s_up = segments.lower_bound(times(t,BETA)); //find the segment that starts just after t. (operator< with double t)
    
    if (s_up == segments.begin()) { //there is no segment that starts before t.
      s_down = segments.end(); s_down--; //let s_down point to the last segment
      if (s_down->t_end() < s_down->t_start()) //that last segment wraps around
        t_down = t - s_down->t_end(); //time difference to last segment
      else
        t_down = t + BETA - s_down->t_end(); //time difference to last segment
    }
    else { //we found a segment that is not the first segment
      s_down = s_up; s_down--; //find the segment just before it
      if (s_down->t_end()>s_down->t_start()) //it does not wrap around
        t_down = t - s_down->t_end(); 
      else //down segment wraps around.
        t_down = t - (BETA+s_down->t_end());
    }
    
    if(s_up == segments.end()) { //up segment is the last segement in the list
      t_up = BETA - t + segments.begin()->t_start();  //wrap around to first segment
    }
    else {
      t_up = s_up->t_start() - t; //just difference to the starting time of this line.
    }
    
  }
  /*  for(int k=0;k<other_segments.size();++k){
   typename S::iterator it=other_segments[k].begin();
   while(it !=other_segments[k].end() && it->t_start()<t){
   if(it->t_end() > t){
   t_down=-1; //enforce overlap - code will think that this position is already occupied.
   }
   }
   typename S::reverse_iterator rit=other_segments[k].rbegin(); //check last segment:
   if(rit!=other_segments[k].rend()){ //not empty
   if(rit->t_end() < rit->t_start()){ //segment wrapping around
   if(rit->t_start() < t || rit->t_end() > t){
   t_down = -1; //enforce overlap
   }
   }
   }
   }*/
}

// compute overlap between a segment and a list of segments
// requires segment with 0<=t_begin<t_end<=BETA
double segment_overlap(times segment, segment_container_t& other_segments, int other_full_line, double BETA) {
  
  double length = (segment.t_start()<segment.t_end() ? segment.t_end()-segment.t_start() : segment.t_end()-segment.t_start()+BETA);
  double t_final = segment.t_start()+length;
  double t = segment.t_start();
  double t_final_segment;    
  double other_length=0;
  if (other_full_line==1)
    other_length=length;
  else if (other_segments.size()>0){
    segment_container_t::iterator it;
    it = other_segments.lower_bound(times(t, BETA));  
    //it = lower_bound(other_segments.begin(), other_segments.end(), t);  
    
    if (it!=other_segments.begin()) {
      it--;
      t_final_segment = (it->t_start()<it->t_end() ? it->t_end() : it->t_end()+BETA);
      if (t<t_final_segment) {
        other_length += (t_final_segment<t_final ? t_final_segment-t : t_final-t);
      }
      it++;
      
    }
    while(it!=other_segments.end() && it->t_start()<t_final) {
      t_final_segment = (it->t_start()<it->t_end() ? it->t_end() : it->t_end()+BETA);
      other_length += (t_final_segment<t_final ? t_final_segment-it->t_start() : t_final-it->t_start());
      it++;
    }
    // check if last segment overlaps
    it=other_segments.end();
    it--;
    if (it->t_end()<it->t_start() && t<it->t_end()) {
      other_length += (t_final<it->t_end() ? t_final-t : it->t_end()-t);
    }
  }   
  return other_length;   
}


double compute_overlap(times segment, segment_container_t& other_segments, int other_full_line, double BETA) {
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
double det_rat_up(times & new_segment, alps_matrix const& M, segment_container_t& segments_old, hybridization_t& F, vector_t& Fs, vector_t& Fe, double BETA, double & det_rat_sign, double & overlap) {
  
  segment_container_t::iterator it=segments_old.begin();
  for (std::size_t i=0; i<segments_old.size(); i++) {
    Fe[i] = interpolate_F(new_segment.t_end()-it->t_start(), BETA, F);
    Fs[i] = interpolate_F(it->t_end()-new_segment.t_start(), BETA, F);
    it++;
  }
  
  double det_rat = interpolate_F(new_segment.t_end()-new_segment.t_start(), BETA, F);

  det_rat -= scalar_product(Fe,M*Fs);
  
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

void compute_M_up(std::size_t k, alps_matrix & M, vector_t& Fs, vector_t &Fe, double det_rat) {
  assert( num_rows(M) == num_cols(M) );
  alps_matrix M_new(num_rows(M)+1,num_cols(M)+1);
  std::size_t i_new, j_new;
  
  // element (k,k)
  M_new(k,k) = 1./det_rat;
  
  // row k and column k
  for (std::size_t i=0; i < num_rows(M); i++) {
    i_new = (i<k ? i : i+1);
    M_new(i_new,k) = 0;
    M_new(k,i_new) = 0;
    
    for (std::size_t n=0; n<num_rows(M); n++) {
      M_new(i_new,k) -= M(i,n)*Fs[n];
      M_new(k,i_new) -= M(n,i)*Fe[n];  
    } 
    M_new(i_new,k) /= det_rat;
    M_new(k,i_new) /= det_rat;
  }
  
  // remaining elements
  for (std::size_t j=0; j<num_cols(M); j++) {
    j_new = (j<k ? j : j+1);
    for (std::size_t i=0; i<num_rows(M); i++) {
      i_new = (i<k ? i : i+1);
      M_new(i_new,j_new) = M(i,j) + det_rat*M_new(i_new,k)*M_new(k,j_new);
    }
  }
  
  swap(M_new, M);
  return;
}  


double det_rat_down(std::size_t k, alps_matrix const& M, segment_container_t& segments_old, double & det_rat_sign) {
  
  double det_rat = M(k,k);
  
  // take care of sign changes produced by segments which "wind around"
  if (k==segments_old.size()-1) {
    segment_container_t::iterator it=segments_old.end(); it--;
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


void compute_M_down(std::size_t k, alps_matrix & M) {
  assert(num_rows(M) == num_cols(M));
  assert(num_rows(M) > 0);
  alps_matrix M_new(num_rows(M)-1, num_cols(M)-1);
  
  for (std::size_t j=0; j<num_cols(M_new); j++) {
    std::size_t j_old = (j<k ? j : j+1);
    for (std::size_t i=0; i<num_rows(M_new); i++) {
      std::size_t i_old = (i<k ? i : i+1);
      M_new(i,j) = M(i_old, j_old)-M(i_old,k)*M(k,j_old)/M(k,k);
    }
  }
  
  swap(M, M_new);
  
}

// move segment without changin its length
double det_rat_move(times & new_segment, std::size_t k, alps_matrix const& M, segment_container_t& segments_old, hybridization_t& F, double BETA, double & det_rat_sign, double & overlap) {
  assert(num_rows(M) == num_cols(M));
  double F_i, F_j;
  segment_container_t::iterator it1, it2;
  
  double det_rat = M(k,k)*interpolate_F(new_segment.t_end()-new_segment.t_start(), BETA, F);
  
  it1=segments_old.begin();
  for (std::size_t i=0; i<num_rows(M); i++) {
    if (i != k) {
      F_i = interpolate_F(new_segment.t_end()-it1->t_start(), BETA, F);
      
      it2=segments_old.begin();
      for (std::size_t j=0; j<num_cols(M); j++) {
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


void compute_M_move(times & new_segment, std::size_t k, alps_matrix & M, segment_container_t & segments_old, hybridization_t& F, double BETA, double det_rat) {
  assert(num_rows(M) == num_cols(M)); 
  alps_matrix M_new(num_rows(M),num_cols(M));
  //double argument;
  
  // row k and column k
  for (std::size_t i=0; i<num_rows(M); i++) {
    if (i!=k) {
      M_new(i,k) = 0;
      M_new(k,i) = 0;
      
      segment_container_t::iterator it=segments_old.begin();
      for (std::size_t n=0; n<num_rows(M); n++) {
        if (n!=k) {
          M_new(i,k) -= 1/det_rat*(M(k,k)*M(i,n)-M(i,k)*M(k,n))*interpolate_F(it->t_end()-new_segment.t_start(), BETA, F);
          M_new(k,i) -= 1/det_rat*(M(k,k)*M(n,i)-M(n,k)*M(k,i))*interpolate_F(new_segment.t_end()-it->t_start(), BETA, F);    
        }
        it++;
      } 
    }
    else {
      M_new(k,k) = M(k,k)/det_rat;
    }
  }
  
  // remaining elements
  for (std::size_t j=0; j<num_cols(M); j++) {
    if (j!=k) {
      for (std::size_t i=0; i<num_rows(M); i++) {
        if (i!=k)
          M_new(i,j) = M(i,j) + (-M(i,k)*M(k,j)+det_rat*M_new(i,k)*M_new(k,j))/M(k,k);
      }
    }
  }
  
  swap(M_new, M);
  return;
}  

// shift end point of segment
double det_rat_shift(times & new_segment, std::size_t k, alps_matrix const& M, segment_container_t& segments_old, hybridization_t& F, double BETA, double & det_rat_sign, double & overlap) {
  
  segment_container_t::iterator it;
  double det_rat = 0;
  
  it=segments_old.begin();
  for (std::size_t i=0; i<num_rows(M); i++) {
    det_rat += interpolate_F(new_segment.t_end()-it->t_start(), BETA, F)*M(i,k);
    it++;
  }
  
  overlap = 1;
  // take care of sign changes produced by segments which "wind around"
  if (k==segments_old.size()-1) {
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


void compute_M_shift(times & new_segment, std::size_t k, alps_matrix & M, segment_container_t & segments_old, hybridization_t& F, double BETA, double det_rat) {
 
  assert(num_rows(M) == num_cols(M));
  std::vector<double> R(num_rows(M),0), M_k(num_rows(M),0), Fe(num_rows(M),0);
  
  segment_container_t::iterator it=segments_old.begin();
  for (std::size_t i=0; i<M_k.size(); i++) {
    M_k[i] = M(i,k);
    Fe[i] = interpolate_F(new_segment.t_end()-it->t_start(), BETA, F);  
    it++;
  }
  
  for (std::size_t i=0; i<R.size(); i++) {
    if (i!=k) {
      for (std::size_t j=0; j<R.size(); j++) 
        R[i] += Fe[j]*M(j,i); 
    }
  }
  
  for (std::size_t m=0; m<num_cols(M); m++) {
    if (m!=k) {
      for (std::size_t n=0; n<num_rows(M); n++) {
        M(n,m) -= M_k[n]*R[m]/det_rat;
      }
    }
    else {
      for (std::size_t n=0; n<num_rows(M); n++) {
        M(n,m) = M_k[n]/det_rat;
      }    
    }
  }
  
  return;
}  


double det_rat_insert_anti(times & anti_segment, alps_matrix const& M, segment_container_t& segments_old, hybridization_t& F, double BETA, double & det_rat_sign, double & overlap, vector_t& R) {
  
  std::vector<double> F_k(R.size());
  
  segment_container_t::iterator it=segments_old.begin();
  for (std::size_t i=0; i<F_k.size(); i++) {
    F_k[i]=interpolate_F(anti_segment.t_start()-it->t_start(), BETA, F);
    it++;
  }
  
  double det_rat = -interpolate_F(anti_segment.t_start()-anti_segment.t_end(), BETA, F);
  
  it=segments_old.begin();
  for (std::size_t i=0; i<R.size(); i++) {
    R[i]=0;
    for (std::size_t l=0; l<R.size(); l++) {  
      R[i] += F_k[l]*M(l,i);
    }
    det_rat += interpolate_F(it->t_end()-anti_segment.t_end(), BETA, F)*R[i];
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

void compute_M_insert_anti(times & anti_segment, int s, int r, alps_matrix & M, segment_container_t& segments_old, hybridization_t& F, double BETA, double det_rat, vector_t& R) {
  assert(num_rows(M) == num_cols(M));
  assert(num_rows(M) == R.size());
  alps_matrix M_new(num_rows(M)+1,num_cols(M)+1);
  alps::numeric::vector<double> F_kp1(R.size()), L(R.size());
  
  segment_container_t::iterator it=segments_old.begin();
  for (std::size_t i=0; i<F_kp1.size(); i++) {
    F_kp1[i]=interpolate_F(it->t_end()-anti_segment.t_end(), BETA, F);
    it++;
  }
  
  L = M*F_kp1;
  
  int i_new, j_new;
  int size=num_rows(M);
  
  // element (k+1,k)
  M_new(r,s) = -1./det_rat;
  
  if (r!=0) { // segments remain in the usual order
    
    // row k+1 and column k
    for (int i=0; i<size; i++) {
      i_new = (i<r ? i : i+1);
      j_new = (i<s ? i : i+1);
      
      M_new(i_new,s) = L[i]/det_rat;
      M_new(r,j_new) = R[i]/det_rat;
    }
    
    // remaining elements
    for (int j=0; j<size; j++) {
      j_new = (j<s ? j : j+1);
      for (int i=0; i<size; i++) {
        i_new = (i<r ? i : i+1);
        M_new(i_new,j_new) = M(i,j) - L[i]*R[j]/det_rat;
      }
    }
  }
  else { // need to permute indices of R, L, M
    
    // row k+1 and column k
    for (int i=0; i<size; i++) {
      i_new = (i<r ? i : i+1);
      j_new = (i<s ? i : i+1);
      
      M_new(i_new,s) = L[i]/det_rat;
      M_new(r,j_new) = R[cycle(i,size)]/det_rat;
    }
    
    // remaining elements
    for (int j=0; j<size; j++) {
      j_new = (j<s ? j : j+1);
      for (int i=0; i<size; i++) {
        i_new = (i<r ? i : i+1);
        M_new(i_new,j_new) = M(i,cycle(j,size)) - L[i]*R[cycle(j,size)]/det_rat;
      }
    }  
  }
  
  swap(M_new, M);
  return;
}

double det_rat_remove_anti(times anti_segment, int r, std::size_t s, alps_matrix const& M, segment_container_t& segments_old, hybridization_t& F, double BETA, double & det_rat_sign) {
  
  // r is the index of the segment which is removed
  // s is the index of the segment which is shifted
  
  segment_container_t::iterator it=segments_old.begin();
  segment_container_t::iterator its(it), itr(it);
  advance(its, s); 
  advance(itr, r);
  
  double inv_det_rat = -interpolate_F(its->t_end()-itr->t_start(), BETA, F);
  
  for (std::size_t i=0; i<segments_old.size(); i++) {
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


void compute_M_remove_anti(alps_matrix & M, int s, int r) {
  assert(num_rows(M) == num_cols(M));
  assert(num_rows(M) > 0); 
  alps_matrix M_new(num_rows(M)-1,num_cols(M)-1);
  
  int size=num_rows(M_new);
  
  if(r!=0) { // order of segments remains unchanged
    for (int j=0; j<size; j++) {
      int j_old = (j<s ? j : j+1);
      for (int i=0; i<size; i++) {
        int i_old = (i<r ? i : i+1);
        M_new(i,j) = M(i_old,j_old) - M(i_old, s)*M(r, j_old)/M(r, s);
      }
    }
  }
  else { // need to permute indices of M
    for (int j=0; j<size; j++) {
      for (int i=0; i<size; i++) {
        M_new(i,cycle(j,size)) = M(i+1,j) - M(i+1, s)*M(r, j)/M(r, s);
      }
    }  
  }
  
  swap(M_new, M);
  return;
}


double get_occupation(segment_container_t &segments, int full_line, double tau, double BETA){
  //get occupation for given flavor and time tau:
  //check if there is a segment at this time
  if(segments.size()==0){
    if(full_line) return 1.0;
    else return 0.0;
  }
  for(segment_container_t::iterator it=segments.begin(); it!=segments.end(); it++){
    if(it->t_end()>it->t_start()){//regular segment
	    if(it->t_start()<=tau && tau <= it->t_end() ) return 1.0;
    }//brackets mandatory
    else//segment winds around the circle
	    if(( tau>=0.0 && tau<=it->t_end()) || (tau>=it->t_start() && tau<=BETA)) return 1.0;
  }//end::for
  return 0.0;
}

