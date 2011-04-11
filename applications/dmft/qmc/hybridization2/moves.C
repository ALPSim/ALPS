#include "impurity.h"
#include "update.h"
#include "moves.h"

void insert_remove_full_line(prng_t& rng, double mu, double BETA, blas_matrix &u, int& full_line, std::vector<segment_container_t>& other_segments, std::vector<int>& other_full_line, int this_flavor) {
  
  int insert = (rng() < 0.5);
  
  if ((insert==1 && full_line==1) || (insert==0 && full_line==0)) return; // insert=1(0) means we want to insert(remove) a full line
  
  int FLAVOR = other_full_line.size();
  
  double otherlength_u=0;
  for (int i=0; i<FLAVOR; i++) {
    if (i==this_flavor) continue;
    
    double other_length=0;
    for (segment_container_t::iterator it=other_segments[i].begin(); it!=other_segments[i].end(); it++)
      other_length += (it->t_end()-it->t_start()>0 ? it->t_end()-it->t_start() : it->t_end()-it->t_start()+BETA);
    
    if (other_full_line[i]==1)
      other_length = BETA;
    
    otherlength_u += other_length*u(i, this_flavor);
    
  }
  
  if (insert) { // try to insert full line
    if (log(rng()) < BETA*mu-otherlength_u)
      full_line = 1;
  }
  else { // try to remove full line
    if (log(rng()) < -BETA*mu+otherlength_u)
      full_line = 0;  
  }
  
}


void insert_remove_segment(prng_t& rng, double t, double BETA, double mu, blas_matrix &u, hybridization_t& F, segment_container_t& segments, blas_matrix & M, double & sign, std::vector<segment_container_t>& other_segments, std::vector<int> other_full_line, int this_flavor) {
  
  double t_up; // distance to next segment up
  double t_down; // distance to next segment down
  segment_container_t::iterator s_up; // iterator of the segment up
  segment_container_t::iterator s_down; // iterator of the segment down
  
  if (rng()<0.5) { // try to insert a segment
    compute_intervals(t, BETA, t_up, t_down, segments,s_up, s_down);
    
    if (t_down>0) { // t does not lie on a segment -> it's possible to insert a new one starting from t
      
      double length = compute_length(rng(), t_up, 0);
      
      times segment_insert;
      segment_insert.set_t_start(t);
      double t_final = t + length;
      if (t_final > BETA)
        segment_insert.set_t_end(t_final-BETA);
      else
        segment_insert.set_t_end(t_final);
      
      double otherlength_u=0;
      int FLAVORS=other_full_line.size();
      for (int i=0; i<FLAVORS; i++) {
        if (i==this_flavor) continue;
        double other_length = compute_overlap(segment_insert, other_segments[i], other_full_line[i], BETA);
        otherlength_u += other_length*u(i, this_flavor);
      }
      double log_prob, overlap, det_rat, det_rat_sign;
      std::vector<double> Fs(segments.size()), Fe(segments.size());
      
      det_rat = det_rat_up(segment_insert, M, segments, F, Fs, Fe, BETA, det_rat_sign, overlap);
      
      log_prob = log(BETA*t_up/(segments.size()+1)*det_rat)+mu*length-otherlength_u;
      
      if (log(rng()) < log_prob) {
        int position=0;
        for (segment_container_t::iterator it=segments.begin(); it!=s_up; it++)
          position++;
        compute_M_up(position, M, Fs, Fe, det_rat*overlap);
        sign *= det_rat_sign;
        segment_container_t::iterator sit=segments.insert(s_up, segment_insert);
        if(sit==segments.end()){std::cerr<<"segment could not be inserted! exiting."<<std::endl;}
      }
    }
  }
  
  else if (segments.size()>0) { // try to remove a segment
    int position = (int)(rng()*segments.size());
    s_down = segments.begin();
    for (int i=0; i<position; i++)
      s_down++;
    s_up=s_down;
    s_up++;
    if (s_up==segments.end())
      s_up = segments.begin();
    
    double length = s_down->t_end()-s_down->t_start();
    if (length < 0) length += BETA;
    
    double t_total = s_up->t_start()-s_down->t_start();
    if (t_total <= 0) t_total += BETA;
    
    times segment_remove = *s_down;
    
    double otherlength_u=0;
    int FLAVORS=other_full_line.size();
    for (int i=0; i<FLAVORS; i++) { 
      if (i==this_flavor) continue;
      double other_length = compute_overlap(segment_remove, other_segments[i], other_full_line[i], BETA);
      otherlength_u += other_length*u(i, this_flavor);
    }
    
    double log_prob, det_rat, det_rat_sign;
    
    det_rat = det_rat_down(position, M, segments, det_rat_sign);    
    
    log_prob = log(BETA*t_total/segments.size()/det_rat)+length*mu-otherlength_u;
    
    if (log(rng()) < -log_prob) {
      compute_M_down(position, M);  
      sign *= det_rat_sign;
      segments.erase(s_down);  
    }
  }
}


void insert_remove_antisegment(prng_t& rng, double t, double BETA, double mu,blas_matrix &u, hybridization_t& F, int& full_line, segment_container_t& segments, blas_matrix& M, double & sign,std::vector<segment_container_t>& other_segments, std::vector<int> other_full_line, int this_flavor) {
  
  double t_up; // distance to next segment up (t_start)
  double t_down; // distance to next segment down (t_end)
  segment_container_t::iterator s_up; // iterator of the segment up
  segment_container_t::iterator s_down; // iterator of the segment down
  
  if (rng()<0.5) { // try to insert an anti-segment
    
    if (full_line==1) {
      t_down = -BETA;
      double length = compute_length(rng(), BETA, 0);
      double t_end = (t+length < BETA ? t+length : t+length-BETA);
      times segment_insert(t_end, t);
      times segment_remove(t,t_end);
      
      double log_prob, overlap, det_rat, det_rat_sign;
      std::vector<double> Fs(segments.size()), Fe(segments.size());
      det_rat = det_rat_up(segment_insert, M, segments, F, Fs, Fe, BETA, det_rat_sign, overlap);
      
      double otherlength_u=0;
      int FLAVORS=other_full_line.size();
      for (int i=0; i<FLAVORS; i++) {
        if (i==this_flavor) continue;
        double other_length = compute_overlap(segment_remove, other_segments[i], other_full_line[i], BETA);
        otherlength_u += other_length*u(i, this_flavor);
      }
      log_prob = log(BETA*BETA*det_rat)-length*mu+otherlength_u;
      
      if (log(rng()) < log_prob) {
        compute_M_up(0, M, Fs, Fe, det_rat*overlap);
        sign *= det_rat_sign;
        //segments.push_back(segment_insert);
        segment_container_t::iterator sit;
        //std::cout<<"segment size: "<<segments.size()<<std::endl;
        sit=segments.insert(segment_insert).first;
        if(sit==segments.end()){std::cerr<<"segment could not be inserted! exiting."<<std::endl;}
        full_line = 0;
        //std::cout<<"insertion completed."<<std::endl;
      }
    }
    
    else {
      compute_intervals(t, BETA, t_up, t_down, segments, s_up, s_down);
      
      if (t_down<0) { // t does lie on a segment -> it's possible to insert an anti-segment starting from t
        
        double length = compute_length(rng(), -t_down, 0);
        
        times segment_shrink(s_down->t_start(),t);
        
        double t_start = t + length;
        if (t_start > BETA)
          t_start-=BETA;
        
        times segment_insert(t_start, s_down->t_end());
        times anti_segment(t,t_start);
        
        double otherlength_u=0;
        int FLAVORS=other_full_line.size();
        for (int i=0; i<FLAVORS; i++) {
          if (i==this_flavor) continue;
          double other_length = compute_overlap(anti_segment, other_segments[i], other_full_line[i], BETA);
          otherlength_u += other_length*u(i, this_flavor);
        }
        double log_prob, overlap, det_rat, det_rat_sign;
        std::vector<double> R(segments.size());
        det_rat = det_rat_insert_anti(anti_segment, M, segments, F, BETA, det_rat_sign, overlap, R);
        
        log_prob = log(BETA*(-t_down)/(segments.size()+1)*det_rat)-length*mu+otherlength_u;
        
        if (log(rng()) < log_prob) {
          
          int s, r; // s is the segment which is shifted, r the segment which is inserted
          s = 0;
          for (segment_container_t::iterator it=segments.begin(); it!=s_down; it++)
            s++;
          if (anti_segment.t_end() > segment_shrink.t_start())
            r = s+1;
          else {
            r = 0;  
            s++;
          }
          
          compute_M_insert_anti(anti_segment, s, r, M, segments, F, BETA, det_rat*overlap, R);    
          //s_down->set_t_end(t);
          times segment_new_endpoint(*s_down);
          segment_container_t::iterator prev_segment=s_down;
          if(s_down !=segments.begin()) prev_segment--;
          else prev_segment=segments.begin();
          segment_new_endpoint.set_t_end(t);
          segments.erase(s_down); //erase old segment (without shifted end
          s_down=segments.insert(segment_new_endpoint).first; //in
          //s_down=segments.insert(prev_segment, segment_new_endpoint); //in
          if(s_down==segments.end()){std::cerr<<"segment could not be inserted! exiting."<<std::endl;}
          segment_container_t::iterator sit=segments.insert(segment_insert).first; //insert  new segment
          //typename S::iterator sit=segments.insert(s_down, segment_insert); //insert  new segment
          if(sit==segments.end()){std::cerr<<"segment could not be inserted! exiting."<<std::endl;}
          if (segment_insert.t_start()>segments.begin()->t_start()) {
            s_down++;
            segment_container_t::iterator sit=segments.insert(s_down, segment_insert);
            if(sit==segments.end()){std::cerr<<"segment could not be inserted! exiting."<<std::endl;}
          }
          else {
            segment_container_t::iterator sit=segments.insert(segments.begin(), segment_insert);    
            if(sit==segments.end()){std::cerr<<"segment could not be inserted! exiting."<<std::endl;}
          }
        }
      }
    }
    
  }
  else if (segments.size()>1) { // try to remove an anti-segment
    int r = (int)(rng()*segments.size());
    s_up = segments.begin();
    for (int i=0; i<r; i++) s_up++;
    
    int s = r-1;
    if (s<0) {
      s=segments.size()-1;
      s_down=segments.end();
      s_down--;    
    }
    else {
      s_down = s_up;
      s_down--;
    }
    
    double length = s_up->t_start() - s_down->t_end();
    if (length < 0) length += BETA;
    
    double t_total = s_up->t_end() - s_down->t_end();
    if (t_total < 0) t_total += BETA;    
    
    times anti_segment(s_down->t_end(),s_up->t_start());
    
    double otherlength_u=0;
    int FLAVORS=other_full_line.size();
    for (int i=0; i<FLAVORS; i++) {
      if (i==this_flavor) continue;
      double other_length = compute_overlap(anti_segment, other_segments[i], other_full_line[i], BETA);
      otherlength_u += other_length*u(i, this_flavor);
    }
    double log_prob, det_rat, det_rat_sign;
    
    det_rat = det_rat_remove_anti(anti_segment, r, s, M, segments, F, BETA, det_rat_sign);
    
    log_prob = log(BETA*t_total/segments.size()/det_rat)-length*mu+otherlength_u;
    
    if (log(rng()) < -log_prob) {
      
      compute_M_remove_anti(M, s, r);
      
      double t_end = s_up->t_end();    
      segments.erase(s_up);
      
      if (r>0) {
        s_up=segments.begin();
        for (int k=0; k<s; k++)
          s_up++;
      }
      else {
        s=segments.size()-1;
        s_up = segments.end();
        s_up--;
      }
      //s_up->set_t_end(t_end);
      times s_up_new(*s_up);
      s_up_new.set_t_end(t_end);
      segments.erase(s_up);
      segment_container_t::iterator sit=segments.insert(s_up_new).first;
      if(sit==segments.end()){std::cerr<<"segment could not be inserted! exiting."<<std::endl;}
    }
  }
  
  else if (segments.size()==1) {
    
    s_down = segments.begin();
    
    double det_rat = fabs(M(0,0));
    double length = s_down->t_start()-s_down->t_end(); 
    if (length<0) length += BETA;
    times anti_segment(s_down->t_end(),s_down->t_start());
    
    double otherlength_u=0;
    int FLAVORS=other_full_line.size();
    for (int i=0; i<FLAVORS; i++) {
      if (i==this_flavor) continue;
      double other_length = compute_overlap(anti_segment, other_segments[i], other_full_line[i], BETA);
      otherlength_u += other_length*u(i, this_flavor);
    }
    double log_prob = log(BETA*BETA/det_rat)-length*mu+otherlength_u;
    
    if (log(rng()) < -log_prob) { 
      full_line=1;
      segments.erase(s_down);
      compute_M_down(0,M); // attention: M.clear() sets elements to zero 
    }  
  }
}

/*
 // move segment without changing its length
 template <class RNG, class S, class G> void move_segment(RNG& rng, S& segments, int N, double BETA, double mu, double u, G& F, blas_matrix& M, double & sign, S& other_segments, int other_full_line) {
 
 int size = segments.size();
 
 if (size < 1) return;
 
 int n = size*rng();
 
 segment_container_t::iterator s=segments.begin();
 double new_t_start, new_t_end, length_segment;
 
 if (size==1) {
 new_t_start=rng()*BETA;
 length_segment = s->t_end()-s->t_start();
 if (length_segment<0)
 length_segment += BETA;    
 }
 else {
 
 double t_down, max_lower_gap; 
 segment_container_t::iterator s_up, s_down;
 for (int i=0; i<n; i++) s++;
 
 s_up = s; s_up++;
 if (s_up == segments.end()) s_up = segments.begin();
 
 if (n == 0) 
 s_down = segments.end();
 else 
 s_down=s;
 s_down--;
 
 if (s_down->t_end()<s->t_start()) 
 t_down = s_down->t_end();
 else
 t_down = 0; // prevent segments from winding around, otherwise the matrix gets messed up
 
 double interval = s_up->t_start() - t_down;
 if (interval < 0)
 interval += BETA;
 
 length_segment = s->t_end()-s->t_start();
 if (length_segment<0)
 length_segment += BETA;
 
 max_lower_gap = interval-length_segment;
 if (t_down+max_lower_gap>BETA)
 max_lower_gap = BETA-t_down; // prevent segments from winding around
 
 new_t_start = t_down + compute_length(rng(), max_lower_gap, 0);
 
 }
 
 new_t_end = new_t_start + length_segment;
 if (new_t_end > BETA)
 new_t_end -= BETA;
 
 times segment_insert(new_t_start, new_t_end);
 times segment_remove=*s;
 double delta = compute_overlap(segment_insert, other_segments, other_full_line, BETA)-compute_overlap(segment_remove, other_segments, other_full_line, BETA);
 
 double det_rat, det_rat_sign, overlap;
 
 det_rat = det_rat_move(segment_insert, n, M, segments, F, BETA, det_rat_sign, overlap);  
 
 if (log(rng()) < log(det_rat)-delta*u) {
 
 compute_M_move(segment_insert, n, M, segments, F, BETA, det_rat*overlap);  
 
 sign *= det_rat_sign;
 s->set_t_start(new_t_start);  
 s->set_t_end(new_t_end);  
 }
 
 }
 */


// shift segment
void shift_segment(prng_t& rng, segment_container_t& segments, double BETA, double mu, blas_matrix &u, hybridization_t& F, blas_matrix & M, double & sign, std::vector<segment_container_t>& other_segments, std::vector<int>& other_full_line, int this_flavor) {
  
  int size = segments.size();
  
  if (size < 1) return;
  
  int n = (int)(size*rng());
  
  segment_container_t::iterator s, s_up;
  s=segments.begin();
  for (int i=0; i<n; i++) s++;
  s_up = s; s_up++;
  if (s_up == segments.end()) s_up = segments.begin();
  
  double interval = s_up->t_start() - s->t_start();
  if (interval <= 0) interval += BETA;
  
  double length = compute_length(rng(), interval, 0);
  double length_old = s->t_end()-s->t_start();
  if (length_old<0)
    length_old += BETA;
  
  double new_t_end = s->t_start() + length;
  if (new_t_end > BETA)
    new_t_end -= BETA;
  
  times segment_insert(s->t_start(), new_t_end);
  times segment_remove=*s;
  
  double otherlength_u=0;
  int FLAVORS=other_full_line.size();
  for (int i=0; i<FLAVORS; i++) {
    if (i==this_flavor) continue;
    double other_length = compute_overlap(segment_insert, other_segments[i], other_full_line[i], BETA)-compute_overlap(segment_remove, other_segments[i], other_full_line[i], BETA);
    otherlength_u += other_length*u(i, this_flavor);
  }
  double det_rat, det_rat_sign, overlap;
  
  det_rat = det_rat_shift(segment_insert, n, M, segments, F, BETA, det_rat_sign, overlap);
  
  if (log(rng()) < log(det_rat)+(length-length_old)*mu-otherlength_u) {
    
    compute_M_shift(segment_insert, n, M, segments, F, BETA, det_rat*overlap);  
    sign *= det_rat_sign;
    //s->set_t_end(new_t_end);  
    times s_new(*s);
    s_new.set_t_end(new_t_end);     
    segments.erase(s);
    segment_container_t::iterator sit;
    sit=segments.insert(s_new).first;
    if(sit==segments.end()){std::cerr<<"segment could not be inserted! exiting."<<std::endl;}
  }
}

// swap segment configurations
void swap_segments(prng_t& rng, double BETA, hybridization_t& F_up, hybridization_t & F_down, segment_container_t& segments_up, segment_container_t& segments_down, int& full_line_up, int& full_line_down, double & sign_up, double & sign_down, blas_matrix & M_up, blas_matrix& M_down) {
  
  blas_matrix M_new_up, M_new_down;
  
  // before swap
  double det_old_up = construct_inverse(M_new_up, segments_up, BETA,  F_up); // here M_new_up is just a dummy
  double det_old_down = construct_inverse(M_new_down, segments_down, BETA,  F_down); // here M_new_down is just a dummy
  // before swap
  double det_new_up = construct_inverse(M_new_up, segments_down, BETA,  F_up); 
  double det_new_down = construct_inverse(M_new_down, segments_up, BETA,  F_down);
  
  double det_rat = (det_new_up/det_old_up)*(det_new_down/det_old_down);
  
  // length of segments, overlap and phonon part are not changed
  if (rng() < fabs(det_rat)) {
    
    //std::cout << "success\n";
    swap(M_new_up, M_up);
    swap(M_new_down, M_down);
    swap(segments_up, segments_down);
    std::swap(full_line_up, full_line_down);
    std::swap(sign_up, sign_down);
  }
}

/*
 // flip segment between up/down
 template <class RNG, class S, class G> void flip_segment(RNG& rng, S& segments, int N, double BETA, blas_matrix& M, double & sign, double & other_sign, G& other_F, blas_matrix& other_M, S& other_segments, int other_full_line) {
 
 int size = segments.size();
 
 if (size < 1) return;
 
 int position = size*rng();
 
 typename S::iterator s;
 s=segments.begin();
 for (int i=0; i<position; i++) s++;
 
 if (compute_overlap(*s, other_segments, other_full_line, BETA)>0) 
 return;
 
 
 // no overlap -> can flip segment (mu and u contributions don't change)
 
 double det_rat_remove, det_rat_remove_sign, det_rat_insert, det_rat_insert_sign, overlap;
 
 det_rat_remove = det_rat_down(position, M, segments, det_rat_remove_sign);    
 
 std::vector<double> other_Fs(other_segments.size()), other_Fe(other_segments.size());  
 det_rat_insert = det_rat_up(*s, other_M, other_segments, other_F, other_Fs, other_Fe, BETA, det_rat_insert_sign, overlap);
 
 if (rng() < det_rat_remove*det_rat_insert*segments.size()/(other_segments.size()+1)) {
 
 typename S::iterator s_up, s_down;
 double t_up, t_down;
 compute_intervals(s->t_start(), BETA, t_up, t_down, other_segments, s_up, s_down);
 
 int n=0;
 for (typename S::iterator it=other_segments.begin(); it!=s_up; it++)
 n++;
 compute_M_up(*s, n, other_M, other_segments, other_F, other_Fs, other_Fe, BETA, det_rat_insert*overlap);
 other_segments.insert(s_up, *s);  
 
 compute_M_down(position, M);  
 segments.erase(s);  
 
 sign *= det_rat_remove_sign;
 other_sign *= det_rat_insert_sign;  
 
 }  
 
 }
 */

void measure_GF_imp(prng_t& rng, int N, double BETA, double mu, blas_matrix &u, hybridization_t& F, segment_container_t& segments, blas_matrix & M, std::vector<segment_container_t>& other_segments, std::vector<int> other_full_line, int this_flavor, std::valarray<double> & G_meas) {
  
  if (segments.size()==0) {
    
    double length = rng()*BETA;
    double t=rng()*BETA;
    double bubble_sign;
    
    times segment_insert;
    segment_insert.set_t_start(t);
    double t_final = t + length;
    if (t_final > BETA)
      segment_insert.set_t_end(t_final-BETA);
    else
      segment_insert.set_t_end(t_final);
    
    times full_segment;
    full_segment.set_t_start(0.);
    full_segment.set_t_end(1.);
    
    
    double otherlength_u=0;
    double full_other_length_u=0;
    int FLAVORS=other_full_line.size();
    for (int i=0; i<FLAVORS; i++) {
      if (i==this_flavor) continue;
      double other_length = compute_overlap(segment_insert, other_segments[i], other_full_line[i], BETA);
      otherlength_u += other_length*u(i, this_flavor);
      double full_other_length = compute_overlap(full_segment, other_segments[i], other_full_line[i], BETA);
      full_other_length_u += full_other_length*u(i, this_flavor);
    }
    double prob, overlap, det_rat, det_rat_sign;
    std::vector<double> Fs(segments.size()), Fe(segments.size());
    
    // note: det_rat>0 by def. 
    det_rat = det_rat_up(segment_insert, M, segments, F, Fs, Fe, BETA, det_rat_sign, overlap);
    
    //prob = BETA*BETA/1.*det_rat*exp(mu*length-otherlength_u);
    prob = BETA*BETA/(1.+exp(mu*BETA-full_other_length_u))*det_rat*exp(mu*length-otherlength_u);
    
    //if (other_full_line[this_flavor]==0) {
    //  prob = 0.5*BETA*BETA/1.*det_rat*exp(mu*length-otherlength_u);
    //}
    //else {
    //  prob = 0.5*BETA*BETA/1.*det_rat*exp(mu*length-otherlength_u-mu*BETA+full_other_length_u); 
    //}      
    if (prob>1) prob=1;
    
    int index = (int)(length/BETA*N+0.5);
    //G_meas[this_flavor*(N+1)+index] += 1./(1+exp(mu*BETA-full_other_length_u))*prob/det_rat; //*det_rat_sign*bubble_sign;        
    G_meas[this_flavor*(N+1)+index] += prob/det_rat;
    
    //G_meas[this_flavor*(N+1)+index] += prob/det_rat;
    
  }
  
  else if (segments.size()==1) {
    
    double length = segments.begin()->t_end()-segments.begin()->t_start();
    double bubble_sign=1;
    if (length<0) { 
      bubble_sign = -1;
      length += BETA;
    }
    
    times full_segment;
    full_segment.set_t_start(0.);
    full_segment.set_t_end(1.);
    
    double otherlength_u=0;
    double full_other_length_u=0;
    int FLAVORS=other_full_line.size();
    for (int i=0; i<FLAVORS; i++) {
      if (i==this_flavor) continue;
      double other_length = compute_overlap(*(segments.begin()), other_segments[i], other_full_line[i], BETA);
      otherlength_u += other_length*u(i,this_flavor);
      double full_other_length = compute_overlap(full_segment, other_segments[i], other_full_line[i], BETA);
      full_other_length_u += full_other_length*u(i, this_flavor);      
    }
    
    //double prob = exp(otherlength_u-mu*length)*M(0,0)*bubble_sign;
    double prob = (1.+exp(mu*BETA-full_other_length_u))/exp(mu*length-otherlength_u)*M(0,0)*bubble_sign;    
    if (prob>1) prob=1;
    int index = (int)(length/BETA*N+0.5);
    G_meas[this_flavor*(N+1)+index] += (1-prob)*M(0,0)*bubble_sign;    
    
  } 
  
  else {
    
    segment_container_t::iterator it1, it2;
    
    for (int i=0; i<M.size1(); i++) {
      (i==0 ? it1 = segments.begin() : it1++);
      for (int k=0; k<M.size1(); k++) {
        (k==0 ? it2 = segments.begin() : it2++);
        if (M(k,i)!=0) {
          double argument = it1->t_end()-it2->t_start();
          double bubble_sign=1;
          if (argument > 0) {
            bubble_sign = 1;
          }
          else {
            bubble_sign = -1;
            argument += BETA;
          }
          
          int index = (int)(argument/BETA*N+0.5);
          G_meas[this_flavor*(N+1)+index] += M(k,i)*bubble_sign;
        }
      }
    }
  }
  
}

