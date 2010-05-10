/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Emanuel Gull <gull@phys.columbia.edu>
 *
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
#ifndef ___MOVES___
#define ___MOVES___

#include "impurity.h"
#include "update.h"

// length of inserted segment (mu>0) or anti-segment (mu<0)
inline double compute_length(double r, double l_max, double mu)
{
  
  if (mu == 0)
    return r * l_max;
  else
    return 1 / mu * log(r * (exp(mu * l_max) - 1) + 1);
}
void print_segments(const std::vector<segment_container_t> &segments);

template < class RNG, class S > void insert_remove_full_line(RNG & rng, double mu, U_matrix & u, double BETA, int &full_line,
                                                             std::vector < S > &other_segments, std::vector < int >&other_full_line,
                                                             int this_flavor)
{
  
  int insert = (rng() < 0.5);
  //std::cout<<"insert remove full line. "<<rng()<<std::endl;
  //print_segments(other_segments);
  if ((insert == 1 && full_line == 1) || (insert == 0 && full_line == 0))
    return;                     // insert=1(0) means we want to insert(remove) a full line
  
  int FLAVOR = other_full_line.size();
  
  double otherlength_u = 0;
  for (int i = 0; i < FLAVOR; i++) {
    if (i == this_flavor)
      continue;
    
    double other_length = 0;
    for (typename S::const_iterator it = other_segments[i].begin(); it != other_segments[i].end(); it++)
      other_length += (it->t_end() - it->t_start() > 0 ? it->t_end() - it->t_start() : it->t_end() - it->t_start() + BETA);
    
    if (other_full_line[i] == 1)
      other_length = BETA;
    
    otherlength_u += other_length * u(i, this_flavor);
    
  }
 
  if (insert) { // try to insert full line
    if (!full_line) {
      if (log(rng()) < BETA*mu-other_length_u) {
        full_line = 1;
      }
    }
  }
  else { // try to remove full line
    if (full_line) {
      if (log(rng()) < -BETA*mu+other_length_u) {
        full_line = 0;  
      }
    }
  }
}


template < class RNG, class S, class G > void insert_remove_segment(RNG & rng, double t, double BETA, double mu, U_matrix & u, G & F,
                                                                    S & segments, blas_matrix & M, double &sign,
                                                                    std::vector < S > &other_segments, std::vector < int >other_full_line,
                                                                    int this_flavor)
{
  
  double t_up;                  // distance to next segment up
  double t_down;                // distance to next segment down
  typename S::iterator s_up;    // iterator of the segment up
  typename S::iterator s_down;  // iterator of the segment down
  
  if (rng() < 0.5) {            // try to insert a segment
    //std::cout<<"insert segments"<<std::endl;
    compute_intervals(t, BETA, t_up, t_down, segments, s_up, s_down);
    
    if (t_down > 0) {           // t does not lie on a segment -> it's possible to insert a new one starting from t
      
      double length = compute_length(rng(), t_up, 0);
      
      times segment_insert;     //the segment we want to insert
      segment_insert.set_t_start(t);
      double t_final = t + length;
      if (t_final > BETA)
        segment_insert.set_t_end(t_final - BETA);
      else
        segment_insert.set_t_end(t_final);
      
      double otherlength_u = 0;
      int FLAVORS = other_full_line.size();
      for (int i = 0; i < FLAVORS; i++) {
        if (i == this_flavor)
          continue;
        double other_length = compute_overlap(segment_insert, other_segments[i], other_full_line[i], BETA);
        otherlength_u += other_length * u(i, this_flavor);
      }
      double log_prob, overlap, det_rat, det_rat_sign;
      blas::vector Fe(segments.size()), Fs(segments.size());
      
      det_rat = det_rat_up(segment_insert, M, segments, F, Fe, Fs, BETA, det_rat_sign, overlap);
      
      log_prob = log(BETA * t_up / (segments.size() + 1) * det_rat) + mu * length - otherlength_u;
      //std::cout<<"insert segments: "<<rng()<<" log prob: "<<log_prob<<std::endl;
      //print_segments(other_segments);
      
      if (log(rng()) < log_prob) {
        int position = 0;
        for (typename S::const_iterator it = segments.begin(); it != s_up; it++)
          position++;
        compute_M_up(position, M, Fs, Fe, det_rat * overlap);
        sign *= det_rat_sign;
        segments.insert(s_up, segment_insert);
      }
    }
  }
  
  else if (segments.size() > 0) {       // try to remove a segment
    int position = (int) (rng() * segments.size());
    s_down = segments.begin();
    for (int i = 0; i < position; i++)
      s_down++;
    s_up = s_down;
    s_up++;
    if (s_up == segments.end())
      s_up = segments.begin();
    
    double length = s_down->t_end() - s_down->t_start();
    if (length < 0)
      length += BETA;
    
    double t_total = s_up->t_start() - s_down->t_start();
    if (t_total <= 0)
      t_total += BETA;
    
    times segment_remove = *s_down;
    
    double otherlength_u = 0;
    int FLAVORS = other_full_line.size();
    for (int i = 0; i < FLAVORS; i++) {
      if (i == this_flavor)
        continue;
      double other_length = compute_overlap(segment_remove, other_segments[i], other_full_line[i], BETA);
      otherlength_u += other_length * u(i, this_flavor);
    }
    
    double log_prob, det_rat, det_rat_sign;
    
    det_rat = det_rat_down(position, M, segments, det_rat_sign);
    
    log_prob = log(BETA * t_total / segments.size() / det_rat) + length * mu - otherlength_u;
    
    if (log(rng()) < -log_prob) {
      compute_M_down(position, M);
      sign *= det_rat_sign;
      segments.erase(s_down);
    }
  }
}


template < class RNG, class S, class G > void insert_remove_antisegment(RNG & rng, double t, double BETA, double mu, U_matrix & u, G & F,
                                                                        int &full_line, S & segments, blas_matrix & M, double &sign,
                                                                        std::vector < S > &other_segments,
                                                                        std::vector < int >&other_full_line, int this_flavor)
{
  
  double t_up;                  // distance to next segment up (t_start)
  double t_down;                // distance to next segment down (t_end)
  typename S::iterator s_up;    // iterator of the segment up
  typename S::iterator s_down;  // iterator of the segment down
  
  if (rng() < 0.5) {            // try to insert an anti-segment
    
    if (full_line == 1) {
      t_down = -BETA;
      double length = compute_length(rng(), BETA, 0);
      double t_end = (t + length < BETA ? t + length : t + length - BETA);
      times segment_insert(t_end, t);
      times segment_remove(t, t_end);
      
      double log_prob, overlap, det_rat, det_rat_sign;
      blas::vector Fe(segments.size()), Fs(segments.size());
      det_rat = det_rat_up(segment_insert, M, segments, F, Fe, Fs, BETA, det_rat_sign, overlap);
      
      double otherlength_u = 0;
      int FLAVORS = other_full_line.size();
      for (int i = 0; i < FLAVORS; i++) {
        if (i == this_flavor)
          continue;
        double other_length = compute_overlap(segment_remove, other_segments[i], other_full_line[i], BETA);
        otherlength_u += other_length * u(i, this_flavor);
      }
      
      log_prob = log(BETA * BETA * det_rat) - length * mu + otherlength_u;
      
      if (log(rng()) < log_prob) {
        compute_M_up(0, M, Fs, Fe, det_rat * overlap);
        sign *= det_rat_sign;
        segments.push_back(segment_insert);
        //segments.insert(segments.end(), segment_insert);
        full_line = 0;
      }
      
    }
    
    else {
      compute_intervals(t, BETA, t_up, t_down, segments, s_up, s_down);
      
      if (t_down < 0) {         // t does lie on a segment -> it's possible to insert an anti-segment starting from t
        
        double length = compute_length(rng(), -t_down, 0);
        
        times segment_shrink(s_down->t_start(), t);
        
        double t_start = t + length;
        if (t_start > BETA)
          t_start -= BETA;
        
        times segment_insert(t_start, s_down->t_end());
        times anti_segment(t, t_start);
        
        double otherlength_u = 0;
        int FLAVORS = other_full_line.size();
        for (int i = 0; i < FLAVORS; i++) {
          if (i == this_flavor)
            continue;
          double other_length = compute_overlap(anti_segment, other_segments[i], other_full_line[i], BETA);
          otherlength_u += other_length * u(i, this_flavor);
        }
        
        double log_prob, overlap, det_rat, det_rat_sign;
        blas::vector R(segments.size());
        //std::vector<double> R(segments.size());
        det_rat = det_rat_insert_anti(anti_segment, M, segments, F, BETA, det_rat_sign, overlap, R);
        
        log_prob = log(BETA * (-t_down) / (segments.size() + 1) * det_rat) - length * mu + otherlength_u;
        
        if (log(rng()) < log_prob) {
          
          int s, r;             // s is the segment which is shifted, r the segment which is inserted
          s = 0;
          for (typename S::const_iterator it = segments.begin(); it != s_down; it++)
            s++;
          if (anti_segment.t_end() > segment_shrink.t_start())
            r = s + 1;
          else {
            r = 0;
            s++;
          }
          
          compute_M_insert_anti(anti_segment, s, r, M, segments, F, BETA, det_rat * overlap, R);
          s_down->set_t_end(t);
          /*                  times segment_new_endpoint(*s_down);
           segment_container_t::iterator prev_segment=s_down;
           prev_segment--;
           segment_new_endpoint.set_t_end(t);
           segments.erase(s_down); //erase old segment (without shifted endpoint) use s_down as an indicator of where we have to insert.
           s_down=segments.insert(prev_segment, segment_new_endpoint); //insert old segment with new endpoint
           segments.insert(s_down, segment_insert); //insert  new segment
           */
          if (segment_insert.t_start() > segments.begin()->t_start()) {
            s_down++;
            segments.insert(s_down, segment_insert);
          } else {
            segments.insert(segments.begin(), segment_insert);
          }
        }
      }
    }
    
  } else if (segments.size() > 1) {     // try to remove an anti-segment
    
    int r = (int) (rng() * segments.size());
    s_up = segments.begin();
    for (int i = 0; i < r; i++)
      s_up++;
    
    int s = r - 1;
    if (s < 0) {
      s = segments.size() - 1;
      s_down = segments.end();
      s_down--;
    } else {
      s_down = s_up;
      s_down--;
    }
    
    double length = s_up->t_start() - s_down->t_end();
    if (length < 0)
      length += BETA;
    
    double t_total = s_up->t_end() - s_down->t_end();
    if (t_total < 0)
      t_total += BETA;
    
    times anti_segment(s_down->t_end(), s_up->t_start());
    
    double otherlength_u = 0;
    int FLAVORS = other_full_line.size();
    for (int i = 0; i < FLAVORS; i++) {
      if (i == this_flavor)
        continue;
      double other_length = compute_overlap(anti_segment, other_segments[i], other_full_line[i], BETA);
      otherlength_u += other_length * u(i, this_flavor);
    }
    
    double log_prob, det_rat, det_rat_sign;
    
    det_rat = det_rat_remove_anti(anti_segment, r, s, M, segments, F, BETA, det_rat_sign);
    
    log_prob = log(BETA * t_total / segments.size() / det_rat) - length * mu + otherlength_u;
    
    if (log(rng()) < -log_prob) {
      
      compute_M_remove_anti(M, s, r);
      
      double t_end = s_up->t_end();
      segments.erase(s_up);
      
      if (r > 0) {
        s_up = segments.begin();
        for (int k = 0; k < s; k++)
          s_up++;
      } else {
        s = segments.size() - 1;
        s_up = segments.end();
        s_up--;
      }
      s_up->set_t_end(t_end);
      /*times s_up_new(*s_up);
       segment_container_t::iterator prev_segment=s_up;
       prev_segment--;
       s_up_new.set_t_end(t_end);
       segments.erase(s_up);
       segments.insert(prev_segment, s_up_new);
       */
    }
  }
  
  else if (segments.size() == 1) {
    
    s_down = segments.begin();
    
    double det_rat = fabs(M(0, 0));
    double length = s_down->t_start() - s_down->t_end();
    if (length < 0)
      length += BETA;
    times anti_segment(s_down->t_end(), s_down->t_start());
    
    double otherlength_u = 0;
    int FLAVORS = other_full_line.size();
    for (int i = 0; i < FLAVORS; i++) {
      if (i == this_flavor)
        continue;
      double other_length = compute_overlap(anti_segment, other_segments[i], other_full_line[i], BETA);
      otherlength_u += other_length * u(i, this_flavor);
    }
    
    double log_prob = log(BETA * BETA / det_rat) - length * mu + otherlength_u;
    
    if (log(rng()) < -log_prob) {
      full_line = 1;
      segments.erase(s_down);
      compute_M_down(0, M);     // attention: M.clear() sets elements to zero 
    }
  }
}

/*
 // move segment without changing its length
 template <class RNG, class S, class G> void move_segment(RNG& rng, S& segments, int N, double BETA, double mu, double u, G& F, blas_matrix& M, double & sign, S& other_segments, int other_full_line) {
 
 int size = segments.size();
 
 if (size < 1) return;
 
 int n = size*rng();
 
 typename S::iterator s=segments.begin();
 double new_t_start, new_t_end, length_segment;
 
 if (size==1) {
 new_t_start=rng()*BETA;
 length_segment = s->t_end()-s->t_start();
 if (length_segment<0)
 length_segment += BETA;      
 }
 else {
 
 double t_down, max_lower_gap; 
 typename S::iterator s_up, s_down;
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
template < class RNG, class S, class G > void shift_segment(RNG & rng, S & segments, double BETA, double mu, U_matrix & u, G & F,
                                                            blas_matrix & M, double &sign, std::vector < S > &other_segments,
                                                            std::vector < int >&other_full_line, int this_flavor)
{
  
  int size = segments.size();
  
  if (size < 1)
    return;
  
  int n = (int) (size * rng());
  
  typename S::iterator s, s_up;
  s = segments.begin();
  for (int i = 0; i < n; i++)
    s++;
  s_up = s;
  s_up++;
  if (s_up == segments.end())
    s_up = segments.begin();
  
  double interval = s_up->t_start() - s->t_start();
  if (interval <= 0)
    interval += BETA;
  
  double length = compute_length(rng(), interval, 0);
  double length_old = s->t_end() - s->t_start();
  if (length_old < 0)
    length_old += BETA;
  
  double new_t_end = s->t_start() + length;
  if (new_t_end > BETA)
    new_t_end -= BETA;
  
  times segment_insert(s->t_start(), new_t_end);
  times segment_remove = *s;
  
  double otherlength_u = 0;
  int FLAVORS = other_full_line.size();
  for (int i = 0; i < FLAVORS; i++) {
    if (i == this_flavor)
      continue;
    double other_length =
    compute_overlap(segment_insert, other_segments[i], other_full_line[i], BETA) - compute_overlap(segment_remove, other_segments[i],
                                                                                                   other_full_line[i], BETA);
    otherlength_u += other_length * u(i, this_flavor);
  }
  
  double det_rat, det_rat_sign, overlap;
  
  det_rat = det_rat_shift(segment_insert, n, M, segments, F, BETA, det_rat_sign, overlap);
  
  if (log(rng()) < log(det_rat) + (length - length_old) * mu - otherlength_u) {
    
    compute_M_shift(segment_insert, n, M, segments, F, BETA, det_rat * overlap);
    sign *= det_rat_sign;
    //s->set_t_end(new_t_end);
    times s_new(*s);
    s_new.set_t_end(new_t_end);
    s->set_t_end(new_t_end);  
    //typename S::iterator sit;
    //sit = segments.insert(s_new).first;
    /*if (sit == segments.end()) {
     std::cerr << "segment could not be inserted! exiting." << std::endl; exit(1);
     }*/
    std::cerr<<"reprogram, this is buggy!"<<std::endl;
    abort();
  }
}
// swap segment configurations
template < class RNG, class S, class G, class MAT > void swap_segments(RNG & rng, double BETA, G & F_up, G & F_down, S & segments_up,
                                                                       S & segments_down, int &full_line_up, int &full_line_down,
                                                                       double &sign_up, double &sign_down, MAT & M_up, MAT & M_down)
{
  
  MAT M_new_up, M_new_down;
  
  // before swap
  double det_old_up = construct_inverse(M_new_up, segments_up, BETA, F_up);     // here M_new_up is just a dummy
  double det_old_down = construct_inverse(M_new_down, segments_down, BETA, F_down);     // here M_new_down is just a dummy
  // before swap
  double det_new_up = construct_inverse(M_new_up, segments_down, BETA, F_up);
  double det_new_down = construct_inverse(M_new_down, segments_up, BETA, F_down);
  
  double det_rat = (det_new_up / det_old_up) * (det_new_down / det_old_down);
  
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

#endif
