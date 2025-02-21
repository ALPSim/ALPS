/* 
 * Copyright Brigitte Surer and Matthias Troyer 2006-2008
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
 */


#ifndef ALPS_RANDOM_PARALLEL_MERSENNE_TWISTER_HPP
#define ALPS_RANDOM_PARALLEL_MERSENNE_TWISTER_HPP

#include <alps/random/mersenne_twister.hpp>
#include <alps/random/parallel/lcg64.hpp>
#include <alps/random/parallel/seed.hpp>

namespace alps { namespace random { namespace parallel {

template<class UIntType, int w, int n, int m, int r, UIntType a, int u,
  int s, UIntType b, int t, UIntType c, int l, UIntType val>
void seed(
    boost::random::mersenne_twister<UIntType,w,n,m,r,a,u,s,b,t,c,l,val>& prng
  , unsigned int num
  , unsigned int total
)
  {
    lcg64a seeder;
    seed(prng, num, total, seeder);
  }

template<class UIntType, int w, int n, int m, int r, UIntType a, int u,
  int s, UIntType b, int t, UIntType c, int l, UIntType val, class SeedType>
void seed(
    boost::random::mersenne_twister<UIntType,w,n,m,r,a,u,s,b,t,c,l,val>& prng
  , unsigned int num
  , unsigned int total
  , SeedType const& the_seed
      )
  {
    lcg64a seeder(the_seed);
    seed(prng, num, total, seeder);
  }

template<class UIntType, int w, int n, int m, int r, UIntType a, int u,
  int s, UIntType b, int t, UIntType c, int l, UIntType val>
void seed(
    boost::random::mersenne_twister<UIntType,w,n,m,r,a,u,s,b,t,c,l,val>& prng
      , unsigned int num
      , unsigned int total
      , lcg64a & engine
      )
  {
     //seeds the seeder, which in turn gives the seedvalue for the mersenne_twister-rng
    typedef boost::uniform_int<unsigned int> dist_t;    
    boost::variate_generator<lcg64a&, dist_t> rng(engine, dist_t(0u, std::numeric_limits<unsigned int>::max()));

    //warm-up to improve decorrelations
    for(unsigned int i = 0; i < 1000; i++)
      rng();
      
    std::vector<UIntType> buffer;
    for (int i = 0; i < 2*n; i++) 
         buffer.push_back(rng());
    
    // seed the generator
    typename std::vector<UIntType>::iterator first=buffer.begin();
    prng.seed(first,buffer.end());
  }

template<class UIntType, int w, int n, int m, int r, UIntType a, int u,
  int s, UIntType b, int t, UIntType c, int l, UIntType val, class Iterator>
void seed(
    boost::random::mersenne_twister<UIntType,w,n,m,r,a,u,s,b,t,c,l,val>& prng
  , unsigned int num
  , unsigned int total
  , Iterator& first
  , Iterator const& last
      )
  {
   first += 2*n*num;
   if(first > last+2*n)
      boost::throw_exception(std::invalid_argument("parallel_seed"));
    prng.seed(first,last);
  }

    

} } } // end namespace random::parallel::alps
  

#endif /*ALPS_RANDOM_PARALLEL_MERSENNE_TWISTERL_HPP*/

