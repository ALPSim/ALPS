// Copyright 2006 Matthias Troyer

 // Permission is hereby granted, free of charge, to any person obtaining
 // a copy of this software and associated documentation files (the “Software”),
 // to deal in the Software without restriction, including without limitation
 // the rights to use, copy, modify, merge, publish, distribute, sublicense,
 // and/or sell copies of the Software, and to permit persons to whom the
 // Software is furnished to do so, subject to the following conditions:
 //
 // The above copyright notice and this permission notice shall be included
 // in all copies or substantial portions of the Software.
 //
 // THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
 // OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 // FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 // AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 // LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 // FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 // DEALINGS IN THE SOFTWARE.

#include <alps/random/parallel/detail/get_prime.hpp>
#include <alps/random/parallel/detail/primelist_64.hpp>
#include <boost/throw_exception.hpp>
#include <boost/assert.hpp>
#include <stdexcept>
#include <string>

// taken from SPRNG implementation

#define MAXPRIME 3037000501U  /* largest odd # < sqrt(2)*2^31+2 */
#define MINPRIME 55108   /* sqrt(MAXPRIME) */
#define MAXPRIMEOFFSET 146138719U /* Total number of available primes */
#define NPRIMES 10000
#define PRIMELISTSIZE1 1000
#define STEP 10000

namespace alps { namespace random { namespace detail {

struct primes_storage
{
  primes_storage()
   : obtained(0)
  {
    int i, j;
    bool isprime;
  
    for(i=3; i < MINPRIME; i += 2) {
      isprime = true;
    
      for(j=0; j < obtained; j++)
        if(i%primes_[j] == 0) {
          isprime = false;
          break;
        }
        else if(primes_[j]*primes_[j] > i)
          break;

      if(isprime) {
        primes_[obtained] = i;
        obtained++;
      }
    }
  }
  
  int operator[](int i) const
  {
    return primes_[i];
  }
  
  int size() const
  {
    return obtained;
  }
  
private:
  int primes_[NPRIMES];
  int obtained;
};

static primes_storage primes;


boost::uint64_t get_prime_64(unsigned int offset)
{
  BOOST_ASSERT(offset <= MAXPRIMEOFFSET);

  if(offset<PRIMELISTSIZE1) 
    return primelist_64[offset];

  unsigned int largest = MAXPRIME;
  
  int index = (unsigned int) ((offset-PRIMELISTSIZE1+1)/STEP) + PRIMELISTSIZE1 -  1;
  largest = primelist_64[index] + 2;
  offset -= (index-PRIMELISTSIZE1+1)*STEP + PRIMELISTSIZE1 - 1;
  
  while(largest > MINPRIME)
  {
    bool isprime = true;
    largest -= 2;
    for(int i=0; i<primes.size(); i++)
      if(largest%primes[i] == 0) {
        isprime = false;
        break;
      }
    
    if(isprime && offset > 0)
      offset--;
    else if(isprime)
      return largest;
  }
  
  // Casting to std::string is a workaround for Fujitsu FCC Compiler
  boost::throw_exception(std::runtime_error(std::string("Insufficient number of primes")));
  return 0; // dummy return

}


} } } // namespace alps::random::detail
