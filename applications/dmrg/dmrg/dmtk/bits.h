#ifndef __DMTK_BITS_H__
#define __DMTK_BITS_H__

namespace dmtk
{

#define IS_EVEN(n) ((n) & 1 == 1 ? false : true)
#define SGN(n) (1 - (((n) & 1) << 1))
#define IEOR(n1,n2) ((n1) ^ (n2))
#define IAND(n1,n2) ((n1) & (1 << (n2)))
#define IOR(n1,n2) ((n1) | (n2))
#define IBITS(n,i) (((n) & (1 << i)) >> i)
#define IBSET(n,i) ((n) | (1 << i))
#define IBCLR(n,i) ((n) ^ (IBITS((n),i) << i))
inline int ISHFTC(int n, int i, int nt)
{
  int mask = (1 << i) - 1;
  int tail = n & mask;
  int r = n >> i;
  r |= tail << (nt - i);
  return r;
}

inline size_t mask(int b1, int b2) { return (1 << (b1-1) | 1 << (b2-1)); }
inline size_t mask(int b) { return (1 << (b-1)); }
inline int mod2(int n) { return (n & 1); }

} // namespace dmtk

#endif // __DMTK_BITS_H__
