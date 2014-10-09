#include<alps/hdf5.hpp>
#include <alps/hdf5/pointer.hpp>
#include <alps/hdf5/complex.hpp>

#include"./vector.h"
#include"./matrix.h"
namespace blas{
alps::hdf5::archive & serialize(alps::hdf5::archive & ar, std::string const & p, blas::vector const &v) {
      ar << alps::make_pvp(p+"/size", v.size());
      ar << alps::make_pvp(p+"/values", &(v(0)), v.size());
      return ar;
}
alps::hdf5::archive & serialize(alps::hdf5::archive & ar, std::string const & p, blas::vector &v) {
      int size;
      ar >> alps::make_pvp(p+"/size", size);
      v.resize(size);
      ar >> alps::make_pvp(p+"/values", &(v(0)), v.size());
      return ar;
}
/*template <> alps::hdf5::archive & serialize(alps::hdf5::archive & ar, std::string
    ** const & p, enum_type & v) {
    int tmp;
    ar >> alps::make_pvp(p, tmp);
    switch (tmp) {
    case 13: v = A; break;
    case 23: v = B; break;
    }
    return ar;
    }*/
alps::hdf5::archive & serialize(alps::hdf5::archive & ar, std::string const & p, blas::matrix const &M) {
      ar << alps::make_pvp(p+"/size", M.size());
      ar << alps::make_pvp(p+"/memory_size", M.memory_size());
      ar << alps::make_pvp(p+"/values", &(M(0,0)), M.memory_size()*M.memory_size());
      return ar;
}
alps::hdf5::archive & serialize(alps::hdf5::archive & ar, std::string const & p, blas::general_matrix<double> const &M) {
      ar << alps::make_pvp(p+"/size1", M.size1());
      ar << alps::make_pvp(p+"/size2", M.size2());
      ar << alps::make_pvp(p+"/values", &(M(0,0)), M.size1()*M.size2());
      return ar;
}
alps::hdf5::archive & serialize(alps::hdf5::archive & ar, std::string const & p, blas::matrix &M) {
      double size;
      double memory_size;
      ar >> alps::make_pvp(p+"/size", size);
      ar >> alps::make_pvp(p+"/memory_size", memory_size);
      std::cout<<"resizing M: "<<size<<" memory size is: "<<memory_size<<std::endl;
      M.resize(memory_size);
      M.memory_size()=memory_size;
      M.size()=size;
      if(size!=0)
        ar >> alps::make_pvp(p+"/values", &(M(0,0)), memory_size*memory_size);
      return ar;
}
alps::hdf5::archive & serialize(alps::hdf5::archive & ar, std::string const & p, blas::general_matrix<double> &M) {
      double size1, size2;
      ar >> alps::make_pvp(p+"/size1", size1);
      ar >> alps::make_pvp(p+"/size2", size2);
      M.resize(size1, size2);
      ar >> alps::make_pvp(p+"/values", &(M(0,0)), size1*size2);
      return ar;
}
alps::hdf5::archive & serialize(alps::hdf5::archive & ar, std::string const & p, blas::general_matrix<std::complex<double> > const &M) {
      ar << alps::make_pvp(p+"/size1", M.size1());
      ar << alps::make_pvp(p+"/size2", M.size2());
      ar << alps::make_pvp(p+"/values", &(M(0,0)), M.size1()*M.size2());
      return ar;
}
}
