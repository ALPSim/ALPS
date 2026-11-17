#include<alps/osiris/dump.h>
#include <alps/alea/obsvalue.h>
#include<vector.h>
namespace alps{
inline alps::IDump& operator>>(alps::IDump &dump, blas::vector &v1){
  return dump;
}
inline alps::ODump& operator<<(alps::ODump &dump, const blas::vector &v1){
  return dump;
}


template<>
struct obs_value_traits<blas::vector >
{
  typedef blas::vector value_type;
  typedef double element_type;
  typedef double count_type;
  typedef unsigned int index_type;
  typedef std::size_t size_type;
  typedef blas::vector time_type;
  typedef double time_element_type;
  BOOST_STATIC_CONSTANT(bool, array_valued = true);
  
  typedef boost::numeric::ublas::matrix<type_traits<double>::average_t> covariance_type;
  typedef blas::vector result_type;
  typedef std::valarray<int> convergence_type;
  BOOST_STATIC_CONSTANT(uint32_t, magic_id = 256+type_traits<double>::type_tag);
  typedef std::vector<std::string> label_type;

  static inline void check_for_max(blas::vector& a,const blas::vector& b) 
  {
/*    for(int32_t i=0;i!=(int32_t)a.size();++i)
      obs_value_traits<double>::check_for_max(a[i],b[i]);*/
    for(int32_t i=0;i!=(int32_t)a.size();++i)
      a[i]=a[i]>b[i]?a[i]:b[i];
  }   

  static inline void check_for_max(blas::vector& a,const double b) 
  {
    for(int32_t i=0;i!=(int32_t)a.size();++i)
      obs_value_traits<double>::check_for_max(a[i],b);
  }   
 
  static inline void check_for_min(blas::vector& a,const blas::vector& b) 
  {
    /*for(int32_t i=0;i!=(int32_t)a.size();++i)
      obs_value_traits<double>::check_for_min(a[i],b[i]);*/
    for(int32_t i=0;i!=(int32_t)a.size();++i)
      a[i]=a[i]<b[i]?a[i]:b[i];
  }    

  static inline void check_for_min(blas::vector& a,const double b) 
  {
    for(int32_t i=0;i!=(int32_t)a.size();++i)
      obs_value_traits<double> ::check_for_min(a[i],b);
  }    

  static void fix_negative(value_type& a) 
  {
    for(int32_t i=0;i!=(int32_t)a.size();++i)
      obs_value_traits<element_type>::fix_negative(a[i]);
  }
  static element_type min() {return obs_value_traits<double>::min();}
  static element_type max() {return obs_value_traits<double>::max();}
  static element_type epsilon() { return obs_value_traits<double>::epsilon();}

  static time_element_type t_max() {return obs_value_traits<double>::max();}

  static inline time_type check_divide(const result_type& a,const result_type& b) 
  {
    time_type retval;
    resize_same_as(retval,b);
    for(int32_t i(0);i<(int32_t)b.size();++i)
      retval[i] = obs_value_traits<element_type>::check_divide(a[i],b[i]);
    return retval;
  }

  /** resize a to given size */
  template <class X, class Y> static void resize_same_as(X& a, const Y& y) {a.resize(y.size());}
  template <class X, class Y> static void copy(X& x,const Y& y) {x.resize(y.size()); for (int i=0;i<(int)y.size();++i) x[i]=y[i]; }
  template <class X> static std::size_t size(const X& a) { return a.size();}
  template <class X> static void resize(X& a, std::size_t s) {a.resize(s);}

  typedef uint32_t slice_iterator;
  static slice_iterator slice_begin(const value_type&) { return 0;}
  static slice_iterator slice_end(const value_type& x) { return x.size();}
  static std::string slice_name(const value_type& ,slice_iterator i) 
    { return boost::lexical_cast<std::string,int>(i); }
  static element_type slice_value(const value_type& x, slice_iterator i) { return x[i];}
  static element_type& slice_value(value_type& x, slice_iterator i) { return  x[i];}
 
  template <class X> static blas::vector convert(const std::valarray<X>& x)  //this is pointless.
  { 
    std::cout<<"converting valarray!"<<std::endl;
    blas::vector res(x.size());
    copy(res,x);
    return res;
  }
  static const blas::vector &convert(const blas::vector& x)  //this is pointless.
  { 
    std::cout<<"this should not be called"<<std::endl;
    return x;
    //return res;
  }

};
template <> const blas::vector &obs_value_cast(const blas::vector& s)
{
    return s;
}

/*template <blas::vector, blas::vector> blas::vector &obs_value_cast(blas::vector& s)
{
  std::cout<<"dude, somebody called obs value cast!"<<std::endl;
    return s;
}
*/
}
