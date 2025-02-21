/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1994-2005 by Matthias Troyer <troyer@comp-phys.org>
*
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
*****************************************************************************/

/* $Id$ */

#include "../diag.h"
#include <alps/numeric/real.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <ietl/interface/ublas.h>
#include <ietl/vectorspace.h>
#include <ietl/lanczos.h>
#include <boost/random.hpp>
#include <boost/regex.hpp> 
#include <complex>

template <class T>
class SparseDiagMatrix : public DiagMatrix<T, boost::numeric::ublas::mapped_vector_of_mapped_vector<T, boost::numeric::ublas::row_major> >
{
public:
  typedef DiagMatrix<T, boost::numeric::ublas::mapped_vector_of_mapped_vector<T, boost::numeric::ublas::row_major> > super_type;
  typedef T value_type;
  typedef typename super_type::magnitude_type magnitude_type;
  typedef typename super_type::matrix_type matrix_type;
  typedef typename super_type::site_iterator site_iterator;
  typedef typename super_type::site_descriptor site_descriptor;
  typedef typename boost::numeric::ublas::vector<value_type> vector_type;
  typedef typename super_type::size_type size_type;
  typedef typename super_type::mag_vector_type mag_vector_type;
  typedef typename super_type::half_integer_type half_integer_type;
  typedef typename super_type::operator_matrix_type operator_matrix_type;
  
  SparseDiagMatrix (const alps::ProcessList& where , const boost::filesystem::path& p);
  void do_subspace();
  void write_xml_body(alps::oxstream&, const boost::filesystem::path&, bool) const;
  void print_eigenvectors(std::ostream& os) const;
private:
  
  std::vector<value_type> calculate(operator_matrix_type const& m) const;
  std::vector<vector_type> eigenvectors;
};

template <class T>
SparseDiagMatrix<T>::SparseDiagMatrix(const alps::ProcessList& where , const boost::filesystem::path& p) 
 : super_type(where,p,true) 
{ 
  this->construct();
}


template <class T>
void SparseDiagMatrix<T>::write_xml_body(alps::oxstream& out, const boost::filesystem::path& p, bool writeallxml) const
{
  if (writeallxml) {   
      // Get minimum energy (over all sectors)
      magnitude_type min_val = std::numeric_limits<magnitude_type>::max();
      magnitude_type exc_val = std::numeric_limits<magnitude_type>::max();
      for (typename std::vector<mag_vector_type>::const_iterator vit=this->eigenvalues_.begin(); vit!=this->eigenvalues_.end(); ++vit)
        if(vit->size()>=1) {
          if ((*vit)[0] < min_val) {
            exc_val = min_val;
            min_val = (*vit)[0];
            if (vit->size()>=2 && (*vit)[1] < exc_val)
              exc_val = (*vit)[1];
          }
          else if ((*vit)[0] < exc_val)
            exc_val = (*vit)[0];
        }

      out << alps::start_tag("AVERAGES");
      if (min_val!=std::numeric_limits<magnitude_type>::max())
        out << alps::start_tag("SCALAR_AVERAGE") << alps::attribute("name","Ground State Energy") << alps::no_linebreak
            << alps::start_tag("MEAN") << min_val << alps::end_tag("MEAN")
            << alps::end_tag("SCALAR_AVERAGE");
      if (exc_val!=std::numeric_limits<magnitude_type>::max())
        out << alps::start_tag("SCALAR_AVERAGE") << alps::attribute("name","Energy Gap") << alps::no_linebreak
            << alps::start_tag("MEAN") << alps::no_linebreak << exc_val-min_val << alps::end_tag("MEAN")
            << alps::end_tag("SCALAR_AVERAGE");
      out << alps::end_tag("AVERAGES");
  }
  super_type::write_xml_body(out,p,writeallxml);
}
   

template <class T>
void SparseDiagMatrix<T>::do_subspace()
{
  typedef ietl::vectorspace<vector_type> vectorspace_type;
  using alps::numeric::real;
  boost::lagged_fibonacci607 generator;
  mag_vector_type ev;
  this->build();
  if (this->dimension()==0)
    return;

  vectorspace_type vec(this->dimension());
  ietl::lanczos<matrix_type,vectorspace_type> lanczos(this->matrix(),vec);
  int n;
  
  if (this->dimension()>1) {
    int max_iter = this->get_parameters().value_or_default("MAX_ITERATIONS",std::min(int(10*this->dimension()),1000));  
    int num_eigenvalues = this->get_parameters().value_or_default("NUMBER_EIGENVALUES",1);
    ietl::lanczos_iteration_nlowest<double> iter(max_iter,num_eigenvalues);
    std::cerr << "Starting Lanczos \n";
    lanczos.calculate_eigenvalues(iter,generator);
    std::cerr << "Finished Lanczos\n";
    n=std::min(num_eigenvalues,int(lanczos.eigenvalues().size()));
    ev.resize(n);
    for (int i=0;i<n;++i) 
      ev[i]=lanczos.eigenvalues()[i];
  }
  else
  {
    ev.resize(1);
    ev[0]=real(value_type(this->matrix()(0,0)));
  }
  if (this->calc_vectors()) {
    if  (this->dimension()>1) {
      // calculate eigen vectors
      ietl::Info<magnitude_type> info; // (m1, m2, ma, eigenvalue, residualm, status).
  
      try {
        eigenvectors.clear();
        lanczos.eigenvectors(lanczos.eigenvalues().begin(),lanczos.eigenvalues().begin()+n,
                             std::back_inserter(eigenvectors),info,generator); 
      }
      catch (std::runtime_error& e) {
        std::cout <<"Exception during eigenvector calculation: " <<  e.what() << "\n";
      }  
    }
    else {
      vector_type v(1);
      v[0]=1.;
      eigenvectors.push_back(v);
    }
  }
  this->perform_measurements();
  eigenvectors.clear();
  std::copy(ev.begin(),ev.end(),std::back_inserter(this->measurements_.rbegin()->average_values["Energy"]));
  this->eigenvalues_.push_back(ev);
}

template <class T>
std::vector<T> SparseDiagMatrix<T>::calculate(operator_matrix_type const& m) const
{
  using namespace boost::numeric::ublas;
  std::vector<value_type> av;
  for(typename std::vector<vector_type>::const_iterator it = eigenvectors.begin();it!=eigenvectors.end();it++)
    av.push_back(inner_prod(boost::numeric::ublas::conj(*it),prod(m,*it)));
  return av;
}

template <class T>
void SparseDiagMatrix<T>::print_eigenvectors(std::ostream& os) const
{
  unsigned int n=0;
  for(typename std::vector<vector_type>::const_iterator it = eigenvectors.begin();it!=eigenvectors.end();it++)
    os << "Eigenvector# " << n++ << ":\n" << *it << "\n";
}
