struct SimpleDCATransformer::Forward {
  template<typename T>
  static void add(green_function<T>& G_new, const green_function<T>& G,
                  const itime_index_t t, const site_t i, const site_t j, const site_t k, const spin_t flavor,
                  double factor) {
    if((G_new.shape()!=diagonal) || (G.shape() != nondiagonal)) throw std::logic_error("problem with shape of GF in add!");
    G_new(t, k, k, flavor) += factor * G(t, i, j, flavor);
  }
  template<typename T>
  static void add(nambu_green_function<T>& G_new, const nambu_green_function<T>& G,
                  const itime_index_t t, const site_t i, const site_t j, const site_t k,
                  double factor) {
    G_new(t, k, k, 0,0) += factor * G(t, i, j, 0,0);
    G_new(t, k, k, 0,1) += factor * G(t, i, j, 0,1);
    G_new(t, k, k, 1,0) += factor * G(t, i, j, 1,0);
    G_new(t, k, k, 1,1) += factor * G(t, i, j, 1,1);
  }
};



struct SimpleDCATransformer::Backward {
  template<typename T>
  static void add(green_function<T>& G_new, const green_function<T>& G,
                  const itime_index_t t, const site_t i, const site_t j, const site_t k, const spin_t flavor,
                  double factor) {
    G_new(t, i, j, flavor) += factor * G(t, k, k, flavor);
  }
  template<typename T>
  static void add(nambu_green_function<T>& G_new, const nambu_green_function<T>& G,
                  const itime_index_t t, const site_t i, const site_t j, const site_t k,
                  double factor) {
    G_new(t, i, j, 0, 0) += factor * G(t, k, k, 0, 0);
    G_new(t, i, j, 0, 1) += factor * G(t, k, k, 0, 1);
    G_new(t, i, j, 1, 0) += factor * G(t, k, k, 1, 0);
    G_new(t, i, j, 1, 1) += factor * G(t, k, k, 1, 1);
  }
};



template<typename T, class D>
green_function<T> SimpleDCATransformer::transform_coordinates(const green_function<T>& G) const
{
  unsigned int nflavor=G.nflavor();
  unsigned int nsite=G.nsite();
  unsigned int nt=G.ntime();
  assert(nsite==n_cluster_sites());
  shape_t new_shape=G.shape()==diagonal?nondiagonal:diagonal;
  green_function<T> G_new(nt, nsite, nflavor, new_shape);
  G_new.clear();
  for (site_t k=0; k<nsite; ++k) {
    for (site_t i=0; i<nsite; ++i) {
      for (site_t j=0; j<nsite; ++j) {
        double phase = 0.;
        for (unsigned int d=0; d<dimension(); ++d)
          phase += cluster_momentum(k,d)*(cluster_coordinate(i,d)-cluster_coordinate(j,d));
        for (itime_index_t t=0; t<nt; ++t) {
          for (spin_t flavor=0; flavor<nflavor; ++flavor) {
            D::add(G_new, G, t, i, j, k, flavor, cos(phase)/double(nsite));
          }
        }
      }
    }
  }
  return G_new;
}


template<class T>
green_function<T> SimpleDCATransformer::spin_symmetrizer(const green_function<T>& G, const bool paramagnet) const
{
  green_function<T> G_new(G);
  if (paramagnet) {
    if(G_new.nflavor()>1 && G_new.nflavor()%2==0) {
      for(spin_t flavor=0; flavor<G_new.nflavor(); flavor+=2) {
        for (site_t k1=0; k1<n_cluster_sites(); ++k1) {
          for (site_t k2=0; k2<n_cluster_sites(); ++k2) {
            if((k1 != k2) && (G.shape()==diagonal)) continue;
            for(frequency_t freq=0; freq<G_new.nfreq(); ++freq) {
              G_new(freq, k1, k2, flavor)= 0.5*(G_new(freq, k1, k2, flavor)+G_new(freq, k1, k2, flavor+1));
              G_new(freq, k1, k2, flavor+1)=G_new(freq, k1, k2, flavor);
            }
          }
        }
      }
    }
  }
  return G_new;
}



inline matsubara_green_function_t SimpleDCATransformer::fourpoint_green_function_inverse(const matsubara_green_function_t& G,const std::size_t nf) const
{
  matsubara_green_function_t G_inv(G);
  for (spin_t f=0; f<G.nflavor(); ++f) {
    Matrix<std::complex<double> >::Type M(nf*G.nsite(), nf*G.nsite());
    for (std::size_t k1=0; k1<G.nsite(); k1++) {
      for (std::size_t k2=0; k2<G.nsite(); k2++) {
        for (std::size_t t1=0; t1<nf; t1++) {
          for (std::size_t t2=0; t2<nf; t2++) {
            M(k1*nf+t1, k2*nf+t2) = G(nf*t1+t2, k1, k2, f);
          }
        }
      }
    }

    Matrix<std::complex<double> >::Type M_inv = Matrix<std::complex<double> >::inverse(M);
    
    for (std::size_t k1=0; k1<G.nsite(); k1++) {
      for (std::size_t k2=0; k2<G.nsite(); k2++) {
        for (std::size_t t1=0; t1<nf; t1++) {
          for (std::size_t t2=0; t2<nf; t2++) {
            G_inv(nf*t1+t2, k1, k2, f) = M_inv(k1*nf+t1, k2*nf+t2);
          }
        }
      }
    }
  }
  return G_inv;
}



template <typename T>
green_function<T> SimpleDCATransformer::fourpoint_green_function_add(const green_function<T>& G1,
                                                                     const green_function<T>& G2) const
{
  green_function<T> G(G1);
  for (time_t t=0; t<G.ntime(); ++t) {
    for (spin_t f=0; f<G.nflavor(); ++f) {
      for (int k1=0; k1<G.nsite(); k1++) {
        for (int k2=0; k2<G.nsite(); k2++) {
          G(t, k1, k2, f) += G2(t, k1, k2, f);
        }
      }
    }
  }
  return G;
}



template <typename T>
green_function<T> SimpleDCATransformer::fourpoint_green_function_subtract(const green_function<T>& G1,
                                                                          const green_function<T>& G2) const
{
  green_function<T> G(G1);
  for (time_t t=0; t<G.ntime(); ++t) {
    for (spin_t f=0; f<G.nflavor(); ++f) {
      for (std::size_t k1=0; k1<G.nsite(); k1++) {
        for (std::size_t k2=0; k2<G.nsite(); k2++) {
          G(t, k1, k2, f) -= G2(t, k1, k2, f);
        }
      }
    }
  }
  return G;
}



template <typename T>
green_function<T> SimpleDCATransformer::fourpoint_green_function_divide(const green_function<T>& G1,
                                                                        const T d) const
{
  green_function<T> G(G1);
  for (time_t t=0; t<G.ntime(); ++t) {
    for (spin_t f=0; f<G.nflavor(); ++f) {
      for (int k1=0; k1<G.nsite(); k1++) {
        for (int k2=0; k2<G.nsite(); k2++) {
          G(t, k1, k2, f) /= d;
        }
      }
    }
  }
  return G;
}

template<typename T, class D>
nambu_green_function<T> SimpleDCATransformer::transform_coordinates(const nambu_green_function<T>& G) const
{
  unsigned int nsite=G.nsite();
  unsigned int nt=G.ntime();
  assert(nsite==n_cluster_sites());
  nambu_green_function<T> G_new(nt, nsite, 2, G.shape()==diagonal?nondiagonal:diagonal);
  G_new.clear();
  //compute G = \sum_{ij} \sum_K exp(i K (R_i-R_j)) G(R_i, R_j)
  for (site_t k=0; k<nsite; ++k) {
    for (site_t i=0; i<nsite; ++i) {
      for (site_t j=0; j<nsite; ++j) {
        double phase = 0.;
        for (unsigned int d=0; d<dimension(); ++d)
          phase += cluster_momentum(k,d)*(cluster_coordinate(i,d)-cluster_coordinate(j,d));
        for (itime_index_t t=0; t<nt; ++t) {
          D::add(G_new, G, t, i, j, k, cos(phase)/double(nsite));
        }
      }
    }
  }
  return G_new;
}

