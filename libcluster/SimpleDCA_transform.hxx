


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
  for (int k=0; k<nsite; ++k) {
    for (int i=0; i<nsite; ++i) {
      for (int j=0; j<nsite; ++j) {
        double phase = 0.;
        for (unsigned int d=0; d<dimension(); ++d)
          phase += cluster_momentum(k,d)*(cluster_coordinate(i,d)-cluster_coordinate(j,d));
        for (int t=0; t<nt; ++t) {
          for (int flavor=0; flavor<nflavor; ++flavor) {
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
      for(int flavor=0; flavor<G_new.nflavor(); flavor+=2) {
        for (int k1=0; k1<n_cluster_sites(); ++k1) {
          for (int k2=0; k2<n_cluster_sites(); ++k2) {
            if((k1 != k2) && (G.shape()==diagonal)) continue;
            for(int freq=0; freq<G_new.nfreq(); ++freq) {
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
  for (int f=0; f<G.nflavor(); ++f) {
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



template<typename T, class D>
nambu_green_function<T> SimpleDCATransformer::transform_coordinates(const nambu_green_function<T>& G) const
{
  unsigned int nsite=G.nsite();
  unsigned int nt=G.ntime();
  assert(nsite==n_cluster_sites());
  nambu_green_function<T> G_new(nt, nsite, 2, G.shape()==diagonal?nondiagonal:diagonal);
  G_new.clear();
  //compute G = \sum_{ij} \sum_K exp(i K (R_i-R_j)) G(R_i, R_j)
  for (int k=0; k<nsite; ++k) {
    for (int i=0; i<nsite; ++i) {
      for (int j=0; j<nsite; ++j) {
        double phase = 0.;
        for (unsigned int d=0; d<dimension(); ++d)
          phase += cluster_momentum(k,d)*(cluster_coordinate(i,d)-cluster_coordinate(j,d));
        for (int t=0; t<nt; ++t) {
          D::add(G_new, G, t, i, j, k, cos(phase)/double(nsite));
        }
      }
    }
  }
  return G_new;
}

