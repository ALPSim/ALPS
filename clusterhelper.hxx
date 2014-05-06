

template<typename T>
inline typename ClusterHelper::Matrix<T>::Type ClusterHelper::Matrix<T>::inverse(typename Matrix<T>::Type M)
{
  typename Matrix<T>::Type M_inv(M.size1(), M.size2());
  for (unsigned int i=0; i<M.size1(); ++i)
    for (unsigned int j=0; j<M.size2(); ++j)
      M_inv(i,j) = (T)(i==j);
  boost::numeric::ublas::vector<fortran_int_t> ipivot(M.size1());
  memset(&(ipivot(0)), 0, sizeof(int)*M.size1());
  boost::numeric::bindings::lapack::gesv(M,ipivot, M_inv);
  return M_inv;
}




template<typename T>
inline T ClusterHelper::interpolate(vector_type k, const std::vector<T>& values, const bool antiferromagnet) const
{
  k = round(matrix_vector_mult(reciprocal_basis_inv_, k));
  std::vector<vector_type> neighbor;
  std::vector<unsigned int> neighbor_index;
  switch (dimension()) {
    case 1:
      neighbor.resize(2);
      neighbor_index.resize(2);
      for (int i=0; i<2; ++i)
        neighbor[i].resize(1);
      neighbor[0][0] = floor(k[0]);
      neighbor[1][0] = ceil(k[0]);
      break;
    case 2:
      neighbor.resize(4);
      neighbor_index.resize(4);
      for (int i=0; i<4; ++i)
        neighbor[i].resize(2);
      for (int x=0; x<=1; ++x) {
        double xt = (x==0) ? floor(k[0]) : ceil(k[0]);
        for (int y=0; y<=1; ++y) {
          double yt = (y==0) ? floor(k[1]) : ceil(k[1]);
          neighbor[y*2+x][0] = xt;
          neighbor[y*2+x][1] = yt;
        }
      }
      break;
    case 3:
      neighbor.resize(8);
      neighbor_index.resize(8);
      for (int i=0; i<8; ++i)
        neighbor[i].resize(3);
      for (int x=0; x<=1; ++x) {
        double xt = (x==0) ? floor(k[0]) : ceil(k[0]);
        for (int y=0; y<=1; ++y) {
          double yt = (y==0) ? floor(k[1]) : ceil(k[1]);
          for (int z=0; z<=1; ++z) {
            double zt = (z==0) ? floor(k[2]) : ceil(k[2]);
            neighbor[z*2*2+y*2+x][0] = xt;
            neighbor[z*2*2+y*2+x][1] = yt;
            neighbor[z*2*2+y*2+x][2] = zt;
          }
        }
      }
      break;
  }
  for (std::size_t i=0; i<neighbor.size(); ++i) {
    vector_type nv = matrix_vector_mult(reciprocal_basis_, neighbor[i]);
    nv = map_back_into_brillouin_zone(nv);
    if (antiferromagnet) {
      if (!inside_reduced_brillouin_zone(nv)) {
        for (std::size_t d=0; d<dimension(); ++d)
          nv[d] -= M_PI;
        nv = map_back_into_brillouin_zone(nv);
      }
      if (!inside_reduced_brillouin_zone(nv)) {
        boost::throw_exception(std::logic_error("something doesn't work in ClusterHelper::interpolate)"));
      }
    }
    neighbor_index[i] = identify_vertex(nv, cluster_momenta_);
  }
  std::vector<T> av(neighbor_index.size());
  for (std::size_t i=0; i<neighbor_index.size(); ++i)
    av[i] = values[neighbor_index[i]];
  for (std::size_t d=0; d<dimension(); ++d) {
    double c = k[d] - floor(k[d]) ;
    std::vector<T> t(av.size()/2);
    for (std::size_t i=0; i<av.size(); i+=2) {
      t[i/2] = (1-c)*av[i] + c*av[i+1];
    }
    av = t;
  }
  return av[0];
}


