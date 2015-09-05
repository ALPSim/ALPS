
template<typename T, class D>
green_function<T> AntiFerroDCATransformer::transform_coordinates(const green_function<T>& G) const
{
  unsigned int nflavor=G.nflavor();
  unsigned int nsite=G.nsite();
  unsigned int nt=G.ntime();
  if(nsite!=n_cluster_sites()) throw std::logic_error("Problem with cluster sites != n_sites.");
  if(G.shape()==diagonal) throw std::logic_error("problem with GF shape!");
  shape_t transformed_shape=G.shape()==blockdiagonal?nondiagonal:blockdiagonal;
  green_function<T> G_new(nt, nsite, nflavor, transformed_shape);
  G_new.clear();
  for (site_t k=0; k<nsite; k+=2) {
    for (site_t q1=0; q1<=1; ++q1) {
      for (site_t q2=0; q2<=1; ++q2) {
        for (site_t i=0; i<nsite; ++i) {
          for (site_t j=0; j<nsite; ++j) {
            double phase = 0.;
            for (unsigned int d=0; d<dimension(); ++d){
              phase += cluster_momentum(k+q1,d)*cluster_coordinate(i,d)-cluster_momentum(k+q2,d)*cluster_coordinate(j,d);
            }for (itime_index_t t=0; t<nt; ++t) {
              for (spin_t flavor=0; flavor<nflavor; ++flavor) {
                D::add(G_new, G, t, i, j, k+q1, k+q2, flavor, cos(phase)/double(nsite));
              }
            }
          }
        }
      }
    }
  }
  return G_new;
}


struct AntiFerroDCATransformer::Forward {
  template<typename T>
  static void add(green_function<T>& G_new, const green_function<T>& G,
                  const itime_index_t t, const site_t i, const site_t j, const site_t k1, const site_t k2, const spin_t flavor,
                  double factor)
  {
    G_new(t, k1, k2, flavor) += factor * G(t, i, j, flavor);
  }
};



struct AntiFerroDCATransformer::Backward {
  template<typename T>
  static void add(green_function<T>& G_new, const green_function<T>& G,
                  const itime_index_t t, const site_t i, const site_t j, const site_t k1, const site_t k2, const spin_t flavor,
                  double factor)
  {
    G_new(t, i, j, flavor) += factor * G(t, k1, k2, flavor);
  }
};

