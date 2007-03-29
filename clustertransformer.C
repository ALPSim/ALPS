/*****************************************************************************
*
* ALPS DMFT Project
*
* Copyright (C) 2007 by Sebastian Fuchs <sfuchs@gwdg.de>
*
* This software is part of the ALPS Applications, published under the ALPS
* Application License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Application License along with
* the ALPS Applications; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/


#include "clustertransformer.h"


ClusterHelper::ClusterHelper(const alps::Parameters& params) :
  alps::graph_helper<alps::graph_helper<>::graph_type>(params)
{
  if (boost::lexical_cast<int>(params["SITES"]) != num_sites() / volume())
    boost::throw_exception(std::logic_error("parameter SITES does not match cluster definition"));
  std::ofstream momstream("LatticeMomenta.dat");
  momentum_iterator itm, itm_end;
  for (boost::tie(itm, itm_end) = momenta(); itm!=itm_end; ++itm) {
    vector_type m=momentum(*itm);
    for (int d=0; d<dimension(); ++d)
      momstream << m[d] << "\t";
    momstream << std::endl;
  } 
  momstream.close();
  std::ofstream clustercoordstream("ClusterCoordinates.dat");
  std::ofstream clustermomstream("ClusterMomenta.dat");
  vertex_iterator itv, itv_end;
  for (boost::tie(itv, itv_end) = cluster_sites(); itv!=itv_end; ++itv) {
    vector_type v(cluster_coordinate(*itv));
    vector_type m = cluster_momentum(*itv);
    for (int d=0; d<dimension(); ++d) {
      clustercoordstream << v[d] << "\t";
      clustermomstream << m[d] << "\t";
    }
    clustermomstream << std::endl;
    clustercoordstream << std::endl;
  }
}


/*template<typename T>
green_function<T> ClusterHelper::transform_into_k_space(const green_function<T>& G) const 
{
  uint nflavor=G.nflavor();
  uint nsite=G.nsite();
  uint nfreq=G.nfreq();
  assert(nsite==n_cluster_sites());
  green_function<T> G_new(nfreq, nsite, nflavor);
  for (spin_t flavor=0; flavor<nflavor; ++flavor) {
    for (site_t i=0; i<nsite; ++i) {
      for (site_t j=0; j<nsite; ++j) {
	for (frequency_t freq=0; freq<nfreq; ++freq) {
	  G_new(freq, i, j, flavor) = 0.;
	}
      }
    }
  }
  vertex_iterator itv_k, itv_k_end;
  for (boost::tie(itv_k, itv_k_end) = cluster_sites(); itv_k!=itv_k_end; ++itv_k) {
    vertex_iterator itv_i, itv_i_end;
    for (boost::tie(itv_i, itv_i_end) = cluster_sites(); itv_i!=itv_i_end; ++itv_i) {
      vertex_iterator itv_j, itv_j_end;
      for (boost::tie(itv_j, itv_j_end) = cluster_sites(); itv_j!=itv_j_end; ++itv_j) {
	double phase = 0.;
	for (int d=0; d<dimension(); ++d) 
	  phase += cluster_momentum(*itv_k)[d]*(cluster_coordinate(*itv_i)[d]-cluster_coordinate(*itv_j)[d]);
	for (frequency_t freq=0; freq<nfreq; ++freq) {
	  for (spin_t flavor=0; flavor<nflavor; ++flavor) {
	    G_new(freq, vertex_index(*itv_k), vertex_index(*itv_k), flavor) 
	      += cos(phase) * G(freq, vertex_index(*itv_i), vertex_index(*itv_j), flavor) / double(nsite);
	  }
	}
      }
    }
  }
  return G_new;
}


template<typename T>
green_function<T> ClusterHelper::transform_into_real_space(const green_function<T>& G)  const
{
  uint nflavor=G.nflavor();
  uint nsite=G.nsite();
  uint nfreq=G.nfreq();
  assert(nsite==n_cluster_sites());
  green_function<T> G_new(nfreq, nsite, nflavor);
  for (spin_t flavor=0; flavor<nflavor; ++flavor) {
    for (site_t i=0; i<nsite; ++i) {
      for (site_t j=0; j<nsite; ++j) {
	for (frequency_t freq=0; freq<nfreq; ++freq) {
	  G_new(freq, i, j, flavor) = 0.;
	}
      }
    }
  }
  vertex_iterator itv_k, itv_k_end;
  for (boost::tie(itv_k, itv_k_end) = cluster_sites(); itv_k!=itv_k_end; ++itv_k) {
    vertex_iterator itv_i, itv_i_end;
    for (boost::tie(itv_i, itv_i_end) = cluster_sites(); itv_i!=itv_i_end; ++itv_i) {
      vertex_iterator itv_j, itv_j_end;
      for (boost::tie(itv_j, itv_j_end) = cluster_sites(); itv_j!=itv_j_end; ++itv_j) {
	double phase = 0.;
	for (int d=0; d<dimension(); ++d) 
	  phase += cluster_momentum(*itv_k)[d]*(cluster_coordinate(*itv_i)[d]-cluster_coordinate(*itv_j)[d]);
	for (frequency_t freq=0; freq<nfreq; ++freq) {
	  for (spin_t flavor=0; flavor<nflavor; ++flavor) {
	    G_new(freq, vertex_index(*itv_i), vertex_index(*itv_j), flavor) 
	      += cos(phase) * G(freq, vertex_index(*itv_k), vertex_index(*itv_k), flavor) / double(nsite);
	  }
	}
      }
    }
  }
  return G_new;
  }*/


/*matsubara_green_function_t ClusterHelper::transform_into_real_space(const matsubara_green_function_t& G) const 
{
  uint nflavor=G.nflavor();
  uint nsite=G.nsite();
  uint nfreq=G.nfreq();
  assert(nsite==n_cluster_sites());
  matsubara_green_function_t G_new(nfreq, nsite, nflavor);
  for (spin_t flavor=0; flavor<nflavor; ++flavor) {
    for (site_t i=0; i<nsite; ++i) {
      for (site_t j=0; j<nsite; ++j) {
	for (frequency_t freq=0; freq<nfreq; ++freq) {
	  G_new(freq, i, j, flavor) = 0.;
	}
      }
    }
  }
  vertex_iterator itv_k, itv_k_end;
  for (boost::tie(itv_k, itv_k_end) = cluster_sites(); itv_k!=itv_k_end; ++itv_k) {
    vertex_iterator itv_i, itv_i_end;
    for (boost::tie(itv_i, itv_i_end) = cluster_sites(); itv_i!=itv_i_end; ++itv_i) {
      vertex_iterator itv_j, itv_j_end;
      for (boost::tie(itv_j, itv_j_end) = cluster_sites(); itv_j!=itv_j_end; ++itv_j) {
	double phase = 0.;
	for (int d=0; d<dimension(); ++d) 
	  phase += cluster_momentum(*itv_k)[d]*(cluster_coordinate(*itv_i)[d]-cluster_coordinate(*itv_j)[d]);
	for (frequency_t freq=0; freq<nfreq; ++freq) {
	  for (spin_t flavor=0; flavor<nflavor; ++flavor) {
	    G_new(freq, vertex_index(*itv_i), vertex_index(*itv_j), flavor).real() 
	      += cos(phase) * G(freq, vertex_index(*itv_k), vertex_index(*itv_k), flavor).real() / double(nsite);
	    G_new(freq, vertex_index(*itv_i), vertex_index(*itv_j), flavor).real() 
	      -= sin(phase) * G(freq, vertex_index(*itv_k), vertex_index(*itv_k), flavor).imag() / double(nsite);
	    G_new(freq, vertex_index(*itv_i), vertex_index(*itv_j), flavor).imag() 
	      += sin(phase) * G(freq, vertex_index(*itv_k), vertex_index(*itv_k), flavor).real() / double(nsite);
	    G_new(freq, vertex_index(*itv_i), vertex_index(*itv_j), flavor).imag() 
	      += cos(phase) * G(freq, vertex_index(*itv_k), vertex_index(*itv_k), flavor).imag() / double(nsite);
	  }
	}
      }
    }
  }
  return G_new;
}*/


itime_green_function_t ClusterHelper::transform_into_real_space(const itime_green_function_t& G) const 
{
  uint nflavor=G.nflavor();
  uint nsite=G.nsite();
  uint nfreq=G.nfreq();
  assert(nsite==n_cluster_sites());
  itime_green_function_t G_new(nfreq, nsite, nflavor);
  for (spin_t flavor=0; flavor<nflavor; ++flavor) {
    for (site_t i=0; i<nsite; ++i) {
      for (site_t j=0; j<nsite; ++j) {
	for (frequency_t freq=0; freq<nfreq; ++freq) {
	  G_new(freq, i, j, flavor) = 0.;
	}
      }
    }
  }
  vertex_iterator itv_k, itv_k_end;
  for (boost::tie(itv_k, itv_k_end) = cluster_sites(); itv_k!=itv_k_end; ++itv_k) {
    vertex_iterator itv_i, itv_i_end;
    for (boost::tie(itv_i, itv_i_end) = cluster_sites(); itv_i!=itv_i_end; ++itv_i) {
      vertex_iterator itv_j, itv_j_end;
      for (boost::tie(itv_j, itv_j_end) = cluster_sites(); itv_j!=itv_j_end; ++itv_j) {
	double phase = 0.;
	for (int d=0; d<dimension(); ++d) 
	  phase += cluster_momentum(*itv_k)[d]*(cluster_coordinate(*itv_i)[d]-cluster_coordinate(*itv_j)[d]);
	for (frequency_t freq=0; freq<nfreq; ++freq) {
	  for (spin_t flavor=0; flavor<nflavor; ++flavor) {
	    G_new(freq, vertex_index(*itv_i), vertex_index(*itv_j), flavor) 
	      += cos(phase) * G(freq, vertex_index(*itv_k), vertex_index(*itv_k), flavor) / double(nsite);
	  }
	}
      }
    }
  }
  return G_new;
}

matsubara_green_function_t ClusterHelper::transform_into_real_space(const matsubara_green_function_t& G) const 
{
  uint nflavor=G.nflavor();
  uint nsite=G.nsite();
  uint nfreq=G.nfreq();
  assert(nsite==n_cluster_sites());
  matsubara_green_function_t G_new(nfreq, nsite, nflavor);
  for (spin_t flavor=0; flavor<nflavor; ++flavor) {
    for (site_t i=0; i<nsite; ++i) {
      for (site_t j=0; j<nsite; ++j) {
	for (frequency_t freq=0; freq<nfreq; ++freq) {
	  G_new(freq, i, j, flavor) = 0.;
	}
      }
    }
  }
  vertex_iterator itv_k, itv_k_end;
  for (boost::tie(itv_k, itv_k_end) = cluster_sites(); itv_k!=itv_k_end; ++itv_k) {
    vertex_iterator itv_i, itv_i_end;
    for (boost::tie(itv_i, itv_i_end) = cluster_sites(); itv_i!=itv_i_end; ++itv_i) {
      vertex_iterator itv_j, itv_j_end;
      for (boost::tie(itv_j, itv_j_end) = cluster_sites(); itv_j!=itv_j_end; ++itv_j) {
	double phase = 0.;
	for (int d=0; d<dimension(); ++d) 
	  phase += cluster_momentum(*itv_k)[d]*(cluster_coordinate(*itv_i)[d]-cluster_coordinate(*itv_j)[d]);
	for (frequency_t freq=0; freq<nfreq; ++freq) {
	  for (spin_t flavor=0; flavor<nflavor; ++flavor) {
	    G_new(freq, vertex_index(*itv_i), vertex_index(*itv_j), flavor) 
	      += cos(phase) * G(freq, vertex_index(*itv_k), vertex_index(*itv_k), flavor) / double(nsite);
	  }
	}
      }
    }
  }
  return G_new;
}

/*matsubara_green_function_t ClusterHelper::transform_into_k_space(const matsubara_green_function_t& G) const 
{
  uint nflavor=G.nflavor();
  uint nsite=G.nsite();
  uint nfreq=G.nfreq();
  assert(nsite==n_cluster_sites());
  matsubara_green_function_t G_new(nfreq, nsite, nflavor);
  for (spin_t flavor=0; flavor<nflavor; ++flavor) {
    for (site_t i=0; i<nsite; ++i) {
      for (site_t j=0; j<nsite; ++j) {
	for (frequency_t freq=0; freq<nfreq; ++freq) {
	  G_new(freq, i, j, flavor) = 0.;
	}
      }
    }
  }
  vertex_iterator itv_k, itv_k_end;
  for (boost::tie(itv_k, itv_k_end) = cluster_sites(); itv_k!=itv_k_end; ++itv_k) {
    vertex_iterator itv_i, itv_i_end;
    for (boost::tie(itv_i, itv_i_end) = cluster_sites(); itv_i!=itv_i_end; ++itv_i) {
      vertex_iterator itv_j, itv_j_end;
      for (boost::tie(itv_j, itv_j_end) = cluster_sites(); itv_j!=itv_j_end; ++itv_j) {
	double phase = 0.;
	for (int d=0; d<dimension(); ++d) 
	  phase += cluster_momentum(*itv_k)[d]*(cluster_coordinate(*itv_i)[d]-cluster_coordinate(*itv_j)[d]);
	for (frequency_t freq=0; freq<nfreq; ++freq) {
	  for (spin_t flavor=0; flavor<nflavor; ++flavor) {
	    G_new(freq, vertex_index(*itv_k), vertex_index(*itv_k), flavor).real() 
	      += cos(phase) * G(freq, vertex_index(*itv_i), vertex_index(*itv_j), flavor).real() / double(nsite);
	    G_new(freq, vertex_index(*itv_k), vertex_index(*itv_k), flavor).real() 
	      += sin(phase) * G(freq, vertex_index(*itv_i), vertex_index(*itv_j), flavor).imag() / double(nsite);
	    G_new(freq, vertex_index(*itv_k), vertex_index(*itv_k), flavor).imag() 
	      -= sin(phase) * G(freq, vertex_index(*itv_i), vertex_index(*itv_j), flavor).real() / double(nsite);
	    G_new(freq, vertex_index(*itv_k), vertex_index(*itv_k), flavor).imag() 
	      += cos(phase) * G(freq, vertex_index(*itv_i), vertex_index(*itv_j), flavor).imag() / double(nsite);
	  }
	}
      }
    }
  }
  return G_new;
}*/


itime_green_function_t ClusterHelper::transform_into_k_space(const itime_green_function_t& G) const 
{
  uint nflavor=G.nflavor();
  uint nsite=G.nsite();
  uint nfreq=G.nfreq();
  assert(nsite==n_cluster_sites());
  itime_green_function_t G_new(nfreq, nsite, nflavor);
  for (spin_t flavor=0; flavor<nflavor; ++flavor) {
    for (site_t i=0; i<nsite; ++i) {
      for (site_t j=0; j<nsite; ++j) {
	for (frequency_t freq=0; freq<nfreq; ++freq) {
	  G_new(freq, i, j, flavor) = 0.;
	}
      }
    }
  }
  vertex_iterator itv_k, itv_k_end;
  for (boost::tie(itv_k, itv_k_end) = cluster_sites(); itv_k!=itv_k_end; ++itv_k) {
    vertex_iterator itv_i, itv_i_end;
    for (boost::tie(itv_i, itv_i_end) = cluster_sites(); itv_i!=itv_i_end; ++itv_i) {
      vertex_iterator itv_j, itv_j_end;
      for (boost::tie(itv_j, itv_j_end) = cluster_sites(); itv_j!=itv_j_end; ++itv_j) {
	double phase = 0.;
	for (int d=0; d<dimension(); ++d) 
	  phase += cluster_momentum(*itv_k)[d]*(cluster_coordinate(*itv_i)[d]-cluster_coordinate(*itv_j)[d]);
	for (frequency_t freq=0; freq<nfreq; ++freq) {
	  for (spin_t flavor=0; flavor<nflavor; ++flavor) {
	    G_new(freq, vertex_index(*itv_k), vertex_index(*itv_k), flavor) 
	      += cos(phase) * G(freq, vertex_index(*itv_i), vertex_index(*itv_j), flavor) / double(nsite);
	  }
	}
      }
    }
  }
  return G_new;
}


matsubara_green_function_t ClusterHelper::transform_into_k_space(const matsubara_green_function_t& G) const 
{
  uint nflavor=G.nflavor();
  uint nsite=G.nsite();
  uint nfreq=G.nfreq();
  assert(nsite==n_cluster_sites());
  matsubara_green_function_t G_new(nfreq, nsite, nflavor);
  for (spin_t flavor=0; flavor<nflavor; ++flavor) {
    for (site_t i=0; i<nsite; ++i) {
      for (site_t j=0; j<nsite; ++j) {
	for (frequency_t freq=0; freq<nfreq; ++freq) {
	  G_new(freq, i, j, flavor) = 0.;
	}
      }
    }
  }
  vertex_iterator itv_k, itv_k_end;
  for (boost::tie(itv_k, itv_k_end) = cluster_sites(); itv_k!=itv_k_end; ++itv_k) {
    vertex_iterator itv_i, itv_i_end;
    for (boost::tie(itv_i, itv_i_end) = cluster_sites(); itv_i!=itv_i_end; ++itv_i) {
      vertex_iterator itv_j, itv_j_end;
      for (boost::tie(itv_j, itv_j_end) = cluster_sites(); itv_j!=itv_j_end; ++itv_j) {
	double phase = 0.;
	for (int d=0; d<dimension(); ++d) 
	  phase += cluster_momentum(*itv_k)[d]*(cluster_coordinate(*itv_i)[d]-cluster_coordinate(*itv_j)[d]);
	for (frequency_t freq=0; freq<nfreq; ++freq) {
	  for (spin_t flavor=0; flavor<nflavor; ++flavor) {
	    G_new(freq, vertex_index(*itv_k), vertex_index(*itv_k), flavor) 
	      += cos(phase) * G(freq, vertex_index(*itv_i), vertex_index(*itv_j), flavor) / double(nsite);
	  }
	}
      }
    }
  }
  return G_new;
}


std::complex<double> *ClusterHelper::exp_iksite(const int sign) const
{
  std::complex<double> *m = new std::complex<double>[n_cluster_sites()*n_cluster_sites()];
  vertex_iterator itv_k, itv_k_end;
  for (boost::tie(itv_k, itv_k_end) = cluster_sites(); itv_k!=itv_k_end; ++itv_k) {
    vertex_iterator itv_i, itv_i_end;
    for (boost::tie(itv_i, itv_i_end) = cluster_sites(); itv_i!=itv_i_end; ++itv_i) {
      double phase = 0.;
      for (int d=0; d<dimension(); ++d) 
	phase += cluster_momentum(*itv_k)[d]*cluster_coordinate(*itv_i)[d];
      phase *= sign;
      m[vertex_index(*itv_k)+n_cluster_sites()*vertex_index(*itv_i)].real() = cos(phase);
      m[vertex_index(*itv_k)+n_cluster_sites()*vertex_index(*itv_i)].imag() = sin(phase);
    }
  }
  return m;
}

matsubara_green_function_t ClusterTransformer::calculate_initial_G0(const alps::Parameters& parms)
{
  int n_matsubara=boost::lexical_cast<int>(parms["NMATSUBARA"]);
  int n_flavor=parms.value_or_default("FLAVORS", 2);
  double beta = static_cast<double>(parms["BETA"]);
  double mu = static_cast<double>(parms["MU"]);
  double h = static_cast<double>(parms["H"]);
  bool paramagnet = static_cast<bool>(parms["PARAMAGNET"]);
  matsubara_green_function_t G_omega(n_matsubara,n_cluster_sites(),n_flavor);
  matsubara_green_function_t G0_omega(n_matsubara,n_cluster_sites(),n_flavor);  
  for (spin_t flavor=0; flavor<n_flavor; ++flavor) {
    for (site_t i=0; i<n_cluster_sites(); ++i) {
      for (site_t j=0; j<n_cluster_sites(); ++j) {
	for (frequency_t freq=0; freq<n_matsubara; ++freq) {
	  G0_omega(freq, i, j, flavor) = (i==j); // zero except for diagonal elements
	}
      }
    }
  }
  G_omega = G0_omega; 
  //evoke transformation with G_omega==G0_omega (selfenergy=0)
  return this->operator()(G_omega, G0_omega, mu, h, beta, paramagnet);
}


matsubara_green_function_t DCATransformer::operator()(matsubara_green_function_t &G_omega, 
						     matsubara_green_function_t &G0_omega, 
						     const double mu, const double h, const double beta, 
						     const bool paramagnet)
{
  assert(h==0);
  uint nflavor=G0_omega.nflavor();
  uint nfreq=G0_omega.nfreq();
  assert(G0_omega.nsite()==n_cluster_sites());
  matsubara_green_function_t selfenergy(G0_omega);
  matsubara_green_function_t G0_omega_new(G0_omega);
  for (spin_t flavor=0; flavor<nflavor; ++flavor) {
    vertex_iterator itv_kc, itv_kc_end;
    for (boost::tie(itv_kc, itv_kc_end) = cluster_sites(); itv_kc!=itv_kc_end; ++itv_kc) {
      for (frequency_t freq=0; freq<nfreq; ++freq) {
	selfenergy(freq, vertex_index(*itv_kc), vertex_index(*itv_kc), flavor) 
	  = 1./G0_omega(freq, vertex_index(*itv_kc), vertex_index(*itv_kc), flavor);
	selfenergy(freq, vertex_index(*itv_kc), vertex_index(*itv_kc), flavor) 
	  -= 1./G_omega(freq, vertex_index(*itv_kc), vertex_index(*itv_kc), flavor);
	std::complex<double> iomega_n(0, (2.*freq+1.)*M_PI/beta);
	G0_omega_new(freq, vertex_index(*itv_kc), vertex_index(*itv_kc), flavor) = 0.;
	momentum_iterator itm_kl, itm_kl_end;
	for (boost::tie(itm_kl, itm_kl_end)=momenta(); itm_kl!=itm_kl_end; ++itm_kl) {
	  vector_type k=momentum(*itm_kl);
	  for (int d=0; d<dimension(); ++d)
	    k[d] += cluster_momentum(*itv_kc)[d];
	  G0_omega_new(freq, vertex_index(*itv_kc), vertex_index(*itv_kc), flavor) 
	    += 1./(iomega_n+mu-epsilon(k)-selfenergy(freq, vertex_index(*itv_kc), vertex_index(*itv_kc), flavor));
	}
	G0_omega_new(freq, vertex_index(*itv_kc), vertex_index(*itv_kc), flavor) /= volume();
	G0_omega_new(freq, vertex_index(*itv_kc), vertex_index(*itv_kc), flavor) = 
	  1. / ( 1./G0_omega_new(freq, vertex_index(*itv_kc), vertex_index(*itv_kc), flavor) 
		 + selfenergy(freq, vertex_index(*itv_kc), vertex_index(*itv_kc), flavor) );
      }
    }
  }
  return G0_omega_new;
}
