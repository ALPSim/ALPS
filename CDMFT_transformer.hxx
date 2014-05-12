template<typename T>
green_function<T> CDMFTTransformer::symmetrize_helper(const green_function<T>& G,
                                                      const std::pair<std::set<pair_class_type>::const_iterator,
                                                      std::set<pair_class_type>::const_iterator> pair_classes) const
{
  unsigned int nflavor=G.nflavor();
  unsigned int nsite=G.nsite();
  unsigned int nfreq=G.nfreq();
  green_function<T> G_new(nfreq, nsite, nflavor, nondiagonal);
  std::set<pair_class_type>::const_iterator itc, itc_end;
  for(boost::tie(itc, itc_end)=pair_classes; itc!=itc_end; itc++) {
    for (spin_t flavor=0; flavor<nflavor; ++flavor) {
      for (frequency_t freq=0; freq<nfreq; ++freq) {
        T av=0;
        for(pair_class_type::const_iterator it=itc->begin(); it!=itc->end(); it++)
          av += G(freq, it->first, it->second, flavor);
        av /= itc->size();
        for(pair_class_type::const_iterator it=itc->begin(); it!=itc->end(); it++)
          G_new(freq, it->first, it->second, flavor) = av;
      }
    }
  }
  return G_new;
}



template<typename T>
nambu_green_function<T> CDMFTTransformer::symmetrize_helper(const nambu_green_function<T>& G,
                                                            const std::pair<std::set<pair_class_type>::const_iterator,
                                                            std::set<pair_class_type>::const_iterator> pair_classes) const
{
  unsigned int nsite=G.nsite();
  unsigned int nfreq=G.nfreq();
  T av;
  nambu_green_function<T> G_new(nfreq, nsite, 2, nondiagonal);
  G_new=G;
  std::set<pair_class_type>::const_iterator itc, itc_end;
  for(boost::tie(itc, itc_end)=pair_classes; itc!=itc_end; itc++) {
    for (spin_t f=0; f<2; ++f) {
      for (frequency_t freq=0; freq<nfreq; ++freq) {
        av=0;
        for(pair_class_type::const_iterator it=itc->begin(); it!=itc->end(); it++){
          av += G(freq, it->first, it->second, f,f);
        }
        av /= itc->size();
        for(pair_class_type::const_iterator it=itc->begin(); it!=itc->end(); it++){
          G_new(freq, it->first, it->second, f,f) = av;
        }
      }
    }
    std::cout<<std::endl;
  }
  
  //transform the dwave factor of K into realspace
  dense_matrix dwave(nsite,nsite);dwave.clear();
  dwave(0,1)=-1;dwave(0,3)=1;dwave(1,0)=-1;dwave(1,2)=1;dwave(2,1)=1;dwave(2,3)=-1;dwave(3,0)=1;dwave(3,2)=-1;
  for(boost::tie(itc, itc_end)=pair_classes; itc!=itc_end; itc++) {
    //std::cout<<"new symmetry: "<<std::endl;
    for (frequency_t freq=0; freq<nfreq; ++freq) {
      T av01=0;
      T av10=0;
      for(pair_class_type::const_iterator it=itc->begin(); it!=itc->end(); it++){
        //if(it->first != it-> second) continue;
        av01+= dwave(it->first,it->second)*G(freq, it->first, it->second, 0,1);
        av10+= dwave(it->first,it->second)*G(freq, it->first, it->second, 1,0);
        //if(freq==0){ std::cout<<"symmetrizing spatial GF: "<<it->first<<" "<<it->second<<" dwave: "<<dwave_prefactor[it->first]<<" entry: "<<dwave_prefactor[it->first]*G(freq, it->first, it->second, 0,1)<<" "<<dwave_prefactor[it->first]*G(freq, it->first, it->second, 1,0)<<std::endl;}
      }
      av01/= itc->size();
      av10/= itc->size();
      //if(freq==0){ std::cout<<"setting average to: "<<av01<<" "<<av10<<std::endl; }
      for(pair_class_type::const_iterator it=itc->begin(); it!=itc->end(); it++){
        //if(it->first != it-> second) continue;
        G_new(freq, it->first, it->second, 0,1) = dwave(it->first,it->second)*av01;
        G_new(freq, it->first, it->second, 1,0) = dwave(it->first,it->second)*av10;
      }
    }
  }
  return G_new;
}

template<class T>
green_function<T> LadderCDMFTTransformer::spin_symmetrizer(const green_function<T>& G, const bool paramagnet) const
{
  green_function<T> G_new(G);
  if (paramagnet) {
    if(G_new.nflavor()>1 && G_new.nflavor()%2==0) {
      for(spin_t flavor=0; flavor<G_new.nflavor(); flavor+=2) {
        for (site_t k1=0; k1<n_cluster_sites(); ++k1) {
          for (site_t k2=0; k2<n_cluster_sites(); ++k2) {
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

