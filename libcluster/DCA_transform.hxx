template<typename T>
green_function<T> DCATransformer::symmetrize_helper(const green_function<T>& G,
                                                    const std::pair<std::set<pair_class_type>::const_iterator,
                                                    std::set<pair_class_type>::const_iterator> pair_classes) const
{
  unsigned int nflavor=G.nflavor();
  unsigned int nsite=G.nsite();
  unsigned int nfreq=G.nfreq();
  green_function<T> G_new(nfreq, nsite, nflavor, G.shape());
  std::set<pair_class_type>::const_iterator itc, itc_end;
  for(boost::tie(itc, itc_end)=pair_classes; itc!=itc_end; itc++) {
    for (spin_t flavor=0; flavor<nflavor; ++flavor) {
      for (frequency_t freq=0; freq<nfreq; ++freq) {
        T av=0;
        for(pair_class_type::const_iterator it=itc->begin(); it!=itc->end(); it++){
          if(G.shape()==diagonal && it->first != it->second) continue;
          av += G(freq, it->first, it->second, flavor);
        }
        av /= itc->size();
        for(pair_class_type::const_iterator it=itc->begin(); it!=itc->end(); it++){
          if(G.shape()==diagonal && it->first != it->second) continue;
          G_new(freq, it->first, it->second, flavor) = av;
        }
      }
    }
  }
  return G_new;
}

template<typename T>
nambu_green_function<T> DCATransformer::symmetrize_helper(const nambu_green_function<T>& G,
                                                          const std::pair<std::set<pair_class_type>::const_iterator,
                                                          std::set<pair_class_type>::const_iterator> pair_classes) const
{
  //  std::cout<<"ignoring symmetrization."<<std::endl;
  //  return G;
  
  unsigned int nsite=G.nsite();
  unsigned int nfreq=G.nfreq();
  nambu_green_function<T> G_new(nfreq, nsite, 2, diagonal);
  G_new=G;
  std::set<pair_class_type>::const_iterator itc, itc_end;
  //std::cout<<"doing spatial symmetries: "<<std::endl;
  for(boost::tie(itc, itc_end)=pair_classes; itc!=itc_end; itc++) {
    for (spin_t f=0; f<2; ++f) {
      //std::cout<<"new symmetry: "<<std::endl;
      for (frequency_t freq=0; freq<nfreq; ++freq) {
        T av=0;
        for(pair_class_type::const_iterator it=itc->begin(); it!=itc->end(); it++){
          //if(freq==0){ std::cout<<"symmetrizing spatial GF: "<<it->first<<" "<<it->second<<" entry: "<<G(freq, it->first, it->second, f,f)<<std::endl;}
          av += G(freq, it->first, it->second, f,f);
        }
        av /= itc->size();
        for(pair_class_type::const_iterator it=itc->begin(); it!=itc->end(); it++)
          G_new(freq, it->first, it->second, f,f) = av;
      }
    }
  }
  
  //this is the k-space symmetrization respecting d-wave symmetry
  std::vector<double> dwave_prefactor(nsite);
  for(unsigned int i=0;i<nsite;++i){
    double dwave=cos(cluster_momentum(i,0))-cos(cluster_momentum(i,1));
    if(dwave < -1.e-5) dwave_prefactor[i]=-1;
    else if(dwave >  1.e-5) dwave_prefactor[i]=1;
    else dwave_prefactor[i]=0.;
  }
  for(boost::tie(itc, itc_end)=pair_classes; itc!=itc_end; itc++) {
    //std::cout<<"new symmetry: "<<std::endl;
    for (frequency_t freq=0; freq<nfreq; ++freq) {
      T av01=0;
      T av10=0;
      for(pair_class_type::const_iterator it=itc->begin(); it!=itc->end(); it++){
        if(it->first != it-> second) continue;
        av01+= dwave_prefactor[it->first]*G(freq, it->first, it->second, 0,1);
        av10+= dwave_prefactor[it->first]*G(freq, it->first, it->second, 1,0);
        //if(freq==0){ std::cout<<"symmetrizing spatial GF: "<<it->first<<" "<<it->second<<" dwave: "<<dwave_prefactor[it->first]<<" entry: "<<dwave_prefactor[it->first]*G(freq, it->first, it->second, 0,1)<<" "<<dwave_prefactor[it->first]*G(freq, it->first, it->second, 1,0)<<std::endl;}
      }
      av01/= itc->size();
      av10/= itc->size();
      //if(freq==0){ std::cout<<"setting average to: "<<av01<<" "<<av10<<std::endl; }
      for(pair_class_type::const_iterator it=itc->begin(); it!=itc->end(); it++){
        if(it->first != it-> second) continue;
        G_new(freq, it->first, it->second, 0,1) = dwave_prefactor[it->first]*av01;
        G_new(freq, it->first, it->second, 1,0) = dwave_prefactor[it->first]*av10;
      }
    }
  }
  return G_new;
}

