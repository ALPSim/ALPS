/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2015 by Emanuel Gull <egull@umich.edu>
 *
 *
 *****************************************************************************/


#include<sstream>
#include<cmath>
#include "alpsgraph.h"

std::vector<std::vector<double> > alps_graph::lookup_basis_from_name(const std::string &cluster_name, const std::string &cluster_library) const{
  std::vector<std::vector<double> > basis;
  //for now we use the ALPS xml file and write a poor man's parser.
  std::string line;
  std::string superlattice_string=cluster_name+" superlattice";
  std::ifstream xml_file(cluster_library.c_str());
  int dimension;
  if(!xml_file.good()) throw std::runtime_error("cluster library: "+cluster_library+" could not be opened.");
  while(getline(xml_file, line)) {
    if (line.find(superlattice_string, 0) != std::string::npos) {
      //on the same line we expect 'dimension'
      if     (line.find("dimension=\"1\"")!=std::string::npos) dimension=1;
      else if(line.find("dimension=\"2\"")!=std::string::npos) dimension=2;
      else if(line.find("dimension=\"3\"")!=std::string::npos) dimension=3;
      else throw std::runtime_error("problem finding dimensionality.");
      getline(xml_file, line); //this should be 'basis'
      for(int i=0;i<dimension;++i){
        getline(xml_file, line); //this should be the first vector
        int svecpos=line.find("<VECTOR>"); if(svecpos==std::string::npos) throw std::runtime_error("could not find <VECTOR> entry in xml line: "+line);
        int evecpos=line.find("</VECTOR>");if(evecpos==std::string::npos) throw std::runtime_error("could not find </VECTOR> entry in xml line: "+line);
        std::stringstream sstr(line.substr(svecpos+8, evecpos-svecpos-8));
        std::vector<double> basis_vector;
        for(int j=0;j<dimension;++j){
          int basis_element;
          sstr>>basis_element;
          basis_vector.push_back(basis_element);
        }
        basis.push_back(basis_vector);
      }
      return basis;
    }
  }
  throw std::logic_error("cluster"+cluster_name+" not found in library: "+cluster_library);
}
double alps_graph::validate_orthogonality() const{
  double deviation_from_orthogonality=0.;
  for(int i=0;i< basis_vectors_.size();++i){
    for(int j=0;j< reciprocal_basis_vectors_.size();++j){
      double inner_product=0.;
      for(int k=0;k<dimension();++k){
        inner_product+= basis_vectors_[i][k]* reciprocal_basis_vectors_[j][k];
      }
      deviation_from_orthogonality+=inner_product-2.*M_PI*(i==j);
    }
  }
  return deviation_from_orthogonality;
}

void alps_graph::compute_unit_cell_volume(){
  if(dimension_==1){ unit_cell_volume_= basis_vectors_[0][0]; }
  else if(dimension_==2){
    unit_cell_volume_=std::abs(basis_vectors_[0][0]* basis_vectors_[1][1]- basis_vectors_[0][1]* basis_vectors_[1][0]);
  }else if(dimension_==3){
    std::vector<double> c=cross(basis_vectors_[0], basis_vectors_[1]);
    unit_cell_volume_=std::abs(basis_vectors_[2][0]*c[0]+ basis_vectors_[2][1]*c[1]+ basis_vectors_[2][2]*c[2]);
  }
}
std::vector<std::vector<double> > alps_graph::compute_reciprocal_from_basis(const std::vector<std::vector<double> > &basis) const{
  std::vector<std::vector<double> > basis_in_3d(3, std::vector<double>(3, 0.));
  std::vector<std::vector<double> > reciprocal_basis_in_3d(3, std::vector<double>(3, 0.));
  std::vector<std::vector<double> > reciprocal_basis(dimension(), std::vector<double>(dimension(), 0.));
  
  for(int i=0;i<dimension();++i){
    for(int j=0;j<dimension();++j){
      basis_in_3d[i][j]=basis[i][j];
    }
  }for(int i=dimension();i<3;++i){
    basis_in_3d[i][i]=1.;
  }
  std::vector<double> c0=cross(basis_in_3d[1], basis_in_3d[2]);
  std::vector<double> c1=cross(basis_in_3d[2], basis_in_3d[0]);
  std::vector<double> c2=cross(basis_in_3d[0], basis_in_3d[1]);
  double numerator=c0[0]*basis_in_3d[0][0]+c0[1]*basis_in_3d[0][1]+c0[2]*basis_in_3d[0][2];

  reciprocal_basis_in_3d[0][0]=2*M_PI/numerator*c0[0];
  reciprocal_basis_in_3d[0][1]=2*M_PI/numerator*c0[1];
  reciprocal_basis_in_3d[0][2]=2*M_PI/numerator*c0[2];
  reciprocal_basis_in_3d[1][0]=2*M_PI/numerator*c1[0];
  reciprocal_basis_in_3d[1][1]=2*M_PI/numerator*c1[1];
  reciprocal_basis_in_3d[1][2]=2*M_PI/numerator*c1[2];
  reciprocal_basis_in_3d[2][0]=2*M_PI/numerator*c2[0];
  reciprocal_basis_in_3d[2][1]=2*M_PI/numerator*c2[1];
  reciprocal_basis_in_3d[2][2]=2*M_PI/numerator*c2[2];
  
  for(int i=0;i<dimension();++i){
    for(int j=0;j<dimension();++j){
      reciprocal_basis[i][j]=reciprocal_basis_in_3d[i][j];
    }
  }
  return reciprocal_basis;
}

bool alps_graph::operator==(const alps_graph &g) const{
  if(volume_!= g.volume_) return false;
  if(unit_cell_volume_!= g.unit_cell_volume_) return false;
  if(dimension_!= g.dimension_) return false;
  
  for(int i=0;i< basis_vectors_.size();++i){
    for(int j=0;j<dimension_;++j){
      if(std::abs(basis_vectors_[i][j]-g.basis_vectors_[i][j])>1.e-8) return false;
    }
  }
  
  for(int i=0;i< reciprocal_basis_vectors_.size();++i){
    for(int j=0;j<dimension_;++j){
      if(std::abs(reciprocal_basis_vectors_[i][j]-g.reciprocal_basis_vectors_[i][j])>1.e-8) return false;
    }
  }
  for(int i=0;i<cell_origins_.size();++i){
    for(int j=0;j<dimension_;++j){
      if(std::abs(cell_origins_[i][j]-g.cell_origins_[i][j])>1.e-8) return false;
    }
  }
  for(int i=0;i<cell_offsets_.size();++i){
    for(int j=0;j<dimension_;++j){
      if(std::abs(cell_offsets_[i][j]-g.cell_offsets_[i][j])>1.e-8) return false;
    }
  }
  return true;
}

void alps_graph::build(const std::vector<std::vector<double> > basis, int L){
  basis_vectors_ =basis;
  dimension_=basis.size(); if(dimension_>3 || dimension_<1) throw std::invalid_argument("clusters work for dimensions 1, 2, and 3.");
  volume_=dimension_==1?L:dimension_==2?L*L:L*L*L;
  compute_unit_cell_volume();
  
  reciprocal_basis_vectors_ =compute_reciprocal_from_basis(basis);
  std::pair<std::vector<std::vector<double> >,std::vector<std::vector<int> > > origins_and_offsets=compute_cell_origins_and_offsets(L);
  cell_origins_=origins_and_offsets.first;
  cell_offsets_=origins_and_offsets.second;
  momenta_=compute_momenta(L);
}

std::vector<std::vector<double> > alps_graph::compute_momenta(int L)const{
  std::vector<std::vector<double> > momenta;
  std::vector<double> momentum(dimension(), 0.);
  if(dimension()==3){
    for(int i=0;i<L;++i){
      for(int j=0;j<L;++j){
        for(int k=0;k<L;++k){
          int offset_x=i>L/2?i-L:i;
          int offset_y=j>L/2?j-L:j;
          int offset_z=k>L/2?k-L:k;
          momentum[0]=
                  reciprocal_basis_vectors_[0][0]*offset_x/(double)L+ reciprocal_basis_vectors_[1][0]*offset_y/(double)L+
                  reciprocal_basis_vectors_[2][0]*offset_z/(double)L;
          momentum[1]=
                  reciprocal_basis_vectors_[0][1]*offset_x/(double)L+ reciprocal_basis_vectors_[1][1]*offset_y/(double)L+
                  reciprocal_basis_vectors_[2][1]*offset_z/(double)L;
          momentum[2]=
                  reciprocal_basis_vectors_[0][2]*offset_x/(double)L+ reciprocal_basis_vectors_[1][2]*offset_y/(double)L+
                  reciprocal_basis_vectors_[2][2]*offset_z/(double)L;
          momenta.push_back(momentum);
        }
      }
    }
  }else if(dimension()==2){
    for(int i=0;i<L;++i){
      for(int j=0;j<L;++j){
        int offset_x=i>L/2?i-L:i;
        int offset_y=j>L/2?j-L:j;
        momentum[0]=
                reciprocal_basis_vectors_[0][0]*offset_x/(double)L+ reciprocal_basis_vectors_[1][0]*offset_y/(double)L;
        momentum[1]=
                reciprocal_basis_vectors_[0][1]*offset_x/(double)L+ reciprocal_basis_vectors_[1][1]*offset_y/(double)L;
        momenta.push_back(momentum);
      }
    }
  }else if(dimension()==1){
    for(int i=0;i<L;++i){
      int offset_x=i>L/2?i-L:i;
      momentum[0]= reciprocal_basis_vectors_[0][0]*offset_x/(double)L;
      momenta.push_back(momentum);
    }
  }
  return momenta;
}
std::pair<std::vector<std::vector<double> >,std::vector<std::vector<int> > > alps_graph::compute_cell_origins_and_offsets(int L){
  std::vector<std::vector<double> > cell_origins;
  std::vector<std::vector<int   > > cell_offsets;
  std::vector<double> cell_origin(dimension(), 0.);
  std::vector<int   > cell_offset(dimension(), 0.);
  if(dimension()==3){
    for(int i=0;i<L;++i){
      for(int j=0;j<L;++j){
        for(int k=0;k<L;++k){
          cell_origin[0] =i* basis_vectors_[0][0]+j* basis_vectors_[1][0]+k* basis_vectors_[2][0];
          cell_origin[1] =i* basis_vectors_[0][1]+j* basis_vectors_[1][1]+k* basis_vectors_[2][1];
          cell_origin[2] =i* basis_vectors_[0][2]+j* basis_vectors_[1][2]+k* basis_vectors_[2][2];
          cell_offset[0] =i; cell_offset[1] =j; cell_offset[2] =k;
          cell_origins.push_back(cell_origin);
          cell_offsets.push_back(cell_offset);
        }
      }
    }
  }else if(dimension()==2){
    for(int i=0;i<L;++i){
      for(int j=0;j<L;++j){
        cell_origin[0] =i* basis_vectors_[0][0]+j* basis_vectors_[1][0];
        cell_origin[1] =i* basis_vectors_[0][1]+j* basis_vectors_[1][1];
        cell_offset[0] =i; cell_offset[1] =j;
        cell_origins.push_back(cell_origin);
        cell_offsets.push_back(cell_offset);
      }
    }
  }else if(dimension()==1){
    for(int i=0;i<L;++i){
      cell_origin[0] =i* basis_vectors_[0][0];
      cell_offset[0] =i;
      cell_origins.push_back(cell_origin);
      cell_offsets.push_back(cell_offset);
    }
  }
  return std::make_pair(cell_origins, cell_offsets);
}

