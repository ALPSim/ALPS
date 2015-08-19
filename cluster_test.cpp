#include "gtest/gtest.h"
#include "alpsgraph.h"

TEST(ALPSGraphTest, FakeXMLParser){
  alps_graph G("Fuchs3D-1000B", "cluster.xml", 20);
}
TEST(ALPSGraphTest, ALPSGraphCreation){
  alps::params p;
  p["LATTICE_LIBRARY"] = "./cluster.xml";
  p["LATTICE"] = "Betts2D-4A";
  p["L"] = "40";
  alps_graph G(p);

  EXPECT_EQ(40*40, G.volume());
  EXPECT_EQ(2, G.dimension());
}
TEST(ALPSGraphTest, BasisVectors){
  alps::params p;
  p["LATTICE_LIBRARY"] = "./cluster.xml";
  p["LATTICE"] = "Betts2D-4A";
  p["L"] = "40";
  alps_graph G(p);
  
  alps_graph::basis_vector_iterator_pair BVIP=G.basis_vector_iterators();
  
  alps_graph::basis_vector basis(G.dimension(), std::vector<double>(2, 0.));
  basis[0][0]=2.;   basis[1][1]=2.;
  
  std::vector<std::vector<double> >::const_iterator my_it=basis.begin();
  for(std::vector<std::vector<double> >::const_iterator it=BVIP.first; it!=BVIP.second;++it){
    for(int i=0;i<G.dimension();++i){
      EXPECT_EQ((*my_it)[i], (*it)[i]);
    }
    ++my_it;
  }
}
TEST(ALPSGraphTest, ReciprocalVectors){
  alps::params p;
  p["LATTICE_LIBRARY"] = "./cluster.xml";
  p["LATTICE"] = "Betts2D-4A";
  p["L"] = "40";
  int size=4;
  alps_graph G(p);
  
  alps_graph::basis_vector_iterator_pair RVIP=G.reciprocal_basis_vector_iterators();
  
  alps_graph::basis_vector basis           (G.dimension(), std::vector<double>(2, 0.));
  alps_graph::basis_vector reciprocal_basis(G.dimension(), std::vector<double>(2, 0.));
  basis[0][0]=2.;   basis[1][1]=2.;
  
  for(int i=0;i<G.dimension();++i){
    for(int j=0;j<G.dimension();++j){
      reciprocal_basis[i][j]=basis[i][j]*2*M_PI/size;
    }
  }

  alps_graph::basis_vector_iterator my_it=reciprocal_basis.begin();
  for(alps_graph::basis_vector_iterator it=RVIP.first; it!=RVIP.second;++it){
    for(int i=0;i<G.dimension();++i){
      EXPECT_EQ((*my_it)[i], (*it)[i]);
    }
    ++my_it;
  }
}

#define TWODCLUSTERS1 "Betts2D-1A","Betts2D-4A","Betts2D-8A","Betts2D-8C","Betts2D-9A","Betts2D-9B","Betts2D-10A","Betts2D-10B","Betts2D-11A","Betts2D-12A","Betts2D-12B","Betts2D-12D","Betts2D-13A","Betts2D-14A","Betts2D-14B","Betts2D-15A","Betts2D-15B","Betts2D-16A","Betts2D-16B","Betts2D-16C","Betts2D-16D","Betts2D-16E","Betts2D-17A","Betts2D-17B","Betts2D-18A","Betts2D-18B","Betts2D-18C","Betts2D-18D","Betts2D-19A"
#define TWODCLUSTERS2 "Betts2D-19B","Betts2D-20A","Betts2D-20B","Betts2D-20C","Betts2D-20E","Betts2D-21A","Betts2D-21B","Betts2D-21C","Betts2D-22A","Betts2D-22B","Betts2D-22D","Betts2D-23A","Betts2D-24A","Betts2D-24B","Betts2D-24C","Betts2D-24D","Betts2D-24E","Betts2D-24F","Betts2D-24G","Betts2D-25A"
#define TWODCLUSTERS3 "Betts2D-25B","Betts2D-26A","Betts2D-26B","Betts2D-26C","Betts2D-27A","Betts2D-27B","Betts2D-28A","Betts2D-28B","Betts2D-28C","Betts2D-28D","Betts2D-29A","Betts2D-29B","Betts2D-30A","Betts2D-30B","Betts2D-30C","Betts2D-30D","Betts2D-30G","Betts2D-31A","Betts2D-31B","Betts2D-32A","Betts2D-32B","Betts2D-32C","Betts2D-32D","Betts2D-32E"
#define THREEDCLUSTERS1 "Kent3D-8A","Kent3D-8B","Kent3D-8C","Kent3D-10A","Kent3D-10B","Kent3D-12A","Kent3D-12B","Kent3D-14A","Kent3D-14B","Kent3D-16C","Kent3D-16D","Kent3D-18A","Kent3D-18B","Kent3D-20A","Kent3D-20B","Kent3D-22A","Kent3D-22B","Kent3D-24C","Kent3D-24D","Kent3D-26A","Kent3D-26B","Kent3D-28C","Kent3D-28D","Kent3D-30A","Kent3D-30B","Kent3D-32C","Kent3D-32D","Kent3D-34A","Kent3D-34B","Kent3D-36A","Kent3D-36B","Kent3D-38A","Kent3D-38B","Kent3D-40A","Kent3D-40B","Kent3D-42A","Kent3D-42B","Kent3D-44A","Kent3D-44B","Kent3D-46A"
#define THREEDCLUSTERS2 "Kent3D-46B","Kent3D-48A","Kent3D-48B","Kent3D-50A","Kent3D-50B","Kent3D-52A","Kent3D-52B","Kent3D-54A","Kent3D-54B","Kent3D-56A","Kent3D-56B","Kent3D-58A","Kent3D-58B","Kent3D-60A","Kent3D-60B","Kent3D-62A","Kent3D-62B","Kent3D-64A","Kent3D-64B","Kent3D-66A","Kent3D-66B","Kent3D-68A","Kent3D-68B","Kent3D-70A","Kent3D-70B","Kent3D-72A","Kent3D-72B","Kent3D-74A","Kent3D-74B","Kent3D-76A","Kent3D-76B","Kent3D-80A","Kent3D-80B","Kent3D-84A","Kent3D-84B","Kent3D-88A","Kent3D-88B","Kent3D-90A","Kent3D-90B","Kent3D-96A","Kent3D-100A","Kent3D-102A"
#define THREEDCLUSTERS3 "Fuchs3D-46B","Fuchs3D-48B","Fuchs3D-50B","Fuchs3D-52B","Fuchs3D-54B","Fuchs3D-56B","Fuchs3D-58B","Fuchs3D-60B","Fuchs3D-62B","Fuchs3D-64B","Fuchs3D-66B","Fuchs3D-68B","Fuchs3D-70B","Fuchs3D-72B","Fuchs3D-74B","Fuchs3D-76B","Fuchs3D-78B","Fuchs3D-80B","Fuchs3D-82B","Fuchs3D-84B","Fuchs3D-86B","Fuchs3D-88B","Fuchs3D-90B","Fuchs3D-92B","Fuchs3D-94B","Fuchs3D-96B","Fuchs3D-98B","Fuchs3D-100B","Fuchs3D-102B","Fuchs3D-104B","Fuchs3D-106B","Fuchs3D-108B","Fuchs3D-110B","Fuchs3D-112B","Fuchs3D-114B"
#define THREEDCLUSTERS4 "Fuchs3D-116B","Fuchs3D-118B","Fuchs3D-120B","Fuchs3D-122B","Fuchs3D-124B","Fuchs3D-126B","Fuchs3D-128B","Fuchs3D-130B","Fuchs3D-132B","Fuchs3D-134B","Fuchs3D-136B","Fuchs3D-138B","Fuchs3D-140B","Fuchs3D-142B","Fuchs3D-144B","Fuchs3D-146B","Fuchs3D-148B","Fuchs3D-150B","Fuchs3D-152B","Fuchs3D-154B","Fuchs3D-156B","Fuchs3D-158B","Fuchs3D-160B","Fuchs3D-162B","Fuchs3D-164B","Fuchs3D-166B","Fuchs3D-168B","Fuchs3D-170B","Fuchs3D-172B","Fuchs3D-174B","Fuchs3D-176B","Fuchs3D-178B","Fuchs3D-180B","Fuchs3D-182B","Fuchs3D-184B","Fuchs3D-186B","Fuchs3D-188B","Fuchs3D-190B","Fuchs3D-192B","Fuchs3D-194B","Fuchs3D-196B","Fuchs3D-198B","Fuchs3D-200B","Fuchs3D-202B","Fuchs3D-204B","Fuchs3D-206B"
#define THREEDCLUSTERS5 "Fuchs3D-208B","Fuchs3D-210B","Fuchs3D-212B","Fuchs3D-214B","Fuchs3D-216B","Fuchs3D-218B","Fuchs3D-220B","Fuchs3D-222B","Fuchs3D-224B","Fuchs3D-226B","Fuchs3D-228B","Fuchs3D-230B","Fuchs3D-232B","Fuchs3D-234B","Fuchs3D-236B","Fuchs3D-238B","Fuchs3D-240B","Fuchs3D-242B","Fuchs3D-244B","Fuchs3D-246B","Fuchs3D-248B","Fuchs3D-250B","Fuchs3D-252B","Fuchs3D-254B","Fuchs3D-256B","Fuchs3D-258B","Fuchs3D-260B","Fuchs3D-262B","Fuchs3D-264B","Fuchs3D-266B","Fuchs3D-268B","Fuchs3D-270B","Fuchs3D-272B","Fuchs3D-274B","Fuchs3D-276B","Fuchs3D-278B","Fuchs3D-280B","Fuchs3D-282B","Fuchs3D-284B","Fuchs3D-286B","Fuchs3D-288B","Fuchs3D-290B","Fuchs3D-292B","Fuchs3D-294B","Fuchs3D-296B","Fuchs3D-298B","Fuchs3D-300B","Fuchs3D-302B","Fuchs3D-304B","Fuchs3D-306B"
#define THREEDCLUSTERS6 "Fuchs3D-308B","Fuchs3D-310B","Fuchs3D-312B","Fuchs3D-314B","Fuchs3D-316B","Fuchs3D-318B","Fuchs3D-320B","Fuchs3D-322B","Fuchs3D-324B","Fuchs3D-326B","Fuchs3D-328B","Fuchs3D-330B","Fuchs3D-332B","Fuchs3D-334B","Fuchs3D-336B","Fuchs3D-338B","Fuchs3D-340B","Fuchs3D-342B","Fuchs3D-344B","Fuchs3D-346B","Fuchs3D-348B","Fuchs3D-350B","Fuchs3D-352B","Fuchs3D-354B","Fuchs3D-356B","Fuchs3D-358B","Fuchs3D-360B","Fuchs3D-362B","Fuchs3D-364B","Fuchs3D-366B","Fuchs3D-368B","Fuchs3D-370B","Fuchs3D-372B","Fuchs3D-374B","Fuchs3D-376B","Fuchs3D-378B","Fuchs3D-380B","Fuchs3D-382B","Fuchs3D-384B","Fuchs3D-386B","Fuchs3D-388B","Fuchs3D-390B","Fuchs3D-392B","Fuchs3D-394B","Fuchs3D-396B","Fuchs3D-398B","Fuchs3D-400B","Fuchs3D-402B","Fuchs3D-404B"
#define THREEDCLUSTERS7 "Fuchs3D-406B","Fuchs3D-408B","Fuchs3D-410B","Fuchs3D-412B","Fuchs3D-414B","Fuchs3D-416B","Fuchs3D-418B","Fuchs3D-420B","Fuchs3D-422B","Fuchs3D-424B","Fuchs3D-426B","Fuchs3D-428B","Fuchs3D-430B","Fuchs3D-432B","Fuchs3D-434B","Fuchs3D-436B","Fuchs3D-438B","Fuchs3D-440B","Fuchs3D-442B","Fuchs3D-444B","Fuchs3D-446B","Fuchs3D-448B","Fuchs3D-450B","Fuchs3D-452B","Fuchs3D-454B","Fuchs3D-456B","Fuchs3D-458B","Fuchs3D-460B","Fuchs3D-462B","Fuchs3D-464B","Fuchs3D-466B","Fuchs3D-468B","Fuchs3D-470B","Fuchs3D-472B","Fuchs3D-474B","Fuchs3D-476B","Fuchs3D-478B","Fuchs3D-480B","Fuchs3D-482B","Fuchs3D-484B","Fuchs3D-486B","Fuchs3D-488B","Fuchs3D-490B","Fuchs3D-492B","Fuchs3D-494B","Fuchs3D-496B","Fuchs3D-498B","Fuchs3D-500B"

#define THREEDCLUSTERS8 "Fuchs3D-502B","Fuchs3D-504B","Fuchs3D-506B","Fuchs3D-508B","Fuchs3D-510B","Fuchs3D-512B","Fuchs3D-514B","Fuchs3D-516B","Fuchs3D-518B","Fuchs3D-520B","Fuchs3D-522B","Fuchs3D-524B","Fuchs3D-526B","Fuchs3D-528B","Fuchs3D-530B","Fuchs3D-532B","Fuchs3D-534B","Fuchs3D-536B","Fuchs3D-538B","Fuchs3D-540B","Fuchs3D-542B","Fuchs3D-544B","Fuchs3D-546B","Fuchs3D-548B","Fuchs3D-550B","Fuchs3D-552B","Fuchs3D-554B","Fuchs3D-556B","Fuchs3D-558B","Fuchs3D-560B","Fuchs3D-562B","Fuchs3D-564B","Fuchs3D-566B","Fuchs3D-568B","Fuchs3D-570B","Fuchs3D-572B","Fuchs3D-574B","Fuchs3D-576B","Fuchs3D-578B","Fuchs3D-580B","Fuchs3D-582B","Fuchs3D-584B","Fuchs3D-586B","Fuchs3D-588B","Fuchs3D-590B","Fuchs3D-592B","Fuchs3D-594B","Fuchs3D-596B","Fuchs3D-598B"
#define THREEDCLUSTERS9 "Fuchs3D-600B","Fuchs3D-602B","Fuchs3D-604B","Fuchs3D-606B","Fuchs3D-608B","Fuchs3D-610B","Fuchs3D-612B","Fuchs3D-614B","Fuchs3D-616B","Fuchs3D-618B","Fuchs3D-620B","Fuchs3D-622B","Fuchs3D-624B","Fuchs3D-626B","Fuchs3D-628B","Fuchs3D-630B","Fuchs3D-632B","Fuchs3D-634B","Fuchs3D-636B","Fuchs3D-638B","Fuchs3D-640B","Fuchs3D-642B","Fuchs3D-644B","Fuchs3D-646B","Fuchs3D-648B","Fuchs3D-650B","Fuchs3D-652B","Fuchs3D-654B","Fuchs3D-656B","Fuchs3D-658B","Fuchs3D-660B","Fuchs3D-662B","Fuchs3D-664B","Fuchs3D-666B","Fuchs3D-668B","Fuchs3D-670B","Fuchs3D-672B","Fuchs3D-674B","Fuchs3D-676B","Fuchs3D-678B","Fuchs3D-680B","Fuchs3D-682B","Fuchs3D-684B","Fuchs3D-686B","Fuchs3D-688B","Fuchs3D-690B","Fuchs3D-692B","Fuchs3D-694B"
#define THREEDCLUSTERS10 "Fuchs3D-696B","Fuchs3D-698B","Fuchs3D-700B","Fuchs3D-702B","Fuchs3D-704B","Fuchs3D-706B","Fuchs3D-708B","Fuchs3D-710B","Fuchs3D-712B","Fuchs3D-714B","Fuchs3D-716B","Fuchs3D-718B","Fuchs3D-720B","Fuchs3D-722B","Fuchs3D-724B","Fuchs3D-726B","Fuchs3D-728B","Fuchs3D-730B","Fuchs3D-732B","Fuchs3D-734B","Fuchs3D-736B","Fuchs3D-738B","Fuchs3D-740B","Fuchs3D-742B","Fuchs3D-744B","Fuchs3D-746B","Fuchs3D-748B","Fuchs3D-750B","Fuchs3D-752B","Fuchs3D-754B","Fuchs3D-756B","Fuchs3D-758B","Fuchs3D-760B","Fuchs3D-762B","Fuchs3D-764B","Fuchs3D-766B","Fuchs3D-768B","Fuchs3D-770B","Fuchs3D-772B","Fuchs3D-774B","Fuchs3D-776B","Fuchs3D-778B","Fuchs3D-780B","Fuchs3D-782B","Fuchs3D-784B","Fuchs3D-786B","Fuchs3D-788B","Fuchs3D-790B","Fuchs3D-792B","Fuchs3D-794B"
#define THREEDCLUSTERS11 "Fuchs3D-796B","Fuchs3D-798B","Fuchs3D-800B","Fuchs3D-802B","Fuchs3D-804B","Fuchs3D-806B","Fuchs3D-808B","Fuchs3D-810B","Fuchs3D-812B","Fuchs3D-814B","Fuchs3D-816B","Fuchs3D-818B","Fuchs3D-820B","Fuchs3D-822B","Fuchs3D-824B","Fuchs3D-826B","Fuchs3D-828B","Fuchs3D-830B","Fuchs3D-832B","Fuchs3D-834B","Fuchs3D-836B","Fuchs3D-838B","Fuchs3D-840B","Fuchs3D-842B","Fuchs3D-844B","Fuchs3D-846B","Fuchs3D-848B","Fuchs3D-850B","Fuchs3D-852B","Fuchs3D-854B","Fuchs3D-856B","Fuchs3D-858B","Fuchs3D-860B","Fuchs3D-862B","Fuchs3D-864B","Fuchs3D-866B","Fuchs3D-868B","Fuchs3D-870B","Fuchs3D-872B","Fuchs3D-874B","Fuchs3D-876B","Fuchs3D-878B","Fuchs3D-880B","Fuchs3D-882B","Fuchs3D-884B","Fuchs3D-886B","Fuchs3D-888B","Fuchs3D-890B","Fuchs3D-892B","Fuchs3D-894B"
#define THREEDCLUSTERS12 "Fuchs3D-896B","Fuchs3D-898B","Fuchs3D-900B","Fuchs3D-902B","Fuchs3D-904B","Fuchs3D-906B","Fuchs3D-908B","Fuchs3D-910B","Fuchs3D-912B","Fuchs3D-914B","Fuchs3D-916B","Fuchs3D-918B","Fuchs3D-920B","Fuchs3D-922B","Fuchs3D-924B","Fuchs3D-926B","Fuchs3D-928B","Fuchs3D-930B","Fuchs3D-932B","Fuchs3D-934B","Fuchs3D-936B","Fuchs3D-938B","Fuchs3D-940B","Fuchs3D-942B","Fuchs3D-944B","Fuchs3D-946B","Fuchs3D-948B","Fuchs3D-950B","Fuchs3D-952B","Fuchs3D-954B","Fuchs3D-956B","Fuchs3D-958B","Fuchs3D-958C","Fuchs3D-960B","Fuchs3D-962B","Fuchs3D-964B","Fuchs3D-966B","Fuchs3D-968B","Fuchs3D-970B","Fuchs3D-972B","Fuchs3D-974B","Fuchs3D-976B","Fuchs3D-978B","Fuchs3D-980B","Fuchs3D-982B","Fuchs3D-984B","Fuchs3D-986B","Fuchs3D-988B","Fuchs3D-990B","Fuchs3D-992B"
#define THREEDCLUSTERS13 "Fuchs3D-994B","Fuchs3D-996B","Fuchs3D-998B","Fuchs3D-1000B"

class ClusterTest : public ::testing::TestWithParam<const char*> {
 protected:
  virtual void SetUp() {
    alps::params p;
    L_=20;
    p["LATTICE_LIBRARY"] = "./cluster.xml";
    p["LATTICE"] = GetParam();
    p["L"] = L_;
    G_=new alps_graph(p);
    extract_unit_cell_volume_from_name(GetParam());
  }
  void extract_unit_cell_volume_from_name(const char* name){
    std::string name_str(name);
    int dashpos=name_str.find('-');
    if(dashpos==std::string::npos) throw std::invalid_argument("could not find location of - in name.");
    name_str.erase(0, dashpos+1);
    std::stringstream n(name_str);
    n>>unit_cell_volume_extracted_;
  }

  virtual void TearDown() { delete G_;}
  int unit_cell_volume_extracted_;
  int L_;
  alps_graph *G_;
};

TEST_P(ClusterTest, TwoDClusterCreation){
  EXPECT_EQ(std::pow((double)L_,G_->dimension()), G_->volume());
  EXPECT_EQ(unit_cell_volume_extracted_, G_->unit_cell_volume());
  EXPECT_NEAR(0., G_->validate_orghogonality(), 1.e-8);

  //test the reciprocal vectors, computed from Ashcroft Mermin and by ALPS
  std::vector<std::vector<double> > reciprocal_from_basis=G_->compute_reciprocal_from_basis(G_->basis_vectors());
  std::vector<std::vector<double> > reciprocal_from_alps=G_->reciprocal_basis_vectors();

  for(int i=0;i<reciprocal_from_basis.size();++i){
    for(int j=0;j<reciprocal_from_basis.size();++j){
      EXPECT_NEAR(reciprocal_from_alps[i][j],reciprocal_from_basis[i][j], 1.e-8);
    }
  }
  
  std::vector<std::vector<double> > reciprocal_from_reciprocal_from_alps=G_->compute_reciprocal_from_basis(reciprocal_from_alps);
  for(int i=0;i<reciprocal_from_basis.size();++i){
    for(int j=0;j<reciprocal_from_basis.size();++j){
      EXPECT_NEAR(G_->basis_vectors()[i][j],reciprocal_from_reciprocal_from_alps[i][j], 1.e-8);
    }
  }
  //test the cell origins as computed by ALPS and by hand
  std::pair<std::vector<std::vector<double> >,std::vector<std::vector<double> > > cell_origins_and_offsets=G_->cell_origins_and_offsets_by_hand(L_);
  std::vector<std::vector<double> > cell_origins=cell_origins_and_offsets.first;
  std::vector<std::vector<double> > cell_offsets=cell_origins_and_offsets.second;

  alps_graph::basis_vector_iterator_pair cell_origins_from_alps=G_->cell_origins();
  std::pair<std::vector<std::vector<int> >::const_iterator,std::vector<std::vector<int> >::const_iterator>  cell_offsets_from_alps=G_->cell_offsets();

  EXPECT_EQ(pow((double)L_, G_->dimension()), cell_origins.size());
  EXPECT_EQ(pow((double)L_, G_->dimension()), cell_offsets.size());

  alps_graph::basis_vector_iterator it=cell_origins_from_alps.first;
  for(int i=0;i<cell_origins.size();++i){
    for(int k=0;k<G_->dimension();++k){
      EXPECT_NEAR((*it)[k],cell_origins[i][k], 1.e-8);
    }
    it++;
  }
  EXPECT_EQ(it,cell_origins_from_alps.second);
  
  //test the cell offsets as computed by ALPS and by hand
  std::vector<std::vector<int> >::const_iterator it2=cell_offsets_from_alps.first;
  for(int i=0;i<cell_offsets.size();++i){
    for(int k=0;k<G_->dimension();++k){
      EXPECT_NEAR((*it2)[k],cell_offsets[i][k], 1.e-8);
    }
    it2++;
  }
  EXPECT_EQ(it2,cell_offsets_from_alps.second);

  //test the momenta as computed by ALPS and by hand
  std::vector<std::vector<double> > momenta_by_hand=G_->compute_momenta_by_hand(L_);
  alps_graph::momenta_iterator_pair momenta_from_alps=G_->cell_momenta();
  std::vector<std::vector<double> >::const_iterator it3=momenta_from_alps.first;
  for(int i=0;i<momenta_by_hand.size();++i){
    for(int k=0;k<G_->dimension();++k){
      EXPECT_NEAR((*it3)[k],momenta_by_hand[i][k], 1.e-8);
    }
    it3++;
  }
  EXPECT_EQ(it3,momenta_from_alps.second);

}
INSTANTIATE_TEST_CASE_P(ClusterCreation2D1 , ClusterTest,::testing::Values(TWODCLUSTERS1));
INSTANTIATE_TEST_CASE_P(ClusterCreation2D2 , ClusterTest,::testing::Values(TWODCLUSTERS2));
INSTANTIATE_TEST_CASE_P(ClusterCreation2D3 , ClusterTest,::testing::Values(TWODCLUSTERS3));
INSTANTIATE_TEST_CASE_P(ClusterCreation3D1 , ClusterTest,::testing::Values(THREEDCLUSTERS1));
INSTANTIATE_TEST_CASE_P(ClusterCreation3D2 , ClusterTest,::testing::Values(THREEDCLUSTERS2));
INSTANTIATE_TEST_CASE_P(ClusterCreation3D3 , ClusterTest,::testing::Values(THREEDCLUSTERS3));
INSTANTIATE_TEST_CASE_P(ClusterCreation3D4 , ClusterTest,::testing::Values(THREEDCLUSTERS4));
INSTANTIATE_TEST_CASE_P(ClusterCreation3D5 , ClusterTest,::testing::Values(THREEDCLUSTERS5));
INSTANTIATE_TEST_CASE_P(ClusterCreation3D6 , ClusterTest,::testing::Values(THREEDCLUSTERS6));
INSTANTIATE_TEST_CASE_P(ClusterCreation3D7 , ClusterTest,::testing::Values(THREEDCLUSTERS7));
INSTANTIATE_TEST_CASE_P(ClusterCreation3D8 , ClusterTest,::testing::Values(THREEDCLUSTERS8));
INSTANTIATE_TEST_CASE_P(ClusterCreation3D9 , ClusterTest,::testing::Values(THREEDCLUSTERS9));
INSTANTIATE_TEST_CASE_P(ClusterCreation3D10, ClusterTest,::testing::Values(THREEDCLUSTERS10));
INSTANTIATE_TEST_CASE_P(ClusterCreation3D11, ClusterTest,::testing::Values(THREEDCLUSTERS11));
INSTANTIATE_TEST_CASE_P(ClusterCreation3D12, ClusterTest,::testing::Values(THREEDCLUSTERS12));
INSTANTIATE_TEST_CASE_P(ClusterCreation3D13, ClusterTest,::testing::Values(THREEDCLUSTERS13));

