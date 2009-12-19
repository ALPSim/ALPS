/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2003-2004 by Salvatore R. Manmana <Salva@theo3.physik.uni-stuttgart.de>,
*                            Reinhard M. Noack <Reinhard.Noack@physik.uni-marburg.de>,
*                            Ian McCulloch <ianmcc@phyics.uq.edu.au>
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

/* $Id$ */

#if !defined(DMRG3_H_)
#define DMRG3_H_

#include <assert.h>
#include <boost/numeric/ublas/matrix.hpp> // matrix library used: uBLAS
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp> // uBLAS input-output 
#include <alps/lattice.h> // interface to ALPS lattice library
#include <alps/lattice/bond_compare.h> // bond comparison utilities
#include <algorithm>
#include <iostream>
#include <iomanip>

// use uBLAS-matrices and vectors with column_major orientation 
typedef boost::numeric::ublas::matrix<
  double, boost::numeric::ublas::column_major> Matrix;
typedef boost::numeric::ublas::vector<double> Vector;

// There is no real reason to make all classes templates on the graph type,
// we make it a typedef instead.

typedef alps::graph_helper<>::graph_type          graph_type;
typedef alps::graph_traits<graph_type>            graph_traits_type;
typedef graph_traits_type::site_descriptor        site_descriptor;
typedef graph_traits_type::site_iterator          site_iterator;
typedef graph_traits_type::bond_descriptor        bond_descriptor;
typedef graph_traits_type::bond_iterator          bond_iterator;
typedef graph_traits_type::neighbor_bond_iterator neighbor_bond_iterator;
typedef int                                       bond_type;

// class 'System': contains all information about the lattice & Hamiltonian terms
// as determined from the parameters.
class System
{
  public:
    // constructs the system from the input parameters
    System(const alps::Parameters& Param);

    // returns the total number of lattice sites in the system
    int num_sites() const;

    // returns the descriptor of the n'th lattice site
    site_descriptor site_at(int n) const;

    // obtain the source and target sites of a given bond of the lattice
    site_descriptor source(bond_descriptor b) const { return helper.source(b); }
    site_descriptor target(bond_descriptor b) const { return helper.target(b); }

    std::pair<site_descriptor, site_descriptor> sites_of_bond(bond_descriptor b) const
      { return std::make_pair(this->source(b), this->target(b)); }

    // returns a [begin, end) pair of all bonds in the lattice
    std::pair<bond_iterator, bond_iterator> bonds() const
    { return helper.bonds(); }

    // returns a [begin, end) pair of the bonds that intersect site s
    std::pair<neighbor_bond_iterator, neighbor_bond_iterator>
      neighbor_bonds(site_descriptor s) const;

    // returns the hopping matrix element of the given bond
    double hopping(bond_descriptor b) const;

    double local_potential(site_descriptor s) const;

    // Displays some information about the system (lattice size, hopping strengths etc)
    void show_info(std::ostream& out);

    void print_debug(std::ostream& out) const;

    graph_type const& get_graph() const { return graph; }

  private:
    typedef int bond_type;  // the 'type' attribute in the <EDGE> element
    typedef int site_type;  // the 'type' attribute in the <VERTEX> element

    typedef alps::property_map<alps::bond_type_t, 
      graph_type, 
      bond_type>::const_type bond_map_type;

    typedef alps::property_map<alps::site_type_t, 
      graph_type, 
      site_type>::const_type site_map_type;

    typedef std::vector<site_descriptor> WormMapType;
    typedef std::map<bond_type, double> HoppingTermsType;
    typedef std::map<site_descriptor, double> LocalPotType;

    alps::graph_helper<graph_type> helper;

    const graph_type&     graph;
    bond_map_type         BondTypeMap;
    site_map_type         SiteTypeMap;
    HoppingTermsType      Hopping;
    LocalPotType          LocalPot;
    WormMapType           WormMap;      // specifies mapping of the lattice sites onto the 1D chain
};

// Class 'Wavefunction' class:
// defines a wavefunction acting on 4 blocks
class Wavefunction 
{
  public:
    // a wave function is a vector with 4 entries.
    // We use an uBLAS-vector because of later convenience, 
    // e.g. when fetching the ground state.

    Wavefunction() 
      : v(4) { }

    Wavefunction(Vector const& V_) : v(V_) { assert(v.size() == 4); }

    double operator[](int i) const { return v[i]; }
    double& operator[](int i) { return v[i]; }

    void print_debug(std::ostream& out) const;
  
  private:
    Vector v;                 
};

// Class 'DensityMatrix': obtain the reduced density matrix of the system block.
// For the single particle problem this reduces to a single eigenvector, which
// is just the (normalized) projection of the wavefunction onto the system block.
class DensityMatrix {                   
  
 public:

  // The entries of the density matrix eigenvector,
  // 'a' corresponds to the system block, 'b' corresponds to the site block.
  double a,b;

  // identify left or right block: 
  enum LR {Left, Right};                  

  // constructor: returns reduced density matrix of 
  // the left or right system block, respectively with the
  // ground state wave function of the actual DMRG-sweep:
  DensityMatrix(const Wavefunction& psi, LR lr) 
  {
    // the reduced density matrix of a non-interacting system is obtained by
    // appropriate 'truncation' of the ground state wave function 
      
    // left block: first two entries of psi
    if(lr == Left)                    
    {
      a = psi[0]; 
      b = psi[1];
    }
    else
    {
      // right block: last two entries of psi
      a = psi[3];
      b = psi[2];
    }

    // normalize
    double Norm = sqrt(a*a+b*b);
    if (Norm == 0)
    {
      a = b = sqrt(0.5);
    }
    else
    {
      a /= Norm;
      b /= Norm;
    }
  }

};

// class 'Block': This holds the wavefunction for a contiguous section of the system,
// and maintains the local Hamiltonian matrix element, and a list of external bonds
// for matrix elements that connect this block to other blocks.
class Block 
{
  private:
    // variable Psi stores the sites that belong with this block
    // together with the wavefunction amplitude at that site
    typedef std::map<int, double> WavefuncType;

  public:
    typedef alps::bond_descriptor_compare_undirected<graph_type> BondDescriptorCompare;

    typedef std::set<bond_descriptor, BondDescriptorCompare> ExternalBondListType;

    typedef WavefuncType::value_type     value_type;
    typedef WavefuncType::const_iterator const_iterator;

    typedef ExternalBondListType::const_iterator external_bond_iterator;

    // constructor: an initial block comprising a single site
    Block(System const& S_, site_descriptor s);

    // constructor: join two blocks together using the reduced density matrix.
    // Typically one of the blocks will be a single site.
    Block(const Block& b1, const Block& b2, const DensityMatrix& dm, bool Debug = false);

    // returns true if the block contains the given site
    bool contains_site(site_descriptor s) const
    {
      return Psi.find(s) != Psi.end();
    }
  
    double local_ham() const { return LocalHam; }

    // returns the value of the wavefunction at the given site
    double wavefunction(site_descriptor site) const
    {
      std::map<int, double>::const_iterator I = Psi.find(site);
      assert(I != Psi.end());
      return I->second;
    }

    // iterates over the block wavefunction
    const_iterator begin() const { return Psi.begin(); }
    const_iterator end() const { return Psi.end(); }

    // Iterate over bonds that have the source site within this block, and the target
    // site outside this block.
    external_bond_iterator begin_external_bonds() const { return ExternalBonds.begin(); }
    external_bond_iterator end_external_bonds() const { return ExternalBonds.end(); }

    void print_debug(std::ostream& out) const;

    System const& system() const { return *Sys; }

  private:
    System const* Sys;

    WavefuncType Psi;

    // The matrix element of the Hamiltonian that represents all terms that
    // are internal to this block
    double LocalHam;

    // All bonds that are external to the block
    ExternalBondListType ExternalBonds;

};

// helper function to find the given site in the blocks provided, 
// returns the block number (0,1,2 or 3), and the wavefunction amplitude
std::pair<int, double>
FindSiteInBlock(site_descriptor s,
                const Block& b1, const Block& b2, 
                const Block& b3, const Block& b4);

// class 'Superblock': contains the Hamiltonian matrix for a full 4-block system.
class Superblock 
{
  public:
    //constructor: build up the superblock hamiltonian
    Superblock(const Block& Left, const Block& LeftSite, 
               const Block& RightSite, const Block& Right, bool Debug = false);
  
    // obtain ground state wave function and energy
    std::pair<double, Wavefunction> GetGroundState(); 

  private:
    Matrix H;
};

// Class 'TotalWavefunction': contains the full wavefunction for the total system,
// used for the final output & calcuating observables etc
class TotalWavefunction
{
  private:
    typedef std::map<site_descriptor, double> WavefuncType;

  public:
    typedef WavefuncType::value_type     value_type;
    typedef WavefuncType::iterator       iterator;
    typedef WavefuncType::const_iterator const_iterator;

    TotalWavefunction() {}

    TotalWavefunction(Block const& b1, Block const& b2, Block const& b3, Block const& b4, 
                      Wavefunction const& Psi);

    iterator begin() { return W.begin(); }
    iterator end() { return W.end(); }

    const_iterator begin() const { return W.begin(); }
    const_iterator end() const { return W.end(); }

    double& operator[](site_descriptor s) { return W[s]; }

    double operator[](site_descriptor s) const 
      { WavefuncType::const_iterator I = W.find(s); return I == W.end() ? 0 : I->first; }

  private:
    WavefuncType W;
};

std::ostream& operator<<(std::ostream& out, TotalWavefunction const& Psi);

#endif // include guard
