/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2003-2005 by Salvatore R. Manmana <Salva@theo3.physik.uni-stuttgart.de>,
*                            Reinhard M. Noack <Reinhard.Noack@physik.uni-marburg.de>,
*                            Ian McCulloch <ianmcc@physics.uq.edu.au>,
*                            Matthias Troyer <troyer@comp-phys.org>
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

#include "dmrg3.h"

#include <alps/config.h> // needed to set up correct bindings
#include <boost/numeric/bindings/lapack/driver/syev.hpp> // interface to LAPACK routines
#include <boost/numeric/bindings/ublas.hpp>
#include <boost/numeric/bindings/upper.hpp>
#include <fstream>
#include <sstream>

//
// System
//

System::System(const alps::Parameters& Param)
  : helper(Param), graph(helper.graph()),
    BondTypeMap(alps::get_or_default(alps::bond_type_t(), graph, 0)),
    SiteTypeMap(alps::get_or_default(alps::site_type_t(), graph, 0))
{

  if (helper.inhomogeneous())
    boost::throw_exception(std::runtime_error("Disordered lattices are not supported by the DMRG program.\n"));

  // iterate over all bonds in the graph and see what bond types
  // we need to get hopping terms for
  std::set<bond_type> BondTypesInGraph;
  bond_iterator bIter, bEnd;
  for (boost::tie(bIter, bEnd) = bonds(); bIter != bEnd; ++bIter)
    BondTypesInGraph.insert(BondTypeMap(*bIter));
    
  // retrieve the hopping terms from the parameters
  for (std::set<bond_type>::const_iterator btIter = BondTypesInGraph.begin();
       btIter != BondTypesInGraph.end(); ++btIter)
  {
    // For bond type 0, either 't0' or 't' will suffice; if neither is present then default to 1.0.
    if (*btIter == 0)
      Hopping[0] = Param.value_or_default("t", Param.value_or_default("t0", 1.0));
    else
      // other parameters must have the form tN where N is the bond type, and these parameters default to 0
      Hopping[*btIter] =
        Param.value_or_default("t"+boost::lexical_cast<std::string>(*btIter), 0.0);

    Hopping[*btIter] = -Hopping[*btIter];  // hopping matrix elements are -t
  }

  // iterate over all sites and determine the local potential.  V applies to all sites,
  // VN applies only to site type N
  site_iterator sIter, sEnd;
  for (boost::tie(sIter, sEnd) = sites(graph); sIter != sEnd; ++sIter)
  {
    alps::Parameters p(Param);
    p << helper.coordinate_as_parameter(*sIter);
    LocalPot[*sIter] = alps::evaluate<double>(
         Param.value_or_default("V"+boost::lexical_cast<std::string>(SiteTypeMap(*sIter)), "0"), 
         p)
      + alps::evaluate<double>(Param.value_or_default("V", "0"), p);
  }

  // Calculate the ordering of the lattice sites.
  // Eventually we want to be able to specify this from the parameters somehow,
  // but for the time being we take advantage of the fact that the lattice
  // library gives us a default ordering.
  // *** Just so that any bugs due to accidentally using a site_descriptor instead
  // of a site number will be more obvious, we reverse the mapping ;-)
  for (int i = this->num_sites()-1; i >= 0; --i)
  {
    WormMap.push_back(site_descriptor(i));
  }

  // Uncomment the following line to slow down the convergence by using
  // a random mapping of the lattice sites onto a 1D chain.
  // This usually makes the interactions very non-local.
  // std::random_shuffle(WormMap.begin(), WormMap.end());
}

void System::show_info(std::ostream& out)
{
  out << "System has " << this->num_sites() << " sites.\n";
  out << "Parameters are:\n";
  for (HoppingTermsType::const_iterator I = Hopping.begin(); I != Hopping.end(); ++I)
  {
    out << "  t";
    if (I->first != 0) out << I->first;
    out << " = " << -I->second << '\n';
  }

  out << "\nLocal potentials are:\n";
  for (LocalPotType::const_iterator I = LocalPot.begin(); I != LocalPot.end(); ++I)
  {
    out << I->first << "  " << I->second << '\n';
  }
  out << std::endl;
}

int System::num_sites() const
{
  return helper.num_sites();
}

site_descriptor System::site_at(int n) const
{
  assert(0 <= n && n < this->num_sites());
  return WormMap[n];
}

std::pair<neighbor_bond_iterator, neighbor_bond_iterator>
System::neighbor_bonds(site_descriptor s) const
{
  return helper.neighbor_bonds(s);
}

double System::hopping(bond_descriptor b) const
{
  HoppingTermsType::const_iterator I = Hopping.find(BondTypeMap[b]);
  assert(I != Hopping.end());  // this should never fail....
  return I->second;
}

double System::local_potential(site_descriptor s) const
{
  LocalPotType::const_iterator I = LocalPot.find(s);
  return I == LocalPot.end() ? 0 : I->second;
}

//
// Wavefunction
//

void Wavefunction::print_debug(std::ostream& out) const
{
  out << "Wavefunction is: " << v << std::endl;
}

//
// TotalWavefunction
//

TotalWavefunction::TotalWavefunction(Block const& b1, Block const& b2, Block const& b3, Block const& b4, 
                                     Wavefunction const& Psi)
{
  // The sign of the wavefunction is arbitrary, but for consistency we make a choice
  // that the sum should be > 0
  double sum = 0;
  for (Block::const_iterator I = b1.begin(); I != b1.end(); ++I)
  {
    W[I->first] = Psi[0] * I->second;
    sum += Psi[0] * I->second;
  }
  for (Block::const_iterator I = b2.begin(); I != b2.end(); ++I)
  {
    W[I->first] = Psi[1] * I->second;
    sum += Psi[1] * I->second;
  }
  for (Block::const_iterator I = b3.begin(); I != b3.end(); ++I)
  {
    W[I->first] = Psi[2] * I->second;
    sum += Psi[2] * I->second;
  }
  for (Block::const_iterator I = b4.begin(); I != b4.end(); ++I)
  {
    W[I->first] = Psi[3] * I->second;
    sum += Psi[3] * I->second;
  }

  // negate the wavefunction if the sum is negative
  if (sum < 0)
  {
    for (iterator I = this->begin(); I != this->end(); ++I)
    {
      I->second = -I->second;
    }
  }
}

std::ostream& operator<<(std::ostream& out, TotalWavefunction const& Psi)
{
  out << "#site          amplitude\n";
  for (TotalWavefunction::const_iterator I = Psi.begin(); I != Psi.end(); ++I)
  {
    out << std::setw(5) << I->first << ' ' << std::setw(18) << I->second << '\n';
  }
  return out;
}

//
// DensityMatrix
//

//
// Block
//

Block::Block(System const& S_, site_descriptor s)
  : Sys(&S_), LocalHam(Sys->local_potential(s)), ExternalBonds(BondDescriptorCompare(&Sys->get_graph()))
{
  // Initialize the external bonds.
  std::pair<neighbor_bond_iterator, neighbor_bond_iterator> Bonds(Sys->neighbor_bonds(s));
  ExternalBonds.insert(Bonds.first, Bonds.second);

  Psi[s] = 1.0;
}

std::pair<int, double>
FindSiteInBlock(site_descriptor s,
                const Block& b1, const Block& b2, 
                const Block& b3, const Block& b4)
{
  if (b1.contains_site(s)) return std::make_pair(0, b1.wavefunction(s));
  if (b2.contains_site(s)) return std::make_pair(1, b2.wavefunction(s));
  if (b3.contains_site(s)) return std::make_pair(2, b3.wavefunction(s));
  if (b4.contains_site(s)) return std::make_pair(3, b4.wavefunction(s));
  // if we got here then the site is not in any block
  return std::make_pair(-1, 0.0);
}

Block::Block(const Block& b1, const Block& b2, const DensityMatrix& dm, bool Debug)
  : Sys(b1.Sys), ExternalBonds(BondDescriptorCompare(&Sys->get_graph()))
{
  if (Debug) std::cout << "Constructing new block\n";

  assert(b1.Sys == b2.Sys);  // sanity check that b1 and b2 refer to the same system

  // our wavefunction is just the combination of b1.Psi and b2.Psi
  // scaled by the appropriate basis elements.
  // Also verify that the sites in b1 and b2 don't overlap.
  for (WavefuncType::const_iterator I = b1.Psi.begin(); I != b1.Psi.end(); ++I)
  {
    Psi[I->first] = dm.a * I->second;
    if (Debug) std::cout << "New block inherits site " << I->first << " from b1\n";
  }

  for (WavefuncType::const_iterator I = b2.Psi.begin(); I != b2.Psi.end(); ++I)
  {
    assert(Psi.count(I->first) == 0);
    Psi[I->first] = dm.b * I->second;
    if (Debug) std::cout << "New block inherits site " << I->first << " from b2\n";
  }

  // Set matrix element for the energy of the block.  This is the sum of the
  // local hamiltonian of each block, plus any cross terms between
  // block b1 and b2.  
  LocalHam = dm.a * dm.a * b1.LocalHam + dm.b * dm.b * b2.LocalHam;

  // The terms that are local to the new block are the set intersection of
  // b1.ExternalBonds and b2.ExternalBonds
  std::list<bond_descriptor> b1b2Interactions;
  std::set_intersection(b1.ExternalBonds.begin(), b1.ExternalBonds.end(), 
                        b2.ExternalBonds.begin(), b2.ExternalBonds.end(),
                        std::back_inserter(b1b2Interactions),
                        BondDescriptorCompare(&Sys->get_graph()));

  for (std::list<bond_descriptor>::const_iterator I = b1b2Interactions.begin();
       I != b1b2Interactions.end(); ++I)
  {
    // Factor of 2 here because it is <source|H|target> + <target|H|source>
    LocalHam += 2 * Sys->hopping(*I) * Psi[Sys->source(*I)] * Psi[Sys->target(*I)];
  }

  // The new set of external bonds is the symmetric set difference
  // of b1.ExternalBonds and b2.ExternalBonds
  std::set_symmetric_difference(b1.ExternalBonds.begin(), b1.ExternalBonds.end(), 
                                b2.ExternalBonds.begin(), b2.ExternalBonds.end(),
                                std::inserter(ExternalBonds, ExternalBonds.end()),
                                BondDescriptorCompare(&Sys->get_graph()));

  if (Debug)
  {
    std::cout << "External bonds of block are:\n";
    for (ExternalBondListType::const_iterator I = ExternalBonds.begin(); 
         I != ExternalBonds.end(); ++I)
    {
      std::cout << "source = " << Sys->source(*I) << ", target = " << Sys->target(*I) << '\n';
    }
  } // if (Debug)
}

void Block::print_debug(std::ostream& out) const
{
  out << "Block contains " << Psi.size() << " sites, wavefunction amplitudes are\n"
      << "site_descriptor         amplitude\n";
  for (WavefuncType::const_iterator I = Psi.begin(); I != Psi.end(); ++I)
  {
    out << std::setw(15) << I->first << std::setw(18) << I->second << '\n';
  }
}

//
// Superblock
//

Superblock::Superblock(const Block& b1, const Block& b2, 
                       const Block& b3, const Block& b4, bool Debug)
  : H(4,4)
{
  System const& S = b1.system();
  assert(&b2.system() == &S && &b3.system() == &S && &b4.system() == &S);
  if (Debug) 
  {
    std::cout << "Constructing superblock\nDetail block 1: ";
    b1.print_debug(std::cout);
    std::cout << "\nDetail block 2: ";
    b2.print_debug(std::cout);
    std::cout << "\nDetail block 3: ";
    b3.print_debug(std::cout);
    std::cout << "\nDetail block 4: ";
    b4.print_debug(std::cout);
  }

  H.clear();
  // effective Hamiltonian for 'a particle in a box': 
  H(0,0) = b1.local_ham();
  H(1,1) = b2.local_ham();
  H(2,2) = b3.local_ham();
  H(3,3) = b4.local_ham();

  // Assemble the cross terms as the union of all of the external bonds of the blocks.
  Block::ExternalBondListType CrossTerms(b1.begin_external_bonds(), b1.end_external_bonds(), &S.get_graph());
  CrossTerms.insert(b2.begin_external_bonds(), b2.end_external_bonds());
  CrossTerms.insert(b3.begin_external_bonds(), b3.end_external_bonds());
  CrossTerms.insert(b4.begin_external_bonds(), b4.end_external_bonds());

  // loop over the cross terms and set the matrix elements
  Block::ExternalBondListType::const_iterator b_end = CrossTerms.end();
  for (Block::ExternalBondListType::const_iterator b = CrossTerms.begin(); b!=b_end; ++b) 
  {
    site_descriptor source = S.source(*b), target = S.target(*b);

    int source_index, target_index;
    double source_wfn, target_wfn;

    boost::tie(source_index, source_wfn) = FindSiteInBlock(source,b1,b2,b3,b4);
    boost::tie(target_index, target_wfn) = FindSiteInBlock(target,b1,b2,b3,b4);

    if (source_index == -1 || target_index == -1) continue;


    if (Debug) std::cout << "adding interaction between sites " << source << " and " << target
                         << ", source block: " << source_index 
                         << ", target block: " << target_index << std::endl;

    H(source_index, target_index) += S.hopping(*b) * source_wfn * target_wfn;
    H(target_index, source_index) += S.hopping(*b) * source_wfn * target_wfn;
  }
  
  if (Debug) std::cout << "H = " << H << std::endl;
}

std::pair<double, Wavefunction> Superblock::GetGroundState()
{
  Vector evalues(4);
  // diagonalize superblock Hamiltonian; H contains eigenvectors on output:
  boost::numeric::bindings::lapack::syev('V',boost::numeric::bindings::upper(H),evalues,boost::numeric::bindings::lapack::optimal_workspace());
  // here: syev() calls LAPACK-routine 'dsyev' 

  // ground state vector: first column of the matrix H
  Wavefunction psi(boost::numeric::ublas::matrix_column<Matrix>(H,0));

  // ground state energy:
  double energy = evalues(0);
  return std::make_pair(energy, psi);
}
