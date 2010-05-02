/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1994-2006 by Matthias Troyer <troyer@comp-phys.org>
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

#include <alps/model.h>
#include <alps/lattice.h>
#include <alps/scheduler.h>
#include <alps/scheduler/measurement_operators.h>

#ifdef ALPS_HAVE_HDF5
#include <alps/hdf5.hpp>
#endif

#define DMRG_VERSION "1.0.0"
#define DMRG_DATE "2006/10/02"

class DMRGTask 
 : public alps::scheduler::Task
 , public alps::graph_helper<>
 , public alps::model_helper<>
 , protected alps::EigenvectorMeasurements<double>
{
public:  
  typedef alps::half_integer<short> half_integer_type;
  DMRGTask(const alps::ProcessList& , const boost::filesystem::path& );
  DMRGTask(const alps::ProcessList& w, const alps::Parameters& p);
  void dostep();
  void write_xml_body(alps::oxstream&, const boost::filesystem::path&, bool) const;

  static void print_copyright(std::ostream& os = std::cout) 
  {
    os << "ALPS/dmrg version " DMRG_VERSION " (" DMRG_DATE ")\n"
       << "  Density Matrix Renormalization Group algorithm\n"
       << "  for low-dimensional interacting systems.\n"
       << "  available from http://alps.comp-phys.org/\n"
       << "  copyright (c) 2006-2013 by Adrian E. Feiguin\n"
       << "  for details see the publication: \n"
       << "  A.F. Albuquerque et al., J. of Magn. and Magn. Materials 310, 1187 (2007).\n\n";
  }

#ifdef ALPS_HAVE_HDF5
  void serialize(alps::hdf5::oarchive &) const;
#endif

private:

  alps::SiteOperator make_site_term(std::string x)
  {
    if (x[x.size()-1]!=')')
      x += "(i)";
    alps::SiteOperator op(x,"i");
    substitute_operators(op,parms);
    return op;
  }
  
  
  void init();
  dmtk::BasicOp<double> create_site_operator(std::string const& name, alps::SiteOperator const& siteop, int type);
  void build_site_operator(alps::SiteOperator const& siteop, int site, dmtk::Hami<double> &);
  void build_bond_operator(alps::BondOperator const& bondop, bond_descriptor const& b, dmtk::Hami<double> &this_hami);  
  void build_2site_operator(std::pair<alps::SiteOperator,alps::SiteOperator> const& siteops, 
                            std::pair<int,int> sites, dmtk::Hami<double> &this_hami);
 
  void save_results(); 
  
  int num_sweeps;
  std::vector<int> num_states;
  std::vector<std::string> quantumnumber_names;
  std::vector<bool> conserved_quantumnumber;
  std::vector<half_integer_type> conserved_quantumnumber_value;
  int qnmask;

  dmtk::System<double> system;
  dmtk::Hami<double> hami;
  dmtk::Lattice lattice;
  std::vector<dmtk::Block<double> > site_block;
  
  int num_eigenvalues;
  int verbose;
  int nwarmup;

  // For resuming a previous run
  int start_sweep;
  int start_dir;
  int start_iter;
  
// measurements during iterations
  std::vector<alps::EigenvectorMeasurements<value_type> > iteration_measurements;
};

