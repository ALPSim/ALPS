/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef APP_ALPS_MODEL_H
#define APP_ALPS_MODEL_H

#include "dmrg/models/model.h"
#include "dmrg/models/alps/lattice.hpp"

#include <alps/parameter.h>
#include <alps/lattice.h>
#include <alps/model.h>

#undef tolower
#undef toupper
#include <boost/tokenizer.hpp>
#include <boost/bind.hpp>
#include <boost/container/flat_map.hpp>

#include "symm_handler.hpp"

namespace detail {
    inline alps::graph_helper<> const& get_graph(Lattice const& lat_)
    {
        alps_lattice const* alattice = static_cast<alps_lattice const*>(lat_.impl().get());
        return alattice->alps_graph();
    }
}

template <class I>
bool safe_is_fermionic(alps::SiteBasisDescriptor<I> const& b, alps::SiteOperator const& op)
{
    using boost::bind;
    std::set<std::string> operator_names = op.operator_names();
    return std::count_if(operator_names.begin(), operator_names.end(), boost::bind(&alps::SiteBasisDescriptor<I>::is_fermionic, b, _1)) % 2;
}


template <class Matrix, class SymmGroup>
class ALPSModel : public model_impl<Matrix, SymmGroup>
{
    typedef model_impl<Matrix, SymmGroup> base;
    
    typedef alps::SiteOperator SiteOperator;
    typedef alps::BondOperator BondOperator;
    
    typedef typename Matrix::value_type value_type;
    typedef typename maquis::traits::scalar_type<Matrix>::type scalar_type;
    typedef boost::multi_array<value_type,2> alps_matrix;
    typedef std::map<std::string, int> qn_map_type;
    
    typedef short I;
    typedef alps::graph_helper<> graph_type;
    
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    
    
public:
    typedef typename base::table_type table_type;
    typedef typename base::table_ptr table_ptr;
    typedef typename base::tag_type tag_type;
    
    typedef typename base::term_descriptor value_term;
    typedef typename alps::expression::Term<value_type> expression_type;
    typedef ::term_descriptor<expression_type> expression_term;
    typedef typename base::terms_type terms_type;
    typedef typename base::op_t op_t;
    typedef typename base::measurement_term_desc_type measurement_term_desc_type;
    typedef typename base::op_prod_type op_prod_type;
    typedef typename base::initializer_ptr initializer_ptr;
    
    typedef typename base::size_t size_t;
    
    typedef std::pair<std::string, int> opkey_type;
    typedef std::map<opkey_type, tag_type> opmap_type;
    typedef typename opmap_type::const_iterator opmap_const_iterator;
    
    typedef typename SymmGroup::charge charge;
    
    
    ALPSModel (Lattice const& lattice_, const alps::Parameters& parms_)
    : parms(parms_)
    , raw_lattice(lattice_)
    , lattice(detail::get_graph(lattice_))
    , model(lattice, parms, true)
    , tag_handler(new table_type())
    {
        locale_shared i;
        
        size_t num_vertex_types = alps::maximum_vertex_type(lattice.graph())+1;
        symm_basis.reserve(num_vertex_types);
        basis_descriptors.reserve(num_vertex_types);
        site_bases.reserve(num_vertex_types);

        /// Parsing conserved quantum numbers
        for (int type=0; type<=alps::maximum_vertex_type(lattice.graph()); ++type) {
            std::set<std::string> type_qn = model.quantum_numbers(type);
            all_qn.insert(type_qn.begin(), type_qn.end());
        }

        if (parms.defined("CONSERVED_QUANTUMNUMBERS")) {
            boost::char_separator<char> sep(" ,");
            std::string qn_string = parms["CONSERVED_QUANTUMNUMBERS"];
            tokenizer qn_tokens(qn_string, sep);
            int n=0;
            for (tokenizer::iterator it=qn_tokens.begin(); it != qn_tokens.end(); it++) {
                if (parms.defined(*it + "_total")) {
                    if (all_qn.find(*it) != all_qn.end())
                        all_conserved_qn.insert( std::make_pair(*it, n++) );
                    else
                        throw std::runtime_error("quantumnumber "+(*it)+" not defined in the model.");
                }
            }
        }
        
        /// Load all possible basis
        for (int type=0; type<=alps::maximum_vertex_type(lattice.graph()); ++type) {
            basis_descriptors.push_back(model.site_basis(type));
            site_bases.push_back(alps::site_basis<I>(basis_descriptors[type]));
            symm_basis.push_back(symmetric_basis_descriptor<SymmGroup>(basis_descriptors[type], all_conserved_qn));
            
            op_t ident, fill;
            for (int i=0; i<symm_basis[type].size(); ++i) {
                charge c = symm_basis[type].charge(i);
                size_t bsize = symm_basis[type].block_size(i);
                // maquis::cout << "Inserting " << c << " for " << site_bases[type][i] << std::endl;
                
                if (!ident.has_block(c, c))
                    ident.insert_block(Matrix::identity_matrix(bsize), c, c);
                
                int sign = (alps::is_fermionic(basis_descriptors[type], site_bases[type][i])) ? -1 : 1;
                if (!fill.has_block(c, c))
                    fill.insert_block(Matrix::identity_matrix(bsize), c, c);
                fill(symm_basis[type].coords(i), symm_basis[type].coords(i)) = sign;
            }
            operators[opkey_type("ident", type)] = tag_handler->register_op(ident, tag_detail::bosonic);
            operators[opkey_type("fill",  type)] = tag_handler->register_op(fill,  tag_detail::bosonic);
        }
        
        
        /// site_term loop with cache to avoid recomputing matrices
        std::vector<std::vector<std::pair<value_type, tag_type> > > site_terms(num_vertex_types);
        for (graph_type::site_iterator it=lattice.sites().first; it!=lattice.sites().second; ++it) {
            int p = lattice.vertex_index(*it);
            int type = lattice.site_type(*it);
            
            if (lattice.inhomogeneous_sites())
                alps::throw_if_xyz_defined(parms,*it); // check whether x, y, or z is set
            alps::expression::ParameterEvaluator<value_type> coords(coordinate_as_parameter(lattice.graph(), *it));
            
            if (site_terms[type].size() == 0) {
                typedef std::vector<boost::tuple<alps::expression::Term<value_type>,alps::SiteOperator> > V;
                V  ops = model.site_term(type).template templated_split<value_type>();
                                                        
                for (int n=0; n<ops.size(); ++n) {
                    SiteOperator op = boost::get<1>(ops[n]);
                    opmap_const_iterator match = operators.find(opkey_type(simplify_name(op), type));
                    if (match == operators.end())
                        match = register_operator(op, type, parms);
                    // site_terms[type].push_back( std::make_pair(boost::get<0>(ops[n]).value(), match->second)  );
                    
                    if (lattice.inhomogeneous_sites())
                        boost::get<0>(ops[n]).partial_evaluate(coords);
                    
                    expression_term term;
                    term.coeff = boost::get<0>(ops[n]);
                    term.is_fermionic = false;
                    term.push_back( boost::make_tuple(p, match->second) );
                    expression_coeff.insert( std::make_pair(term.coeff, value_type()) );
                    expression_terms.push_back(term);
                }
            }

            // All site terms summed into one
//            if (site_terms[type].size() > 0) {
//                opmap_const_iterator match = operators.find(opkey_type("site_terms", type));
//                if (match == operators.end()) {
//                    op_t op_matrix;
//                    for (int n=0; n<site_terms[type].size(); ++n)
//                        op_matrix += site_terms[type][n].first * tag_handler->get_op(site_terms[type][n].second);
//                    tag_type mytag = tag_handler->register_op(op_matrix, tag_detail::bosonic);
//                    boost::tie(match, boost::tuples::ignore) = operators.insert( std::make_pair(opkey_type("site_terms", type), mytag) );
//                }
//
//                term_descriptor term;
//                term.coeff = 1.;
//                term.is_fermionic = false;
//                term.push_back( boost::make_tuple(p, match->second) );
//                this->terms_.push_back(term);
//            }
            
            
        }
        
        /// bond terms loop
        for (graph_type::bond_iterator it=lattice.bonds().first; it!=lattice.bonds().second; ++it) {
            int p_s = lattice.source(*it);
            int p_t = lattice.target(*it);
            int type = lattice.bond_type(*it);
            int type_s = lattice.site_type(lattice.source(*it));
            int type_t = lattice.site_type(lattice.target(*it));
            
            bool wrap_pbc = boost::get(alps::boundary_crossing_t(), lattice.graph(), *it);
            
            BondOperator bondop = model.bond_term(type);
            
            typedef std::vector<boost::tuple<alps::expression::Term<value_type>,alps::SiteOperator,alps::SiteOperator > > V;
            alps::SiteBasisDescriptor<I> const& b1 = basis_descriptors[type_s];
            alps::SiteBasisDescriptor<I> const& b2 = basis_descriptors[type_t];
            
            if (lattice.inhomogeneous_bonds())
                alps::throw_if_xyz_defined(parms, lattice.graph()); // check whether x, y, or z is set
            alps::expression::ParameterEvaluator<value_type> coords(coordinate_as_parameter(lattice.graph(), *it));
            
            V  ops = bondop.template templated_split<value_type>(b1,b2);
            for (typename V::iterator tit=ops.begin(); tit!=ops.end();++tit) {
                SiteOperator op1 = boost::get<1>(*tit);
                SiteOperator op2 = boost::get<2>(*tit);
                
                opmap_const_iterator match1 = operators.find(opkey_type(simplify_name(op1), type_s));
                if (match1 == operators.end())
                    match1 = register_operator(op1, type_s, parms);
                opmap_const_iterator match2 = operators.find(opkey_type(simplify_name(op2), type_t));
                if (match2 == operators.end())
                    match2 = register_operator(op2, type_t, parms);
                
                bool with_sign = fermionic(b1, op1, b2, op2);
                
                if (lattice.inhomogeneous_bonds())
                    boost::get<0>(*tit).partial_evaluate(coords);
                
                expression_term term;
                term.coeff = boost::get<0>(*tit);
                term.is_fermionic = with_sign;
                
                {
                    tag_type mytag = match1->second;
                    if (with_sign && !wrap_pbc) {
                        // Note inverse notation because of notation in operator.
                        std::pair<tag_type, value_type> ptag = tag_handler->get_product_tag(operators[opkey_type("fill",type_s)],
                                                                                            mytag);
                        mytag = ptag.first;
                        term.coeff *= ptag.second;
                    }
                    if (with_sign && wrap_pbc)
                        term.coeff *= value_type(-1.);
                    term.push_back( boost::make_tuple(p_s, mytag) );
                }
                {
                    tag_type mytag = match2->second;
                    if (with_sign && wrap_pbc) {
                        // Note inverse notation because of notation in operator.
                        std::pair<tag_type, value_type> ptag = tag_handler->get_product_tag(operators[opkey_type("fill",type_t)],
                                                                                            mytag);
                        mytag = ptag.first;
                        term.coeff *= ptag.second;
                    }
                    term.push_back( boost::make_tuple(p_t, mytag) );
                }
                
                expression_coeff.insert( std::make_pair(term.coeff, value_type()) );
                expression_terms.push_back(term);
            }
        }
        
        generate_terms();
    }
    
    void update(BaseParameters const& p)
    {
        parms << p;
        generate_terms();
    }
    
    Index<SymmGroup> const& phys_dim(size_t type) const
    {
        return symm_basis[type].phys_dim();
    }
    
    typename SymmGroup::charge total_quantum_numbers(BaseParameters& parms_) const
    {
        return init_charge<SymmGroup>(parms_, all_conserved_qn);
    }
    
    tag_type identity_matrix_tag(size_t type) const
    {
        return operators[opkey_type("ident", type)]; // TODO: avoid using map here
    }

    tag_type filling_matrix_tag(size_t type) const
    {
        return operators[opkey_type("fill", type)]; // TODO: avoid using map here
    }
    
    bool has_operator(std::string const & name, size_t type) const
    {
        alps::SiteBasisDescriptor<I> const& b = basis_descriptors[type];
        return b.has_operator(name);
    }

    bool has_operator(op_prod_type const & name, size_t type) const
    {
        alps::SiteBasisDescriptor<I> const& b = basis_descriptors[type];
        bool ret = true;
        for (typename op_prod_type::const_iterator it = name.begin(); it != name.end(); ++it)
            ret = ret && b.has_operator(*it);
        return ret;
    }

    tag_type get_operator_tag(std::string const & name, size_t type) const
    {
        if (name == "id" || name == "ident" || name == "identity") {
            return operators[opkey_type("ident", type)];
        } else {
            opmap_const_iterator match = operators.find(opkey_type(name, type));
            if (match == operators.end()) {
                SiteOperator op = make_site_term(name, parms);
                match = register_operator(op, type, parms);
            }
            return match->second;
        }
    }

    tag_type get_operator_tag(op_prod_type const & name, size_t type) const
    {
        if (name.size() == 0) throw std::runtime_error("name.size() must be > 0");
        if (name.size() == 1)
            return get_operator_tag(name[0], type);
        
        /// Compute internal product name
        std::string prod_name = name[0];
        for (typename op_prod_type::const_iterator it = name.begin()+1; it != name.end(); ++it)
            prod_name += "__times__" + (*it);
        
        opmap_const_iterator match = operators.find(opkey_type(prod_name, type));
        if (match == operators.end()) {
            /// Compute product
            tag_type prod_tag = get_operator_tag(*(name.begin()), type);
            std::string prod_name = *(name.begin());
            for (typename op_prod_type::const_iterator it = name.begin()+1; it != name.end(); ++it)
            {
                /// Look if partial product is cached
                prod_name += "__times__" + (*it);
                match = operators.find(opkey_type(prod_name, type));
                if (match != operators.end()) {
                    prod_tag = match->second;
                } else {
                    /// Calculate partial product
                    typename Matrix::value_type coeff;
                    // Note: inverse multiplication to keep consistency with Op notation
                    boost::tie(prod_tag, coeff) = tag_handler->get_product_tag(get_operator_tag(*it, type), prod_tag);
                    if (coeff != 1.) { // TODO: find good test candidate
                        /// force register op with trivial coeff: get, multiply and register
                        tag_detail::operator_kind kind = tag_handler->is_fermionic(prod_tag) ? tag_detail::fermionic : tag_detail::bosonic;
                        op_t m = tag_handler->get_op(prod_tag);
                        m *= coeff;
                        prod_tag = tag_handler->register_op(m, kind);
                    }
                    boost::tie(match, boost::tuples::ignore) = operators.insert( std::make_pair(opkey_type(prod_name, type), prod_tag) );
                }
            }
        }

        return match->second;
    }

    table_ptr operators_table() const
    {
        return tag_handler;
    }

    initializer_ptr initializer(Lattice const& lat, BaseParameters & p_) const;
    std::vector<measurement_term_desc_type> unpack_measurement_terms(std::string const & name) const;

private:
    
    template <class SiteOp>
    std::string simplify_name(const SiteOp &op) const
    {
        std::string term = op.term();
        std::string arg = "("+op.site()+")";
        boost::algorithm::replace_all(term,arg,"");
        return term;
    }
    
    bool fermionic (alps::SiteBasisDescriptor<I> const& b1, SiteOperator const& op1,
                    alps::SiteBasisDescriptor<I> const& b2, SiteOperator const& op2) const
    {  
        bool is_ferm1 = safe_is_fermionic(b1, op1);
        bool is_ferm2 = safe_is_fermionic(b2, op2);
        return is_ferm1 || is_ferm2;
    }
    
    inline op_t convert_matrix (const alps_matrix& m, int type) const
    {
        op_t newm;
        for (int i=0; i<m.shape()[0]; ++i) {
            for (int j=0; j<m.shape()[1]; ++j) {
                if (m[i][j] != 0.) {
                    charge c_i = symm_basis[type].charge(i);
                    size_t bsize_i = symm_basis[type].block_size(i);
                    charge c_j = symm_basis[type].charge(j);
                    size_t bsize_j = symm_basis[type].block_size(j);
                    
                    if (!newm.has_block(c_i, c_j))
                        newm.insert_block(Matrix(bsize_i, bsize_j, 0), c_i, c_j);
                    // Notation: going from state i to state j
                    newm(symm_basis[type].coords(i), symm_basis[type].coords(j)) = m[i][j];
                }
            }
        }
        return newm;
    }
    
    alps::SiteOperator make_site_term(std::string x, alps::Parameters const & parms) const
    {
        if (x[x.size()-1]!=')')
            x += "(i)";
        alps::SiteOperator op(x,"i");
        model.substitute_operators(op, parms);
        return op;
    }
    
    opmap_const_iterator register_operator(SiteOperator const& op, int type, alps::Parameters const& p) const
    {
        alps::SiteBasisDescriptor<I> const& b = basis_descriptors[type];
        alps_matrix m = alps::get_matrix(value_type(), op, b, p, true);
        tag_detail::operator_kind kind = safe_is_fermionic(b, op) ? tag_detail::fermionic : tag_detail::bosonic;
        tag_type mytag = tag_handler->register_op(convert_matrix(m, type), kind);
        
        opmap_const_iterator match;
        boost::tie(match, boost::tuples::ignore) = operators.insert( std::make_pair(opkey_type(simplify_name(op), type), mytag) );
        return match;
    }
    
    std::vector<std::string> split_term_factors(SiteOperator const& op) const
    {
        std::vector<std::string> term_factors;
        alps::expression::Expression<value_type> op_expression(op.term());
        if (!op_expression.is_single_term())
            throw std::runtime_error("Something went wront with splitting factors in " + op.term());
        alps::expression::Term<value_type> const& term = *(op_expression.terms().first);
        for (typename alps::expression::Term<value_type>::factor_iterator fit = term.factors().first;
             fit != term.factors().second; ++fit)
        {
            // Clenup factor name, i.e. remove "(x)" or "(i)"
            std::string factor_name = boost::lexical_cast<std::string>(*fit);
            std::string arg = "("+op.site()+")";
            boost::algorithm::replace_all(factor_name,arg,"");
            term_factors.push_back(factor_name);
        }
        return term_factors;
    }
    
    void generate_terms()
    {
        this->terms_.clear();
        this->terms_.reserve(expression_terms.size());
        
        alps::Parameters parms_with_defaults(parms);
        parms_with_defaults.copy_undefined(model.model().default_parameters());
        
        // typedef typename boost::container::flat_map<expression_type, value_type>::iterator coeff_iterator;
        // for(coeff_iterator it = expression_coeff.begin(); it != expression_coeff.end(); ++it)
        //     it->second = alps::partial_evaluate<value_type>(it->first, parms_with_defaults);
        
        typedef typename std::vector<expression_term>::const_iterator terms_iterator;
        for(terms_iterator it = expression_terms.begin(); it != expression_terms.end(); ++it) {
            
            value_type val = 0.;
            if (lattice.inhomogeneous_sites() && it->size() == 1) {
                alps::Parameters p(parms_with_defaults);
                alps::throw_if_xyz_defined(p, lattice.graph()); // check whether x, y, or z is set
                p << coordinate_as_parameter(lattice.graph(), lattice.site(it->position(0)));
                
                val = alps::evaluate<value_type>(it->coeff, p);
            } else if(lattice.inhomogeneous_bonds() && it->size() == 2) {
                alps::Parameters p(parms_with_defaults);
                alps::throw_if_xyz_defined(p, lattice.graph()); // check whether x, y, or z is set
                p << coordinate_as_parameter(lattice.graph(), lattice.site(it->position(0)), lattice.site(it->position(1)));
                
                val = alps::evaluate<value_type>(it->coeff, p);
            } else {
                val = alps::evaluate<value_type>(it->coeff, parms_with_defaults);
            }
            
            // value_type const& val = expression_coeff[it->coeff];
            if ( alps::numeric::is_nonzero(val) ) {
                value_term term;
                term.is_fermionic = it->is_fermionic;
                term.insert(term.end(), it->begin(), it->end());
                term.coeff = val;
                this->terms_.push_back(term);
            }
        }
    }
    
    
    alps::Parameters parms;
    Lattice raw_lattice;
    graph_type const& lattice;
    alps::model_helper<I> model;
    mutable table_ptr tag_handler;

    std::set<std::string> all_qn;
    qn_map_type all_conserved_qn;
    std::vector<symmetric_basis_descriptor<SymmGroup> > symm_basis;
    std::vector<alps::SiteBasisDescriptor<I> > basis_descriptors;
    std::vector<alps::site_basis<I> > site_bases;
    
    mutable opmap_type operators; // key=<name,type>
    std::vector<expression_term> expression_terms;
    boost::container::flat_map<expression_type, value_type> expression_coeff;
};

// Initial states
template <class Matrix, class SymmGroup>
typename ALPSModel<Matrix, SymmGroup>::initializer_ptr ALPSModel<Matrix, SymmGroup>::initializer(Lattice const& lat, BaseParameters & p_) const
{
    if ( p_["init_state"] == "local_quantumnumbers" ) {
        int max_site_type = 0;
        std::vector<int> site_types(lat.size(), 0);
        for (int p = 0; p < lat.size(); ++p) {
            site_types[p] = lat.get_prop<int>("type", p);
            max_site_type = std::max(site_types[p], max_site_type);
        }
        
        std::cout << "site_types: ";
        std::copy(site_types.begin(), site_types.end(), std::ostream_iterator<int>(std::cout, " "));
        std::cout << std::endl;
        
        std::vector<Index<SymmGroup> > phys_bases(symm_basis.size());
        for (int type = 0; type < phys_bases.size(); ++type) {
            phys_bases[type] = symm_basis[type].phys_dim();
            maquis::cout << "phys["<< type <<"]: " << phys_bases[type] << std::endl;
        }
        
        // TODO: avoid QN of size=1
        std::map<std::string, std::vector<double> > initial_local_charges;
        for(std::set<std::string>::const_iterator it = all_qn.begin(); it != all_qn.end(); ++it) {
            const std::string pname = "initial_local_" + *it;
            if (!p_.defined(pname))
                throw std::runtime_error(pname + " required for local_quantumnumbers initial state.");
            initial_local_charges[*it] = p_[pname].as<std::vector<double> >();
            if (initial_local_charges[*it].size() != lat.size())
                throw std::runtime_error(pname + " does not match the lattice size.");
        }
        
        std::vector<boost::tuple<charge, size_t> > state(lat.size());
        for (size_t p=0; p<lat.size(); ++p) {
            const int type = site_types[p];
            alps::SiteBasisDescriptor<I> const& b = basis_descriptors[type];
            alps::site_state<I> local_state;
            for (size_t i=0; i<b.size(); ++i)
                local_state.push_back( initial_local_charges[b[i].name()][p] );
            state[p] = symm_basis[type].coords(site_bases[type].index(local_state));
        }
        
        return initializer_ptr(new basis_mps_init_generic<Matrix, SymmGroup>(state, phys_bases, this->total_quantum_numbers(p_), site_types));
    } else {
        return base::initializer(lat, p_);
    }
}


// Unpack measurements
template <class Matrix, class SymmGroup>
std::vector<typename ALPSModel<Matrix, SymmGroup>::measurement_term_desc_type>
ALPSModel<Matrix, SymmGroup>::unpack_measurement_terms(std::string const & name) const
{
    std::vector<measurement_term_desc_type> measurement_terms;
    
    std::string packed_name = boost::trim_copy(name);
    int ntypes = alps::maximum_vertex_type(lattice.graph())+1;

    if (model.has_bond_operator(packed_name)) {
        BondOperator bondop = model.get_bond_operator(packed_name);
        typedef std::vector<boost::tuple<alps::expression::Term<value_type>,alps::SiteOperator,alps::SiteOperator > > V;
        /// Assumption: operator names are the same for all type1-type2 combinations
        
        alps::SiteBasisDescriptor<I> const& b1 = basis_descriptors[0];
        alps::SiteBasisDescriptor<I> const& b2 = basis_descriptors[0];
        
        V  bond_terms = bondop.template templated_split<value_type>(b1,b2);
        for (typename V::iterator tit=bond_terms.begin(); tit!=bond_terms.end();++tit) {
            SiteOperator op1 = boost::get<1>(*tit);
            SiteOperator op2 = boost::get<2>(*tit);
            
            model.substitute_operators(op1, parms);
            model.substitute_operators(op2, parms);

            measurement_term_desc_type desc;
            desc.coeff = boost::get<0>(*tit).value();
            desc.op_names.push_back(split_term_factors(op1));
            desc.op_names.push_back(split_term_factors(op2));
            measurement_terms.push_back(desc);
        }
    } else {
        bool site_term_found = false;
        for (int type=0; type<ntypes; ++type) {
            alps::SiteBasisDescriptor<I> const& b = basis_descriptors[type];
            if (b.has_operator(packed_name)) {
                measurement_term_desc_type desc;
                desc.coeff = 1.;
                desc.op_names.push_back(std::vector<std::string>(1,packed_name));
                measurement_terms.push_back(desc);
                site_term_found = true;
                break;
            }
        }
        if (!site_term_found && model.has_site_operator(packed_name)) {
            SiteOperator op = model.get_site_operator(packed_name);
            model.substitute_operators(op, parms);
            
            // Split operator into terms
            typedef std::vector<boost::tuple<alps::expression::Term<value_type>,SiteOperator> > V;
            V site_terms = op.template templated_split<value_type>();
            for (typename V::const_iterator tit=site_terms.begin(); tit!=site_terms.end();++tit) {
                measurement_term_desc_type desc;
                desc.coeff = boost::get<0>(*tit).value();
                desc.op_names.push_back(split_term_factors(boost::get<1>(*tit)));
                measurement_terms.push_back(desc);
            }
        }
    }
    if (measurement_terms.size() == 0)
        throw std::runtime_error("Operator "+packed_name+" not found.");
    return measurement_terms;
}


#endif
