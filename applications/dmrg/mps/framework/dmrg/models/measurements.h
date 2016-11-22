/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef MEASUREMENTS_H
#define MEASUREMENTS_H

#include "dmrg/models/measurements/average.h"
#include "dmrg/models/measurements/local.h"
#include "dmrg/models/measurements/local_at.h"
#include "dmrg/models/measurements/correlations.h"
#include "dmrg/models/measurements/custom.h"
#include "dmrg/models/measurements/overlap.h"
#include "dmrg/models/measurements/entanglement.h"

#include "dmrg/models/model.h"
#include "dmrg/models/lattice.h"
#include "dmrg/models/measurement_utils.hpp"

#undef tolower
#undef toupper
#include <boost/tokenizer.hpp>
#include <boost/regex.hpp>


namespace measurements {
    namespace detail {
        /// Build matrices needed for some measurements
        template <class Matrix, class SymmGroup>
        std::pair<std::vector<std::pair<std::vector<typename Model<Matrix, SymmGroup>::op_t>, bool> >, short>
        operators_for_meas(std::string const& ops,
                           Lattice const& lattice,
                           Model<Matrix, SymmGroup> const& model,
                           bool repeat_one=false)
        {
            typedef boost::ptr_vector<measurement<Matrix, SymmGroup> > measurements_container;
            typedef typename Model<Matrix, SymmGroup>::op_t op_t;
            typedef std::vector<op_t> op_vec;
            typedef std::vector<std::pair<op_vec, bool> > meas_operators_type;
            typedef typename Model<Matrix, SymmGroup>::table_ptr table_ptr;
            typedef typename Model<Matrix, SymmGroup>::tag_type tag_type;
            typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
            typedef typename Model<Matrix, SymmGroup>::measurement_term_desc_type measurement_term_desc_type;
            
            meas_operators_type ret;
            
            table_ptr const & tag_handler = model.operators_table();
            int ntypes = lattice.maximum_vertex_type()+1;
            int f_ops = 0;
            short op_types = 0;
            
            boost::char_separator<char> sep(":");
            tokenizer corr_tokens(ops, sep);
            for (tokenizer::iterator it=corr_tokens.begin(); it != corr_tokens.end(); it++)
            {
                std::string name = boost::trim_copy(*it);
                
                std::vector<measurement_term_desc_type> meas_terms_vec = model.unpack_measurement_terms(name);
                if (meas_terms_vec.size() != 1)
                    throw std::runtime_error("Cannot measure sum of terms in LOCAL_AT and CORRELATIONS");
                
                measurement_term_desc_type const& meas_term_desc = meas_terms_vec[0];
                
                for (int i=0; i<meas_term_desc.op_names.size(); ++i) {
                    enum {uknown, bosonic, fermionic} kind = uknown;
                    op_vec tops(ntypes);
                    for (int type=0; type<ntypes; ++type) {
                        if (model.has_operator(meas_term_desc.op_names[i], type)) {
                            tag_type tag = model.get_operator_tag(meas_term_desc.op_names[i], type);
                            
                            tops[type] = tag_handler->get_op(tag);
                            if (i == 0)
                                tops[type] *= meas_term_desc.coeff;
                            
                            bool is_ferm = tag_handler->is_fermionic(tag);
                            if (kind == uknown)
                                kind = is_ferm ? fermionic : bosonic;
                            else if ((is_ferm && kind==bosonic) || (!is_ferm && kind==fermionic))
                                throw std::runtime_error("Model is inconsitent. On some site the operator " + name + "fermionic, on others is bosonic.");
                        }
                    }
                    ret.push_back(std::make_pair( tops , (kind == fermionic) ));
                    if (kind == fermionic) ++f_ops;
                }
                if (meas_term_desc.op_names.size() == 2)
                    op_types |= 1<<1; // flag that at least one bond term is in the list
                else
                    op_types |= 1<<0; // flag that at least one site term is in the list
            }
            
            if (op_types == 3) throw std::runtime_error("Can only have either site operators or bond operators.");
            
            /// repeat last site term in case only one in the input
            if (repeat_one && op_types == 1 && ret.size() == 1) {
                ret.push_back(ret[0]);
                if (ret[1].second) ++f_ops;
            }
            /// repeat last bond term (=two site terms) in case only one in the input
            if (repeat_one && op_types == 2 && ret.size() == 2) {
                ret.push_back(ret[0]);
                ret.push_back(ret[1]);
                if (ret[2].second) ++f_ops;
                if (ret[3].second) ++f_ops;
            }
            
            if (f_ops % 2 != 0)
                throw std::runtime_error("Number of fermionic operators has to be even.");
            
            return std::make_pair(ret, op_types);
        }
        
    }
    
    
    /// Parse and create vector of all measurements from parameters
    template <class Matrix, class SymmGroup>
    boost::ptr_vector<measurement<Matrix, SymmGroup> >
    parse_and_create(Lattice const& lattice,
                     Model<Matrix, SymmGroup> const& model,
                     BaseParameters const& parms)
    {
        typedef boost::ptr_vector<measurement<Matrix, SymmGroup> > measurements_container;
        typedef typename Model<Matrix, SymmGroup>::op_t op_t;
        typedef std::vector<op_t> op_vec;
        typedef std::vector<std::pair<op_vec, bool> > meas_operators_type;
        typedef typename Model<Matrix, SymmGroup>::table_ptr table_ptr;
        typedef typename Model<Matrix, SymmGroup>::tag_type tag_type;
        typedef typename Model<Matrix, SymmGroup>::measurement_term_desc_type measurement_term_desc_type;
        typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

        measurements_container meas;
        
        int ntypes = lattice.maximum_vertex_type()+1;
        table_ptr const & tag_handler = model.operators_table();
        
        std::vector<op_t> identitities(ntypes), fillings(ntypes);
        for (int type=0; type<ntypes; ++type) {
            identitities[type] = model.identity_matrix(type);
            fillings[type]     = model.filling_matrix(type);
        }
        
        std::vector<MeasurementSpecs> meas_specs = parse_measurements(parms);
        for (typename std::vector<MeasurementSpecs>::const_iterator itm = meas_specs.begin();
             itm != meas_specs.end(); ++itm) {
            
            if (itm->meas_class == "LOCAL" || itm->meas_class == "AVERAGE") {
                
                std::vector<measurement_term_desc_type> meas_term_desc = model.unpack_measurement_terms(itm->args);
                assert(meas_term_desc.size() > 0); // this is a post-condition of unpack_measurement_terms
                
                typedef std::vector<std::pair<op_vec, bool> > bond_element;
                std::vector<bond_element> operators;
                
                if (meas_term_desc[0].op_names.size() == 1) {
                    /// Site measurement
                    
                    operators.resize(meas_term_desc.size(), bond_element(1, std::make_pair(op_vec(ntypes), false)));
                    for (int i=0; i<meas_term_desc.size(); ++i) {
                        for (int type=0; type<ntypes; ++type) {
                            if (model.has_operator(meas_term_desc[i].op_names[0], type))
                                operators[i][0].first[type] = model.get_operator(meas_term_desc[i].op_names[0], type);
                        }
                    }
                    
                } else {
                    /// Bond measurement
                    
                    operators.resize(meas_term_desc.size(), bond_element(2, std::make_pair(op_vec(ntypes), false)));
                    for (int i=0; i<meas_term_desc.size(); ++i) {
                        for (int type=0; type<ntypes; ++type) {
                            // Op1
                            if (model.has_operator(meas_term_desc[i].op_names[0], type)) {
                                tag_type tag = model.get_operator_tag(meas_term_desc[i].op_names[0], type);
                                operators[i][0].second = tag_handler->is_fermionic(tag);
                                operators[i][0].first[type] = meas_term_desc[i].coeff * tag_handler->get_op(tag);
                            }
                            // Op2
                            if (model.has_operator(meas_term_desc[i].op_names[0], type)) {
                                tag_type tag = model.get_operator_tag(meas_term_desc[i].op_names[1], type);
                                operators[i][1].second = tag_handler->is_fermionic(tag);
                                operators[i][1].first[type] = tag_handler->get_op(tag);
                            }
                        }
                    }
                }
                
                if (itm->meas_class == "LOCAL")
                    meas.push_back( new measurements::local<Matrix, SymmGroup>(itm->name, lattice, identitities, fillings, operators) );
                else
                    meas.push_back( new measurements::average<Matrix, SymmGroup>(itm->name, lattice, identitities, fillings, operators) );
                
                
            } else if (itm->meas_class == "LOCAL_AT") {
                // Example: MEASURE_LOCAL_AT[Custom correlation] = "bdag:b|(1,2),(3,4),(5,6)"
                
                boost::char_separator<char> part_sep("|");
                tokenizer part_tokens(itm->args, part_sep);
                std::vector<std::string> parts;
                std::copy( part_tokens.begin(), part_tokens.end(), std::back_inserter(parts) );
                
                if (parts.size() != 2)
                    throw std::runtime_error("MEASURE_LOCAL_AT must contain a `|` delimiter.");
                
                /// parse operators
                meas_operators_type operators;
                boost::tie(operators, boost::tuples::ignore) = detail::operators_for_meas(parts[0], lattice, model, false);
                
                /// parse positions
                std::vector<std::vector<std::size_t> > positions;
                boost::regex pos_re("\\(([^(^)]*)\\)");
                boost::sregex_token_iterator it_pos(parts[1].begin(), parts[1].end(), pos_re, 1);
                boost::sregex_token_iterator it_pos_end;
                for (; it_pos != it_pos_end; ++it_pos)
                {
                    boost::char_separator<char> int_sep(", ");
                    std::string raw = *it_pos;
                    tokenizer int_tokens(raw, int_sep);
                    
                    std::vector<std::size_t> pos;
                    BOOST_FOREACH(std::string t, int_tokens) {
                        pos.push_back(boost::lexical_cast<std::size_t, std::string>(t));
                    }
                    positions.push_back(pos);
                }
                
                meas.push_back( new measurements::local_at<Matrix, SymmGroup>(itm->name, lattice, positions, identitities, fillings, operators) );
                
                
            } else if (itm->meas_class.find("CORRELATIONS") != std::string::npos) {
                bool half_only=true, nearest_neighbors_only=false;
                if (itm->meas_class == "CORRELATIONS") {
                    half_only = false;
                    nearest_neighbors_only = false;
                }
                if (itm->meas_class == "HALF_CORRELATIONS") {
                    half_only = true;
                    nearest_neighbors_only = false;
                }
                if (itm->meas_class == "NN_CORRELATIONS") {
                    half_only = false;
                    nearest_neighbors_only = true;
                }
                if (itm->meas_class == "HALF_NN_CORRELATIONS") {
                    half_only = true;
                    nearest_neighbors_only = true;
                }
                
                /// split op1:op2:...@p1,p2,p3,... into {op1:op2:...}, {p1,p2,p3,...}
                std::vector<std::string> value_split;
                boost::split( value_split, itm->args, boost::is_any_of("@"));
                
                meas_operators_type operators;
                short ops_type;
                boost::tie(operators, ops_type) = detail::operators_for_meas(value_split[0], lattice, model, true);
                
                if (ops_type == 2) nearest_neighbors_only = true;
                
                /// parse positions p1,p2,p3,... (or `space`)
                std::vector<std::size_t> positions;
                if (value_split.size() > 1) {
                    boost::char_separator<char> pos_sep(", ");
                    tokenizer pos_tokens(value_split[1], pos_sep);
                    BOOST_FOREACH(std::string t, pos_tokens) {
                        positions.push_back(boost::lexical_cast<std::size_t, std::string>(t));
                    }
                }
                
                meas.push_back( new measurements::correlations<Matrix, SymmGroup>(itm->name, lattice, identitities, fillings, operators,
                                                                                  half_only, nearest_neighbors_only,
                                                                                  positions) );
            }
        }
        
        return meas;
    }
    

    
}




///////////////////////////////////////////////////
/// UTILITIES
///////////////////////////////////////////////////

template <class Matrix, class SymmGroup>
class measure_and_save {
public:
    measure_and_save(std::string const& rfile_, std::string const& archive_path_,
                     MPS<Matrix, SymmGroup> const& mps_, int eigenstate_=0)
    : rfile(rfile_)
    , archive_path(archive_path_)
    , eigenstate(eigenstate_)
    , mps(mps_)
    , rmps(mps)
    { }
    
    void operator()(measurement<Matrix, SymmGroup> & meas) const
    {
        maquis::cout << "Measuring " << meas.name() << std::endl;
        meas.eigenstate_index() = eigenstate;
        meas.evaluate(mps, rmps);
        storage::archive ar(rfile, "w");
        ar[archive_path] << meas;
    }
    
private:
    std::string rfile, archive_path;
    int eigenstate;
    MPS<Matrix, SymmGroup> const& mps;
    reduced_mps<Matrix, SymmGroup> rmps;
};


namespace detail {
    class name_not_in_list {
    public:
        name_not_in_list(std::vector<std::string> const& list_)
        : list(list_)
        { }
        
        template <class Matrix, class SymmGroup>
        bool operator() (measurement<Matrix, SymmGroup> const& term) const
        {
            return std::find(list.begin(), list.end(), term.name()) == list.end();
        }
        
    private:
        std::vector<std::string> const& list;
    };
}

template <class Matrix, class SymmGroup>
boost::ptr_vector<measurement<Matrix, SymmGroup> > &
operator<< (boost::ptr_vector<measurement<Matrix, SymmGroup> > & lhs,
            boost::ptr_vector<measurement<Matrix, SymmGroup> > const& rhs)
{
    lhs.insert(lhs.end(), rhs.begin(), rhs.end());
    return lhs;
}

template <class Matrix, class SymmGroup>
boost::ptr_vector<measurement<Matrix, SymmGroup> >
meas_sublist(boost::ptr_vector<measurement<Matrix, SymmGroup> > const& m,
             std::vector<std::string> const& meas_list)
{
    boost::ptr_vector<measurement<Matrix, SymmGroup> > sublist(m.clone());
    sublist.erase_if( ::detail::name_not_in_list(meas_list) );
    return sublist;
}

//template <class Matrix, class SymmGroup>
//class DMOverlapMeasurement : public Measurement_Term<Matrix, SymmGroup> {
//public:
//    typedef Measurement_Term<Matrix, SymmGroup> base;
//    
//    MPS<Matrix, SymmGroup> mps_ident;
//    std::vector<MPS<Matrix, SymmGroup> > overlaps_mps;
//    std::vector<std::string> labels;
//    
//protected:
//    virtual Measurement_Term<Matrix, SymmGroup> * do_clone() const
//    {
//        return new DMOverlapMeasurement<Matrix, SymmGroup>(*this);
//    }
//};


template <class Matrix, class SymmGroup>
boost::ptr_vector<measurement<Matrix, SymmGroup> >
overlap_measurements(BaseParameters const & parms, boost::optional<size_t> sweep = boost::none)
{
    /* Syntax for MEASURE_OVERLAP:
     *  (1) MEASURE_OVERLAP[obsname] = "/path/to/ckp.h5"
     *  (2) MEASURE_OVERLAP[obsname(sweep)] = "/path/to/ckp.h5"
     *
     * `obsname` is the name that will be given in archive output.
     * if `sweep` is prensent, the overlap will only be computed when the sweep number
     * matches the given one.
     */
    boost::ptr_vector<measurement<Matrix, SymmGroup> > meas;
    boost::regex expression("^MEASURE_OVERLAP\\[([a-zA-Z]+)(\\(([0-9]+)\\))?\\]$");
    boost::smatch what;
    for (BaseParameters::const_iterator it=parms.begin();it != parms.end();++it) {
        std::string lhs = it->key();
        if (boost::regex_match(lhs, what, expression)) {
            if (sweep && !what[3].matched) continue;
            if (sweep && what[3].matched && boost::lexical_cast<long>(what.str(3)) != sweep.get()) continue;
            
            std::string name = what.str(1), bra_chkp = it->value();
            meas.push_back( new measurements::overlap<Matrix, SymmGroup>(name, bra_chkp) );
        }
    }
    return meas;
}

#endif
