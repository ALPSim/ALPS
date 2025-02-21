/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 * 
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#ifndef MAQUIS_DMRG_MODELS_MODELS_H
#define MAQUIS_DMRG_MODELS_MODELS_H

#include "dmrg/utils/BaseParameters.h"
#include "dmrg/models/term_descriptor.h"
#include "dmrg/models/measurement_term_desc.hpp"

#include "dmrg/models/lattice.h"
#include "dmrg/models/op_handler.h"
#include "dmrg/mp_tensors/mps_initializers.h"
#include "dmrg/block_matrix/block_matrix.h"


#include <boost/shared_ptr.hpp>


/// base type for all models
template <class Matrix, class SymmGroup>
class model_impl {
public:
    typedef boost::shared_ptr<mps_initializer<Matrix, SymmGroup> > initializer_ptr;

    typedef TagHandler<Matrix, SymmGroup> table_type;
    typedef boost::shared_ptr<table_type> table_ptr;
    typedef typename table_type::tag_type tag_type;

    typedef ::term_descriptor<typename Matrix::value_type> term_descriptor;
    typedef typename std::vector<term_descriptor> terms_type;
    typedef block_matrix<Matrix, SymmGroup> op_t;
    typedef MeasurementTermDescriptor<typename Matrix::value_type> measurement_term_desc_type;
    typedef typename measurement_term_desc_type::op_prod_type op_prod_type;
    
    typedef std::size_t size_t;
    
    virtual ~model_impl() {}
    
    virtual void update(BaseParameters const& p) =0;

    virtual Index<SymmGroup> const& phys_dim(size_t type) const=0;
    virtual op_t const& identity_matrix(size_t type) const { return operators_table()->get_op( identity_matrix_tag(type) ); }
    virtual tag_type identity_matrix_tag(size_t type) const=0;
    virtual op_t const& filling_matrix(size_t type) const { return operators_table()->get_op( filling_matrix_tag(type) ); }
    virtual tag_type filling_matrix_tag(size_t type) const=0;

    virtual typename SymmGroup::charge total_quantum_numbers(BaseParameters & parms) const=0;

    virtual terms_type const& hamiltonian_terms() const { return terms_; }
    
    virtual bool has_operator(std::string const & name, size_t type) const=0;
    virtual op_t const& get_operator(std::string const & name, size_t type) const { return operators_table()->get_op( get_operator_tag(name, type) ); }
    virtual tag_type get_operator_tag(std::string const & name, size_t type) const=0;

    virtual bool has_operator(op_prod_type const & name, size_t type) const
    {
        if (name.size() != 1)
            throw std::runtime_error("Default fallback on has_operator(name, type) is only valid for trivial op_prod names.");
        return has_operator(name[0], type);
    }
    virtual op_t const& get_operator(op_prod_type const & name, size_t type) const { return operators_table()->get_op( get_operator_tag(name, type) ); }
    virtual tag_type get_operator_tag(op_prod_type const & name, size_t type) const
    {
        if (name.size() != 1)
            throw std::runtime_error("Default fallback on get_operator_tag(name, type) is only valid for trivial op_prod names.");
        return get_operator_tag(name[0], type);
    }

    virtual std::vector<measurement_term_desc_type> unpack_measurement_terms(std::string const & name) const=0;

    virtual table_ptr operators_table() const=0;
    
    virtual initializer_ptr initializer(Lattice const& lat, BaseParameters & parms) const;
    
protected:
    terms_type terms_;
};

/// model factory
template <class Matrix, class SymmGroup>
boost::shared_ptr<model_impl<Matrix, SymmGroup> >
model_factory(Lattice const& lattice, BaseParameters & parms);


/// pimpl for Model
template <class Matrix, class SymmGroup>
class Model {
    typedef model_impl<Matrix, SymmGroup> impl_type;
    typedef boost::shared_ptr<impl_type> impl_ptr;
public:
    typedef typename impl_type::initializer_ptr initializer_ptr;
    
    typedef typename impl_type::table_type table_type;
    typedef typename impl_type::table_ptr table_ptr;
    typedef typename impl_type::tag_type tag_type;
    
    typedef typename impl_type::term_descriptor term_descriptor;
    typedef typename impl_type::terms_type terms_type;
    typedef typename impl_type::op_t op_t;
    typedef typename impl_type::measurement_term_desc_type measurement_term_desc_type;
    typedef typename impl_type::op_prod_type op_prod_type;
    
    typedef typename impl_type::size_t size_t;
    
    Model() { }
    
    Model(Lattice const& lattice, BaseParameters & parms)
    : impl_(model_factory<Matrix, SymmGroup>(lattice, parms))
    { }
    
    Model(impl_ptr impl) : impl_(impl) { }
    
    void update(BaseParameters const& p) { return impl_->update(p); }
    
    Index<SymmGroup> const& phys_dim(size_t type=0) const { return impl_->phys_dim(type); }
    op_t const& identity_matrix(size_t type=0) const { return impl_->identity_matrix(type); }
    tag_type identity_matrix_tag(size_t type=0) const { return impl_->identity_matrix_tag(type); }
    op_t const& filling_matrix(size_t type=0) const { return impl_->filling_matrix(type); }
    tag_type filling_matrix_tag(size_t type=0) const { return impl_->filling_matrix_tag(type); }
    
    typename SymmGroup::charge total_quantum_numbers(BaseParameters & parms) const { return impl_->total_quantum_numbers(parms); }
    
    terms_type const& hamiltonian_terms() const { return impl_->hamiltonian_terms(); }
    
    bool has_operator(std::string const & name, size_t type) const { return impl_->has_operator(name, type); }
    op_t const& get_operator(std::string const & name, size_t type=0) const { return impl_->get_operator(name, type); }
    tag_type get_operator_tag(std::string const & name, size_t type=0) const { return impl_->get_operator_tag(name, type); }

    bool has_operator(op_prod_type const & name, size_t type) const { return impl_->has_operator(name, type); }
    op_t const& get_operator(op_prod_type const & name, size_t type=0) const { return impl_->get_operator(name, type); }
    tag_type get_operator_tag(op_prod_type const & name, size_t type=0) const { return impl_->get_operator_tag(name, type); }

    std::vector<measurement_term_desc_type> unpack_measurement_terms(std::string const & name) const { return impl_->unpack_measurement_terms(name); }
    
    table_ptr operators_table() const { return impl_->operators_table(); }
    
    initializer_ptr initializer(Lattice const& lat, BaseParameters & parms) const { return impl_->initializer(lat, parms); }

private:
    impl_ptr impl_;
};


#endif
