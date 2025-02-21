/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2016 Institute for Theoretical Physics, ETH Zurich
 *               2011-2016 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef DMRG_EVOLE_TROTTER_ITERATION_H
#define DMRG_EVOLE_TROTTER_ITERATION_H

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <cmath>

#include "utils/io.hpp"


class trotter_steps_iterator;


class trotter_decomposer {
    
public:
    enum tevol_order_tag {order_unknown, first_order, second_order, fourth_order};
    
    typedef std::size_t size_type;
    typedef std::vector<std::size_t> sequence_type;
    typedef std::pair<std::size_t, double> term_type;
    typedef trotter_steps_iterator steps_iterator;

    inline trotter_decomposer(size_type, std::string const&, bool);
    
    bool is_optimized() const
    {
        return te_optim;
    }
    
    size_type size() const
    {
        return gates_coeff.size();
    }
    
    term_type const& trotter_term(size_type i) const
    {
        return gates_coeff[i];
    }
    
    sequence_type const& simple_terms_sequence() const
    {
        return Useq;
    }
    sequence_type const& initial_terms_sequence() const
    {
        return Useq_initial;
    }
    sequence_type const& double_terms_sequence() const
    {
        return Useq_double;
    }
    sequence_type const& final_terms_sequence() const
    {
        return Useq_final;
    }
    
    inline steps_iterator steps_begin(size_type nsteps) const;
    
    std::string description() const
    {
        std::stringstream ss;
        print_trotter_order(ss, trotter_order);
        ss << (is_optimized() ? " with " : " without ");
        ss << "optimization";
        return ss.str();
    }
private:
    tevol_order_tag parse_trotter_order (std::string const & trotter_param)
    {
        if (trotter_param == "first")
            return first_order;
        else if (trotter_param == "second")
            return second_order;
        else if (trotter_param == "fourth")
            return fourth_order;
        else {
            throw std::runtime_error("Don't know this Trotter decomposition");
            return order_unknown;
        }
    }

    inline std::ostream& print_trotter_order (std::ostream& os, trotter_decomposer::tevol_order_tag const& o) const
    {
        switch (o)
        {
            case trotter_decomposer::first_order:
                os << "First order Trotter decomposition";
                break;
            case trotter_decomposer::second_order:
                os << "Second order Trotter decomposition";
                break;
            case trotter_decomposer::fourth_order:
                os << "Fourth order Trotter decomposition";
                break;
            default:
                os << "uknown Trotter decomposition";
        }
        return os;
    }

    
    tevol_order_tag trotter_order;
    bool te_optim;
    std::vector<term_type> gates_coeff;
    sequence_type Useq; // trivial sequence
    sequence_type Useq_double, Useq_initial, Useq_final; // sequence with two sweeps; meas before; meas after
};




/// Sequences iterator
class trotter_steps_iterator {
public:
    typedef std::size_t size_type;
    typedef trotter_decomposer::sequence_type value_type;
    typedef value_type const& reference_type;
    
    trotter_steps_iterator() : end_(true) {}
    
    trotter_steps_iterator(trotter_decomposer const& decomp, size_type nsteps)
    : end_(false)
    , i_(0)
    , decomp_(&decomp)
    , nsteps_(nsteps)
    { }
    
    inline trotter_steps_iterator& operator++()
    {
        ++i_;
        return *this;
    }
    
    inline bool operator!=(trotter_steps_iterator const& rhs) const
    {
//        std::cout << " (!=) i_=" << i_ << ", end_=" << end_ << ", rhs.i_=" << rhs.i_ << ", rhs.i_=" << rhs.i_ << std::endl;
        if (rhs.end_) return !(end_ || i_ == nsteps_);
        else if (end_) return !(rhs.end_ || rhs.i_ == rhs.nsteps_);
        else return (i_ != rhs.i_);
    }
    
    inline reference_type operator*() const
    {
        if (nsteps_ < 2 || !decomp_->is_optimized()) {
            return decomp_->simple_terms_sequence();
        } else {
            if (i_ == 0)
                return decomp_->initial_terms_sequence();
            else if (i_ == nsteps_-1)
                return decomp_->final_terms_sequence();
            else
                return decomp_->double_terms_sequence();
        }
    }
    
    inline size_type index() const
    {
        return i_;
    }
    
private:
    bool end_;
    size_type i_;
    const trotter_decomposer * decomp_;
    size_type nsteps_;
};



/// Definitions of decomposer

trotter_decomposer::trotter_decomposer(trotter_decomposer::size_type nterms,
                                       std::string const& decomposition,
                                       bool is_optimized_)
: trotter_order(parse_trotter_order(decomposition))
, te_optim(is_optimized_)
{
    /// alpha coeffiecients and set sequence of Uterms according to trotter order
    switch (trotter_order){
        case first_order:
        {
            /// All terms with alpha=1
            for (size_type i=0; i<nterms; ++i)
                gates_coeff.push_back(std::make_pair(i,1.));
            
            /// Insert all terms to the sequence
            for (size_type k=0; k<size(); ++k)
                Useq.push_back(k);
            
            if (is_optimized())
            {
                Useq_initial = Useq;
                Useq_double = Useq;
                Useq_final = Useq;
            }
            
            maquis::cout << "Sequence initialized with " << Useq.size() << " terms." << std::endl;
            break;
        }
        case second_order:
        {
            /// h_1(dt/2) * h_2(dt/2) ... h_l(dt/2) * h_l(dt/2) * h_l(dt/2) * h_l-1(dt/2) ... h_1(dt/2)
            /// -->
            /// h_1(dt/2) * h_2(dt/2) ... h_l-1(dt/2) * h_l(dt) * h_l-1(dt/2) ... h_1(dt/2)
            
            /// All terms with alpha=0.5
            for (size_type i=0; i<nterms-1; ++i)
                gates_coeff.push_back(std::make_pair(i,.5));
            /// Last term with alpha=1.
            gates_coeff.push_back(std::make_pair(nterms-1,1.)); // ix: nterms-1
            
            /// Storage is:
            /// 0: (0,0.5), 1: (1,0.5), ..., nterms-2: (nterms-2,0.5), nterms-1: (nterms-1,1.)
            
            
            /// Insert 1 --> l-1 with alpha=0.5
            for (size_type k=0; k<nterms-1; ++k)
                Useq.push_back(k);
            /// Insert l with alpha=1
            Useq.push_back(nterms-1);
            /// Insert l-1 --> 1 with alpha=0.5
            for (size_type kk=0; kk<nterms-1; ++kk)
                Useq.push_back(nterms-2-kk);
            
            
            if (is_optimized())
            {
                /// First term with alpha=1.
                gates_coeff.push_back(std::make_pair(0,1.)); // ix: nterms
                
                /// Insert 1 --> l-1 with alpha=0.5
                for (size_type k=0; k<nterms-1; ++k)
                    Useq_initial.push_back(k);
                /// Insert l with alpha=1
                Useq_initial.push_back(nterms-1);
                /// Insert l-1 --> 2 with alpha=0.5
                for (size_type kk=0; kk<nterms-2; ++kk)
                    Useq_initial.push_back(nterms-2-kk);
                /// Insert 0 with alpha=1
                Useq_initial.push_back(nterms);
                
                
                /// Insert 2 --> l-1 with alpha=0.5
                for (size_type k=1; k<nterms-1; ++k)
                    Useq_double.push_back(k);
                /// Insert l with alpha=1
                Useq_double.push_back(nterms-1);
                /// Insert l-1 --> 2 with alpha=0.5
                for (size_type kk=0; kk<nterms-2; ++kk)
                    Useq_double.push_back(nterms-2-kk);
                /// Insert 0 with alpha=1
                Useq_double.push_back(nterms);

                /// Insert 2 --> l-1 with alpha=0.5
                for (size_type k=1; k<nterms-1; ++k)
                    Useq_final.push_back(k);
                /// Insert l with alpha=1
                Useq_final.push_back(nterms-1);
                /// Insert l-1 --> 1 with alpha=0.5
                for (size_type kk=0; kk<nterms-1; ++kk)
                    Useq_final.push_back(nterms-2-kk);
            }
            
            maquis::cout << "Sequence initialized with " << Useq.size() << " terms." << std::endl;
            break;
        }
        case fourth_order:
        {
            double alpha_1=1./(4.0-pow(4.0,0.33333));
            double alpha_3=1.-4.0*alpha_1;
            
            /// All terms with alpha_1/2
            for (size_type i=0; i<nterms-1; ++i)
                gates_coeff.push_back(std::make_pair(i,alpha_1*0.5));
            /// Last term with alpha_1
            gates_coeff.push_back(std::make_pair(nterms-1,alpha_1)); // ix: nterms-1
            /// All terms with alpha_3/2
            for (size_type i=0; i<nterms-1; ++i)
                gates_coeff.push_back(std::make_pair(i,alpha_3*0.5));
            /// Last term with alpha_3
            gates_coeff.push_back(std::make_pair(nterms-1,alpha_3)); // ix: 2*nterms-2
            
            /// Storage is:
            /// 0: (0,0.5*alpha_1), 1: (1,0.5*alpha_1), ..., nterms-2: (nterms-2,0.5*alpha_1), nterms-1: (nterms-1,alpha_1)
            /// nterms: (0,0.5*alpha_3), nterms+1: (1,0.5*alpha_3), ..., 2*nterms-2: (nterms-2,0.5*alpha_3), 2*nterms-1: (nterms-1,alpha_3)
            
            
            /// Insert 1 --> l-1 with alpha=0.5*alpha_1 and l with alpha=alpha_1
            for (size_type k=0; k<nterms; ++k)
                Useq.push_back(k);
            /// Insert l-1 --> 1 with alpha=0.5*alpha_1
            for (size_type kk=0; kk<nterms-1; ++kk)
                Useq.push_back(nterms-2-kk);

            /// Insert 1 --> l-1 with alpha=0.5*alpha_1 and l with alpha=alpha_1
            for (size_type k=0; k<nterms; ++k)
                Useq.push_back(k);
            /// Insert l-1 --> 1 with alpha=0.5*alpha_1
            for (size_type kk=0; kk<nterms-1; ++kk)
                Useq.push_back(nterms-2-kk);

            /// Insert 1 --> l-1 with alpha=0.5*alpha_3 and l with alpha=alpha_3
            for (size_type k=0; k<nterms; ++k)
                Useq.push_back(nterms+k);
            /// Insert l-1 --> 1 with alpha=0.5*alpha_3
            for (size_type kk=0; kk<nterms-1; ++kk)
                Useq.push_back(2*nterms-2-kk);
            
            /// Insert 1 --> l-1 with alpha=0.5*alpha_1 and l with alpha=alpha_1
            for (size_type k=0; k<nterms; ++k)
                Useq.push_back(k);
            /// Insert l-1 --> 1 with alpha=0.5*alpha_1
            for (size_type kk=0; kk<nterms-1; ++kk)
                Useq.push_back(nterms-2-kk);

            /// Insert 1 --> l-1 with alpha=0.5*alpha_1 and l with alpha=alpha_1
            for (size_type k=0; k<nterms; ++k)
                Useq.push_back(k);
            /// Insert l-1 --> 1 with alpha=0.5*alpha_1
            for (size_type kk=0; kk<nterms-1; ++kk)
                Useq.push_back(nterms-2-kk);

            
            if (is_optimized())
            {
                gates_coeff.push_back(std::make_pair(0,alpha_1));                 // ix: 2*nterms
                gates_coeff.push_back(std::make_pair(0,alpha_3*0.5+alpha_1*0.5)); // ix: 2*nterms+1
                
                /////////////////////////////////////////////////////////////////////
                
                /// Insert 1 --> l-1 with 0.5*alpha_1 and l with alpha_1
                for (size_type k=0; k<nterms; ++k)
                    Useq_initial.push_back(k);
                /// Insert l-1 --> 2 with 0.5*alpha_1
                for (size_type kk=0; kk<nterms-2; ++kk)
                    Useq_initial.push_back(nterms-2-kk);
                /// Insert 1 with alpha_1
                Useq_initial.push_back(2*nterms);

                /// Insert 2 --> l-1 with 0.5*alpha_1 and l with alpha_1
                for (size_type k=1; k<nterms; ++k)
                    Useq_initial.push_back(k);
                /// Insert l-1 --> 2 with 0.5*alpha_1
                for (size_type kk=0; kk<nterms-2; ++kk)
                    Useq_initial.push_back(nterms-2-kk);

                /// Insert 1 with alpha_3*0.5+alpha_1*0.5
                Useq_initial.push_back(2*nterms+1);

                /// Insert 2 --> l-1 with 0.5*alpha_3 and l with alpha_3
                for (size_type k=1; k<nterms; ++k)
                    Useq_initial.push_back(nterms+k);
                /// Insert l-1 --> 2 with 0.5*alpha_3
                for (size_type kk=0; kk<nterms-2; ++kk)
                    Useq_initial.push_back(2*nterms-2-kk);
                
                /// Insert 1 with alpha_3*0.5+alpha_1*0.5
                Useq_initial.push_back(2*nterms+1);

                /// Insert 2 --> l-1 with 0.5*alpha_1 and l with alpha_1
                for (size_type k=1; k<nterms; ++k)
                    Useq_initial.push_back(k);
                /// Insert l-1 --> 2 with alpha=0.5*alpha_1
                for (size_type kk=0; kk<nterms-2; ++kk)
                    Useq_initial.push_back(nterms-2-kk);
                /// Insert 1 with alpha_1
                Useq_initial.push_back(2*nterms);

                /// Insert 2 --> l-1 with 0.5*alpha_1 and l with alpha_1
                for (size_type k=1; k<nterms; ++k)
                    Useq_initial.push_back(k);
                /// Insert l-1 --> 2 with 0.5*alpha_1
                for (size_type kk=0; kk<nterms-2; ++kk)
                    Useq_initial.push_back(nterms-2-kk);
                /// Insert 1 with alpha_1
                Useq_initial.push_back(2*nterms);

                /////////////////////////////////////////////////////////////////////
                
                std::copy(Useq_initial.begin()+1, Useq_initial.end(), std::back_inserter(Useq_double));
                
                /////////////////////////////////////////////////////////////////////
                
                Useq_final = Useq_double;
                Useq_final[Useq_final.size()-1] = 0;
                
                /////////////////////////////////////////////////////////////////////
            }
            maquis::cout << "Sequence initialized with " << Useq.size() << " terms." << std::endl;
            
            break;
        }
        default:
        {
            throw std::runtime_error("uknown Trotter decomposition");
            break;
        }
    }
    
}


trotter_decomposer::steps_iterator trotter_decomposer::steps_begin(size_type nsteps) const
{
    return steps_iterator(*this, nsteps);
}



#endif
