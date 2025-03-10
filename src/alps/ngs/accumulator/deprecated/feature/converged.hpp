/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2011 - 2012 by Mario Koenz <mkoenz@ethz.ch>                       *
 *                                                                                 *
 * Permission is hereby granted, free of charge, to any person obtaining           *
 * a copy of this software and associated documentation files (the “Software”),    *
 * to deal in the Software without restriction, including without limitation       *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,        *
 * and/or sell copies of the Software, and to permit persons to whom the           *
 * Software is furnished to do so, subject to the following conditions:            *
 *                                                                                 *
 * The above copyright notice and this permission notice shall be included         *
 * in all copies or substantial portions of the Software.                          *
 *                                                                                 *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS         *
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,     *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE     *
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER          *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING         *
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER             *
 * DEALINGS IN THE SOFTWARE.                                                       *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
 
#ifndef ALPS_NGS_ALEA_DETAIL_CONVERGED_HPP
#define ALPS_NGS_ALEA_DETAIL_CONVERGED_HPP

#include <alps/ngs/alea/feature/mean.hpp>
#include <alps/ngs/alea/feature/feature_traits.hpp>
#include <alps/ngs/alea/feature/generate_property.hpp>

namespace alps
{
    namespace accumulator
    {
        enum error_convergence
        {
              yes
            , no
            , maybe
        };
        //=================== converged proxy ===================
        //=================== converged trait ===================
        template <typename T>
        struct converged_type
        {
            typedef error_convergence type;
        };
        //=================== converged implementation ===================
        namespace detail
        {

        //set up the dependencies for the tag::autocorrelation-Implementation
            template<> 
            struct Dependencies<tag::detail::converged> 
            {
                typedef MakeList<>::type type;
            };

            template<typename base_type> 
            class AccumulatorImplementation<tag::detail::converged, base_type> : public base_type 
            {
                typedef typename base_type::value_type value_type_loc;
                typedef typename converged_type<value_type_loc>::type converged_type;
                typedef typename mean_type<value_type_loc>::type mean_type;
                typedef AccumulatorImplementation<tag::detail::converged, base_type> ThisType;
                
                public:
                    AccumulatorImplementation<tag::detail::converged, base_type>(ThisType const & arg): base_type(arg)
                    
                    {}
                    template<typename ArgumentPack>
                    AccumulatorImplementation<tag::detail::converged, base_type>(ArgumentPack const & args
                                                 , typename boost::disable_if<
                                                                              boost::is_base_of<ThisType, ArgumentPack>
                                                                            , int
                                                                            >::type = 0
                                             ): base_type(args)
                    {}
                    
                    inline converged_type const converged() const 
                    {
                        //~ //Simplebinning.h Zeile 300
                        //~ template <class T>
                        //~ typename SimpleBinning<T>::convergence_type SimpleBinning<T>::converged_errors() const
                        //~ {
                          //~ convergence_type conv;
                          //~ result_type err=error();
                          //~ resize_same_as(conv,err);
                          //~ const unsigned int range=4;
                          //~ typename slice_index<convergence_type>::type it;
                          //~ if (binning_depth()<range) {
                            //~ for (it= slices(conv).first; it!= slices(conv).second; ++it)
                              //~ slice_value(conv,it) = MAYBE_CONVERGED;
                          //~ }
                          //~ else {
                            //~ for (it= slices(conv).first; it!= slices(conv).second; ++it)
                              //~ slice_value(conv,it) = CONVERGED;
                        //~ 
                            //~ for (unsigned int i=binning_depth()-range;i<binning_depth()-1;++i) {
                              //~ result_type this_err(error(i));
                              //~ for (it= slices(conv).first; it!= slices(conv).second; ++it)
                                //~ if (std::abs(slice_value(this_err,it)) >= std::abs(slice_value(err,it)))
                                  //~ slice_value(conv,it)=CONVERGED;
                                //~ else if (std::abs(slice_value(this_err,it)) < 0.824 * std::abs(slice_value(err,it)))
                                  //~ slice_value(conv,it)=NOT_CONVERGED;
                                //~ else if (std::abs(slice_value(this_err,it)) <0.9* std::abs(slice_value(err,it))  &&
                                    //~ slice_value(conv,it)!=NOT_CONVERGED)
                                  //~ slice_value(conv,it)=MAYBE_CONVERGED;
                            //~ }
                          //~ }
                          //~ return conv;
                        //~ }
                        
                        //TODO: implement
                        return maybe;
                    }
                    
                    inline void operator()(value_type_loc const & val) 
                    {
                        base_type::operator()(val);
                    }
                    inline ThisType& operator<<(value_type_loc const & val) 
                    {
                        (*this)(val);
                        return (*this);
                    }
                    
                    template<typename Stream> 
                    inline void print(Stream & os) 
                    {
                        base_type::print(os);
                        os << "tag::detail::converged: " << std::endl;
                    }
                    inline void reset()
                    {
                        base_type::reset();
                    }
                private:
            };

            template<typename base_type> class ResultImplementation<tag::detail::converged, base_type> : public base_type  {

                public:

                    template<typename Accumulator> ResultImplementation(Accumulator const & accum)
                        : base_type(accum)
                    {}
// TODO: implement!
            };
        }

        //=================== call GENERATE_PROPERTY macro ===================
        GENERATE_PROPERTY(converged, tag::detail::converged)

    }
}
#endif
