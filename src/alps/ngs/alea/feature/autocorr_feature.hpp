/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2011 - 2012 by Mario Koenz <mkoenz@ethz.ch>                       *
 *                                                                                 *
 * This software is part of the ALPS libraries, published under the ALPS           *
 * Library License; you can use, redistribute it and/or modify it under            *
 * the terms of the license, either version 1 or (at your option) any later        *
 * version.                                                                        *
 *                                                                                 *
 * You should have received a copy of the ALPS Library License along with          *
 * the ALPS Libraries; see the file LICENSE.txt. If not, the license is also       *
 * available from http://alps.comp-phys.org/.                                      *
 *                                                                                 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,        *
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT       *
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE       *
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,     *
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER     *
 * DEALINGS IN THE SOFTWARE.                                                       *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


#ifndef ALPS_NGS_ALEA_DETAIL_AUTOCORR_IMPLEMENTATION_HEADER
#define ALPS_NGS_ALEA_DETAIL_AUTOCORR_IMPLEMENTATION_HEADER

#include <alps/ngs/alea/accumulator/accumulator_impl.hpp>
#include <alps/ngs/alea/global_enum.hpp>
#include <alps/ngs/alea/features.hpp>

#include <boost/cstdint.hpp>

#include <ostream>
#include <algorithm>
#include <vector>
#include <cmath>

namespace alps
{
    namespace alea
    {
        //=================== autocorr proxy ===================
        template<typename value_type>
        class autocorr_proxy_type
        {
            typedef typename mean_type<value_type>::type mean_type;
            static std::vector<value_type> unused;
        public:
            autocorr_proxy_type(): bin2_(unused)
                                , bin1_(unused)
                                , count_(0) 
            {}
            autocorr_proxy_type(  std::vector<value_type> const & bin2
                                , std::vector<value_type> const & bin1
                                , boost::uint64_t const & count):
                                                                  bin2_(bin2)
                                                                , bin1_(bin1)
                                                                , count_(count)
            {}
            
            inline std::vector<value_type> const & bins() const 
            {
                return bin2_;
            }
            inline std::vector<value_type> const & sum() const 
            {
                return bin1_;
            }
            
            inline mean_type const error(boost::uint64_t level = -1) const
            {
                if(level == -1)
                    level = std::max(bin2_.size() - 5, typename std::vector<value_type>::size_type(0));
                using std::sqrt;
                return sqrt((bin2_[level] - bin1_[level] * bin1_[level]) / (count_ - 1));
            }
            
            template<typename T>
            friend std::ostream & operator<<(std::ostream & os, autocorr_proxy_type<T> const & arg);
        private:
            std::vector<value_type> const & bin2_;
            std::vector<value_type> const & bin1_;
            boost::uint64_t const count_;
        };
        
        template<typename T>
        inline std::ostream & operator<<(std::ostream & os, autocorr_proxy_type<T> const & arg)
        {
            os << "autocorr_proxy" << std::endl;
            return os;
            
        };
        //=================== autocorr trait ===================
        template <typename T>
        struct autocorr_type
        {
            typedef autocorr_proxy_type<T> type;
        };
        //=================== autocorr implementation ===================
        namespace detail
        {

        //set up the dependencies for the tag::autocorrelation-Implementation
            template<> 
            struct Dependencies<tag::autocorrelation> 
            {
                typedef MakeList<tag::mean, tag::error, tag::detail::tau, tag::detail::converged>::type type;
            };

            template<typename base_type> 
            class Implementation<tag::autocorrelation, base_type> : public base_type 
            {
                typedef typename base_type::value_type value_type_loc;
                typedef typename autocorr_type<value_type_loc>::type autocorr_type;
                typedef typename std::vector<value_type_loc>::size_type size_type;
                typedef typename mean_type<value_type_loc>::type mean_type;
                typedef Implementation<tag::autocorrelation, base_type> ThisType;
                
                public:
                    Implementation<tag::autocorrelation, base_type>(ThisType const & arg): base_type(arg)
                                                                        , bin_(arg.bin_)
                                                                        , partial_(arg.partial_)
                                                                        , bin_size_now_(arg.bin_size_now_)
                    {}
                    
                    template<typename ArgumentPack>
                    Implementation<tag::autocorrelation, base_type>(ArgumentPack const & args
                                                 , typename boost::disable_if<
                                                                              boost::is_base_of<ThisType, ArgumentPack>
                                                                            , int
                                                                            >::type = 0
                                             ): base_type(args)
                                              , bin_size_now_(0)
                    {}
                    
                    inline autocorr_type const autocorr() const 
                    { 
                        return autocorr_proxy_type<value_type_loc>(bin_, partial_, base_type::count());
                    }
                    
                    inline ThisType& operator <<(value_type_loc const &  val) 
                    {
                        base_type::operator<<(val);
                        
                        if(base_type::count() == (1 << bin_size_now_))
                        {
                            bin_.push_back(value_type_loc());
                            partial_.push_back(value_type_loc());
                            ++bin_size_now_;
                        }
                        for (unsigned i = 0; i < bin_size_now_; ++i)
                        {
                            if(base_type::count() % (1lu<<i) == 0)
                            {
                                partial_[i] = base_type::sum_ - partial_[i];
                                bin_[i] += partial_[i]*partial_[i];
                                partial_[i] = base_type::sum_;
                            }
                        }
                        
                        return *this;
                    }
                    
                    template<typename Stream> 
                    inline void print(Stream & os) 
                    {
                        base_type::print(os);
                        os << "tag::autocorrelation: " << std::endl;
                        
                        //~ os << std::endl;
                        //~ for (unsigned int i = 0; i < bin_.size(); ++i)
                        //~ {
                            //~ os << "bin[" << i << "] = " << bin_[i] << std::endl;
                        //~ }
                    }
                    
                private:
                    std::vector<value_type_loc> bin_;
                    std::vector<value_type_loc> partial_;
                    size_type bin_size_now_;
            };
            
        } // end namespace detail
    }//end alea namespace 
}//end alps namespace
#endif // ALPS_NGS_ALEA_DETAIL_AUTOCORR_IMPLEMENTATION