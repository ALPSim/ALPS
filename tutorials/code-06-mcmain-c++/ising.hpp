/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2013 by Lukas Gamper <gamperl@gmail.com>                   *
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

#ifndef ALPS_TUTORIAL_ISING_HPP
#define ALPS_TUTORIAL_ISING_HPP

#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>

#include <alps/ngs/params.hpp>
#include <alps/ngs/accumulator.hpp>

#include <boost/function.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>

#include <vector>
#include <string>

class ALPS_DECL ising_sim {

    typedef alps::accumulator::accumulator_set accumulators_type;

    public:

        typedef alps::params parameters_type;
        typedef std::vector<std::string> result_names_type;
        typedef alps::accumulator::result_set results_type;

        ising_sim(parameters_type const & params);

        void update();
        void measure();
        double fraction_completed() const;
        bool run(boost::function<bool ()> const & stop_callback);

        result_names_type result_names() const;
        result_names_type unsaved_result_names() const;
        results_type collect_results() const;
        results_type collect_results(result_names_type const & names) const;

        void save(boost::filesystem::path const & filename) const;
        void load(boost::filesystem::path const & filename);
        void save(alps::hdf5::archive & ar) const;
        void load(alps::hdf5::archive & ar);

    protected:

        parameters_type parameters;
        boost::variate_generator<boost::mt19937, boost::uniform_real<> > random;
        accumulators_type measurements;

    private:
        
        int length;
        int sweeps;
        int thermalization_sweeps;
        int total_sweeps;
        double beta;
        std::vector<int> spins;
};

#endif
