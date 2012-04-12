/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2012 by Lukas Gamper <gamperl@gmail.com>                   *
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

#include <alps/ngs.hpp>

#include <boost/lambda/lambda.hpp>

template<typename Base> class ising_sim : public Base {

	public:

		ising_sim(typename Base::parameters_type const & params, std::size_t seed_offset = 42)
			: Base(params, seed_offset)
			, length(params["L"])
			, thermalization_sweeps(int(params["THERMALIZATION"]))
			, total_sweeps(int(params["SWEEPS"]))
			, beta(1. / double(params["T"]))
			, spins(length)
		{
			init();
		}

		template <typename Arg> ising_sim(typename Base::parameters_type const & params, Arg comm)
			: Base(params, comm)
			, length(params["L"])
			, thermalization_sweeps(int(params["THERMALIZATION"]))
			, total_sweeps(int(params["SWEEPS"]))
			, beta(1. / double(params["T"]))
			, spins(length)
		{
			init();
		}

// TODO: rename to update()
		void do_update() {
			for (int j = 0; j < length; ++j) {
				using std::exp;
				int i = int(double(length) * this->random());
				int right = ( i + 1 < length ? i + 1 : 0 );
				int left = ( i - 1 < 0 ? length - 1 : i - 1 );
				double p = exp( 2. * beta * spins[i] * ( spins[right] + spins[left] ));
				if ( p >= 1. || this->random() < p )
					spins[i] = -spins[i];
			}
		};

// TODO: rename to measure() ... use verbs
		void do_measurements() {
			sweeps++;
			if (sweeps > thermalization_sweeps) {
				double tmag = 0;
				double ten = 0;
				double sign = 1;
				std::vector<double> corr(length);
				for (int i = 0; i < length; ++i) {
					tmag += spins[i];
					sign *= spins[i];
					ten += -spins[i] * spins[ i + 1 < length ? i + 1 : 0 ];
					for (int d = 0; d < length; ++d)
						corr[d] += spins[i] * spins[( i + d ) % length ];
				}
				std::transform(corr.begin(), corr.end(), corr.begin(), boost::lambda::_1 / double(length));
				ten /= length;
				tmag /= length;
				this->measurements["Energy"] << ten;
				this->measurements["Magnetization"] << tmag;
				this->measurements["Magnetization^2"] << tmag * tmag;
				this->measurements["Magnetization^4"] << tmag * tmag * tmag * tmag;
				this->measurements["Correlations"] << corr;
			}
		};

		double fraction_completed() const {
			return (sweeps < thermalization_sweeps ? 0. : ( sweeps - thermalization_sweeps ) / double(total_sweeps));
		}

	private:

		void init() {
			for(int i = 0; i < length; ++i)
				spins[i] = (this->random() < 0.5 ? 1 : -1);
			this->measurements
				<< alps::ngs::RealObservable("Energy")
				<< alps::ngs::RealObservable("Magnetization")
				<< alps::ngs::RealObservable("Magnetization^2")
				<< alps::ngs::RealObservable("Magnetization^4")
				<< alps::ngs::RealVectorObservable("Correlations")
			;
		}

		int length;
		int sweeps;
		int thermalization_sweeps;
		int total_sweeps;
		double beta;
		std::vector<int> spins;
};
