//  (C) Copyright 2010 Lukas Gamper <gamperl -at- gmail.com>
//  Use, modification, and distribution are subject to the Boost Software 
//  License, Version 1.0. (See at <http://www.boost.org/LICENSE_1_0.txt>.)

#include "ngs.hpp"

#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>

class ising_impl : public alps::mcbase {
	public:
		ising_impl(parameters_type const & params)
			: alps::mcbase(params)
			, length(params["L"])
			, beta(1. / double(params["T"]))
			, sweeps(0)
			, thermalization_sweeps(int(params["THERMALIZATION"]))
			, total_sweeps(int(params["SWEEPS"]))
			, spins(length)
			, random(boost::mt19937(), boost::uniform_real<>())
		{
			for(int i = 0; i < length; ++i)
				spins[i] = (random() < 0.5 ? 1 : -1);
			results << alps::RealObservable("Energy");
			results << alps::RealObservable("Magnetization");
			results << alps::RealObservable("Magnetization^2");
			results << alps::RealObservable("Magnetization^4");
			results << alps::RealVectorObservable("Correlations");
		}
		void do_update() {
			for (int j = 0; j < length; ++j) {
				int i = int(double(length) * random());
				int right = ( i + 1 < length ? i + 1 : 0 );
				int left = ( i - 1 < 0 ? length - 1 : i - 1 );
				double p = exp( 2. * beta * spins[i] * ( spins[right] + spins[left] ));
				if ( p >= 1. || random() < p )
					spins[i] =- spins[i];
			}
		};
		void do_measurements() {
			if (sweeps > thermalization_sweeps) {
				tmag = 0;
				ten = 0;
				corr.resize(length, 0.);
				for (int i = 0; i < length; ++i) {
					tmag += spins[i];
					ten += -spins[i] * spins[ i + 1 < length ? i + 1 : 0 ];
					for (int d = 0; d < length; ++d)
						corr[d] += spins[i] * spins[( i + d ) % length ];
				}
				corr /= double(length);
				ten /= length;
				tmag /= length;
				results["Energy"] << ten;
				results["Magnetization"] << tmag;
				results["Magnetization^2"] << tmag * tmag;
				results["Magnetization^4"] << tmag * tmag * tmag * tmag;
				results["Correlations"] << corr;
			}
			sweeps++;
		};
		double fraction_completed() const {
			return (sweeps < thermalization_sweeps ? 0. : ( sweeps - thermalization_sweeps ) / double(total_sweeps));
		}
	private:
		int length;
		int sweeps;
		int thermalization_sweeps;
		int total_sweeps;
		double beta;
		double tmag;
		double ten;
		std::vector<int> spins;
		std::valarray<double> corr;
		boost::variate_generator<boost::mt19937, boost::uniform_real<> > random;
};
typedef alps::mcrun<ising_impl> simple_sim;
#ifdef ALPS_HAVE_MPI
	typedef alps::mcmpirun<ising_impl> parallel_sim;
#endif
