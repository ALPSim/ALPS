//  (C) Copyright 2010 Lukas Gamper <gamperl -at- gmail.com>
//  Use, modification, and distribution are subject to the Boost Software 
//  License, Version 1.0. (See at <http://www.boost.org/LICENSE_1_0.txt>.)

#include "isingng.hpp"

template<typename S> void run_sim(alps::mcoptions const & options, int argc, char *argv[]) {
	typename S::parameters_type params(options);
	S s(params, argc, argv);
	s.run(boost::bind(&S::stop, boost::ref(s)));
	while (!s.stop())
		usleep(500000);
}
int main(int argc, char *argv[]) {
	if (argc < 2)
		throw std::invalid_argument("parameter file missing " + std::string(argv[0]) + " param.h5");
	alps::mcoptions options(argc, argv);
#ifdef ALPS_HAVE_MPI
	if (options.is_valid() && options.use_mpi())
		run_sim<parallel_sim>(options, argc, argv);
	else 
#endif
	if (options.is_valid())
		run_sim<simple_sim>(options, argc, argv);
}
