//  (C) Copyright 2010 Lukas Gamper <gamperl -at- gmail.com>
//  Use, modification, and distribution are subject to the Boost Software 
//  License, Version 1.0. (See at <http://www.boost.org/LICENSE_1_0.txt>.)

/* $Id: abstract_task.C 3822 2010-01-30 22:02:39Z troyer $ */

#include <alps/hdf5.hpp>
#include <alps/parameter.h>

#include <iostream>

int main(int argc, char **argv) {
	if (argc < 2)
		throw std::invalid_argument("no name passed");
	alps::Parameters parms;
	std::cin >> parms;
	alps::hdf5::oarchive ar(argv[1]);
	ar << make_pvp("/parameters", parms);
}