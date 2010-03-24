//  (C) Copyright 2010 Lukas Gamper <gamperl -at- gmail.com>
//  Use, modification, and distribution are subject to the Boost Software 
//  License, Version 1.0. (See at <http://www.boost.org/LICENSE_1_0.txt>.)

/* $Id: abstract_task.C 3822 2010-01-30 22:02:39Z troyer $ */

#include "isingng.hpp"

#include <iostream>

bool callback() {
	return false;
}
int main(int argc, char **argv) {
	alps::Parameters parms;
	std::cin >> parms;
	ising_sim s(parms);
	s.run(&callback);
	s.save("sim.h5");
}