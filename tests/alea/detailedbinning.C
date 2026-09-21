/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2007 by Matthias Troyer <troyer@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <iostream>
#include <alps/alea.h>
#include <alps/parameter.h>
#include <boost/random.hpp>

int main()
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  //DEFINE RANDOM NUMBER GENERATOR
  //------------------------------
  typedef boost::minstd_rand0 random_base_type;
  typedef boost::uniform_01<random_base_type> random_type;
  random_base_type random_int;
  random_type random(random_int);

  //DEFINE OBSERVABLES
  //------------------
  alps::RealObservable obs_a("observable a");
  alps::RealObservable obs_b("observable b");
  alps::SimpleRealObservable obs_c("observable c");
  alps::RealObservable obs_d("observable d");

  //READ PARAMETERS
  //---------------
  alps::Parameters parms(std::cin);
  unsigned thermalization_steps=parms.value_or_default("THERMALIZATION",1000);
  unsigned number_of_steps=parms.value_or_default("STEPS",10000);

  //THERMALIZATION STEPS, NO MEASUREMENTS
  //-----------------------------------
  for(unsigned i = 0; i < thermalization_steps; ++i){
    random();
    random();
  }


  //ADD MEASUREMENTS TO THE OBSERVABLES
  //-----------------------------------
  for(unsigned i = 0; i < number_of_steps; ++i){
    obs_a << random();
    obs_b << random()+1;
    obs_c << 1.2;
    obs_d << 1.2;
  }

  //OUTPUT OBSERVABLES
  //---------------------
  std::cout << obs_a << "\n"
            << obs_b << "\n"
            << obs_c << "\n"
            << obs_d << "\n";

  //JACKKNIVE ANALYSIS
  //------------------
  alps::RealObsevaluator obseval_a(obs_a);
  alps::RealObsevaluator obseval_b(obs_b);
  alps::RealObsevaluator obseval_c;
  obseval_c = obseval_b / obseval_a;
  std::cout << obseval_c;

#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& exc) {
  std::cerr << exc.what() << "\n";
  return -1;
}
catch (...) {
  std::cerr << "Fatal Error: Unknown Exception!\n";
  return -2;
}
#endif
  return 0;
}
