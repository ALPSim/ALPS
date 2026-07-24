# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2009-2010 by Matthias Troyer <troyer@phys.ethz.ch> 
#                            Jan Gukelberger
#                            Brigitte Surer
# 
# SPDX-License-Identifier: MIT
# 
# ****************************************************************************

import ising

L = 4       # Linear lattice size
N = 5000    # of simulation steps

print('# L:', L, 'N:', N)

# Scan beta range [0,1] in steps of 0.1
for beta in [0.,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.]:
    print('-----------')
    print('beta =', beta)
    sim = ising.Simulation(beta,L)
    sim.run(N/2,N)
