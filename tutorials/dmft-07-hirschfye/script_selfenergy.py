# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2012-2013 by Jakub Imriska  <jimriska@phys.ethz.ch>
# 
# This software is part of the ALPS libraries, published under the ALPS
# Library License; you can use, redistribute it and/or modify it under
# the terms of the license, either version 1 or (at your option) any later
# version.
#  
# You should have received a copy of the ALPS Library License along with
# the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
# available from http://alps.comp-phys.org/.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
# SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
# FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
# DEALINGS IN THE SOFTWARE.
# 
# ****************************************************************************

import pyalps
import numpy as np
import matplotlib.pyplot as plt
import pyalps.plot
from math import pi
import sys

if '-?' in sys.argv or '--?' in sys.argv:
    print "Script for plotting of the selfenergy in Matsubara representation, Sigma_{flavor=0}(i w_n)."
    print "Selfenergy is computed employing the Dyson equation from G0_omega and G_omega."
    print "Optional parameters:"
    print "  -min_it value: to show only iterations>=value"
    print "  -real_part: to plot the real part of Sigma_{flavor=0}(i w_n) [default: imaginary part]"
    sys.exit()

min_it_=-1
if '-min_it' in sys.argv:
    min_it_=int(sys.argv[sys.argv.index('-min_it')+1])

real_part = ('-real_part' in sys.argv)

listobs=['0']   # we look at convergence of a single flavor (=0) 

## load all results
data_G = pyalps.loadDMFTIterations(pyalps.getResultFiles(pattern='parm_*.h5'), observable="G_omega", measurements=listobs, verbose=True)
data_G0 = pyalps.loadDMFTIterations(pyalps.getResultFiles(pattern='parm_*.h5'), observable="G0_omega", measurements=listobs, verbose=True)

## create a figure for each BETA
grouped_G = pyalps.groupSets(pyalps.flatten(data_G), ['BETA'])
for sim in grouped_G:
    common_props = pyalps.dict_intersect([ d.props for d in sim ])
    sim_ = [s for s in pyalps.flatten(sim) if int(s.props['iteration'])>=min_it_]
    
    ## compute selfenergy, rescale x-axis and set label
    for d_G in sim_:
        d_G0 = [s for s in pyalps.flatten(data_G0) if s.props['iteration']==d_G.props['iteration'] and s.props['BETA']==common_props['BETA']][0]
        d_G.x = np.array([(2.*n+1)*pi/common_props['BETA'] for n in d_G.x])
        Sigma = np.array([1./d_G0.y[w] - 1./d_G.y[w] for w in range(len(d_G.y))])
        if real_part:
            d_G.y = np.array(Sigma.real)
        else:
            d_G.y = np.array(Sigma.imag)
        d_G.props['label'] = "it"+d_G.props['iteration']
        d_G.props['line']="scatter"
        d_G.props['fillmarkers'] = False
    
    ## plot all iterations for this BETA
    plt.figure()
    plt.xlabel(r'$i\omega_n$')
    if real_part:
        plt.ylabel(r'$Re\ \Sigma_{flavor=0}(i\omega_n)$')
    else:
        plt.ylabel(r'$Im\ \Sigma_{flavor=0}(i\omega_n)$')
    plt.title('Simulation at ' + r'$\beta = %.4s$' % common_props['BETA'])
    pyalps.plot.plot(sim_)
    plt.legend()

plt.show()
