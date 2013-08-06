# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2010 by Michael L. Wall <mwall@mines.edu> 
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
import matplotlib.pyplot as plt
import pyalps.plot
import numpy as np
import math
import scipy.special

#prepare the input parameters
parms=[]
count=0
for nsteps in [100, 250, 500, 750, 1000]:
        count+=1
        parms.append({ 
                  'L'                         : 50,
                  'MODEL'                     : 'spin',
                  'local_S'                   : 0.5,
                  'CONSERVED_QUANTUMNUMBERS'  : 'Sz_total',
                  'Jxy'                         : 1,
                  'INITIAL_STATE' : 'kink',
                  'CHI_LIMIT' : 20,
                  'TRUNC_LIMIT' : 1E-12,
                  'NUM_THREADS' : 1,
                  'TAUS' : [20.0],
                  'POWS' : [0.0],
                  'GS' : ['H'],
                  'GIS' : [0.0],
                  'GFS' : [0.0],
                  'NUMSTEPS' : [nsteps],
                  'STEPSFORSTORE' : [int(math.floor(nsteps/100))],
                  'SIMID': count
                })


baseName='tutorial_2b_'
nmlnameList=pyalps.writeTEBDfiles(parms, baseName)
res=pyalps.runTEBD(nmlnameList)

#Get magnetization data
Magdata=pyalps.load.loadTimeEvolution( pyalps.getResultFiles(prefix='tutorial_2b'), measurements=['Local Magnetization'])

#Postprocessing-get the exact result for comparison
for q in Magdata:
        syssize=q[0].props['L']
        #Get the exact result of M(1,t)=-(1/2)*(j_0(t)^2), where j_0(t) is the 0^{th} order
        # bessel function and M(1,t) is the magnetization one site to the right of the chain center
        loc=-0.5*scipy.special.jn(0,q[0].props['Time'])*scipy.special.jn(0,q[0].props['Time'])
        #Get the difference between the computed and exact results
        q[0].y=[abs(q[0].y[syssize/2+1-1]-loc)]



#Plot the Error in the magnetization one site to the right of the chain center
Mag=pyalps.collectXY(Magdata, x='Time', y='Local Magnetization', foreach=['SIMID'])
for q in Mag:
    dt=round(q.props['TAUS']/q.props['NUMSTEPS'],3)
    q.props['label']='dt='+str(dt)

plt.figure()
pyalps.plot.plot(Mag)
plt.xlabel('Time $t$')
plt.yscale('log')
plt.ylabel('Magnetization Error')
plt.title('Error in the magnetization vs. time')
plt.legend(loc='lower left')
plt.show()


