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
import math
import scipy.special

#prepare the input parameters
parms=[]
count=0
for z in [0.0, 0.3, 0.9, 1.0, 1.1, 1.5]:
    count+=1
    parms.append({ 
              'L'                         : 50,
              'MODEL'                     : 'spin',
              'local_S'                   : 0.5,
              'CONSERVED_QUANTUMNUMBERS'  : 'Sz_total',
              'Jxy'                         : 1,
              'Jz'                         : z,
          'INITIAL_STATE' : 'kink',
          'CHI_LIMIT' : 40,
          'TRUNC_LIMIT' : 1E-12,
          'NUM_THREADS' : 1,
          'TAUS' : [20.0],
          'POWS' : [0.0],
          'GS' : ['H'],
          'GIS' : [0.0],
          'GFS' : [0.0],
          'NUMSTEPS' : [500],
          'STEPSFORSTORE' : [5],
          'SIMID': count
            })

baseName='tutorial_2d'
nmlnameList=pyalps.writeTEBDfiles(parms, baseName)
res=pyalps.runTEBD(nmlnameList)

#Get magnetization data
Magdata=pyalps.load.loadTimeEvolution( pyalps.getResultFiles(prefix='tutorial_2d'), measurements=['Local Magnetization'])

#Compute the integrated magnetization across the center
for q in Magdata:
    syssize=q[0].props['L']
    #Compute the integrated flow of magnetization through the center \Delta M=\sum_{n>L/2}^{L} (<S_n^z(t)>+1/2)
    #\Delta M= L/4
    loc=0.5*(syssize/2)
    #\Delta M-=<S_n^z(t)> from n=L/2 to L
    q[0].y=[0.5*(syssize/2)+sum(q[0].y[syssize/2:syssize])]


#Plot the integrated magnetization
Mag=pyalps.collectXY(Magdata, x='Time', y='Local Magnetization', foreach=['Jz'])

plt.figure()
pyalps.plot.plot(Mag)
plt.xlabel('Time $t$')
plt.ylabel('Integrated Magnetization $\Delta M(t)$')
plt.title('Integrated Magnetization vs. Time')
plt.legend(loc='upper left')
plt.show()


