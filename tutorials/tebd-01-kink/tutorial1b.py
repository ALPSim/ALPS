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
import pyalps.pyplot
import numpy as np
import math
import scipy.special

#prepare the input parameters
parms=[]
count=0
for nsteps in [100, 250, 500, 750, 1000]:
	count+=1
	parms.append([{ 
	          'L'                         : 50,
	          'MODEL'                     : 'spin',
	          'local_S'                   : 0.5,
	          'CONSERVED_QUANTUMNUMBERS'  : 'Sz',
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
	        }])

baseName='tutorial_1b_'
for p in parms:
	nmlname=pyalps.write_TEBD_files(p, baseName+str(p[0]['SIMID']))
	res=pyalps.run_TEBD(nmlname)

syssize=parms[0][0]['L']
Magdata=[]
#Get magnetization data
for p in parms:
	Magdata.extend(pyalps.load.loadTimeEvolution(p, baseName+str(p[0]['SIMID']), measurements=['Local Magnetization']))

#Postprocessing
postData=[]
for q in Magdata:
	loc=-0.5*scipy.special.jn(0,q[0].props['Time'])*scipy.special.jn(0,q[0].props['Time'])
	q[0].y=[abs(q[0].y[syssize/2+1-1]-loc)]


Mag=pyalps.collectXY(Magdata, x='Time', y='Local Magnetization', foreach=['SIMID'])
plt.figure()
pyalps.pyplot.plot(Mag)
plt.xlabel('Time $t$')
plt.yscale('log')
plt.ylabel('Magnetization Error')
plt.show()


