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
parms = [{ 
          'L'                         : 50,
          'MODEL'                     : 'spin',
          'local_S'                   : 0.5,
          'CONSERVED_QUANTUMNUMBERS'  : 'true',
          'Jxy'                         : 1,
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
	  'STEPSFORSTORE' : [2]
        }]

baseName='tutorial_1a'
nmlname=pyalps.write_TEBD_files(parms, baseName)
res=pyalps.run_TEBD(nmlname)

ll=pyalps.load.Hdf5Loader()

#Get magnetization data
stepper=parms[0]['NUMSTEPS']
counter=0
syssize=parms[0]['L']
for d in parms[0]['STEPSFORSTORE']:
	stepper[counter]/=d
	counter+=1
stepper=[i+1 for i in range(sum(stepper))]
Magdata=[]
for d in stepper:
	for i in range(1,5):
		data=ll.ReadMeasurementFromFile(['./'+baseName+'.h5'],proppath='/timesteps/'+str(d).rjust(8)+'/Local Props', \
		respath='/timesteps/'+str(d).rjust(8)+'/results', measurements=['Local Magnetization'])
		for q in data:
			q[0].props['Distance']=i
			if i%2==0:
				loc=0.0
				for n in range(1-i+1,i-1):
					loc-=0.5*scipy.special.jn(n,q[0].props['Time'])*scipy.special.jn(n,q[0].props['Time'])
				q[0].y=[loc]
			else :
				q[0].y=[q[0].y[syssize/2+i-1] ]
		Magdata.extend(data)

Mag=pyalps.collectXY(Magdata, x='Time', y='Local Magnetization', foreach=['Distance'])
plt.figure()
pyalps.pyplot.plot(Mag)
plt.xlabel('Time $t$')
plt.ylabel('Magnetization')


#Get Scaling data
Scaldata=[]
for d in stepper:
	for i in range(1,10,2):
		data=ll.ReadMeasurementFromFile(['./'+baseName+'.h5'],proppath='/timesteps/'+str(d).rjust(8)+'/Local Props', \
		respath='/timesteps/'+str(d).rjust(8)+'/results', measurements=['Local Magnetization'])
		for q in data:
			if i==1:
				q[0].props['Time']=i/q[0].props['Time']
				q[0].y=[-(1.0/3.1415926)*math.asin(min(q[0].props['Time'],1.0))]
				q[0].props['Distance']=i
			else:
				q[0].props['Time']=i/q[0].props['Time']
				q[0].props['Distance']=i
				q[0].y=[q[0].y[syssize/2+i-1] ]
		Scaldata.extend(data)


Scal=pyalps.collectXY(Scaldata, x='Time', y='Local Magnetization', foreach=['Distance'])
plt.figure()
pyalps.pyplot.plot(Scal)
plt.xlabel('Scaling variable $n/t$')
plt.ylabel('Magnetization$(n,t)$')
plt.xlim(0,1.5)
plt.show()


