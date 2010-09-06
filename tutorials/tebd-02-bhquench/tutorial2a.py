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

#prepare the input parameters
parms=[]
count=0
for A in [5.0]:#, 10.0, 15.0, 25.0, 50.0]:
	count+=1
	parms.append([{ 
	          'L'                         : 10,
	          'MODEL'                     : 'boson Hubbard',
	          'Nmax'                   : 5,
	          'CONSERVED_QUANTUMNUMBERS'  : 'true',
	          'VERBOSE'  : 'true',
		  'totQ' : 10,
	          'ITP_t'                         : 1.0,
	          'ITP_U'                         : 10.0,
	          't'                         : 1.0,
	          'U'                         : 10.0,
		  'ITP_CHIS' : [20, 30, 35],
		  'ITP_DTS' : [0.05, 0.05,0.025],
		  'ITP_CONVS' : [1E-8, 1E-8, 1E-9],
		  'INITIAL_STATE' : 'ground',
		  'CHI_LIMIT' : 40,
		  'TRUNC_LIMIT' : 1E-12,
		  'NUM_THREADS' : 1,
		  'TAUS' : [A,  A],
		  'POWS' : [1.0, 1.0],
		  'GS' : ['U',  'U'],
		  'GIS' : [10.0,  1.0],
		  'GFS' : [1.0,  10.0],
		  'NUMSTEPS' : [1000,  1000],
		  'STEPSFORSTORE' : [10, 5],
		  'SIMID' : count
	        }])
		

baseName='tutorial_2a_tau'
for p in parms:
	nmlname=pyalps.write_TEBD_files(p, baseName+str(p[0]['TAUS'][0]))
	res=pyalps.run_TEBD(nmlname)

ll=pyalps.load.Hdf5Loader()
Udata=[]
LEdata=[]

for p in parms:
	stepper=p[0]['NUMSTEPS']
	counter=0
	for d in p[0]['STEPSFORSTORE']:
		stepper[counter]/=d
		counter+=1
	stepper=[i+1 for i in range(sum(stepper))]
	for d in stepper:
		data=ll.ReadMeasurementFromFile(['./'+baseName+str(p[0]['TAUS'][0])+'.h5'],proppath='/timesteps/'+str(d).rjust(8)+'/Local Props', respath='/timesteps/'+str(d).rjust(8)+'/results', measurements=['U'])
		Udata.extend(data)
		data=ll.ReadMeasurementFromFile(['./'+baseName+str(p[0]['TAUS'][0])+'.h5'],proppath='/timesteps/'+str(d).rjust(8)+'/Local Props', respath='/timesteps/'+str(d).rjust(8)+'/results', measurements=['Loschmidt Echo'])
		LEdata.extend(data)


LE=pyalps.collectXY(LEdata, x='Time', y='Loschmidt Echo',foreach=['SIMID'])
plt.figure()
pyalps.pyplot.plot(LE)
plt.xlabel('Scaled Time $t/T$')
plt.ylabel('Loschmidt Echo $|< \psi(0)|\psi(t) > |^2$')

Ufig=pyalps.collectXY(Udata, x='Time', y='U',foreach=['SIMID'])
plt.figure()
pyalps.pyplot.plot(Ufig)
plt.xlabel('Time $t$')
plt.ylabel('U')
plt.show()







