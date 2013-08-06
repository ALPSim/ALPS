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
import copy
import math
import scipy.special

#prepare the input parameters
parms = [{ 
          'L'                         : 50,
          'MODEL'                     : 'spin',
          'local_S'                   : 0.5,
          'CONSERVED_QUANTUMNUMBERS'  : 'Sz_total',
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


baseName='tutorial_2a'
nmlname=pyalps.writeTEBDfiles(parms, baseName)
res=pyalps.runTEBD(nmlname)

#Get the results of the simulation
Data=pyalps.load.loadTimeEvolution(pyalps.getResultFiles(prefix='tutorial_2a'), measurements=['Local Magnetization'])
#define a dataset numericalSolution to contain the numerical result
numericalResult=[]
#define a dataset exactSolution to contain the exact solution
exactResult=[]
#define a dataset scalingForm to contain the scaling form
scalingForm=[]

#Compute the exact result M(n,t)=<S_n^z>=-(1/2)*sum_{i=1-n}^{n-1} j_i(t)^2, where
# j_i(t) is the Bessel function of order i and compare to the numerically obtained result
for q in Data:
        syssize=q[0].props['L']
        #Assign a label 'Distance' denoting the distance from the center n (only do the first two sites
        #to avoid cluttering the plot)
        for n in range(1,3):
                #Create copies of the data for postprocessing
                numericalCopy=copy.deepcopy(q)
                exactCopy=copy.deepcopy(q)
                
                numericalCopy[0].props['Distance']=n
                numericalCopy[0].props['SIMID']='Numerical at n='+str(n)
                exactCopy[0].props['Distance']=n
                exactCopy[0].props['SIMID']='Exact at n='+str(n)

                #compute the exact result of the manetization n sites from the center
                loc=0.0
                for i in range(1-n,n):
                        loc-=0.5*scipy.special.jn(i,q[0].props['Time'])*scipy.special.jn(i,q[0].props['Time'])                        
                exactCopy[0].y=[loc]
                #add to the the exact dataset
                exactResult.extend(exactCopy)

                #get the numerical result of the magnetization n sites from the center
                numericalCopy[0].y=[q[0].y[syssize/2+n-1]]
                #add to the the numerical dataset
                numericalResult.extend(numericalCopy)

#compute the scaling form
# \phi(n/t)=-(1/pi)*arcsin(n/t) that M(n,t) approaches as n->infinity and t->infinity
# and compare it with the numerically computed values of M(n/t)
for q in Data:
        syssize=q[0].props['L']
        #Assign a label 'Distance' denoting the distance from the center n (only do the first few sites
        #to avoid cluttering the plot)
        for n in range(0,5):
                #Create a copy of the data for postprocessing
                scalingCopy=copy.deepcopy(q)
                scalingCopy[0].props['Distance']=n

                #The first distance contains the exact scaling form \phi(n/t)=-(1/pi)*arcsin(n/t)
                if n==0:
                        scalingCopy[0].props['Time']=1.0/scalingCopy[0].props['Time']
                        scalingCopy[0].y=[-(1.0/3.1415926)*math.asin(min(scalingCopy[0].props['Time'],1.0))]
                        scalingCopy[0].props['SIMID']='Exact'

                #The other distances contain the numerical data as a function of the scaling variable M(n/t)
                else:
                        scalingCopy[0].props['Time']=n/scalingCopy[0].props['Time']
                        scalingCopy[0].y=[scalingCopy[0].y[syssize/2+n-1] ]
                        scalingCopy[0].props['SIMID']='Numerical at n='+str(n)
                #add to the scaling dataset
                scalingForm.extend(scalingCopy)



#Plot the numerical and exact magnetization for comparison
exactMag=pyalps.collectXY(exactResult, x='Time', y='Local Magnetization',foreach=['SIMID'])
for q in exactMag:
    q.props['label']=q.props['SIMID']
numericalMag=pyalps.collectXY(numericalResult, x='Time', y='Local Magnetization',foreach=['SIMID'])
for q in numericalMag:
    q.props['label']=q.props['SIMID']

plt.figure()
pyalps.plot.plot([exactMag, numericalMag])
plt.xlabel('Time $t$')
plt.ylabel('Magnetization')
plt.legend(loc='lower right')
plt.title('Magnetization vs. time')

#Plot the scaling form with the numerical data for comparison
Scal=pyalps.collectXY(scalingForm, x='Time', y='Local Magnetization', foreach=['SIMID'])
for q in Scal:
    q.props['label']=q.props['SIMID']

plt.figure()
pyalps.plot.plot(Scal)
plt.xlabel('Scaling variable $n/t$')
plt.ylabel('Magnetization$(n,t)$')
plt.legend()
plt.xlim(0,1.5)
plt.title('Magnetization scaling function; numerical and exact results')
plt.show()





