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

#prepare the input parameters
parms=[]
count=0
for A in [5.0, 10.0, 15.0, 25.0, 50.0]:
        count+=1
        parms.append({ 
                  'L'                         : 10,
                  'MODEL'                     : 'hardcore boson',
                  'CONSERVED_QUANTUMNUMBERS'  : 'N_total',
                  'N_total' : 5,
                  't'                         : 1.0,
                  'V'                         : 10.0,
                  'ITP_CHIS' : [20, 30, 35],
                  'ITP_DTS' : [0.05, 0.05,0.025],
                  'ITP_CONVS' : [1E-8, 1E-8, 1E-9],
                  'INITIAL_STATE' : 'ground',
                  'CHI_LIMIT' : 40, 
                  'TRUNC_LIMIT' : 1E-12,
                  'NUM_THREADS' : 1,
                  'TAUS' : [A,  A],
                  'POWS' : [ 1.0, 1.0],
                  'GS' : ['V',  'V'],
                  'GIS' : [10.0,  0.0],
                  'GFS' : [0.0,  10.0],
                  'NUMSTEPS' : [500,  500],
                  'STEPSFORSTORE' : [5, 3],
                  'SIMID' : count
                })
                

baseName='tutorial_1a'
#write output files
nmlnameList=pyalps.writeTEBDfiles(parms, baseName)
#run the application
res=pyalps.runTEBD(nmlnameList)

#Load the loschmidt echo and U
LEdata=pyalps.load.loadTimeEvolution(pyalps.getResultFiles(prefix='tutorial_1a'), measurements=['Loschmidt Echo', 'V'])

LE=pyalps.collectXY(LEdata, x='Time', y='Loschmidt Echo',foreach=['SIMID'])
for q in LE:
    q.props['label']=r'$\tau=$'+str(q.props['TAUS'][0])

plt.figure()
pyalps.plot.plot(LE)
plt.xlabel('Time $t$')
plt.ylabel('Loschmidt Echo $|< \psi(0)|\psi(t) > |^2$')
plt.title('Loschmidt Echo vs. Time')
plt.legend(loc='lower right')

Ufig=pyalps.collectXY(LEdata, x='Time', y='V',foreach=['SIMID'])
for q in Ufig:
    q.props['label']=r'$\tau=$'+str(q.props['TAUS'][0])

plt.figure()
pyalps.plot.plot(Ufig)
plt.xlabel('Time $t$')
plt.ylabel('V')
plt.title('Interaction parameter $V$ vs. Time')
plt.legend(loc='lower right')
plt.show()









