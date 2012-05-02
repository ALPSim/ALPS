# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2009-2010 by Brigitte Surer <surerb@phys.ethz.ch> 
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
import pyalps.fit_wrapper as fw

#prepare the input parameters
parms = []
for l in [32,48,64]:
    for t in [2.24, 2.25, 2.26, 2.27, 2.28, 2.29, 2.30, 2.31, 2.32, 2.33, 2.34, 2.35]:
        parms.append(
            { 
              'LATTICE'        : "square lattice", 
              'T'              : t,
              'J'              : 1 ,
              'THERMALIZATION' : 5000,
              'SWEEPS'         : 150000,
              'UPDATE'         : "cluster",
              'MODEL'          : "Ising",
              'L'              : l
            }
    )

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm7b',parms)
pyalps.runApplication('spinmc',input_file,Tmin=5)
# use the following instead if you have MPI
#pyalps.runApplication('spinmc',input_file,Tmin=5,MPI=4)

pyalps.evaluateSpinMC(pyalps.getResultFiles(prefix='parm7b'))

#load the susceptibility and collect it as function of temperature T
# data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm7b'),['|Magnetization|', 'Connected Susceptibility', 'Specific Heat', 'Binder Cumulant', 'Binder Cumulant U2'])
# magnetization_abs = pyalps.collectXY(data,x='T',y='|Magnetization|',foreach=['L'])
# connected_susc = pyalps.collectXY(data,x='T',y='Connected Susceptibility',foreach=['L'])
# spec_heat = pyalps.collectXY(data,x='T',y='Specific Heat',foreach=['L'])
# binder_u4 = pyalps.collectXY(data,x='T',y='Binder Cumulant',foreach=['L'])
# binder_u2 = pyalps.collectXY(data,x='T',y='Binder Cumulant U2',foreach=['L'])
# 
# #make a plot of the Binder cumulant:
# plt.figure()
# pyalps.plot.plot(binder_u4)
# plt.xlabel('Temperature $T$')
# plt.ylabel('Binder Cumulant U4 $g$')
# plt.title('2D Ising model')
# plt.show()
# 
# #perform a data collapse of the Binder cumulant: 
# Tc=2.269 #your estimate
# a=1  #your estimate
# 
# for d in binder_u4:
#     d.x -= Tc
#     d.x = d.x/Tc
#     l = d.props['L']
#     d.x = d.x * pow(float(l),a)
#     
# plt.figure()
# pyalps.plot.plot(binder_u4)
# plt.xlabel('Rescaled Temperature $(T-T_c)/T_c L^{1/\nu}$')
# plt.ylabel('Binder Cumulant U4 $g$')
# plt.title('2D Ising model')
# plt.show()
#     
# #make a plot of the specific heat and connected susceptibility:
# plt.figure()
# pyalps.plot.plot(connected_susc)
# plt.xlabel('Temperature $T$')
# plt.ylabel('Connected Susceptibility $\chi_c$')
# plt.title('2D Ising model')
# 
# plt.figure()
# pyalps.plot.plot(spec_heat)
# plt.xlabel('Temperature $T$')
# plt.ylabel('Specific Heat $c_v$')
# plt.title('2D Ising model')
# plt.show()
# 
# #make a fit of the connected susceptibility as a function of L:
# cs_mean=[]
# for q in connected_susc:
#     cs_mean.append(np.array([d.mean for d in q.y]))
# 
# peak_cs = pyalps.DataSet()
# peak_cs.props = pyalps.dict_intersect([q.props for q in connected_susc])
# peak_cs.y = np.array([np.max(q) for q in cs_mean])
# peak_cs.x = np.array([q.props['L'] for q in connected_susc])
# 
# sel = np.argsort(peak_cs.x)
# peak_cs.y = peak_cs.y[sel]
# peak_cs.x = peak_cs.x[sel]
# 
# pars = [fw.Parameter(1), fw.Parameter(1)]
# f = lambda self, x, pars: pars[0]()*np.power(x,pars[1]())
# fw.fit(None, f, pars, peak_cs.y, peak_cs.x)
# prefactor = pars[0].get()
# gamma_nu = pars[1].get()
# 
# plt.figure()
# plt.plot(peak_cs.x, f(None, peak_cs.x, pars))
# pyalps.plot.plot(peak_cs)
# plt.xlabel('System Size $L$')
# plt.ylabel('Connected Susceptibility $\chi_c(T_c)$')
# plt.title('2D Ising model, $\gamma$ is %.4s' % gamma_nu)
# plt.show()
# 
# #make a fit of the specific heat as a function of L:
# sh_mean=[]
# for q in spec_heat:
#     sh_mean.append(np.array([d.mean for d in q.y]))
#     
# peak_sh = pyalps.DataSet()
# peak_sh.props = pyalps.dict_intersect([q.props for q in spec_heat])
# peak_sh.y = np.array([np.max(q) for q in sh_mean])
# peak_sh.x = np.array([q.props['L'] for q in spec_heat])
# 
# sel = np.argsort(peak_sh.x)
# peak_sh.y = peak_sh.y[sel]
# peak_sh.x = peak_sh.x[sel]
# 
# pars = [fw.Parameter(1), fw.Parameter(1)]
# f = lambda self, x, pars: pars[0]()*np.power(x,pars[1]())
# fw.fit(None, f, pars, peak_sh.y, peak_sh.x)
# prefactor = pars[0].get()
# alpha_nu = pars[1].get()
# 
# plt.figure()
# plt.plot(peak_sh.x, f(None, peak_sh.x, pars))
# pyalps.plot.plot(peak_cs)
# plt.xlabel('System Size $L$')
# plt.ylabel('Specific Heat $c_v(T_c)$')
# plt.title(r'2D Ising model, $\alpha$ is %.4s' % alpha_nu)
# plt.show()
# 
# #make a data collapse of the connected susceptibility as a function of (T-Tc)/Tc:
# for d in connected_susc:
#     d.x -= Tc
#     d.x = d.x/Tc
#     l = d.props['L']
#     d.x = d.x * pow(float(l),a)
# 
# two_minus_eta=1.75 #your estimate
# for d in connected_susc:
#     l = d.props['L']
#     d.y = d.y/pow(float(l),two_minus_eta)
# 
# plt.figure()
# pyalps.plot.plot(connected_susc)
# plt.xlabel('Temperature $T$')
# plt.ylabel('Connected Susceptibility $\chi_c$')
# plt.title('2D Ising model')
# plt.show()
# 
# #make a data collapse of the |magnetization| as a function of (T-Tc)/Tc
# for d in magnetization_abs:
#     d.x -= Tc
#     d.x = d.x/Tc
#     l = d.props['L']
#     d.x = d.x * pow(float(l),a)
# beta_over_nu=... #your estimate    
# for d in magnetization_abs:
#     l = d.props['L']
#     d.y = d.y / pow(float(l),-beta_over_nu)
#     
# plt.figure()
# pyalps.plot.plot(magnetization_abs)
# plt.xlabel('Temperature $T$')
# plt.ylabel('Magnetization $|m|$')
# plt.title('2D Ising model')
# plt.show()
