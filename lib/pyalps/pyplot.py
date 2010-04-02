# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2009-2010 by Bela Bauer <bauerb@phys.ethz.ch>
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

import matplotlib.pyplot as plt
import numpy as np
from hlist import flatten
from dataset import DataSet

colors = ['k','b','g','m','c','y']
markers = ['s', 'o', '^', '>', 'v', '<', 'd', 'p', 'h', '8', '+', 'x']


def plot(data):
    lines = []
    icolor = 0
    imarker = 0
    if isinstance(data,DataSet):
      s = [data]
    else:
      s = data
    for q in flatten(s):
        try:
            xmeans = np.array([xx.mean for xx in q.x])
            xerrors = np.array([xx.error for xx in q.x])
        except AttributeError:
            xmeans = [float(vvv) for vvv in q.x]
            xerrors = None
        
        try:
            ymeans = np.array([xx.mean for xx in q.y])
            yerrors = np.array([xx.error for xx in q.y])
        except AttributeError:
            ymeans = [float(vvv) for vvv in q.y]
            yerrors = None

        if 'label' in q.props and q.props['label'] != 'none':
            lab = q.props['label']
        elif 'filename' in q.props:
            lab = q.props['filename']
        else:
            lab = None

        if 'xlabel' in q.props:
            plt.xlabel(q.props['xlabel'])

        if 'ylabel' in q.props:
            plt.ylabel(q.props['ylabel'])

        if 'title' in q.props:
            plt.title(q.props['title'])

        if 'line' in q.props and q.props['line'] == 'scatter':
            plt.scatter(xmeans, ymeans, c=colors[icolor], marker=markers[imarker], label=lab)
            imarker = (imarker+1)%len(markers)
        else:
            line_props = colors[icolor]
            if 'line' in q.props:
                line_props += q.props['line']
            
            plt.errorbar(xmeans,ymeans,yerr=yerrors,xerr=xerrors,fmt=line_props,label=lab)
                    
        
        icolor = (icolor+1)%len(colors)


