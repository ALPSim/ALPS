# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 1994-2010 by Bela Bauer <bauerb@phys.ethz.ch>
#                            Brigitte Surer <surerb@phys.ethz.ch>
#
# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
# 
# ****************************************************************************

import numpy as np
from xml.etree import ElementTree
from .dataset import DataSet
from .floatwitherror import FloatWithError as fwe
from .floatwitherror import get_mean
from .hlist import flatten

def read_xml(filename):
    root = ElementTree.parse(filename).getroot()
    data = DataSet()

    data.props['xlabel'] = root.find('xaxis').attrib['label']
    data.props['ylabel'] = root.find('yaxis').attrib['label']

    x = []
    y = []
    dx = []
    dy = []
    for point in list(root.find('set')):
        x.append(float(point.find('x').text))
        y.append(float(point.find('y').text))
        if point.find('dx') is not None:
            x[-1] = fwe(x[-1],float(point.find('dx').text))
        if point.find('dy') is not None:
            y[-1] = fwe(y[-1],float(point.find('dy').text))

    data.x = np.array(x)
    data.y = np.array(y)
    
    parameters = root.find('PARAMETERS')
    for par in list(parameters):
        data.props[par.attrib['name']] = par.text
    
    return data

def Axis(label=None,mmin=None,mmax=None,log=False):
    d = {}
    if label is not None:
        d['label'] = label
    if mmin is not None:
        d['min'] = min
    if mmax is not None:
        d['max'] = max
    if log is not None:
        d['log'] = log
    return d

def Legend(location=None):
    d = {}
    if location is not None:
        d['location'] = location
    return d

def Plot(data,xaxis=None,yaxis=None,legend=None):
    d = {'data':data}
    if xaxis is not None:
        d['xaxis'] = xaxis
    if yaxis is not None:
        d['yaxis'] = yaxis
    if legend is not None:
        d['legend'] = legend

def convertToText(data,title=None,xaxis=None,yaxis=None):
    output = ''
    if  title is not None:
        output += title + '\n'
    
    if xaxis is not None:
        output += '# X'
        if 'label' in xaxis:
            output += ': ' + xaxis['label']
        if 'min' in xaxis and 'max' in xaxis:
            output += ': ' + str(xaxis['min']) + ' to ' + str(xaxis['max'])
        output+='\n'

    if yaxis is not None:
        output += '# Y'
        if 'label' in yaxis:
            output += ': ' + yaxis['label']
        if 'min' in yaxis and 'max' in yaxis:
            output += ': ' + str(yaxis['min']) + ' to ' + str(yaxis['max'])
        output+='\n\n'
            
    
    for q in flatten(data):
        if 'label' in q.props and q.props['label'] != 'none':
            output += '# ' + q.props['label']
        elif 'filename' in q.props:
            output += '# ' + q.props['filename']
        output += '\n'
        if 'xlabel' in q.props:
            output += '# X: ' + q.props['xlabel'] + '\n'
        if 'ylabel' in q.props:
            output += '# Y: ' + q.props['ylabel'] + '\n'

        for i in range(len(q.x)):
            output += str(q.x[i]) + '\t' + str(q.y[i]) + '\n'
        output+='\n\n'                
    return output
        
def convert_to_text(desc):
    """ converts a plot descriptor to a text string """
    if 'title' in desc: t = desc['title'] 
    else: t = None
    if 'xaxis' in desc: x = desc['xaxis'] 
    else: x = None
    if 'yaxis' in desc: y = desc['yaxis'] 
    else: y = None
    return convertToText(desc['data'],title=t,xaxis=x,yaxis=y)    

            
def makeGracePlot(data,title=None,xaxis=None,yaxis=None,legend=None):
        output =  '# Grace project file\n'
        output += '#\n@    g0 on\n@    with g0\n'
        output += '@     frame linewidth 2.0\n'
        output += '@    page background fill off\n'

        xrange = [0,1]
        yrange = [0,1]
        if xaxis is not None:
          if 'min' in xaxis and 'max' in xaxis: 
            xrange = [ xaxis['min'],xaxis['max']]
        if yaxis is not None:
          if 'min' in yaxis and 'max' in yaxis:
            yrange = [ yaxis['min'],yaxis['max']]

        output += '@    world ' + str(xrange[0])+', ' + str (yrange[0]) + ','
        output +=                 str(xrange[1])+', ' + str (yrange[1]) + '\n'

        if title is not None:
            output += '@    title "'+ title + '"\n'           
            output += '@    title size 1.500000\n'

        xlog = False
        ylog = False
        if xaxis is not None:
          if  'logarithmic' in xaxis: 
            xlog = xaxis['logarithmic']
        if yaxis is not None:
          if 'logarithmic' in yaxis and yaxis['logarithmic'] is not None:
            ylog = yaxis['logarithmic']
            
        if xlog:
            output += '@    xaxes scale Logarithmic\n'
        else:
            output += '@    xaxes scale Normal\n'

        if ylog:
            output += '@    yaxes scale Logarithmic\n'
        else:
            output += '@    yaxes scale Normal\n'

        if xaxis is not None:
          if 'label' in xaxis:
            output += '@    xaxis  label "' + xaxis['label'] +'"\n'
            output += '@    xaxis  label char size 1.500000\n'
        else:
            for ds in data:
                if 'xlabel' in ds.props:
                    output += '@    xaxis label "%s"\n' % ds.props['xlabel']
                    output += '@    xaxis  label char size 1.500000\n'
                    break
        output += '@    xaxis  ticklabel char size 1.250000\n'
        output += '@    xaxis  tick minor ticks 4\n'

        if yaxis is not None:
          if 'label' in yaxis:
            output += '@    yaxis  label "' + yaxis['label'] +'"\n'
            output += '@    yaxis  label char size 1.500000\n'
        else:
            for ds in data:
                if 'ylabel' in ds.props:
                    output += '@    yaxis  label "%s"\n' % ds.props['ylabel']
                    output += '@    yaxis  label char size 1.500000\n'
                    break
        output += '@    yaxis  ticklabel char size 1.250000\n'
        output += '@    yaxis  tick minor ticks 4\n'
        
        if legend is not None and legend != False:
            output += '@    legend on\n'
            output += '@    legend loctype view\n'
            output += '@    legend 0.85, 0.8\n'
        
        num = 0
        symnum = 0
        for q in flatten(data):
            output += '@target G0.S'+str(num)+'\n'
            output += '@    s'+str(num)+' symbol ' + str(num+1) +'\n'
            output += '@    s'+str(num)+' symbol size 0.500000\n'
            if 'line' in q.props and q.props['line'] == 'scatter':
              symnum += 1
              output += '@    s'+str(num)+' line type 0\n'
              output += '@    s'+str(num)+' symbol ' + str(symnum) + '\n'
              output += '@    s'+str(num)+' symbol size 1.000000\n'
            else:
              output += '@    s'+str(num)+' line type 1\n'
            if 'label' in q.props and q.props['label'] != 'none':
                output += '@    s'+str(num)+' legend "' + q.props['label'] + '"\n'
            elif 'filename' in q.props:
                output += '@    s'+str(num)+' legend "' + q.props['filename'] + '"\n'
            output += '\n'

            if len(q.y):
                try:
                    xerrors = np.array([xx.error for xx in q.x])
                except AttributeError:
                    xerrors = None
                
                try:
                    yerrors = np.array([xx.error for xx in q.y])
                except AttributeError:
                    yerrors = None
                    
                if xerrors is None and yerrors is None:
                    output += '@type xy\n'
                    for i in range(len(q.x)):
                        output += str(q.x[i]) + '\t' + str(q.y[i]) + '\n'
                if xerrors is None and yerrors is not None:
                    output += '@type xydy\n'
                    for i in range(len(q.x)):
                        output += str(q.x[i]) + '\t' + str(q.y[i].mean) + '\t' + str(q.y[i].error) + '\n'
                if xerrors is not None and yerrors is None:
                    output += '@type xydx\n'
                    for i in range(len(q.x)):
                        output += str(q.x[i]) + '\t' + str(q.y[i].mean) + '\t' + str(q.x[i].error) + '\n'
                if xerrors is not None and yerrors is not None:
                    output += '@type xydxdy\n'
                    for i in range(len(q.x)):
                        output += str(q.x[i]) + '\t' + str(q.y[i].mean) + '\t' + str(q.x[i].error) + '\t' + str(q.x[i].error) + '\n'
                output += '&\n'
                num+=1
                     
        return output

def convert_to_grace(desc):
    """ converts a plot descriptor to a grace plot string """
    if 'title' in desc: t = desc['title'] 
    else: t = None
    if 'xaxis' in desc: x = desc['xaxis'] 
    else: x = None
    if 'yaxis' in desc: y = desc['yaxis'] 
    else: y = None
    if 'legend' in desc: l = desc['legend'] 
    else: l = None
    return makeGracePlot(desc['data'],title=t,xaxis=x,yaxis=y,legend=l)    

  
def makeGnuplotPlot(data,title=None,xaxis=None,yaxis=None,legend=None, outfile=None, terminal=None, fontsize=24):
    output =  '# Gnuplot project file\n'
    if outfile is not None:
        output += 'set output "' + outfile + '"\n'
        if terminal is None:
          if outfile[-3:]=='eps':
            terminal='postscript color eps enhanced ' + str(fontsize)
          if outfile[-3:]=='pdf':
            terminal='pdf color enhanced'
    if terminal is not None:
        output += 'set terminal ' + str(terminal) +'\n'
    if xaxis is not None:
      if 'min' in xaxis and 'max' in xaxis: 
        xrange = [ xaxis['min'],xaxis['max']]
        output += 'set xrange [' + str(xrange[0])+': ' + str (xrange[1]) + ']\n'
    if yaxis is not None:
      if 'min' in yaxis and 'max' in yaxis:
        yrange = [ yaxis['min'],yaxis['max']]
        output += 'set yrange [' + str(yrange[0])+': ' + str (yrange[1]) + ']\n'
    
    if title is not None:
        output += 'set title "'+ title + '"\n'           

    xlog = False
    ylog = False
    if xaxis is not None:
      if 'logarithmic' in xaxis:
        xlog = xaxis['logarithmic']
        output += 'set xlogscale \n'
    else:
        output += '# no xlogscale \n'

    if yaxis is not None:
      if 'logarithmic' in yaxis:
        ylog = yaxis['logarithmic']
        output += 'set ylogscale\n'
    else:
        output += '# no ylogscale\n'
        

    if xaxis is not None:
      if 'label' in xaxis:
        output += 'set xlabel "' + xaxis['label'] +'"\n'

    if yaxis is not None:
      if 'label' in yaxis:
        output += 'set ylabel "' + yaxis['label'] +'"\n'
                    
    if legend is not None and legend != False:
        output += 'set key top right\n'
        
    num = 0
    output += 'plot '
    
    for q in flatten(data):
        if len(q.y):
            try:
                xerrors = np.array([xx.error for xx in q.x])
            except AttributeError:
                xerrors = None
                
            try:
                yerrors = np.array([xx.error for xx in q.y])
            except AttributeError:
                yerrors = None
        if 'line' in q.props and q.props['line'] == 'scatter':
            if 'label' in q.props:
                if xerrors is None and yerrors is None:
                    output += ' "-" using 1:2 title "' + q.props['label'] + '",'
                if xerrors is None and yerrors is not None:
                    output += ' "-" using 1:2:3 w yerrorbars  title "' + q.props['label'] + '",'
                if xerrors is not None and yerrors is None:
                    output += ' "-" using 1:2:3 w xerrorbars  title "' + q.props['label'] + '",'
                if xerrors is not None and yerrors is not None:
                    output += ' "-" using 1:2:3:4 w xyerrorbars  title "' + q.props['label'] + '",'
            else:
                if xerrors is None and yerrors is None:
                    output += ' "-" using 1:2 notitle ,"' 
                if xerrors is None and yerrors is not None:
                    output += ' "-" using 1:2:3 w yerrorbars  notitle ,' 
                if xerrors is not None and yerrors is None:
                    output += ' "-" using 1:2:3 w xerrorbars  notitle ,' 
                if xerrors is not None and yerrors is not None:
                    output += ' "-" using 1:2:3:4 w xyerrorbars  notitle ,'
        else:
            if 'label' in q.props:
                if xerrors is None and yerrors is None:
                    output += ' "-" using 1:2 title "' + q.props['label'] + '",'
                if xerrors is None and yerrors is not None:
                    output += ' "-" using 1:2:3 w yerrorline  title "' + q.props['label'] + '",'
                if xerrors is not None and yerrors is None:
                    output += ' "-" using 1:2:3 w xerrorline  title "' + q.props['label'] + '",'
                if xerrors is not None and yerrors is not None:
                    output += ' "-" using 1:2:3:4 w xyerrorline  title "' + q.props['label'] + '",'
            else:
                if xerrors is None and yerrors is None:
                    output += ' "-" using 1:2 notitle ,"' 
                if xerrors is None and yerrors is not None:
                    output += ' "-" using 1:2:3 w yerrorline  notitle ,' 
                if xerrors is not None and yerrors is None:
                    output += ' "-" using 1:2:3 w xerrorline  notitle ,' 
                if xerrors is not None and yerrors is not None:
                    output += ' "-" using 1:2:3:4 w xyerrorline  notitle ,'
    output=output[:-1]
    output+='\n'
    
    for q in flatten(data):    
            if xerrors is None and yerrors is None:
                output += '# X Y \n'
                for i in range(len(q.x)):
                    output += str(q.x[i]) + '\t' + str(q.y[i]) + '\n'
                output += 'end \n'
            if xerrors is None and yerrors is not None:
                output += '# X Y DY \n'
                for i in range(len(q.x)):
                    output += str(q.x[i]) + '\t' + str(q.y[i].mean) + '\t' + str(q.y[i].error) + '\n'
                output += 'end \n'
            if xerrors is not None and yerrors is None:
                output += '# X Y DX \n'
                for i in range(len(q.x)):
                    output += str(q.x[i].mean) + '\t' + str(q.y[i]) + '\t' + str(q.x[i].error) + '\n'
                output += 'end \n'
            if xerrors is not None and yerrors is not None:
                output += '# X Y DXY \n'
                for i in range(len(q.x)):
                    output += str(q.x[i].mean) + '\t' + str(q.y[i].mean) + '\t' + str(q.x[i].error) + '\t' + str(q.y[i].error) + '\n'
                output += 'end \n'
            output += '\n'
            num+=1
                     
    return output

def convert_to_gnuplot(desc, outfile=None, terminal=None, fontsize=24):
    """ converts a plot descriptor to a gnuplot string """
    if 'title' in desc: t = desc['title'] 
    else: t = None
    if 'xaxis' in desc: x = desc['xaxis'] 
    else: x = None
    if 'yaxis' in desc: y = desc['yaxis'] 
    else: y = None
    if 'legend' in desc: l = desc['legend'] 
    else: l = None
    return makeGnuplotPlot(desc['data'],title=t,xaxis=x,yaxis=y,legend=l,outfile=outfile, terminal=terminal,fontsize=fontsize)    





