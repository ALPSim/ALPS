import matplotlib.pyplot as plt
import numpy as np

def Axis(label=None,mmin=None,mmax=None,log=False):
	d = {}
	if label != None:
		d['label'] = label
	if mmin != None:
		d['min'] = min
	if mmax != None:
		d['max'] = max
	if log != None:
		d['log'] = log
	return d

def Legend(location=None):
	d = {}
	if location != None:
		d['location'] = location
	return d

def Plot(data,xaxis=None,yaxis=None,legend=None):
	d = {'data':data}
	if xaxis != None:
		d['xaxis'] = xaxis
	if yaxis != None:
		d['yaxis'] = yaxis
	if legend != None:
		d['legend'] = legend

class MplXYPlot_core:
	colors = ['k','b','g','m','c','y']
	
	def __init__(self):
		self.icolor = 0
		self.output = ''
	
	def draw_lines(self):
		self.lines = []
		self.icolor = 0
		
		xlog = False
		ylog = False
		if 'xaxis' in self.plt and 'logarithmic' in self.plt['xaxis']:
			xlog = self.plt['xaxis']['logarithmic']
		if 'yaxis' in self.plt and 'logarithmic' in self.plt['yaxis']:
			ylog = self.plt['yaxis']['logarithmic']
		
		print self.plt	
		for q in self.plt['data']:
			if xlog and ylog:
				self.acc('self.lines.append(plt.loglog())')
				self.lines.append(plt.loglog(q.x,q.y,self.colors[self.icolor]))
			elif xlog:
				self.lines.append(plt.semilogx(q.x,q.y,self.colors[self.icolor]))
			elif ylog:
				self.lines.append(plt.semilogy(q.x,q.y,self.colors[self.icolor]))
			else:
				self.lines.append(plt.plot(q.x,q.y,self.colors[self.icolor]))
			
			if 'label' in q.props and q.props['label'] != 'none':
				self.lines[-1][0].set_label(q.props['label'])
			elif 'filename' in q.props:
				self.lines[-1][0].set_label(q.props['filename'])
			
			self.icolor = (self.icolor+1)%len(self.colors)
			
	def __call__(self, desc):
		self.plt = desc
		
		self.draw_lines()
		
		if 'xaxis' in self.plt:
			if 'label' in self.plt['xaxis']:
				plt.xlabel(self.plt['xaxis']['label'])
			if 'min' in self.plt['xaxis'] and 'max' in self.plt['xaxis']:
				plt.xlim(self.plt['xaxis']['min'],self.plt['xaxis']['max'])
				
		if 'yaxis' in self.plt:
			if 'label' in self.plt['yaxis']:
				plt.ylabel(self.plt['yaxis']['label'])
			if 'min' in self.plt['yaxis'] and 'max' in self.plt['yaxis']:
				plt.ylim(self.plt['yaxis']['min'],self.plt['yaxis']['max'])
		
		if 'legend' in self.plt:
			if 'location' in self.plt['legend']:
				plt.legend(loc=self.plt['legend']['location'])
			else:
				plt.legend()
		
		if 'title' in self.plt:
			plt.title(self.plt['title'])