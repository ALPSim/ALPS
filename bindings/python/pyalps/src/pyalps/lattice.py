#  Copyright Bela Bauer 2010-2011.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

import pyalps

import matplotlib.pyplot as plt
from xml.etree import ElementTree

import sys

colors = ['black', 'green', 'blue', 'red']

# parse file describing the lattice, and return a tuple of vertices and edges
# both are returned as dictionaries, with the ID as key
# and the value:
# - a coordinate tuple in the case of vertices
# - a dict in the case of edges, containing source, target, type and all other attributes from the XML
# A third dict is returned, which maps vertex IDs to vertex types
def parse(fn):
    root = ElementTree.parse(fn).getroot()
    vertices = {}
    vtype = {}
   
    dimension = int(root.get('dimension'))
    for vertex in root.findall('VERTEX'):
        vid = int(vertex.get('id'))
        vpos = tuple([float(x) for x in vertex.find('COORDINATE').text.split()])
        if(len(vpos) != dimension):
            raise RuntimeError('Dimension of the coordinates does not match the dimension of the lattice.')
        vertices[vid] = vpos
        vtype[vid] = vertex.get('type')
    
    edges = {}
    for edge in root.findall('EDGE'):
        eid = int(edge.get('id'))
        edges[eid] = dict([(k,edge.get(k)) for k in ['source', 'target', 'type', 'vector']])
        for c in ['source', 'target', 'type']:
            edges[eid][c] = int(edges[eid][c])
    
    return (dimension,vertices,edges,vtype)

def showgraph(graph):
    dimension = graph[0]
    vertices = graph[1]
    edges = graph[2]
    vtypes = graph[3]
    
    if(dimension > 2):
        raise RuntimeError('This function only supports 1 and 2 dimensional lattices.')

    if(len(vertices.values()[0]) == 1):
        vertices = {k: (v[0],0) for k, v in vertices.items()}
    
    x = [v[0] for v in vertices.values()]
    y = [v[1] for v in vertices.values()]
    plt.scatter(x, y)
    for k, v in vertices.items():
        plt.annotate('%s (%s)' % (k, vtypes[k]), v)
    
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    for edge in edges.values():
        s = edge['source']
        t = edge['target']
        
        p0 = vertices[s]
        p1 = vertices[t]
        
        c = colors[edge['type'] % len(colors)] if 'type' in edge else None
        plt.plot([p0[0], p1[0]], [p0[1], p1[1]], color=c)

