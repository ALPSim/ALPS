#  Copyright Bela Bauer 2010-2011.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

import pyalps

import matplotlib.pyplot as plt
from xml.etree import ElementTree

import sys

# parse file describing the lattice, and return a tuple of vertices and edges
# both are returned as dictionaries, with the ID as key
# and the value:
# - a coordinate tuple in the case of vertices
# - a dict in the case of edges, containing source, target, type and all other attributes from the XML
def parse(fn):
    root = ElementTree.parse(fn).getroot()
    vertices = {}
    
    for vertex in root.findall('VERTEX'):
        vid = int(vertex.get('id'))
        vpos = tuple([float(x) for x in vertex.find('COORDINATE').text.split()])
        vertices[vid] = vpos
    
    edges = {}
    for edge in root.findall('EDGE'):
        eid = int(edge.get('id'))
        edges[eid] = dict([(k,edge.get(k)) for k in ['source', 'target', 'type', 'vector']])
        for c in ['source', 'target', 'type']:
            edges[eid][c] = int(edges[eid][c])
    
    return (vertices,edges)

def showgraph(graph):
    vertices = graph[0]
    edges = graph[1]
    
    x = [v[0] for v in vertices.values()]
    y = [v[1] for v in vertices.values()]
    plt.scatter(x, y)
    for k, v in vertices.items():
        plt.annotate(k, v)
    
    for edge in edges.values():
        s = edge['source']
        t = edge['target']
        
        p0 = vertices[s]
        p1 = vertices[t]
        
        plt.plot([p0[0], p1[0]], [p0[1], p1[1]])