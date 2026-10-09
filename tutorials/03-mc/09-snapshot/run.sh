#!/bin/sh

parameter2xml parm9a
simplemc parm9a.in.xml
python plot9a.py

parameter2xml parm9b
simplemc parm9b.in.xml
snap2vtk parm9b.*.snap

parameter2xml parm9c
simplemc parm9c.in.xml
snap2vtk parm9c.*.snap

parameter2xml parm9d
simplemc parm9d.in.xml
snap2vtk parm9d.*.snap
