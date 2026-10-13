#  Copyright Bela Bauer 2010-2011.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

import pyalps
from pyalps.lattice import parse, showgraph

import matplotlib.pyplot as plt

import sys

if __name__ == '__main__':
    graph = parse(sys.argv[1])
    showgraph(graph)
    plt.show()
