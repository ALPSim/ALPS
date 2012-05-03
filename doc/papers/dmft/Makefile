#  Copyright Matthias Troyer 2009 - 2011.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)



all:
				python strip.py
				latex loop.tex
				dvips -E -o loop.eps loop.dvi
				latex codes.tex
				bibtex codes
				latex codes.tex
				latex codes.tex
				dvipdf codes
				open codes.pdf
