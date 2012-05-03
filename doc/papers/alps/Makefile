#!/bin/sh
#  Copyright Matthias Troyer 2009 - 2011.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)


all:
				pdflatex --shell-escape alps20.tex
				bibtex alps20
				pdflatex --shell-escape alps20.tex
				pdflatex --shell-escape alps20.tex
				svn add cached vistrails_images
clean:
				\rm *.aux *.bbl *.log *.cmdline *.cmdline.* *.out *.blg
clear_cache:
				svn remove cached vistrails_images
