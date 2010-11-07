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
