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
