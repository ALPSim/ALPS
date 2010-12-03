all:
				pdflatex --shell-escape alps20.tex
				bibtex alps20
				pdflatex --shell-escape alps20.tex
				pdflatex --shell-escape alps20.tex
clean:
				\rm *.aux *.bbl *.log *.cmdline *.cmdline.* *.out *.blg
clear_cache:
				\rm -rf cached/*.tex vistrails_images/
