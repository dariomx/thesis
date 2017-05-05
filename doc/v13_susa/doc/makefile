doc = thesis

all:
	pdflatex ${doc}
	bibtex ${doc}.aux
	makeglossaries ${doc}
	pdflatex ${doc}
	pdflatex ${doc}

clean:
	rm -f *.glg *.glo *.gls *.lot *.lof *~ \
              *.xdy *.pdf *.log *.out *.toc *.aux *.bbl *.blg


