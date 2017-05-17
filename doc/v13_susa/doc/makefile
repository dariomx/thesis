doc = thesis
cover = images/portada.pdf

all:
	pdflatex ${doc}
	bibtex ${doc}.aux
	makeglossaries ${doc}
	pdflatex ${doc}
	pdflatex ${doc}
	pdftk $(cover) $(doc).pdf cat output $(doc)-with-cover.pdf

clean:
	rm -f *.glg *.glo *.gls *.lot *.lof *~ \
              *.xdy *.pdf *.log *.out *.toc *.aux *.bbl *.blg


