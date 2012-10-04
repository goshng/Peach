OUTEXTS=aux bbl blg log pdf dvi lof lot toc tex-temp a.tex RData tex

DOCS = tomato

all: ${DOCS:%=%.pdf}
	cp tomato.pdf ~/Google\ Drive/Shares/tomato/

tomato.pdf: tomato.tex tomato.R
	cp /Users/goshng/Dropbox/Documents/Projects/Peach/notes/latex/bst/jbactas.bst .
	cp /Users/goshng/Dropbox/Documents/Projects/Peach/notes/latex/sty/Sweave.sty .
	cp /Users/goshng/Dropbox/Documents/Projects/Peach/notes/latex/bib/burne.bib .
	pdflatex $<
	bibtex tomato
	cp $< tomato.a.tex
	pdflatex $<
	pdflatex $<

tomato.tex: tomato.Rnw
	R CMD Sweave $<

tomato.R: tomato.Rnw
	R CMD Stangle $<

clean:
	rm -f $(foreach d,${DOCS},${OUTEXTS:%=$d.%})
	rm -f ${DOCS:%=%*.pdf}
	rm -f ${DOCS:%=%*.tex}
	rm -f jbactas.bst Sweave.sty burne.bib
