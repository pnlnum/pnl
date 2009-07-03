do: manual

TEXFILES=$(wildcard *.tex)

manual: pnl-manual.tex
	latex -src $<
	makeindex -s index.isty ${<:%.tex=%}
	bibtex ${<:%.tex=%}
	latex -src $<
	dvipdf ${<:%.tex=%.dvi}

clean:
	for i in aux idx ilg ind log out dvi toc;\
	do \
		$(RM) ${TEXFILES:%.tex=%.$$i} ;\
	done

