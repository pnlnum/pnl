all: pdf

pdf: manual

TEXFILES=$(wildcard *.tex)

manual: pnl-manual.pdf

pnl-manual.pdf: deterministic_methods.tex  pnl-manual.tex types.tex \
		linear_algebra.tex probabilistic_methods.tex mathematical_functions.tex \
		some_bindings.tex pnl-manual.tex
	pdflatex pnl-manual
	makeindex -s index.isty pnl-manual
	pdflatex pnl-manual

clean:
	for i in aux idx ilg ind log out dvi toc;\
	do \
		$(RM) ${TEXFILES:%.tex=%.$$i} ;\
	done

html:
	(cd html; make)

cleanhtml:
	(cd html; make clean)

proto:
	grep extern `ls ../src/include/pnl/*.h | grep -v pnl_f2c` | awk '{ printf $$3 "\n"; }' | grep -v '^Pnl' | grep -v '^{' | sed 's/(.*$$//g'  | grep -v '\(hmat\|mat\|vect\)_\(complex\|int\|uint\)' | grep -v '^\*' | grep -v C2F | grep -v '.*;$$' > protos.h
	for f in `cat protos.h`; \
		do \
		`egrep -q $$f *.tex` || `egrep -q $$f exclude-protos.txt` || echo $$f; \
		done

