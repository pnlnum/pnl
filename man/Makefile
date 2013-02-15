all: pdf html

pdf: manual

TEXFILES=$(wildcard *.tex)

manual: pnl-manual.pdf

pnl-manual.pdf: bases.tex cdf.tex complex.tex fft.tex finance.tex integ.tex\
		laplace.tex  linear_algebra.tex mathematical_functions.tex \
		mathtools.tex ode.tex  optim.tex pnl-manual.tex rng.tex roots.tex \
		some_bindings.tex specfun.tex  types.tex 
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
	grep extern `ls ../src/include/pnl/*.h | grep -v pnl_f2c` | awk '{ printf $$3 "\n"; }' | grep -v '^Pnl' | grep -v '^{' | sed 's/(.*$$//g'  | grep -v '\(hmat\|mat\|vect\)_\(complex\|int\|uint\)' | grep -v '^\*' | grep -v C2F | grep -v '.*;$$'

missing-proto:
	grep extern `ls ../src/include/pnl/*.h | grep -v pnl_f2c` | awk '{ printf $$3 "\n"; }' | grep -v '^Pnl' | grep -v '^{' | sed 's/(.*$$//g'  | grep -v '\(hmat\|mat\|vect\)_\(complex\|int\|uint\)' | grep -v '^\*' | grep -v C2F | grep -v '.*;$$' > protos.h
	for f in `cat protos.h`; \
		do \
		`egrep -q $$f *.tex` || echo $$f; \
		done

.PHONY: html pdf all manual
