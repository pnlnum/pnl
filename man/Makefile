do: manual

TEXFILES=$(wildcard *.tex)

manual: pnl-manual.pdf

pnl-manual.pdf: deterministic_methods.tex  pnl-manual.tex types.tex linear_algebra.tex probabilistic_methods.tex mathematical_functions.tex some_bindings.tex pnl-manual.tex
	latex -src pnl-manual
	makeindex -s index.isty pnl-manual
	latex -src pnl-manual
	dvipdf pnl-manual.dvi

clean:
	for i in aux idx ilg ind log out dvi toc;\
	do \
		$(RM) ${TEXFILES:%.tex=%.$$i} ;\
	done

html:
	(cd pnl-man.html; make)

cleanhtml:
	(cd pnl-man.html; make clean)

proto:
	grep extern ../src/include/pnl/*.h | awk '{ printf $$3 "\n"; }' | grep -v '^Pnl' | grep -v '^{' | sed 's/(.*$$//g'  | grep -v '\(hmat\|mat\|vect\)_\(complex\|int\|uint\)' | grep -v sparse | grep -v morse | grep -v '^\*' | grep -v '.*;$$' > protos.h
	for f in `cat protos.h`; \
		do \
		`egrep -q $$f *.tex` || echo $$f; \
		done

