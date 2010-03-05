do: manual

TEXFILES=$(wildcard *.tex)

manual: pnl-manual.tex
	latex -src $<
	makeindex -s index.isty ${<:%.tex=%}
	latex -src $<
	dvipdf ${<:%.tex=%.dvi}

clean:
	for i in aux idx ilg ind log out dvi toc;\
	do \
		$(RM) ${TEXFILES:%.tex=%.$$i} ;\
	done

proto:
	grep extern ../src/include/*.h | awk '{ printf $$3 "\n"; }' | grep -v '^Pnl' | grep -v '^{' | sed 's/(.*$$//g'  | grep -v '\(hmat\|mat\|vect\)_\(complex\|int\|uint\)' | grep -v sparse | grep -v morse | grep -v '^\*' | grep -v '.*;$$' > protos.h
	for f in `cat protos.h`; \
		do \
		`egrep -q $$f *.tex` || echo $$f; \
		done
		

