all: manual index manual finalize

manual:
	@echo "running htlatex"
	@export TEXINPUTS=`pwd`/..//:: ; \
	htlatex ../pnl-manual 'xhtml,uni-html4,2,subsection+,next' 


index:
	tex '\def\filename{{pnl-manual}{idx}{4dx}{ind}} \input  idxmake.4ht'
	makeindex -o pnl-manual.ind pnl-manual.4dx

finalize: 
	cp pnl-manual.css pnl-manual.css.bak
	cat ../pnl.css >> pnl-manual.css

clean:
	rm -f pnl-manual* idxmake.dvi  idxmake.log

