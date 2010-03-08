all: manual index manual

manual:
	@echo "running htlatex"
	@export TEXINPUTS=`pwd`/..//:: ; \
	htlatex ../pnl-manual 'xhtml,uni-html4,2,section+,next' 
	cp pnl-manual.css pnl-manual.css.bak
	cat ../pnl.css >> pnl-manual.css


	#  cp ../pnl.css pnl-manual.css

index:
	tex '\def\filename{{pnl-manual}{idx}{4dx}{ind}} \input  idxmake.4ht'
	makeindex -o pnl-manual.ind pnl-manual.4dx

clean:
	rm -f pnl-manual* idxmake.dvi  idxmake.log

