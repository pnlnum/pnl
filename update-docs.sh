#!/bin/bash

DOCS=docs
MANDIR=man

echo "Updating docs/pnl-manual.pdf"
cp -a $MANDIR/pnl-manual.pdf $DOCS/
echo "Updating docs/manual-html"
cp -a $MANDIR/html/pnl-manual*.{html,png,css} $DOCS/manual-html/
