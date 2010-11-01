" This file adds detection of Pnl types for syntax highlighting
" Jerome Lelong <jerome.lelong@imag.fr>
"
" Add the content of this file to ~/.vim/after/syntax/c.vim
" or create it if it does not exist

syn keyword cType PnlRng PnlType  PnlVect PnlVectInt PnlVectComplex PnlMat PnlMatInt PnlMatComplex PnlTridiagMat PnlBandMat PnlHmat PnlHmatInt PnlHmatComplex PnlRng PnlBasis PnlList PnlObject
syn match cType 'Pnl[a-zA-Z]*Object'
