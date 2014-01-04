" This file adds detection of Pnl types for syntax highlighting
" Jerome Lelong <jerome.lelong@imag.fr>
"
" Add the content of this file to ~/.vim/after/syntax/c.vim
" or create it if it does not exist

syn keyword cType PnlRng PnlType  PnlVect PnlVectInt PnlVectComplex 
            \ PnlMat PnlMatInt PnlMatComplex PnlTridiagMat PnlTridiagMatLU 
            \ PnlBandMat PnlHmat PnlHmatInt PnlHmatComplex PnlRng PnlBasis
            \ PnlList PnlObject dcomplex PnlArray PnlPermutation PnlCell 
            \ PnlSpMat PnlSpMatInt PnlSpMatComplex
syn match cType 'Pnl[a-zA-Z]*Object'
syn match cType 'Pnl[0-9a-zA-Z]*Func[0-9a-zA-Z]*'
syn keyword cConstant CONE CZERO
