" To load the function into Vim : source this file
" This function only operates on a region. To run this function, 
" first select a region and run :call JLExtractHeader()

function! PnlDocHeader ()
    let l:line = getline(".")
    let l:line = substitute(l:line, 'extern \+', '', 'eI')
    let l:line = substitute(l:line, '\([a-zA-Z0-9_ \*]\+\) \+\([a-zA-Z0-9_]\+\) *(\(.*\));', '\\item \\describefun{\1}{\2}{\3}', 'eI')
    let l:line = substitute(l:line, '*', '\\ptr ', 'ge')
    let l:line = substitute(l:line, '\(Pnl[a-zA-Z]\+\)', '\\\1', 'geI')
    call setline(".",l:line)
endfunction
    

