set (libinline_SRCS
    build-inline.c)

install(FILES 
    pnl_array.h 
    pnl_band_matrix.h 
    pnl_basis.h 
    pnl_cdf.h 
    pnl_complex.h 
    pnl_deprecated.h
    pnl_fft.h 
    pnl_finance.h
    pnl_integration.h 
    pnl_interpolation.h 
    pnl_laplace.h 
    pnl_linalgsolver.h 
    pnl_list.h 
    pnl_mathtools.h 
    pnl_matrix.h 
    pnl_matvect.h 
    pnl_mpi.h 
    pnl_object.h 
    pnl_optim.h 
    pnl_perm.h 
    pnl_random.h 
    pnl_root.h 
    pnl_sp_matrix.h 
    pnl_specfun.h
    pnl_templates_off.h 
    pnl_templates_on.h 
    pnl_tridiag_matrix.h 
    pnl_vector.h 
    DESTINATION include/pnl
)

