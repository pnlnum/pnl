
add_library(libblas OBJECT 
    dasum.c daxpy.c dcabs1.c dcopy.c ddot.c dgbmv.c dgemm.c dgemv.c dger.c
    dnrm2.c drot.c drotg.c drotm.c drotmg.c dsbmv.c dscal.c dsdot.c dspmv.c
    dspr.c dspr2.c dswap.c dsymm.c dsymv.c dsyr.c dsyr2.c dsyr2k.c dsyrk.c
    dtbmv.c dtbsv.c dtpmv.c dtpsv.c dtrmm.c dtrmv.c dtrsm.c dtrsv.c dzasum.c
    dznrm2.c icamax.c idamax.c isamax.c izamax.c scabs1 
    zaxpy.c zcopy.c
    zdotc.c zdotu.c zdscal.c zgbmv.c zgemm.c zgemv.c zgerc.c zgeru.c zhbmv.c
    zhemm.c zhemv.c zher.c zher2.c zher2k.c zherk.c zhpmv.c zhpr.c zhpr2.c
    zrotg.c zscal.c zswap.c zsymm.c zsyr2k.c zsyrk.c ztbmv.c ztbsv.c ztpmv.c
    ztpsv.c ztrmm.c ztrmv.c ztrsm.c ztrsv.c csrot.c zdrot.c)

### The following files must not be compiled because they
### are already included in liblapack
### lsame.c xerbla.c xerbla_array.c 
