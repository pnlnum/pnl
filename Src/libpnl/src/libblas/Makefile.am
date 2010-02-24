SUBDIRS = 

include $(top_srcdir)/Make.incl
AM_CPPFLAGS = $(PNL_INCLUDES)

noinst_LTLIBRARIES = libblas.la
libblas_la_SOURCES = \
caxpy.c ccopy.c cdotc.c cdotu.c cgbmv.c cgemm.c cgemv.c cgerc.c cgeru.c chbmv.c\
chemm.c chemv.c cher.c cher2.c cher2k.c cherk.c chpmv.c chpr.c chpr2.c crotg.c\
cscal.c csscal.c cswap.c csymm.c csyr2k.c csyrk.c ctbmv.c ctbsv.c ctpmv.c\
ctpsv.c ctrmm.c ctrmv.c ctrsm.c ctrsv.c dasum.c daxpy.c dcabs1.c dcopy.c ddot.c\
dgbmv.c dgemm.c dgemv.c dger.c dnrm2.c drot.c drotg.c dsbmv.c dscal.c dspmv.c\
dspr.c dspr2.c dswap.c dsymm.c dsymv.c dsyr.c dsyr2.c dsyr2k.c dsyrk.c dtbmv.c\
dtbsv.c dtpmv.c dtpsv.c dtrmm.c dtrmv.c dtrsm.c dtrsv.c dzasum.c dznrm2.c\
icamax.c idamax.c isamax.c izamax.c  sasum.c saxpy.c scasum.c scnrm2.c\
scopy.c sdot.c sgbmv.c sgemm.c sgemv.c sger.c snrm2.c srot.c srotg.c ssbmv.c\
sscal.c sspmv.c sspr.c sspr2.c sswap.c ssymm.c ssymv.c ssyr.c ssyr2.c ssyr2k.c\
ssyrk.c stbmv.c stbsv.c stpmv.c stpsv.c strmm.c strmv.c strsm.c strsv.c \
zaxpy.c zcopy.c zdotc.c zdotu.c zdscal.c zgbmv.c zgemm.c zgemv.c zgerc.c zgeru.c\
zhbmv.c zhemm.c zhemv.c zher.c zher2.c zher2k.c zherk.c zhpmv.c zhpr.c zhpr2.c\
zrotg.c zscal.c zswap.c zsymm.c zsyr2k.c zsyrk.c ztbmv.c ztbsv.c ztpmv.c ztpsv.c\
ztrmm.c ztrmv.c ztrsm.c ztrsv.c

## The following files must not be compiled because they
## are already included in liblapck
## lsame.c xerbla.c

libblas_la_LDFLAGS=-avoid-version

