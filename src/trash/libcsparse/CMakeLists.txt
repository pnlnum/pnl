
########### next target ###############

SET(csparse_STAT_SRCS
   cs_add.c
   cs_dupl.c
   cs_load.c
   cs_post.c
   cs_sqr.c
   cs_amd.c
   cs_entry.c
   cs_lsolve.c
   cs_print.c
   cs_symperm.c
   cs_chol.c
   cs_ereach.c
   cs_ltsolve.c
   cs_pvec.c
   cs_tdfs.c
   cs_cholsol.c
   cs_etree.c
   cs_lu.c
   cs_qr.c
   cs_transpose.c
   cs_compress.c
   cs_fkeep.c
   cs_lusol.c
   cs_qrsol.c
   cs_updown.c
   cs_counts.c
   cs_gaxpy.c
   cs_malloc.c
   cs_randperm.c
   cs_usolve.c
   cs_cumsum.c
   cs.h
   cs_maxtrans.c
   cs_reach.c
   cs_util.c
   cs_dfs.c
   cs_happly.c
   cs_multiply.c
   cs_scatter.c
   cs_utsolve.c
   cs_dmperm.c
   cs_house.c
   cs_norm.c
   cs_scc.c
   cs_droptol.c
   cs_ipvec.c
   cs_permute.c
   cs_schol.c
   cs_dropzeros.c
   cs_leaf.c
   cs_pinv.c
   cs_spsolve.c
)

add_library(csparse STATIC ${csparse_STAT_SRCS})


########### install files ###############




#original Makefile.am contents follow:

#SUBDIRS =
#
#include $(top_srcdir)/Make.incl
#AM_CPPFLAGS = $(PNL_INCLUDES)
#
#noinst_LTLIBRARIES = libcsparse.la
#libcsparse_la_SOURCES = \
#cs_add.c        cs_dupl.c    cs_load.c      cs_post.c      cs_sqr.c\
#cs_amd.c        cs_entry.c   cs_lsolve.c    cs_print.c     cs_symperm.c\
#cs_chol.c       cs_ereach.c  cs_ltsolve.c   cs_pvec.c      cs_tdfs.c\
#cs_cholsol.c    cs_etree.c   cs_lu.c        cs_qr.c        cs_transpose.c\
#cs_compress.c   cs_fkeep.c   cs_lusol.c     cs_qrsol.c     cs_updown.c\
#cs_counts.c     cs_gaxpy.c   cs_malloc.c    cs_randperm.c  cs_usolve.c\
#cs_cumsum.c     cs.h         cs_maxtrans.c  cs_reach.c     cs_util.c\
#cs_dfs.c        cs_happly.c  cs_multiply.c  cs_scatter.c   cs_utsolve.c\
#cs_dmperm.c     cs_house.c   cs_norm.c      cs_scc.c       cs_droptol.c\
#cs_ipvec.c   cs_permute.c   cs_schol.c      cs_dropzeros.c  cs_leaf.c  \
#cs_pinv.c      cs_spsolve.c
#libcsparse_la_LDFLAGS=-avoid-version	
