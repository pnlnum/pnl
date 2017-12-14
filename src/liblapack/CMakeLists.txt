set(ALLAUX  maxloc.c ilaenv.c ieeeck.c lsamen.c xerbla.c xerbla_array.c iparmq.c	
   ilaprec.c ilatrans.c ilauplo.c iladiag.c chla_transtype.c  
	slamch.c lsame.c dlamch.c)

set (DZLAUX  
  dbdsdc.c 
  dbdsqr.c ddisna.c dlabad.c dlacpy.c dladiv.c dlae2.c  dlaebz.c 
  dlaed0.c dlaed1.c dlaed2.c dlaed3.c dlaed4.c dlaed5.c dlaed6.c 
  dlaed7.c dlaed8.c dlaed9.c dlaeda.c dlaev2.c dlagtf.c 
  dlagts.c dlamrg.c dlanst.c 
  dlapy2.c dlapy3.c dlarnv.c 
  dlarra.c dlarrb.c dlarrc.c dlarrd.c dlarre.c dlarrf.c dlarrj.c 
  dlarrk.c dlarrr.c dlaneg.c 
  dlartg.c dlaruv.c dlas2.c  dlascl.c 
  dlasd0.c dlasd1.c dlasd2.c dlasd3.c dlasd4.c dlasd5.c dlasd6.c 
  dlasd7.c dlasd8.c dlasda.c dlasdq.c dlasdt.c 
  dlaset.c dlasq1.c dlasq2.c dlasq3.c dlasq4.c dlasq5.c dlasq6.c 
  dlasr.c  dlasrt.c dlassq.c dlasv2.c dpttrf.c dstebz.c dstedc.c 
  dsteqr.c dsterf.c dlaisnan.c disnan.c )

set (DLASRC_MINIMAL  
  dgbbrd.c dgbcon.c dgbequ.c dgbrfs.c dgbsv.c  
  dgbsvx.c dgbtf2.c dgbtrf.c dgbtrs.c dgebak.c dgebal.c dgebd2.c 
  dgebrd.c dgecon.c dgeequ.c dgees.c  dgeesx.c dgeev.c  dgeevx.c 
  dgegs.c  dgegv.c  dgehd2.c dgehrd.c dgelq2.c dgelqf.c 
  dgels.c  dgelsd.c dgelss.c dgelsx.c dgelsy.c dgeql2.c dgeqlf.c 
  dgeqp3.c dgeqpf.c dgeqr2.c dgeqrf.c dgerfs.c dgerq2.c dgerqf.c 
  dgesc2.c dgesdd.c dgesv.c  dgesvd.c dgesvx.c dgetc2.c dgetf2.c 
  dgetrf.c dgetri.c 
  dgetrs.c dggbak.c dggbal.c dgges.c  dggesx.c dggev.c  dggevx.c 
  dggglm.c dgghrd.c dgglse.c dggqrf.c 
  dggrqf.c dggsvd.c dggsvp.c dgtcon.c dgtrfs.c dgtsv.c  
  dgtsvx.c dgttrf.c dgttrs.c dgtts2.c dhgeqz.c 
  dhsein.c dhseqr.c dlabrd.c dlacon.c dlacn2.c 
  dlaein.c dlaexc.c dlag2.c  dlags2.c dlagtm.c dlagv2.c dlahqr.c 
  dlahrd.c dlahr2.c dlaic1.c dlaln2.c dlals0.c dlalsa.c dlalsd.c 
  dlangb.c dlange.c dlangt.c dlanhs.c dlansb.c dlansp.c 
  dlansy.c dlantb.c dlantp.c dlantr.c dlanv2.c 
  dlapll.c dlapmt.c 
  dlaqgb.c dlaqge.c dlaqp2.c dlaqps.c dlaqsb.c dlaqsp.c dlaqsy.c 
  dlaqr0.c dlaqr1.c dlaqr2.c dlaqr3.c dlaqr4.c dlaqr5.c 
  dlaqtr.c dlar1v.c dlar2v.c iladlr.c iladlc.c 
  dlarf.c  dlarfb.c dlarfg.c dlarft.c dlarfx.c dlargv.c 
  dlarrv.c dlartv.c dlarfp.c 
  dlarz.c  dlarzb.c dlarzt.c dlaswp.c dlasy2.c dlasyf.c 
  dlatbs.c dlatdf.c dlatps.c dlatrd.c dlatrs.c dlatrz.c dlatzm.c dlauu2.c 
  dlauum.c dopgtr.c dopmtr.c dorg2l.c dorg2r.c 
  dorgbr.c dorghr.c dorgl2.c dorglq.c dorgql.c dorgqr.c dorgr2.c 
  dorgrq.c dorgtr.c dorm2l.c dorm2r.c 
  dormbr.c dormhr.c dorml2.c dormlq.c dormql.c dormqr.c dormr2.c 
  dormr3.c dormrq.c dormrz.c dormtr.c dpbcon.c dpbequ.c dpbrfs.c 
  dpbstf.c dpbsv.c  dpbsvx.c 
  dpbtf2.c dpbtrf.c dpbtrs.c dpocon.c dpoequ.c dporfs.c dposv.c  
  dposvx.c dpotf2.c dpotrf.c dpotri.c dpotrs.c dpstrf.c dpstf2.c 
  dppcon.c dppequ.c 
  dpprfs.c dppsv.c  dppsvx.c dpptrf.c dpptri.c dpptrs.c dptcon.c 
  dpteqr.c dptrfs.c dptsv.c  dptsvx.c dpttrs.c dptts2.c drscl.c  
  dsbev.c  dsbevd.c dsbevx.c dsbgst.c dsbgv.c  dsbgvd.c dsbgvx.c 
  dsbtrd.c  dspcon.c dspev.c  dspevd.c dspevx.c dspgst.c 
  dspgv.c  dspgvd.c dspgvx.c dsprfs.c dspsv.c  dspsvx.c dsptrd.c 
  dsptrf.c dsptri.c dsptrs.c dstegr.c dstein.c dstev.c  dstevd.c dstevr.c 
  dstevx.c dsycon.c dsyev.c  dsyevd.c dsyevr.c 
  dsyevx.c dsygs2.c dsygst.c dsygv.c  dsygvd.c dsygvx.c dsyrfs.c 
  dsysv.c  dsysvx.c 
  dsytd2.c dsytf2.c dsytrd.c dsytrf.c dsytri.c dsytrs.c dtbcon.c 
  dtbrfs.c dtbtrs.c dtgevc.c dtgex2.c dtgexc.c dtgsen.c 
  dtgsja.c dtgsna.c dtgsy2.c dtgsyl.c dtpcon.c dtprfs.c dtptri.c 
  dtptrs.c 
  dtrcon.c dtrevc.c dtrexc.c dtrrfs.c dtrsen.c dtrsna.c dtrsyl.c 
  dtrti2.c dtrtri.c dtrtrs.c dtzrqf.c dtzrzf.c dstemr.c 
  dlag2s.c slag2d.c dlat2s.c 
  dlansf.c dpftrf.c dpftri.c dpftrs.c dsfrk.c dtfsm.c dtftri.c dtfttp.c 
  dtfttr.c dtpttf.c dtpttr.c dtrttf.c dtrttp.c 
  dgejsv.c  dgesvj.c  dgsvj0.c  dgsvj1.c 
  dgeequb.c dsyequb.c dpoequb.c dgbequb.c )

set (ZLASRC_MINIMAL  
  zbdsqr.c zgbbrd.c zgbcon.c zgbequ.c zgbrfs.c zgbsv.c  zgbsvx.c 
  zgbtf2.c zgbtrf.c zgbtrs.c zgebak.c zgebal.c zgebd2.c zgebrd.c 
  zgecon.c zgeequ.c zgees.c  zgeesx.c zgeev.c  zgeevx.c 
  zgegs.c  zgegv.c  zgehd2.c zgehrd.c zgelq2.c zgelqf.c 
  zgels.c  zgelsd.c zgelss.c zgelsx.c zgelsy.c zgeql2.c zgeqlf.c zgeqp3.c 
  zgeqpf.c zgeqr2.c zgeqrf.c zgerfs.c zgerq2.c zgerqf.c 
  zgesc2.c zgesdd.c zgesv.c  zgesvd.c zgesvx.c zgetc2.c zgetf2.c zgetrf.c 
  zgetri.c zgetrs.c 
  zggbak.c zggbal.c zgges.c  zggesx.c zggev.c  zggevx.c zggglm.c 
  zgghrd.c zgglse.c zggqrf.c zggrqf.c 
  zggsvd.c zggsvp.c 
  zgtcon.c zgtrfs.c zgtsv.c  zgtsvx.c zgttrf.c zgttrs.c zgtts2.c zhbev.c  
  zhbevd.c zhbevx.c zhbgst.c zhbgv.c  zhbgvd.c zhbgvx.c zhbtrd.c 
  zhecon.c zheev.c  zheevd.c zheevr.c zheevx.c zhegs2.c zhegst.c 
  zhegv.c  zhegvd.c zhegvx.c zherfs.c zhesv.c  zhesvx.c zhetd2.c 
  zhetf2.c zhetrd.c 
  zhetrf.c zhetri.c zhetrs.c zhgeqz.c zhpcon.c zhpev.c  zhpevd.c 
  zhpevx.c zhpgst.c zhpgv.c  zhpgvd.c zhpgvx.c zhprfs.c zhpsv.c  
  zhpsvx.c 
  zhptrd.c zhptrf.c zhptri.c zhptrs.c zhsein.c zhseqr.c zlabrd.c 
  zlacgv.c zlacon.c zlacn2.c zlacp2.c zlacpy.c zlacrm.c zlacrt.c zladiv.c 
  zlaed0.c zlaed7.c zlaed8.c 
  zlaein.c zlaesy.c zlaev2.c zlags2.c zlagtm.c 
  zlahef.c zlahqr.c 
  zlahrd.c zlahr2.c zlaic1.c zlals0.c zlalsa.c zlalsd.c zlangb.c zlange.c 
  zlangt.c zlanhb.c 
  zlanhe.c 
  zlanhp.c zlanhs.c zlanht.c zlansb.c zlansp.c zlansy.c zlantb.c 
  zlantp.c zlantr.c zlapll.c zlapmt.c zlaqgb.c zlaqge.c 
  zlaqhb.c zlaqhe.c zlaqhp.c zlaqp2.c zlaqps.c zlaqsb.c 
  zlaqr0.c zlaqr1.c zlaqr2.c zlaqr3.c zlaqr4.c zlaqr5.c 
  zlaqsp.c zlaqsy.c zlar1v.c zlar2v.c ilazlr.c ilazlc.c 
  zlarcm.c zlarf.c  zlarfb.c 
  zlarfg.c zlarft.c zlarfp.c 
  zlarfx.c zlargv.c zlarnv.c zlarrv.c zlartg.c zlartv.c 
  zlarz.c  zlarzb.c zlarzt.c zlascl.c zlaset.c zlasr.c  
  zlassq.c zlaswp.c zlasyf.c 
  zlatbs.c zlatdf.c zlatps.c zlatrd.c zlatrs.c zlatrz.c zlatzm.c zlauu2.c 
  zlauum.c zpbcon.c zpbequ.c zpbrfs.c zpbstf.c zpbsv.c  
  zpbsvx.c zpbtf2.c zpbtrf.c zpbtrs.c zpocon.c zpoequ.c zporfs.c 
  zposv.c  zposvx.c zpotf2.c zpotrf.c zpotri.c zpotrs.c zpstrf.c zpstf2.c 
  zppcon.c zppequ.c zpprfs.c zppsv.c  zppsvx.c zpptrf.c zpptri.c zpptrs.c 
  zptcon.c zpteqr.c zptrfs.c zptsv.c  zptsvx.c zpttrf.c zpttrs.c zptts2.c 
  zrot.c   zspcon.c zspmv.c  zspr.c   zsprfs.c zspsv.c  
  zspsvx.c zsptrf.c zsptri.c zsptrs.c zdrscl.c zstedc.c 
  zstegr.c zstein.c zsteqr.c zsycon.c zsymv.c  
  zsyr.c   zsyrfs.c zsysv.c  zsysvx.c zsytf2.c zsytrf.c zsytri.c 
  zsytrs.c ztgevc.c ztgex2.c 
  ztgexc.c ztgsen.c ztgsja.c ztgsna.c ztgsyl.c ztgsy2.c 
  ztptri.c 
  ztrevc.c ztrexc.c ztrsen.c ztrsna.c 
  ztrsyl.c ztrti2.c ztrtri.c ztrtrs.c ztzrqf.c ztzrzf.c zung2l.c 
  zung2r.c zungbr.c zunghr.c zungl2.c zunglq.c zungql.c zungqr.c zungr2.c 
  zungrq.c zungtr.c zunm2l.c zunm2r.c zunmbr.c zunmhr.c zunml2.c 
  zunmlq.c zunmql.c zunmqr.c zunmr2.c zunmr3.c zunmrq.c zunmrz.c 
  zunmtr.c zupgtr.c 
  zupmtr.c izmax1.c dzsum1.c zstemr.c) 

add_library(liblapack OBJECT
   ${DLASRC_MINIMAL}
   ${ZLASRC_MINIMAL}
   ${DZLAUX}
   ${ALLAUX}
)

