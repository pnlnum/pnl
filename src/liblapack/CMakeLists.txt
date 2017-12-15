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
  dgbsv.c  
  dgbtf2.c dgbtrf.c dgbtrs.c dgebak.c dgebal.c
  dgecon.c dgeev.c 
  dgehd2.c dgehrd.c
  dgelsy.c 
  dgeqp3.c dgeqpf.c dgeqr2.c dgeqrf.c dgerfs.c dgerq2.c dgerqf.c 
  dgesc2.c dgetc2.c dgetf2.c 
  dgetrf.c dgetri.c 
  dgetrs.c dggbak.c dggbal.c 
  dgghrd.c
  dggrqf.c dgtsv.c  
  dgttrf.c dgttrs.c dgtts2.c dhgeqz.c 
  dhseqr.c dlabrd.c dlacn2.c 
  dlaein.c dlaexc.c dlag2.c  dlags2.c dlagtm.c dlagv2.c dlahqr.c 
  dlahr2.c dlaic1.c dlaln2.c dlals0.c dlalsa.c dlalsd.c 
  dlangb.c dlange.c dlangt.c dlanhs.c 
  dlansy.c dlantb.c dlantr.c dlanv2.c 
  dlapll.c dlapmt.c 
  dlaqgb.c dlaqge.c dlaqp2.c dlaqps.c 
  dlaqr0.c dlaqr1.c dlaqr2.c dlaqr3.c dlaqr4.c dlaqr5.c
  dlaqtr.c dlar1v.c dlar2v.c iladlr.c iladlc.c 
  dlarf.c  dlarfb.c dlarfg.c dlarft.c dlarfx.c dlargv.c 
  dlarrv.c dlartv.c dlarfp.c 
  dlarz.c  dlarzb.c dlarzt.c dlaswp.c dlasy2.c dlasyf.c 
  dlatbs.c dlatdf.c dlatrd.c dlatrs.c dlatrz.c dlatzm.c dlauu2.c 
  dlauum.c dorg2l.c dorg2r.c 
  dorgbr.c dorghr.c dorgl2.c dorglq.c dorgql.c dorgqr.c dorgr2.c 
  dorgtr.c dorm2r.c 
  dormbr.c dormhr.c  dorml2.c dormlq.c dormqr.c dormr2.c 
  dormr3.c dormrq.c dormrz.c
  dpbtf2.c dpbtrf.c dpbtrs.c  
  dpotf2.c dpotrf.c dpotri.c dpotrs.c dpstrf.c dpstf2.c 
  drscl.c  
  dsbtrd.c 
  dsptrd.c 
  dsptrf.c dsptrs.c dstein.c 
  dsyev.c
  dsytd2.c dsytf2.c dsytrd.c dsytrf.c dsytri.c dsytrs.c 
  dtgevc.c dtgex2.c dtgexc.c dtgsen.c 
  dtgsja.c dtgsna.c dtgsy2.c dtgsyl.c dtptri.c 
  dtrevc.c dtrexc.c dtrsen.c dtrsna.c dtrsyl.c 
  dtrti2.c dtrtri.c dtrtrs.c dtzrqf.c dtzrzf.c dstemr.c 
  )

set (ZLASRC_MINIMAL  
  zgebak.c zgebal.c 
  zgeev.c 
  zgehd2.c zgehrd.c 
  zgetrf.c zgetf2.c 
  zgetri.c zgetrs.c 
  zheev.c 
  zhetd2.c 
  zhetrd.c 
  zhseqr.c 
  zlacgv.c zlacn2.c zlacpy.c zladiv.c 
  zlahef.c zlahqr.c 
  zlahr2.c zlange.c 
  zlanhe.c 
  zlanhs.c 
  zlaqr0.c zlaqr1.c zlaqr2.c zlaqr3.c zlaqr4.c zlaqr5.c 
  zlar2v.c ilazlr.c ilazlc.c 
  zlarf.c  zlarfb.c 
  zlarfg.c zlarft.c zlarfp.c 
  zlargv.c zlartg.c zlartv.c 
  zlascl.c zlaset.c zlasr.c  
  zlassq.c zlaswp.c  
  zlatrd.c zlatrs.c zlauu2.c 
  zlauum.c
  zpotf2.c zpotrf.c zpotri.c zpotrs.c zpstrf.c zpstf2.c 
  zrot.c
  zdrscl.c 
  zsteqr.c 
  ztrevc.c ztrexc.c 
  ztrsyl.c ztrti2.c ztrtri.c ztrtrs.c zung2l.c 
  zung2r.c zunghr.c zungql.c zungqr.c 
  zungtr.c zunm2r.c zunmhr.c zunml2.c 
  zunmlq.c zunmqr.c
  izmax1.c dzsum1.c
  ) 

add_library(liblapack OBJECT
   ${DLASRC_MINIMAL}
   ${ZLASRC_MINIMAL}
   ${DZLAUX}
   ${ALLAUX}
)

