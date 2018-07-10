ARFLAGS=rvU
DIRLIBS  = -L.

OBJS = \
 libsvd.a(DISNAN.o)\
 libsvd.a(DLAISNAN.o)\
 libsvd.a(DLANGE.o)\
 libsvd.a(DLAPY2.o)\
 libsvd.a(DLAMCH.o)\
 libsvd.a(DLAMC3.o)\
 libsvd.a(IEEECK.o)\
 libsvd.a(ILADLC.o)\
 libsvd.a(ILADLR.o)\
 libsvd.a(ILAENV.o)\
 libsvd.a(IPARMQ.o)\
 libsvd.a(LSAME.o)\
 libsvd.a(DNRM2.o)\
 libsvd.a(SVD.o)\
 libsvd.a(DBDSQR.o)\
 libsvd.a(DGEBD2.o)\
 libsvd.a(DGEBRD.o)\
 libsvd.a(DGELQ2.o)\
 libsvd.a(DGELQF.o)\
 libsvd.a(DGEQR2.o)\
 libsvd.a(DGEQRF.o)\
 libsvd.a(DGESVD.o)\
 libsvd.a(DLABRD.o)\
 libsvd.a(DLACPY.o)\
 libsvd.a(DLARF.o)\
 libsvd.a(DLARFB.o)\
 libsvd.a(DLARFG.o)\
 libsvd.a(DLARFT.o)\
 libsvd.a(DLARTG.o)\
 libsvd.a(DLAS2.o)\
 libsvd.a(DLASCL.o)\
 libsvd.a(DLASET.o)\
 libsvd.a(DLASQ1.o)\
 libsvd.a(DLASQ2.o)\
 libsvd.a(DLASQ3.o)\
 libsvd.a(DLASQ4.o)\
 libsvd.a(DLASQ5.o)\
 libsvd.a(DLASQ6.o)\
 libsvd.a(DLASR.o)\
 libsvd.a(DLASRT.o)\
 libsvd.a(DLASSQ.o)\
 libsvd.a(DLASV2.o)\
 libsvd.a(DORG2R.o)\
 libsvd.a(DORGBR.o)\
 libsvd.a(DORGL2.o)\
 libsvd.a(DORGLQ.o)\
 libsvd.a(DORGQR.o)\
 libsvd.a(DORM2R.o)\
 libsvd.a(DORMBR.o)\
 libsvd.a(DORML2.o)\
 libsvd.a(DORMLQ.o)\
 libsvd.a(DORMQR.o)\
 libsvd.a(XERBLA.o)\
 libsvd.a(DROT.o)\
 libsvd.a(DSCAL.o)\
 libsvd.a(DGEMM.o)\
 libsvd.a(DGEMV.o)\
 libsvd.a(DGER.o)\
 libsvd.a(DCOPY.o)\
 libsvd.a(DSWAP.o)\
 libsvd.a(DTRMM.o)\
 libsvd.a(DTRMV.o)

.POSIX:

all : libsvd.a

libsvd.a : $(OBJS)
	ar $(ARFLAGS) $@ $?

clean :
	rm -rf *.o *.a *~
