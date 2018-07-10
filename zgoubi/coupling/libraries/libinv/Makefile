# libinv/makefile

.SUFFIXES:
.PHONY: clean

CC = gfortran
CFLAGS = -W
OBJ =\
cher2.o  csyr2k.o  dcopy.o   dlaswp.o  dsymm.o   dtrsv.o   izamax.o   sgemm.o   sswap.o   strsm.o   zgeru.o   zswap.o\
caxpy.o  cher2k.o  csyrk.o   ddot.o    dnrm2.o   dsymv.o   dtrti2.o   lsame.o   sgemv.o   ssymm.o   strsv.o   zhbmv.o\
ccopy.o  cher.o    ctbmv.o   drot.o    dsyr2.o   dtrtri.o  sger.o     ssymv.o   xerbla.o  zhemm.o   zsyr2k.o  zsymm.o\
cdotc.o  cherk.o   ctbsv.o   drotg.o   dsyr2k.o  dzasum.o  snrm2.o    ssyr2.o   zaxpy.o   zhemv.o   zsyrk.o   ztrsv.o\
cdotu.o  chpmv.o   ctpmv.o   drotm.o   dsyr.o    dznrm2.o  sasum.o    srot.o    ssyr2k.o  zcopy.o   zher2.o   ztbmv.o\
cgbmv.o  chpr2.o   ctpsv.o   dgbmv.o   drotmg.o  dsyrk.o   icamax.o   saxpy.o   srotg.o   ssyr.o    zdotc.o   zher2k.o  ztbsv.o\
cgemm.o  chpr.o    ctrmm.o   dgemm.o   dsbmv.o   dtbmv.o   idamax.o   scabs1.o  srotm.o   ssyrk.o   zdotu.o   zher.o    ztpmv.o\
cgemv.o  crotg.o   ctrmv.o   dgemv.o   dscal.o   dtbsv.o   ieeeck.o   scasum.o  srotmg.o  stbmv.o   zdrot.o   zherk.o   ztpsv.o\
cgerc.o  cscal.o   ctrsm.o   dger.o    dsdot.o   dtpmv.o   ilaenv.o   scnrm2.o  ssbmv.o   stbsv.o   zdscal.o  zhpmv.o   ztrmm.o\
cgeru.o  csrot.o   ctrsv.o   dgetf2.o  dspmv.o   dtpsv.o   scopy.o    sscal.o   stpmv.o   zgbmv.o   zhpr2.o   ztrmv.o\
chbmv.o  csscal.o  dasum.o   dgetrf.o  dspr2.o   dtrmm.o   sdot.o     sspmv.o   stpsv.o   zgemm.o   zhpr.o    ztrsm.o\
chemm.o  cswap.o   daxpy.o   dgetri.o  dspr.o    dtrmv.o   iparmq.o   sdsdot.o  sspr2.o   strmm.o   zgemv.o   zrotg.o\
chemv.o  csymm.o   dcabs1.o  dlamch.o  dswap.o   dtrsm.o   isamax.o   sgbmv.o   sspr.o    strmv.o   zgerc.o   zscal.o\

all: $(OBJ)

%o: %f
	$(CC) -c $^ -o $@ $(CFLAGS)

clean:
	rm -rf *.o *~ *#