FC=g77
INCLUDE=-I../include
FFLAGS=$(INCLUDE) -O -Wall -fno-automatic -finit-local-zero 

.PHONY : all clean

OBJS=\
 libzg.a(date2.o)\
 libzg.a(empty.o)\
 libzg.a(iraz.o)\
 libzg.a(okkp.o)\
 libzg.a(raz.o)\
 libzg.a(rffag.o)\
 libzg.a(rhist.o)\
 libzg.a(rndm.o)\
 libzg.a(strcon.o)\
 libzg.a(strget.o)\
 libzg.a(time2.o)


all : $(OBJS)

clean:
	$(RM) libzg.a *.o *~

libsg.a(okkp.o) : ../include/OKKP.H
libsg.a(rffag.o) : ../include/MXLD.H
