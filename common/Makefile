INC=../include

OBJS=\
 libzg.a(date2.o)\
 libzg.a(empty.o)\
 libzg.a(gttext.o)\
 libzg.a(iraz.o)\
 libzg.a(okkp.o)\
 libzg.a(raz.o)\
 libzg.a(rffag.o)\
 libzg.a(rhist.o)\
 libzg.a(rndm.o)\
 libzg.a(strcon.o)\
 libzg.a(strget.o)\
 libzg.a(time2.o)\
 libzg.a(vscal.o)


.POSIX:

.f.o:
	$(FC) $(FFLAGS) -I$(INC) -c $<

all : $(OBJS)

clean:
	$(RM) libzg.a *.o *~

libzg.a(okkp.o) : $(INC)/OKKP.H
libzg.a(rffag.o) : $(INC)/MXLD.H
