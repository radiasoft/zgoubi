INC=../include

OBJS=\
 libzg.a(date2.o)\
 libzg.a(empty.o)\
 libzg.a(endjob.o)\
 libzg.a(gttext.o)\
 libzg.a(iraz.o)\
 libzg.a(raz.o)\
 libzg.a(rffag.o)\
 libzg.a(rhist.o)\
 libzg.a(rndm.o)\
 libzg.a(strcon.o)\
 libzg.a(strget.o)\
 libzg.a(time2.o)\
 libzg.a(vscal.o)\
 libzg.a(vvect.o)


.POSIX:

.f.o:
	$(FC) $(FFLAGS) -I$(INC) -c $<

all : $(OBJS) 

clean:
	$(RM) libzg.a *.o *~

libzg.a(rffag.o) : $(INC)/MXLD.H
