INC=../include

OBJS=\
 libzg.a(date2.o)\
 libzg.a(empty.o)\
 libzg.a(iraz.o)\
 libzg.a(okkp.o)\
 libzg.a(raz.o)\
 libzg.a(rffags.o)\
 libzg.a(rhist.o)\
 libzg.a(rndm.o)\
 libzg.a(strcon.o)\
 libzg.a(strget.o)\
 libzg.a(time2.o)


.POSIX:

.f.o:
	$(FC) $(FFLAGS) -I$(INC) -c $<

all : $(OBJS)

clean:
	$(RM) libzg.a *.o *~

libsg.a(okkp.o) : $(INC)/OKKP.H
libsg.a(rffags.o) : $(INC)/MXLD.H
