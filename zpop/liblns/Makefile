LIBS     = -ldl -lcurses

# Ne rien modifier au delà de cette ligne
DIRLIBS  = -L.

FMINIGRAF = fminigraf.o tfalco.o
RPSLD   = bdrpsld.o inkey.o  rpsld.o  scrloc.o vt200.o ifmes.o wrstat.o

.POSIX:

all : librpsld.a libfminigraf.a

docs :
	cd doc ; $(MAKE)

ttraxe  : ttraxe.o libfminigraf.a librpsld.a
	$(FC) $(FFLAGS) $(LDFLAGS) ttraxe.o \
        -o $@ libfminigraf.a librpsld.a \
        $(DIRLIBS) $(LIBS)

sttraxe : ttraxe.o
	$(F77) ttraxe.o -o $@   \
        $(DIRLIBS) -lfminigraf -lrpsld $(LIBS)

librpsld.a : $(RPSLD)
	ar rv $@ $?

libfminigraf.a : $(FMINIGRAF)
	ar rv $@ $?

clean :
	rm -f *.o *.a ttraxe ftraxe sftraxe sttraxe *~
	cd doc ; $(MAKE) clean
