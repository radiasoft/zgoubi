LIBS     = -ldl -lcurses
LDSHFLAGS= -s -shared
SHCFLAGS =  -O -fPIC -DLINUX
CFLAGS   = -g -DLINUX
F77      = g77
FFLAGS	 = -g
SHFFLAGS =  -O -fPIC -DLINUX
SL       = so
ENDO     =

# Ne rien modifier au delà de cette ligne
DIRLIBS  = -L.
INCLUDES =

.PHONY : all doc

.c.o:
	$(CC)  $(CFLAGS) $(INCLUDES) -c $<
	cd sh ; $(CC)  $(SHCFLAGS) $(INCLUDES) -c ../$<
.f.o:
	$(F77)  $(FFLAGS) -c $<
	cd sh ; $(F77) $(SHFFLAGS) -c ../$<


FMINIGRAF = fminigraf.o tfalco.o
RPSLD   = bdrpsld.o inkey.o  rpsld.o  scrloc.o vt200.o ifmes.o wrstat.o
SFMINIGRAF = sh/fminigraf.o sh/tfalco.o
SRPSLD   = sh/bdrpsld.o sh/inkey.o sh/rpsld.o sh/scrloc.o sh/vt200.o \
           sh/ifmes.o sh/wrstat.o

all : librpsld.a libfminigraf.a

doc :
	$(MAKE) -C doc

ttraxe  : ttraxe.o libfminigraf.a librpsld.a
	$(F77)   ttraxe.o $(ENDO) \
        -o $@ libfminigraf.a librpsld.a \
        $(DIRLIBS) $(LIBS)

sttraxe : ttraxe.o
	$(F77) ttraxe.o $(ENDO) -o $@   \
        $(DIRLIBS) -lfminigraf -lrpsld $(LIBS)

librpsld.a : $(RPSLD)
	ar rv $@ $?

libfminigraf.a : $(FMINIGRAF)
	ar rv $@ $?

librpsld.$(SL) : $(SRPSLD)
	$(F77) -shared $(SRPSLD) -o $@

libfminigraf.$(SL) : $(SFMINIGRAF)
	$(F77) -shared $(SFMINIGRAF) -o $@

clean :
	rm -f *.o *.a *.$(SL) sh/*.o ttraxe ftraxe sftraxe sttraxe *~
	$(MAKE) -C doc clean

install:
	cp *.a *.$(SL) /usr/local/lib
