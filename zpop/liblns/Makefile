# Ne rien modifier au delà de cette ligne
DIRLIBS  = -L.

FMINIGRAF = fminigraf.o

.POSIX:

all : libfminigraf.a

docs :
	cd doc ; $(MAKE)

libfminigraf.a : $(FMINIGRAF)
	ar rv $@ $?

clean :
	rm -f *.o *.a *~
	cd doc ; $(MAKE) clean
