CFLAGS=-g -DGFORTRAN4
FC=gfortran
FFLAGS=-O -Wall -fno-automatic

.POSIX:

all :
	cd common ; $(MAKE) CFLAGS="$(CFLAGS)" FC="$(FC)" FFLAGS="$(FFLAGS)"
	cd zgoubi ; $(MAKE) CFLAGS="$(CFLAGS)" FC="$(FC)" FFLAGS="$(FFLAGS)"
	cd zpop ; $(MAKE) CFLAGS="$(CFLAGS)" FC="$(FC)" FFLAGS="$(FFLAGS)"

clean :
	$(RM) *~
	cd common ; $(MAKE) clean
	cd zgoubi ; $(MAKE) clean
	cd zpop ; $(MAKE) clean
	cd doc ; $(MAKE) clean

docs :
	cd doc ; $(MAKE)
	cd zpop/liblns/doc ; $(MAKE)
