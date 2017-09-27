CFLAGS=-g -DGFORTRAN4
FC=gfortran
#FFLAGS=-O -g -fno-automatic -pg
FFLAGS=-O9 -Wall -fno-automatic -pedantic -cpp
#FFLAGS=-O9 -Wall -fno-automatic -pedantic -mtune=native
#FFLAGS=-O4 -Wall -fno-automatic -pedantic -fbacktrace -cpp
#FFLAGS=-Og -O4 -Wall -fno-automatic -pedantic -cpp
#FFLAGS=-O4 -Wall                -pedantic

.POSIX:

all :
	cd modules ;  make CFLAGS="$(CFLAGS)" FC="$(FC)" FFLAGS="$(FFLAGS)"
	cd common ; rm libzg.a || true  ; make CFLAGS="$(CFLAGS)" FC="$(FC)" FFLAGS="$(FFLAGS)"
	cd zgoubi ; make CFLAGS="$(CFLAGS)" FC="$(FC)" FFLAGS="$(FFLAGS)"
	cd zpop ; make CFLAGS="$(CFLAGS)" FC="$(FC)" FFLAGS="$(FFLAGS)"
#	cd guide ; make

clean :
	$(RM) *~
	cd modules ;  make clean
	cd common ;  make clean
	cd zgoubi/coupling ; make clean
	cd zgoubi ;  make clean
	cd zpop/liblns ; $(MAKE) clean
	cd zpop ; $(MAKE) clean
#	cd guide ; $(MAKE) clean



# Comments : 
# -fno-automatic      Treat each program unit (except those marked as RECURSIVE) as if the SAVE statement were specified for every local variable and array referenced in it. Does not affect common blocks. (Some Fortran compilers provide this option under the name -static or -save.) The default, which is -fautomatic, uses the stack for local variables smaller than the value given by -fmax-stack-var-size. 
