# Main Makefile
.SUFFIXES:
.PHONY:

LD = gfortran
LDFLAGS = -W
OBJ = libraries/libeis/*.o libraries/libinv/*.o libraries/libdet/*.o tunesc/*.o 

all:
	cd tunesc && make
	cd libraries/libeis && make
	cd libraries/libinv && make
	cd libraries/libdet && make
	$(LD) $(OBJ) -o ETparam

clean:
	cd tunesc && make clean
	cd libraries/libeis && make clean
	cd libraries/libinv && make clean
	cd libraries/libdet && make clean
	rm -rf *.o *~ *#
	rm -rf ETparam