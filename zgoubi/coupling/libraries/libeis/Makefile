# libeis/makefile

.SUFFIXES:
.PHONY: clean

CC = gfortran
CFLAGS = -W
OBJ = balanc.o balbak.o cdiv.o elmhes.o eltran.o hqr2.o hqr.o rg.o

all: $(OBJ)

%o: %f
	$(CC) -c $^ -o $@ $(CFLAGS)

clean:
	rm -rf *.o *~ *#