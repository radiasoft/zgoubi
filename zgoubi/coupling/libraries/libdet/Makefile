# libdet/makefile

.SUFFIXES:
.PHONY: clean

CC = gfortran
CFLAGS = -W
OBJ = caxpy.o  cgeco.o  cgedi.o  cgefa.o  cscal.o  cswap.o getdet.o 

all: $(OBJ)

%o: %f
	$(CC) -c $^ -o $@ $(CFLAGS)

clean:
	rm -rf *.o *~ *#