# libdet/makefile

.SUFFIXES:
.PHONY: clean

CC = gfortran
CFLAGS = -W
OBJ = getdet.o 

all: $(OBJ)

%o: %f
	$(CC) -c $^ -o $@ $(CFLAGS)

clean:
	rm -rf *.o *~ *#