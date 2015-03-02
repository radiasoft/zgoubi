# tunesc/Makefile
.SUFFIXES:
.PHONY: clean

CC = gfortran
CFLAGS = -W
OBJ = debstr.o extrac.o finstr.o hamilt.o positi.o propag.o strcon.o

all: $(OBJ)

%o: %f
	$(CC) -c $^ -o $@ $(CFLAGS)

clean:
	rm -rf *.o *~ *#