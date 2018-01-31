# E1 for obtaining non-standard features, e.g. the DATE routine
# -g for debbug, -O for optimizer, +O3 for higher optim, 
#    -K for statiq mem alloc.
CFLAGS 	=  
FFLAGS  =  -ansi -g -fno-automatic -finit-local-zero
MODULEC = 
MODULEF = \
	madzg.f
# obtient la liste des .o a partir des .f :
OBJECTF = $(MODULEF:.f=.o)
OBJECTC = $(MODULEC:.c=.o)
OBJECTS = $(OBJECTF) $(OBJECTC)
madzg  :$(OBJECTS)
	g77  $(OBJECTS) -lm -o $@
