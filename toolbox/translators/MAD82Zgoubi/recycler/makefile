#
CFLAGS 	=  
FFLAGS	= 
MODULEC = 
MODULEF = \
	/home/room1/meot/mad/tools/trmadzgoubi/recycler/madzg.f
# obtient la liste des .o a partir des .f :
OBJECTF = $(MODULEF:.f=.o)
OBJECTC = $(MODULEC:.c=.o)
OBJECTS = $(OBJECTF) $(OBJECTC)
madzg  :$(OBJECTS)
	f77  $(OBJECTS) -o $@
