.PHONY : all clean doc

all :
	$(MAKE) -C common
	$(MAKE) -C zgoubi
	$(MAKE) -C zpop

clean :
	$(RM) *~
	$(MAKE) -C common clean
	$(MAKE) -C zgoubi clean
	$(MAKE) -C zpop clean
	$(MAKE) -C doc clean

doc :
	$(MAKE) -C doc
	$(MAKE) -C zpop/liblns/doc
