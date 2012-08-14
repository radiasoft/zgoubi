
.POSIX:

all :
	$(MAKE) -f Makefile_zgoubi
	$(MAKE) -f Makefile_zpop

clean :
	$(RM) *~
	cd common ; $(MAKE) clean
	cd zgoubi ; $(MAKE) clean
	cd zpop ; $(MAKE) clean
	cd doc ; $(MAKE) clean

docs :
	cd doc ; $(MAKE)
	cd zpop/liblns/doc ; $(MAKE)
