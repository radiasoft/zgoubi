INCLUDE = -I../include
FFLAGS	= $(INCLUDE) -O -fno-automatic -finit-local-zero
LDFLAGS=-L../common -lzg

MODULEF =\
 amini.f\
 asinh.f\
 begvcf.f\
 binarf.f\
 bin.f\
 block.f\
 caltra.f\
 caltri.f\
 class.f\
 clorbi.f\
 clscr.f\
 debstr.f\
 dfring.f\
 doree.f\
 drawspi.f\
 drctry.f\
 dsherm.f\
 endvcf.f\
 erf.f\
 filcoo.f\
 finstr.f\
 flush2.f\
 fmag.f\
 fmatch.f\
 fmtcgr.f\
 fringe.f\
 funk.f\
 gk1.f\
 grad.f\
 graph1.f\
 graph2.f\
 header.f\
 help.f\
 hermit.f\
 hmatch.f\
 hmtcgr.f\
 homclr.f\
 idluni.f\
 impv.f\
 inigr.f\
 inpech.f\
 inpvar.f\
 inrang.f\
 insy.f\
 isofld.f\
 libspa.f\
 lintyp.f\
 lips.f\
 lpscnt.f\
 lpsdrw.f\
 lpsfit.f\
 match.f\
 menvcf.f\
 minmax.f\
 mnfich.f\
 mnfile.f\
 mnzgra.f\
 okkl.f\
 okkt.f\
 opndef.f\
 opngnu.f\
 opnmn.f\
 opnmnl.f\
 opnwrn.f\
 plhist.f\
 plinit.f\
 pllips.f\
 ploter.f\
 pltemi.f\
 plthis.f\
 pltlps.f\
 pltopt.f\
 ppar.f\
 prarc.f\
 prdata.f\
 prdip.f\
 prffag.f\
 prtime.f\
 raiman.f\
 range.f\
 rautor.f\
 ravorb.f\
 razi.f\
 rbend.f\
 rbinar.f\
 rcarte.f\
 rcavit.f\
 rcible.f\
 rcoils.f\
 rcolli.f\
 rcsr.f\
 rdips.f\
 readco.f\
 rebmul.f\
 rel2tu.f\
 relcmi.f\
 relcyl.f\
 relmir.f\
 rewin2.f\
 rffags.f\
 rfit.f\
 rfitgt.f\
 rfmtch.f\
 rgasca.f\
 rhelix.f\
 rhmtch.f\
 rmcdes.f\
 rmcobj.f\
 rmulti.f\
 rnorma.f\
 robjet.f\
 robjta.f\
 rqsod.f\
 rscale.f\
 rsepar.f\
 rsimb.f\
 rsolen.f\
 rspn.f\
 rsrlos.f\
 rtrans.f\
 rundul.f\
 runipo.f\
 rwienf.f\
 savplt.f\
 shermi.f\
 spctra.f\
 speana.f\
 spedrw.f\
 spegr.f\
 speimp.f\
 spepr.f\
 splcli.f\
 square.f\
 srax.f\
 srcul.f\
 srdnu.f\
 srdwc.f\
 srdw.f\
 srdwm.f\
 srdwpl.f\
 srdwsc.f\
 srdwst.f\
 sref.f\
 srefm.f\
 srinc.f\
 srinm.f\
 srmdl1.f\
 srmen.f\
 srmodl.f\
 srpli.f\
 srplt.f\
 srplv.f\
 storco.f\
 storig.f\
 superp.f\
 synaxe.f\
 synbox.f\
 tracsy.f\
 traxes.f\
 trkvar.f\
 trtxt.f\
 tylab.f\
 typtra.f\
 vectpl.f\
 xyzbr.f\
 zgoubi_plot.f

# obtient la liste des .o a partir des .f :
OBJECTF = $(MODULEF:.f=.o)
OBJECTS = $(OBJECTF)

.PHONY : clean

zpop  : $(OBJECTS)
	$(MAKE) -C liblns
	$(FC) $(LDFLAGS) $(OBJECTS) liblns/libfminigraf.a liblns/librpsld.a -L../common -lzg -lm -lcurses -o $@

clean :
	$(MAKE) -C liblns clean
	$(RM) *.o zpop *~

bin.o: MXVAR.H
block.o: ../include/MAXTRA.H
caltra.o: ../include/MXLD.H ../include/MXFS.H ../include/LSTKEY.H
filcoo.o: MAXNTR.H
graph1.o: MXVAR.H MAXNTR.H
graph2.o: MXVAR.H
inigr.o: MXVAR.H MAXNTR.H
inpvar.o: MXVAR.H
isofld.o: MXVAR.H
lips.o: MAXNTR.H FILFAI.H
lpscnt.o: MAXNTR.H
lpsdrw.o: MXVAR.H MAXNTR.H
lpsfit.o: MAXNTR.H ../include/MAXCOO.H
okkl.o: OKKL.H
okkt.o: MAXNTR.H ../include/OKKT.H
opnmn.o: ../include/FILPLT.H FILFAI.H
opnmnl.o: ../include/FILPLT.H FILFAI.H
pllips.o: MXVAR.H MAXNTR.H
ploter.o: MXVAR.H ../include/MAXTRA.H
pltemi.o: MXVAR.H ../include/MXLD.H FILFAI.H
plthis.o: MXVAR.H
pltlps.o: ../include/MAXTRA.H MXVAR.H
pltopt.o: MXVAR.H
prdata.o: ../include/MXLD.H
rautor.o: ../include/MXLD.H
rbend.o: ../include/MXLD.H
rbinar.o: ../include/MXLD.H
rcarte.o: ../include/MXLD.H
rcavit.o: ../include/MXLD.H
rcible.o: ../include/MXLD.H
rcolli.o: ../include/MXLD.H
rcsr.o: ../include/MXLD.H
rdips.o: ../include/MXLD.H
readco.o: MXVAR.H ../include/MXLD.H MAXNTR.H ../include/FRMFAI.H ../include/FRMPLT.H
rebmul.o: ../include/MXLD.H
rel2tu.o: ../include/MXLD.H
relcmi.o: ../include/MXLD.H
relcyl.o: ../include/MXLD.H
relmir.o: ../include/MXLD.H
rffags.o: ../include/MXLD.H
rfit.o: ../include/MXLD.H
rfitgt.o: ../include/MXLD.H
rgasca.o: ../include/MXLD.H
rhelix.o: ../include/MXLD.H
rmcdes.o: ../include/MXLD.H
rmcobj.o: ../include/MXLD.H
rnorma.o: ../include/MXLD.H
robjet.o: ../include/MXLD.H
robjta.o: ../include/MXLD.H
rscale.o: ../include/MXLD.H ../include/MXFS.H
rsepar.o: ../include/MXLD.H
rsimb.o: ../include/MXLD.H
rsolen.o: ../include/MXLD.H
rspn.o: ../include/MXLD.H ../include/MAXTRA.H
rsrlos.o: ../include/MXLD.H
rtrans.o: ../include/MXLD.H
rundul.o: ../include/MXLD.H
runipo.o: ../include/MXLD.H
rwienf.o: ../include/MXLD.H
spctra.o: MAXNTR.H FILFAI.H
speana.o: MAXNTR.H
speimp.o: MAXNTR.H
srdwc.o: ../include/MXSTEP.H
srdw.o: ../include/MXSTEP.H
srdwm.o: ../include/FILPLT.H
srdwpl.o: MXVAR.H
sref.o: ../include/MXSTEP.H MXVAR.H
srefm.o: ../include/FILPLT.H
srinc.o: ../include/MXSTEP.H
srinm.o: MXVAR.H ../include/FILPLT.H
srmen.o: ../include/FILPLT.H
storco.o: MAXNTR.H MXVAR.H
storig.o: ../include/MXLD.H
synaxe.o: PARCST.H
synbox.o: PARCST.H
xyzbr.o: MXVAR.H ../include/MXSTEP.H
zpop_diagQ.o: LMKDIM.H
