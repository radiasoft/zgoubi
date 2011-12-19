INC=../include

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
 graph7.f\
 graph8.f\
 header.f\
 help.f\
 hermit.f\
 hmatch.f\
 hmtcgr.f\
 homclr.f\
 idluni.f\
 idlg.f\
 impv.f\
 inigr.f\
 inpech.f\
 inpvar.f\
 inrang.f\
 insy.f\
 isofld.f\
 lastp.f\
 lintyp.f\
 lips.f\
 lpscnt.f\
 lpsdrw.f\
 lpsfit.f\
 match.f\
 menvcf.f\
 minmax.f\
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
 rmatrx.f\
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
 square.f\
 srax.f\
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
 unitr.f\
 vectpl.f\
 xyzbr.f\
 zgoubi_plot.f

# obtient la liste des .o a partir des .f :
OBJECTF = $(MODULEF:.f=.o)
OBJECTS = $(OBJECTF)
LNSLIBS = liblns/libfminigraf.a

.POSIX:

.f.o:
	$(FC) $(FFLAGS) -I$(INC) -c $<

zpop  : $(OBJECTS) build_liblns
	$(FC) $(LDFLAGS) $(OBJECTS) $(LNSLIBS) -L../common -lzg -lm -o $@

build_liblns :
	cd liblns ; $(MAKE)

clean :
	cd liblns ; $(MAKE) clean
	$(RM) *.o zpop *~

bin.o: MXVAR.H
block.o: $(INC)/MAXTRA.H
caltra.o: $(INC)/MXLD.H $(INC)/MXFS.H $(INC)/LSTKEY.H
filcoo.o: MAXNTR.H
graph7.o: MXVAR.H MAXNTR.H
graph8.o: MXVAR.H
inigr.o: MXVAR.H MAXNTR.H
inpvar.o: MXVAR.H
isofld.o: MXVAR.H
lastp.o: MXVAR.H
lips.o: MAXNTR.H FILFAI.H
lpscnt.o: MAXNTR.H
lpsdrw.o: MXVAR.H MAXNTR.H
lpsfit.o: MAXNTR.H $(INC)/MAXCOO.H
okkl.o: OKKL.H
okkt.o: MAXNTR.H $(INC)/OKKT.H
opnmn.o: $(INC)/FILPLT.H FILFAI.H
opnmnl.o: $(INC)/FILPLT.H FILFAI.H
pllips.o: MXVAR.H MAXNTR.H
ploter.o: MXVAR.H $(INC)/MAXTRA.H
pltemi.o: MXVAR.H $(INC)/MXLD.H FILFAI.H
plthis.o: MXVAR.H
pltlps.o: $(INC)/MAXTRA.H MXVAR.H
pltopt.o: MXVAR.H
prdata.o: $(INC)/MXLD.H
rautor.o: $(INC)/MXLD.H
rbend.o: $(INC)/MXLD.H
rbinar.o: $(INC)/MXLD.H
rcarte.o: $(INC)/MXLD.H
rcavit.o: $(INC)/MXLD.H
rcible.o: $(INC)/MXLD.H
rcolli.o: $(INC)/MXLD.H
rcsr.o: $(INC)/MXLD.H
rdips.o: $(INC)/MXLD.H
readco.o: MXVAR.H $(INC)/MXLD.H MAXNTR.H $(INC)/FRMFAI.H $(INC)/FRMPLT.H $(INC)/FRMSPN.H
rebmul.o: $(INC)/MXLD.H
rel2tu.o: $(INC)/MXLD.H
relcmi.o: $(INC)/MXLD.H
relcyl.o: $(INC)/MXLD.H
relmir.o: $(INC)/MXLD.H
rffags.o: $(INC)/MXLD.H
rfit.o: $(INC)/MXLD.H
rfitgt.o: $(INC)/MXLD.H
rgasca.o: $(INC)/MXLD.H
rhelix.o: $(INC)/MXLD.H
rmcdes.o: $(INC)/MXLD.H
rmatrx.o: $(INC)/MXLD.H
rmcobj.o: $(INC)/MXLD.H
rnorma.o: $(INC)/MXLD.H
robjet.o: $(INC)/MXLD.H
robjta.o: $(INC)/MXLD.H
rscale.o: $(INC)/MXLD.H $(INC)/MXFS.H
rsepar.o: $(INC)/MXLD.H
rsimb.o: $(INC)/MXLD.H
rsolen.o: $(INC)/MXLD.H
rspn.o: $(INC)/MXLD.H $(INC)/MAXTRA.H
rsrlos.o: $(INC)/MXLD.H
rtrans.o: $(INC)/MXLD.H
rundul.o: $(INC)/MXLD.H
runipo.o: $(INC)/MXLD.H
rwienf.o: $(INC)/MXLD.H
spctra.o: MAXNTR.H FILFAI.H
speana.o: MAXNTR.H
speimp.o: MAXNTR.H
srdwc.o: $(INC)/MXSTEP.H
srdw.o: $(INC)/MXSTEP.H
srdwm.o: $(INC)/FILPLT.H
srdwpl.o: MXVAR.H
sref.o: $(INC)/MXSTEP.H MXVAR.H
srefm.o: $(INC)/FILPLT.H
srinc.o: $(INC)/MXSTEP.H
srinm.o: MXVAR.H $(INC)/FILPLT.H
srmen.o: $(INC)/FILPLT.H
storco.o: MAXNTR.H MXVAR.H
storig.o: $(INC)/MXLD.H
synaxe.o: PARCST.H
synbox.o: PARCST.H
xyzbr.o: MXVAR.H $(INC)/MXSTEP.H
zpop_diagQ.o: LMKDIM.H
