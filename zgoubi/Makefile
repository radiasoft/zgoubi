FC=g77
INCLUDE=-I../include
FFLAGS=$(INCLUDE) -O -Wall -fno-automatic -finit-local-zero 
LDFLAGS=

.PHONY : clean

OBJ=\
 accen.o\
 accep.o\
 adpol.o\
 aimant.o\
 apdpo1.o\
 apdpo.o\
 apherf.o\
 autorf.o\
 avorpr.o\
 bar2xy.o\
 baxbxr.o\
 beaimp.o\
 beamat.o\
 bendf.o\
 bendi.o\
 binari.o\
 binary.o\
 block.o\
 bonul.o\
 boost.o\
 butee.o\
 bxrxyz.o\
 carla.o\
 carshi.o\
 cavite.o\
 chamb.o\
 chamc.o\
 chamk.o\
 charef.o\
 chmbre.o\
 chmbr.o\
 chrefe.o\
 chreg.o\
 chxc.o\
 chxp.o\
 cible.o\
 cinema.o\
 cntinl.o\
 cntpar.o\
 coeffs.o\
 cofin.o\
 coils.o\
 coilsf.o\
 collim.o\
 cptfct.o\
 csri.o\
 csrint.o\
 ctrlb.o\
 ctrlc.o\
 cubspl.o\
 dbdxyz.o\
 debstr.o\
 dejaca.o\
 delta.o\
 depla.o\
 desbin.o\
 devtra.o\
 dipi.o\
 dipolm.o\
 dipsi.o\
 dist2.o\
 dpkick.o\
 drtent.o\
 drvg.o\
 dstefb.o\
 eax3tu.o\
 eaxial.o\
 el2tub.o\
 el2tu.o\
 elcmid.o\
 elcmif.o\
 elcmii.o\
 elcyld.o\
 elcyl.o\
 ellip.o\
 elmirf.o\
 elmiri.o\
 end.o\
 endfit.o\
 endjob.o\
 enrgy.o\
 ererz.o\
 erexyz.o\
 erzxyz.o\
 esl.o\
 eslvir.o\
 event.o\
 fborne.o\
 ffagi.o\
 ff.o\
 ffgspi.o\
 ffspd.o\
 fgf.o\
 finstr.o\
 fitarr.o\
 fitest.o\
 fitexp.o\
 fitgtv.o\
 fitlav.o\
 fitlma.o\
 fitnu.o\
 fitset.o\
 fitsta.o\
 flush2.o\
 fmapw.o\
 focal1.o\
 focale.o\
 funcd.o\
 gammln.o\
 gang.o\
 gascat.o\
 gasesl.o\
 gasini.o\
 geneta.o\
 gg.o\
 header.o\
 helix.o\
 helixf.o\
 histab.o\
 hist.o\
 icentr.o\
 idlg.o\
 idluni.o\
 ifmes.o\
 impaju.o\
 impcsr.o\
 impctr.o\
 impdev.o\
 impfai.o\
 impplt.o\
 imptra.o\
 impvar.o\
 inicin.o\
 inidrt.o\
 initra.o\
 integr.o\
 interm.o\
 interv.o\
 intplf.o\
 itavan.o\
 iter.o\
 itsens.o\
 jcentr.o\
 krtax.o\
 kstop.o\
 lign.o\
 lingua.o\
 lnul.o\
 lpsfit.o\
 lstins.o\
 lstset.o\
 majtra.o\
 maplim.o\
 mat1.o\
 mat2.o\
 matimp.o\
 matric.o\
 mcdes.o\
 mcdesi.o\
 mcdesl.o\
 mcobj1.o\
 mcobj3.o\
 mcobj.o\
 mino1.o\
 minonm.o\
 mix.o\
 mksa.o\
 modord.o\
 multip.o\
 multpo.o\
 nmfonc.o\
 nmmin1.o\
 nmmin.o\
 normal.o\
 normb.o\
 obj3.o\
 obj5.o\
 obj6.o\
 obj8.o\
 objerr.o\
 objeta.o\
 objets.o\
 objtrk.o\
 okkt.o\
 open2.o\
 partic.o\
 path.o\
 pavel.o\
 pckup.o\
 pickup.o\
 pltdat.o\
 poidev.o\
 polmes.o\
 prdata.o\
 pscal.o\
 pvect.o\
 quasex.o\
 raiman.o\
 rautor.o\
 raysyn.o\
 razdrv.o\
 rbend.o\
 rbinar.o\
 rcarte.o\
 rcavit.o\
 rcible.o\
 rcoils.o\
 rcolli.o\
 rcsr.o\
 rdip.o\
 rdips.o\
 rdpois.o\
 rebel.o\
 rebmul.o\
 refer.o\
 rel2tu.o\
 relcmi.o\
 relcyl.o\
 relmir.o\
 rempli.o\
 reset.o\
 rfaist.o\
 rffags.o\
 rfit.o\
 rfitgt.o\
 rgasca.o\
 rhelix.o\
 rmcdes.o\
 rmcobj.o\
 rmulti.o\
 rnorma.o\
 robjet.o\
 robjta.o\
 rotx.o\
 rparti.o\
 rpckup.o\
 rqsod.o\
 rscal.o\
 rsepar.o\
 rsimb.o\
 rsolen.o\
 rspn.o\
 rsrlos.o\
 rtnewt.o\
 rtrans.o\
 rundul.o\
 runipo.o\
 rwienf.o\
 scal0.o\
 scaler.o\
 scalin.o\
 scumul.o\
 separa.o\
 serv2.o\
 serv32.o\
 serv34.o\
 serv4.o\
 serv5.o\
 settim.o\
 solenf.o\
 soleno.o\
 solv.o\
 splint.o\
 spn.o\
 spnprn.o\
 spnprt.o\
 spnrot.o\
 spntrk.o\
 sref.o\
 srloss.o\
 srprn.o\
 stpsiz.o\
 straco.o\
 symmed.o\
 sympl.o\
 synpar.o\
 syntrk.o\
 toroid.o\
 toscac.o\
 toscap.o\
 trace.o\
 track.o\
 transf.o\
 transm.o\
 traobj.o\
 trarot.o\
 trbeam.o\
 tridi.o\
 trobj.o\
 trrote.o\
 tunes.o\
 twiss.o\
 undulfa.o\
 undulf.o\
 undulia.o\
 unduli.o\
 unipo.o\
 unipot.o\
 valqv.o\
 vequal.o\
 vrble.o\
 vscal.o\
 vzero.o\
 waist.o\
 wedgki.o\
 wienfi.o\
 xnorm.o\
 xrotb.o\
 ymoiny.o\
 zerof.o\
 zgoubi.o\
 zgoubi_main.o

zgoubi : $(OBJ)
	$(FC) $(LDFLAGS) $(OBJ) -L../common -lzg -o $@

clean :
	$(RM) *.o zgoubi *~

accen.o: ../include/MXLD.H ../include/MAXCOO.H ../include/MAXTRA.H
accep.o: ../include/MXLD.H ../include/MAXCOO.H ../include/MAXTRA.H
aimant.o: ../include/MXLD.H ../include/MAXCOO.H
autorf.o: ../include/MXLD.H ../include/MAXCOO.H ../include/MAXTRA.H
bendi.o: ../include/MXLD.H
binary.o: ../include/MXLD.H
block.o: ../include/MAXCOO.H ../include/MAXTRA.H ../include/MXFS.H
bonul.o: ../include/MAXCOO.H ../include/MAXTRA.H
carla.o: PARIZ.H ../include/MAXTRA.H ../include/MXLD.H
carshi.o: PARIZ.H ../include/MXLD.H
cavite.o: ../include/MXLD.H ../include/MAXCOO.H ../include/MAXTRA.H ../include/MXFS.H
chamb.o: ../include/MAXTRA.H ../include/MXLD.H ../include/MAXCOO.H
chamc.o: ../include/MAXTRA.H
chamk.o: PARIZ.H ../include/MAXTRA.H
charef.o: ../include/MAXTRA.H ../include/MAXCOO.H
chmbre.o: ../include/MAXTRA.H
chmbr.o: ../include/MAXTRA.H ../include/MAXCOO.H
chrefe.o: ../include/MAXTRA.H ../include/MXSTEP.H CSR.H ../include/MXLD.H ../include/MAXCOO.H
chxc.o: PARIZ.H ../include/MAXTRA.H ../include/MXLD.H ../include/MXFS.H ../include/FILPLT.H
chxp.o: PARIZ.H ../include/MXLD.H ../include/MXFS.H ../include/FILPLT.H
cible.o: ../include/MAXCOO.H ../include/MAXTRA.H ../include/MXLD.H
cntinl.o: ../include/MAXCOO.H ../include/MAXTRA.H
coeffs.o: ../include/MXLD.H ../include/MAXCOO.H ../include/MAXTRA.H
cofin.o: ../include/MAXTRA.H ../include/MXSTEP.H CSR.H ../include/MAXCOO.H
coils.o: ../include/MXLD.H
collim.o: ../include/MAXTRA.H ../include/MXLD.H ../include/MAXCOO.H
csri.o: ../include/MXSTEP.H ../include/MAXTRA.H CSR.H ../include/MXLD.H ../include/MAXCOO.H
csrint.o: ../include/MXSTEP.H ../include/MAXTRA.H CSR.H ../include/MAXCOO.H
dejaca.o: ../include/MAXTRA.H ../include/MAXCOO.H
dipi.o: ../include/MXLD.H
dipolm.o: PARIZ.H ../include/MAXTRA.H ../include/MXLD.H
dipsi.o: ../include/MXLD.H
dist2.o: ../include/MAXTRA.H ../include/MAXCOO.H
dpkick.o: ../include/MAXTRA.H ../include/MXLD.H ../include/MAXCOO.H ../include/MXFS.H
el2tub.o: ../include/MXLD.H
elcmii.o: ../include/MXLD.H
elcyld.o: ../include/MXLD.H
elmiri.o: ../include/MXLD.H
end.o: ../include/MXSTEP.H ../include/MAXTRA.H CSR.H ../include/MAXCOO.H
endfit.o: ../include/MXLD.H
esl.o: ../include/MAXTRA.H ../include/MXLD.H ../include/MAXCOO.H
eslvir.o: ../include/MAXTRA.H ../include/MXLD.H ../include/MAXCOO.H
event.o: ../include/MAXTRA.H
fborne.o: ../include/MXLD.H
ffagi.o: ../include/MXLD.H
ff.o: ../include/MAXTRA.H ../include/MAXCOO.H
ffgspi.o: ../include/MXLD.H
fitarr.o: ../include/MXLD.H
fitest.o: ../include/MAXCOO.H ../include/MAXTRA.H ../include/MXLD.H
fitgtv.o: ../include/MXLD.H ../include/MAXCOO.H ../include/MAXTRA.H
fitlav.o: ../include/MAXCOO.H ../include/MAXTRA.H
fitlma.o: ../include/MAXCOO.H ../include/MAXTRA.H
fitnu.o: ../include/MXLD.H
fitset.o: ../include/MXLD.H
fmapw.o: PARIZ.H
focal1.o: ../include/MAXTRA.H ../include/MAXCOO.H
focale.o: ../include/MXLD.H ../include/MAXTRA.H ../include/MAXCOO.H
gascat.o: ../include/MAXTRA.H ../include/MAXCOO.H
gasesl.o: ../include/MAXTRA.H ../include/MAXCOO.H
gasini.o: ../include/MXLD.H
helix.o: ../include/MXLD.H
hist.o: ../include/MAXTRA.H ../include/MXLD.H ../include/MAXCOO.H
icentr.o: PARIZ.H
impcsr.o: ../include/MAXTRA.H ../include/MXSTEP.H CSR.H ../include/MXLD.H ../include/MAXCOO.H ../include/MXFS.H
impdev.o: ../include/MAXTRA.H
impfai.o: ../include/MAXTRA.H ../include/MAXCOO.H ../include/FRMFAI.H
impplt.o: ../include/MAXTRA.H ../include/MXLD.H ../include/MAXCOO.H ../include/MXFS.H ../include/FRMPLT.H
imptra.o: ../include/MAXTRA.H ../include/MAXCOO.H
impvar.o: ../include/MXLD.H
inidrt.o: ../include/MXLD.H
initra.o: ../include/MAXTRA.H ../include/MAXCOO.H
integr.o: ../include/MXSTEP.H
intplf.o: ../include/MAXTRA.H
jcentr.o: PARIZ.H
krtax.o: PARIZ.H
kstop.o: ../include/MAXTRA.H ../include/MXLD.H ../include/MAXCOO.H ../include/MXFS.H
lign.o: ../include/MAXTRA.H ../include/MAXCOO.H
lpsfit.o: ../include/MAXCOO.H ../include/MAXTRA.H
majtra.o: ../include/MAXCOO.H ../include/MAXTRA.H
mat1.o: ../include/MAXTRA.H ../include/MAXCOO.H
mat2.o: ../include/MXLD.H ../include/MAXTRA.H ../include/MAXCOO.H
matimp.o: ../include/MAXCOO.H ../include/MAXTRA.H
matric.o: ../include/MXLD.H ../include/MAXTRA.H ../include/MAXCOO.H
mcdes.o: ../include/MAXTRA.H ../include/MAXCOO.H
mcdesi.o: ../include/MAXTRA.H ../include/MXLD.H ../include/MAXCOO.H
mcdesl.o: ../include/MAXTRA.H ../include/MAXCOO.H
mcobj1.o: ../include/MXLD.H ../include/MAXCOO.H ../include/MAXTRA.H
mcobj3.o: ../include/MXLD.H ../include/MAXCOO.H ../include/MAXTRA.H
mcobj.o: ../include/MXLD.H ../include/MAXCOO.H ../include/MAXTRA.H
mix.o: ../include/MAXCOO.H
mksa.o: ../include/MAXCOO.H
modord.o: ../include/MXLD.H
multip.o: ../include/MAXTRA.H
multpo.o: ../include/MAXTRA.H ../include/MXLD.H ../include/MXFS.H
normal.o: ../include/MXLD.H ../include/MAXTRA.H ../include/MAXCOO.H
obj3.o: ../include/MAXTRA.H ../include/MXLD.H ../include/MAXCOO.H ../include/FRMFAI.H
obj5.o: ../include/MXLD.H ../include/MAXCOO.H ../include/MAXTRA.H
obj6.o: ../include/MXLD.H ../include/MAXCOO.H ../include/MAXTRA.H
obj8.o: ../include/MXLD.H ../include/MAXCOO.H ../include/MAXTRA.H
objeta.o: ../include/MXLD.H ../include/MAXCOO.H ../include/MAXTRA.H
objets.o: ../include/MXLD.H ../include/MAXCOO.H ../include/MAXTRA.H
objtrk.o: ../include/MAXTRA.H ../include/MXLD.H ../include/MAXCOO.H
okkt.o: ../include/OKKT.H
partic.o: ../include/MXLD.H ../include/MAXTRA.H ../include/MAXCOO.H
path.o: ../include/MAXTRA.H ../include/MXLD.H ../include/MAXCOO.H
pavel.o: PARIZ.H
pckup.o: ../include/MAXTRA.H ../include/MAXCOO.H
pltdat.o: ../include/MXLD.H ../include/MAXTRA.H ../include/MAXCOO.H
polmes.o: PARIZ.H ../include/MXLD.H
prdata.o: ../include/MXLD.H
quasex.o: ../include/MXLD.H ../include/MAXCOO.H
rautor.o: ../include/MXLD.H
raysyn.o: ../include/MXLD.H ../include/MXFS.H ../include/MAXTRA.H
rbend.o: ../include/MXLD.H
rbinar.o: ../include/MXLD.H
rcarte.o: ../include/MXLD.H
rcavit.o: ../include/MXLD.H
rcible.o: ../include/MXLD.H
rcolli.o: ../include/MXLD.H
rcsr.o: ../include/MXLD.H
rdip.o: ../include/MXLD.H
rdips.o: ../include/MXLD.H
rebel.o: ../include/MAXTRA.H ../include/MXLD.H ../include/MAXCOO.H
rebmul.o: ../include/MXLD.H
refer.o: ../include/MAXCOO.H ../include/MAXTRA.H
rel2tu.o: ../include/MXLD.H
relcmi.o: ../include/MXLD.H        
relcyl.o: ../include/MXLD.H
relmir.o: ../include/MXLD.H        
rempli.o: ../include/MXLD.H
reset.o: ../include/MAXTRA.H ../include/MAXCOO.H ../include/MXFS.H
rfaist.o: ../include/MXLD.H
rffags.o: ../include/MXLD.H
rfit.o: ../include/MXLD.H
rfitgt.o: ../include/MXLD.H
rgasca.o: ../include/MXLD.H
rhelix.o: ../include/MXLD.H
rmcdes.o: ../include/MXLD.H
rmcobj.o: ../include/MXLD.H ../include/MAXCOO.H ../include/MAXTRA.H
rnorma.o: ../include/MXLD.H
robjet.o: ../include/MXLD.H ../include/MAXCOO.H ../include/MAXTRA.H
robjta.o: ../include/MXLD.H ../include/MAXTRA.H ../include/MAXCOO.H
rparti.o: ../include/MXLD.H
rscal.o: ../include/MXLD.H ../include/MXFS.H
rsepar.o: ../include/MXLD.H
rsimb.o: ../include/MXLD.H
rspn.o: ../include/MXLD.H ../include/MAXTRA.H ../include/MAXCOO.H
rsrlos.o: ../include/MXLD.H
rtrans.o: ../include/MXLD.H
rundul.o: ../include/MXLD.H
runipo.o: ../include/MXLD.H
rwienf.o: ../include/MXLD.H        
scaler.o: ../include/MXFS.H ../include/MXLD.H
scalin.o: ../include/MXLD.H ../include/MXFS.H
scumul.o: ../include/MAXCOO.H
separa.o: ../include/MAXTRA.H ../include/MXLD.H ../include/MAXCOO.H
serv2.o: PARIZ.H
serv32.o: PARIZ.H
serv34.o: PARIZ.H
serv4.o: PARIZ.H
serv5.o: PARIZ.H
settim.o: ../include/MAXCOO.H ../include/MAXTRA.H
soleno.o: ../include/MXLD.H
spn.o: ../include/MXLD.H ../include/MAXTRA.H ../include/MAXCOO.H
spnprn.o: ../include/MXLD.H ../include/MAXTRA.H ../include/MAXCOO.H
spnprt.o: ../include/MAXTRA.H ../include/MAXCOO.H
spnrot.o: ../include/MAXTRA.H ../include/MAXCOO.H
spntrk.o: ../include/MAXTRA.H
sref.o: ../include/MXSTEP.H ../include/MAXTRA.H CSR.H ../include/MAXCOO.H
srloss.o: ../include/MXLD.H
srprn.o: ../include/MXLD.H
stpsiz.o: ../include/MXLD.H
synpar.o: ../include/MXLD.H
toroid.o: ../include/MAXTRA.H
toscac.o: PARIZ.H ../include/MAXTRA.H ../include/MXLD.H
toscap.o: PARIZ.H ../include/MAXTRA.H ../include/MXLD.H
trace.o: ../include/MAXTRA.H ../include/MAXCOO.H
track.o: ../include/MXLD.H
transf.o: ../include/MAXTRA.H ../include/MXSTEP.H CSR.H ../include/MXLD.H ../include/MAXCOO.H ../include/MXFS.H
transm.o: ../include/MXLD.H ../include/MAXCOO.H ../include/MAXTRA.H
traobj.o: ../include/MAXTRA.H ../include/MXLD.H ../include/MAXCOO.H
trarot.o: ../include/MAXTRA.H ../include/MXSTEP.H CSR.H ../include/MXLD.H ../include/MAXCOO.H
trbeam.o: ../include/MXLD.H ../include/MXFS.H
trobj.o: ../include/MAXTRA.H ../include/MXLD.H ../include/MAXCOO.H
trrote.o: ../include/MAXTRA.H ../include/MAXCOO.H
twiss.o: ../include/MXLD.H ../include/MAXCOO.H ../include/MAXTRA.H
undulia.o: ../include/MXLD.H
unduli.o: ../include/MXLD.H
unipot.o: ../include/MXLD.H
waist.o: ../include/MAXTRA.H ../include/MAXCOO.H
wedgki.o: ../include/MAXCOO.H
wienfi.o: ../include/MXLD.H
ymoiny.o: ../include/MAXTRA.H ../include/MAXCOO.H
zgoubi.o: ../include/MXLD.H ../include/MAXTRA.H ../include/MAXCOO.H ../include/MXFS.H PARIZ.H ../include/FILPLT.H ../include/LSTKEY.H
zgoubi_main.o: ../include/MXLD.H
