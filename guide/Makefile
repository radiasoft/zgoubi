ZGOUBITEX=\
 Zgoubi.tex\
 ZGlossTitl.tex\
 ZFigCap.tex\
 Zgoucover.tex\
 ZGlossaires.tex\
 ZOptElem.tex\
 Zlst1.tex\
 ZGlossaires2.tex\
 ZOptElem2.tex\
 Zlst2.tex\
 Zlst3.tex\
 Zbiblio.tex

ZGOUBIFIG=\
 ffagTriplet.eps\
 Fig10.ps\
 Fig11.ps\
 Fig12.ps\
 Fig13.ps\
 Fig14.eps\
 Fig15.ps\
 Fig17.ps\
 Fig18.ps\
 Fig19.ps\
 Fig1.ps\
 Fig20.ps\
 Fig21bis.ps\
 Fig21.ps\
 Fig22.ps\
 Fig23.ps\
 Fig24.ps\
 Fig25.ps\
 Fig26.ps\
 Fig27.ps\
 Fig28.ps\
 Fig29.ps\
 Fig2.ps\
 Fig30.ps\
 Fig31.ps\
 Fig32.ps\
 Fig33.ps\
 Fig34A.eps\
 Fig34.eps\
 Fig35.eps\
 Fig36.ps\
 Fig3.ps\
 Fig4.ps\
 Fig5.ps\
 Fig6.eps\
 Fig7.ps\
 Fig8.ps\
 Fig9a.eps\
 Fig9b.eps\
 FigC1-1.eps\
 FigC2-1.ps\
 FigC3-1.eps\
 FigC5-1.ps\
 FigC6-2a.ps\
 FigC6-2b.ps\
 FigC6-2c.ps\
 FigC6-2d.ps\
 FigC6-4.ps\
 FigC7-1.ps\
 FigC7-2.ps\
 FigCHAREFa.eps\
 FigCHAREFb.eps\
 FigCover1.ps\
 FigCover2.ps\
 FigCover3.ps\
 FigCover4.ps\
 FigELMIRC.eps\
 FigELMIR.eps\
 FigKPOS3.eps\
 G.eps\
 grid.eps\
 NewFig5.ps

.POSIX:

Zgoubi.pdf : Zgoubi.ps
	ps2pdf14 Zgoubi.ps

Zgoubi.ps : Zgoubi.dvi
	dvips -Ppdf -j0 -t a4 Zgoubi -o Zgoubi.ps

Zgoubi.dvi : Zgoubi.aux Zgoubi.ind
	latex Zgoubi

Zgoubi.ind : Zgoubi.idx
	makeindex Zgoubi

Zgoubi.aux Zgoubi.idx : $(ZGOUBITEX)
	latex Zgoubi

clean :
	$(RM) *.aux *.toc Zgoubi.dvi Zgoubi.idx Zgoubi.ind Zgoubi.ilg Zgoubi.log Zgoubi.ps Zgoubi.pdf
