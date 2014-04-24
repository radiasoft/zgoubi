
1/ execute ./scanSpinResonances/scanSpinResonances_launch
This will create the zgoubi_xxx.dat files and run these
Note : this will use /home/meot/zgoubi/struct/tools/spin/xing_geneZgDat/geneZGDat4Xing_fromCalcStrength,
so it requires geneZGDat4Xing.data and scanSpinResonances.In

2/ execute ./xing_dataTreatment/dataTreatment
This will compute resonance strengths etc. from the b_zgoubi.fai files

3/ execute ./xing_gnuplots/gnuplot.cmd
This will create the plots :  Sz vs. energy, xx' zz' dpPhase etc. for monitoring

4/ execute xing_geneTexLog
This will create part of log.tex for these particular xing experiments

5/ possibly, merge these log.tex files into a single one
