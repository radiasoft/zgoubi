C  ZGOUBI, a program for computing the trajectories of charged particles
C  in electric and magnetic fields
C  Copyright (C) 1988-2007  François Méot
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.
C
C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.
C
C  You should have received a copy of the GNU General Public License
C  along with this program; if not, write to the Free Software
C  Foundation, Inc., 51 Franklin Street, Fifth Floor,
C  Boston, MA  02110-1301  USA
C
C  François Méot <fmeot@bnl.gov>
C  Brookhaven National Laboratory                    és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE RAYSYN(DS,IT,IMAX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CHAVE/ B(5,3),V(5,3),E(5,3)
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      COMMON/PTICUL/ AM,Q,G,TO
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      INCLUDE 'MXLD.H'
      INCLUDE 'MXFS.H'
      PARAMETER (LBLSIZ=8)
      PARAMETER (KSIZ=10)
      CHARACTER FAM*(KSIZ),LBF*(LBLSIZ),KLEY*(KSIZ),LABEL*(LBLSIZ)
      COMMON/SCALT/ FAM(MXF),LBF(MXF,2),KLEY,LABEL(MXL,2)
      COMMON/TRAJ/ YT,T,Z,PT,XT,SAR,TAR,KEX,ITT,AMT,QT

      INCLUDE "MAXTRA.H"
      DIMENSION TPHOT(MXT),TLOSS(MXT)

      PARAMETER (NSPLIN=43)
      DIMENSION X(NSPLIN), Y(NSPLIN)
      CHARACTER TYPMAG*8,TYPMAGI*8
      PARAMETER(GAM13=2.6789385347D0,UT=1.D0/3.D0)
      
      SAVE TPHOT, TLOSS, TL2, NSTEP, ECMEAN
      SAVE X, Y
      SAVE CQ,AM2,UNIT,UNITE,IRA,DXSPLI,FEC,FAC,FACG
      SAVE TYPMAG
      DATA TYPMAG / 'ALL' /
      DATA X / .00123D0, .0123D0, .0265D0, .0571D0, .1228D0, .1544D0, 
C Y                     .01    .0136   .02   .0292    
     >.194D0, .221D0, .2619D0, .2893D0, .327D0,.36914D0, 
C Y      .04    .052     .065    .09    .13       .17      .2    
     >.4074D0,.4417D0,.47273D0,.5209D0,.5791D0,.62385D0,.6517D0,
C Y     .25     .29           .35           .4
     >  .6905D0,  .7165D0,     .7494D0,     .7726D0, 
C Y     .5      .6       .68        .8     .85     .88      .9
     >.8104D0,.8401D0, .85954D0, .8835D0,.892D0,.89665D0,.8997D0,
C Y      .93       .95       1.      1.2       1.5       1.7
     >  .90397D0,.9067D0,.9132D0, .93444D0, .9561D0,  .9661D0,
C Y      2.       2.5       3.       4.       5.5        7        10.
     >.9768D0, .98734D0, .9930D0, .99776D0, .999574D0, .9999156D0, 1.D0/
      DATA Y / 1.D-9, 1.D-6, 1.D-5, 1.D-4, 1.D-3, 2.D-3, 4.D-3, 6.D-3, 
     >.01D0, .0136D0, 
     >.02D0, .0292D0, .04D0,  .052D0, .065D0, .09D0, .13D0,.17D0,
     > .2D0,.25D0,  .29D0, .35D0, 
     >.4D0, .5D0, .6D0, .68D0, .8D0, .85D0, .88D0, .9D0,.93D0,
     > .95D0, 1.D0, 1.2D0, 1.5D0,1.7D0, 
     > 2.D0, 2.5D0, 3.D0, 4.D0,5.5D0,7.D0,1.D1/
      DATA UNIT,UNITE / 1.D-2, 1.D-6/

      IF(TYPMAG.NE.'ALL') THEN
          IF(KLEY.NE.TYPMAG) RETURN
      ENDIF

C      N is the maximum value of k with non-zero proba of Poisson law
C----- Y(X) : data to be splined (i.e., the integral of the K_5/3 sum)

C----- Curvature (/m)
      CURV=SQRT(B(1,1)*B(1,1)+B(1,2)*B(1,2)+B(1,3)*B(1,3)) / UNIT
C----- Momentum (MeV/c), energy (MeV) 
C         BR=Brho (kG.cm)
      P = QBR*CL9
      QBFIELD=(ABS(QBR)*1.D-3)*CURV
      EN2 = P*P+AM2
      EN=SQRT(EN2)
      BTA=P/EN
      G3=EN*EN2/(AM*AM2)
C----- Parameter of Poisson law (step DS is in cm)
      A=FAC*BTA*BTA*QBFIELD*(DS*UNIT)
C         K= INT(APDPO(A))
         K= INT(POIDEV(A))
         TPHOT(IT)=TPHOT(IT)+K
C----- Photon energies (ELOSS*EC) at current step
C----- Cumulated energy loss TLOSS per particle (over all integr. steps) 
      EC=FEC*G3*CURV
      NSTEP=NSTEP+1
      ECMEAN=ECMEAN+EC
      ELOSS=0.D0
      DO 19 JJ=1,K
        R1=RNDM()
C-------- 0.26 is a (somewhat arbitrary) limit within which 
C    G(E/Ec)~(E/Ec)^(1/3), giving (E/Ec) with better than 1% precision  
        IF(R1.LT.0.26D0) THEN
          DELOSS=(R1/FACG)**3
        ELSE
C          X2=X(1)+DXSPLI*R1
          X2=DXSPLI*R1
          DELOSS=SPLINT(X,Y,NSPLIN,X2)
        ENDIF
        ELOSS=ELOSS+DELOSS
 19   CONTINUE
      DTI=ELOSS*EC
      TLOSS(IT)=TLOSS(IT)+DTI
      TL2=TL2+DTI*DTI
C----- Correction to particle rigidity
      EN = EN-DTI
      QBR = SQRT(EN*EN-AM2)/CL9
      BRI = QT/QBR

      RETURN

      ENTRY RAYSY1(IMAX,IRAI)
      IRA=IRAI
C      DXSPLI=(X(NSPLIN)-X(1))
      DXSPLI=X(NSPLIN)
      CLQE=CL/QE
      CQ = CLQE*1.D-9*Q
C------ AM=rest mass (MeV)
      AM2 = AM*AM
      HBAR=6.6260755D-34/(2.D0*PI)
      R0=QE*1.D-7*CL*CL/(AM*1.D6)
C------- Working unit for energies is MeV
      FEC=1.5D0*HBAR*CLQE   * UNITE
      FAC=2.5D0*R0*QE/(SQRT(3.D0)*HBAR) 
      FACG=12.D0*SQRT(3.D0)/(5.D0*2.D0**UT*GAM13)
      CALL RAZ(TPHOT,IMAX)
      CALL RAZ(TLOSS,IMAX)
      TL2=0.D0
      NSTEP=0
      ECMEAN=0.D0
      RETURN

      ENTRY RAYSY2(IMAX,LUN)
      IF(LUN.GT.0) WRITE(LUN,FMT='(/,
     >    '' Total #pass * #part,  #photons,  total loss (MeV)'')')
      TTPHOT=0.D0
      TTLOSS=0.D0
      DO 55 I=1,IMAX
        TTPHOT=TTPHOT+TPHOT(I)
        TTLOSS=TTLOSS+TLOSS(I)
        IF(LUN.GT.0) WRITE(LUN,FMT='(I6,'' *'',I6,2G15.7)') 
     >                                IPASS,I,TPHOT(I),TLOSS(I)
 55   CONTINUE
      IF(NRES.GT.0) THEN
C        PP = BORO*CL*1.D-9*Q/QE
C        PP = BORO*CL*1.D-9*Q
        PP = BORO*CL9*Q
        EE = SQRT(PP*PP+AM*AM)
        WRITE(NRES,FMT='(/,10X,
     >  '' S.R. statistics, from beginning of structure,'',/,15X,
     >  '' on a total of '',1P,G12.4,'' integration steps :'',/)') NSTEP
        XEVNT=DBLE(IMAX*IPASS)
        XSTEP= NSTEP
        WRITE(NRES,FMT='(5X,'' Mean energy loss per particle :'',1P,
     >  T50,G12.4,'' keV.       Relative to initial energy :'',G12.4)') 
     >  TTLOSS/XEVNT *1.D3,TTLOSS/(XEVNT*EE)
        WRITE(NRES,FMT='(5X,'' Critical energy of photons (mean) :'',1P,
     >  T50,G12.4,'' keV'')') ECMEAN/XSTEP *1.D3
        WRITE(NRES,FMT='(5X,'' Mean energy of radiated photons :'',1P,
     >  T50,G12.4,'' keV'')') TTLOSS/TTPHOT *1.D3
        WRITE(NRES,FMT='(5X,'' rms energy of radiated photons :'',1P,
     >  T50,G12.4,'' keV'')') 
     >      SQRT(TL2/TTPHOT-(TTLOSS/TTPHOT)**2) *1.D3
        WRITE(NRES,FMT='(5X,'' Number of photons radiated - Total :'',
     >  1P,T65,G12.4)') TTPHOT
        WRITE(NRES,FMT='(5X,''                            - per'',
     >  '' particle :'',1P,T65,G12.4)') TTPHOT/XEVNT
        WRITE(NRES,FMT='(5X,''                            - per'',
     >  '' particle, per step :'',1P,T65,G12.4)') TTPHOT/XSTEP
      ENDIF
      RETURN

      ENTRY RAYSY3(TYPMAGI)
      TYPMAG=TYPMAGI
      RETURN
      END
