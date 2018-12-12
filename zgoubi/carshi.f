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
      SUBROUTINE CARSHI(NBSHIM)
      USE dynhc
      use pariz_namelist_interface, only : ID
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ------------------------------------------------
C     CORRECTION DE LA CARTE DE Champ CALCULEE PAR CARLA. ON INTRODUIT
C     DES SHIMS (OU DES ILOTS) D'AMPLITUDE DB/B0.
C     LES LIMITES DU SHIM SONT TMIN-MAX (DEG) ET RMIN-MAX (CM).
C     ------------------------------------------------
      INCLUDE "C.AIM_3.H"     ! COMMON/AIM/ ATO,AT,ATOS,RM,XI,XF,EN,EB1,EB2,EG1,EG2
      INCLUDE "C.CARSH.H"     ! COMMON/CARSH/ ATS,RMS
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
 
      DOUBLE PRECISION LAMBDA
 
      ATS=ATS*RAD
      ISHIM=0
 
10    CONTINUE
      RMIN  = A(NOEL,13)
      AI     = A(NOEL,62+ISHIM*9+1)
      AF     = A(NOEL,62+ISHIM*9+2)
      RI     = A(NOEL,62+ISHIM*9+3)
      RF     = A(NOEL,62+ISHIM*9+4)
      LAMBDA = A(NOEL,62+ISHIM*9+5)
      GAM = A(NOEL,62+ISHIM*9+6)
      ALP = A(NOEL,62+ISHIM*9+7)
      AMU = A(NOEL,62+ISHIM*9+8)
      BET = A(NOEL,62+ISHIM*9+9)
 
      ISHIM=ISHIM+1
      LAMBDA = 0D0
      DBSB=AMU
      IF(NRES.GT.0)
     > WRITE(NRES,100) ISHIM,AI,AF,RI,RF,LAMBDA,GAM,ALP,AMU,BET
100   FORMAT(/,20X,' CADRE DU SHIM ',I2,' : THETA =',F7.3,' DEG A ',F7.3
     1,' DEG, R =',F7.3,' CM A ',F7.3,' CM',/,20X,'LAMBDA TJRS  = ',F7.3
     1,/,20X,' PARAMETRES DU SHIM : GAMA =',F10.4,'     ALPHA =',F10.4
     1,'     MU (DB/B0) =',F10.4,'     BETA =',F10.4)
      AI=AI*RAD
      AF=AF*RAD
C     ****ON PLACE L'ORIGINE DES COORDONNEES DU SHIM AU CENTRE DU SHIM
      AFAI=AF-AI
      IRI=1+INT((RI-RMIN)/RMS)
      RFRI=RF-RI
      RMSHI=(RF+RI)*.5D0
      IRF=IRI+INT( RFRI /RMS)
C     ***RPLAT(APLAT)=LIMITES RADIALE (ANGULAIRE) DU PLATEAU CENTRAL DU SHIM
      APLAT= AI* RMSHI-2.D0*LAMBDA
      RI=RFRI/2.D0
      RPLAT= RI -2.D0*LAMBDA
C     WRITE(NRES,102) IAI,IAF,IRI,IRF,AI,RI,APLAT,RPLAT
C102  FORMAT(5X,4I5,4G12.4)
      IF(RPLAT.LE.0.D0) RPLAT=1.D-6
      IF(NRES.GT.0) WRITE(NRES,103) APLAT,RPLAT
 103  FORMAT(20X,' DEMI-DIMENSIONS DU SHIM, AU CENTRE : ANGULAIRE  +/-',
     1F10.4,' CM ;  RADIALE  +/-',F10.4,' CM')
      FACSHI=(GAM+ALP/AMU)*BET/ RMSHI / RMSHI
C
C     *** LAMBDA EST TOUJOURS NUL DANS CETTE VERSION , DONC :
      FA = 1D0
      FR = 1D0
C
      DO 2 IR=IRI,IRF
         RAY = (IR-IRI)*RMS-RI
         RAY2= RAY*RAY
C
C        ABSRAY=ABS(RAY)
C        IF(RAY       .LT.-RPLAT) THEN
C           FR=EXP(- ((RAY +RPLAT)/LAMBDA) **4)
C        ELSEIF(ABSRAY.LT. RPLAT) THEN
C           FR=1.
C        ELSEIF(RAY   .GT. RPLAT) THEN
C           FR=EXP(- ((RAY -RPLAT)/LAMBDA) **4)
C        ENDIF
C
         AIP=AI+FACSHI*RAY2
         AFP=AF-FACSHI*RAY2
         IAI=1  +INT( AIP /ATS )
         AFAI= AFP-AIP
         IAF=IAI+INT( AFAI/ATS)
         APLAT= AFAI* RMSHI *.5D0-2.D0*LAMBDA
         IF(NRES.GT.0) WRITE(NRES,105) RAY,IAI,IAF,AIP,AFP,APLAT
105      FORMAT(' RAY,IAI,IAF,AIP,AFP,APLAT :',F10.2,2I5,3G12.4)
         IF(APLAT.LE.0.D0) APLAT=1.D-6
         DO 1 IA=IAI,IAF
C
C           ANG= (IA-IAI)*ATS-AI
C           RANG=ABSRAY*SIN(ANG)
C           IF(RANG         .LT.-APLAT) THEN
C              FA=EXP(- ((RANG+APLAT)/LAMBDA) **4)
C           ELSEIF(ABS(RANG).LT. APLAT) THEN
C              FA=1.
C           ELSEIF(RANG     .GT. APLAT) THEN
C              FA=EXP(- ((RANG-APLAT)/LAMBDA) **4)
C           ENDIF
C
            HC(ID,IA,IR,1,1)=HC(ID,IA,IR,1,1)*(1.D0+FA*FR*DBSB)
C
1        CONTINUE
2     CONTINUE
C
      IF(ISHIM.LT.NBSHIM) GOTO 10
C
      RETURN
      END
