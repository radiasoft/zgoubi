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
C  François Méot <meot@lpsc.in2p3.fr>
C  Service Accélerateurs
C  LPSC Grenoble
C  53 Avenue des Martyrs
C  38026 Grenoble Cedex
C  France
C1       2         3         4         5         6         7
C2345678901234567890123456789012345678901234567890123456789012
      BLOCK DATA BLOCK
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

C----- NUMERO DES UNITES LOGIQUES D'ENTREES-SORTIE
      COMMON/CDF/ 
     >      IES,IORDRE,LCHA,LIST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN
 
C----- CONSTANTES
      COMMON/CONST/ CL,PI,DPI,RAD,DEG,QE,AH

      PARAMETER (MRD=9)
      COMMON/DROITE/ AM(MRD),BM(MRD),CM(MRD),IDRT

      COMMON/EFBS/ AFB(MRD), BFB(MRD), CFB(MRD), IFB

      CHARACTER  KAR(41)
      COMMON/KAR/ KAR
 
      INCLUDE "MAXTRA.H"
      COMMON/OBJET/ FO(6,1),KOBJ,IDMAX,IMAXT
 
      LOGICAL ZSYM
      COMMON/OPTION/ KORD,KFLD,MG,LC,ML,ZSYM
 
      COMMON/PTICUL/ AAM,Q,G,TO
 
      COMMON/REBELO/ NRBLT,IPASS,KWRI,NNDES,STDVM
 
      CHARACTER FAM*8,KLEY*10
      PARAMETER (MXF=30,MXC=10) 
      COMMON/SCAL/
     >  SCL(MXF,MXC),TIM(MXF,MXC),FAM(MXF),KTI(MXF),KSCL,KLEY
 
C----- CONVERSION COORD. (CM,MRD) -> (M,RD)
      PARAMETER (MXJ=7)
      COMMON/UNITS/ UNIT(MXJ-1)

      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB

      DATA NDAT,NRES,NPLT,NFAI,NMAP,NSPN/3,4,1,2,8,9/
 
      DATA CL, PI, DPI, RAD, DEG, QE, AH /
     >2.99792458D8 , 3.141592653589D0, 6.283185307178D0,
     > .01745329252D0 , 57.29577951D0, 1.60217733D-19, 6.626075D-34 /
 
      DATA IDRT / 0 /

      DATA IFB / 0 /

      DATA (KAR(I),I=1,41) /
     > 'O','A','B','C','D','E','F','G','H','I','J','K','L','M','N'
     >,'P','Q','R','T','U','V','W','X','Y','Z','2','3','4','5','6'
     >,'7','8','(',')','+','-','/','=','"','0','*'/
 
      DATA KOBJ /0/
 
      DATA MG,LC,ML,ZSYM/ 1,2,3,.TRUE./
 
      DATA Q / 1.60217733D-19  /
 
      DATA NRBLT,IPASS/0, 1/
 
      DATA (FAM(I),I=1,MXF)/
     > 'AIMANT' , 'QUADRUPO', 'SEXTUPOL', 'QUADISEX' , 'SEXQUAD'
     >,'TOSCA3D', 'OCTUPOLE', 'DECAPOLE', 'DODECAPO'
     >, 'TOSCA' , 'MULTIPOL' , 'DIPOLE'
     >, 'BEND'    , 'SOLENOID' , 'CAVITE'
     >,'POISSON', 14*' ' /
 
C                    Y     T     Z       P     X,S   dp/p
      DATA UNIT / .01D0,.001D0,.01D0, .001D0, .01D0, 1.D0 /

C      DATA KX,KY,IAX,LIS,NB /6, 2, 1, 1, 100 /
      DATA KX,KY,IAX,LIS,NB /2, 3, 1, 1, 100 /

      END
