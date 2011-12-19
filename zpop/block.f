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
C  Brookhaven National Laboratory                                                               és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
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

      LOGICAL OKECH, OKVAR, OKBIN
      COMMON/ECHL/OKECH, OKVAR, OKBIN

      CHARACTER  KAR(41)
      COMMON/KAR/ KAR
 
      INCLUDE 'MXVAR.H'
      CHARACTER KVAR(MXVAR)*7, KPOL(2)*9, KDIM(MXVAR)*7
      COMMON/INPVR/ KVAR, KPOL, KDIM

      INCLUDE 'MAXNTR.H'
      COMMON/OBJET/ FO(6,1),KOBJ,IDMAX,IMAXT
 
      LOGICAL ZSYM
      COMMON/OPTION/ KORD,KFLD,MG,LC,ML,ZSYM
 
      COMMON/PTICUL/ AAM,Q,G,TO
 
      COMMON/REBELO/ NRBLT,IPASS,KWRI,NNDES,STDVM
 
      INCLUDE 'MXFS.H'
      COMMON/SCAL/SCL(MXF,MXS),TIM(MXF,MXS),KTI(MXF),KSCL
 
      PARAMETER (NCANAL=2500)
      COMMON/SPEDF/BORNE(6),SPEC(NCANAL,3),PMAX(3),NC0(3)

      COMMON/TRACKM/COOR(NTRMAX,9),NPTS,NPTR

C----- CONVERSION COORD. (CM,MRD) -> (M,RD)
      PARAMETER (MXJ=7)
      COMMON/UNITS/ UNIT(MXJ-1)

      COMMON/VEN/Y1,S1,D,DQ,DS,DM,DSS,DPU,TETA1,
     >           SMIN,SMAX,YMIN,YMAX,SB0,YB0

      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB

      DATA SB0, YB0, TETA1 / 0.D0, 0.D0, 0.D0 / 

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
 
C                    Y     T     Z       P     X,S   dp/p
      DATA UNIT / .01D0,.001D0,.01D0, .001D0, .01D0, 1.D0 /

C      DATA KX,KY,IAX,LIS,NB /6, 2, 1, 1, 100 /
      DATA KX,KY,IAX,LIS,NB /2, 3, 1, 1, 100 /

      DATA DS / 0.D0 /

      DATA  OKECH ,OKVAR, OKBIN / .FALSE., .TRUE., .FALSE.  /

      DATA KVAR/
     >' dp/p ','   Y  ','  Y'' ','   Z  ','  Z'' ','   s  ',' Time ',
     >'   X  ',' Step ','   r  ',
     >'dp/p|o',' Y_o  ',' Y''_o ',' Z_o  ',' Z''_o  ',' s_o  ','Time_0',
     >'Phase ',' dp/p ','KinEnr',
     >'  Sx  ','  Sy  ','  Sz  ',' <S>  ',
     >' <Sx> ',' <Sy> ',' <Sz> ','COUNT ','      ',
     >'  Bx  ','  By  ','  Bz  ','  Br  ',
     >'  Ex  ','  Ey  ','  Ez  ','  Er  ',
     >' S_out',' Pass#'  ,2*'      ',
     >' Y_Lab','      ',' Z_Lab','      ','      ','      ',' X_Lab',
     >8*' ',
     >' lmnt#' ,
     >13*' '
     >/
      DATA KPOL/ 'cartesian' , 'cylindr.' /

C      DATA BORNE/ .01D0, .99D0, .01D0, .99D0, .001D0, .999D0 /
      DATA BORNE/ .65D0, .7D0, .65D0, .7D0, .001D0, .999D0 /
      DATA NC0/ 2000, 2000, 2000 /

      DATA NPTS / NTRMAX /

      DATA KDIM/
     >'       ','  (m)  ',' (rad) ','  (m)  ',' (rad) ','  (m)  ',
     >'(mu_s) ','  (m)  ','  (m)  ','  (m)  '                  ,
     >'       ','  (m)  ',' (rad) ','  (m)  ',' (rad) ','  (m)  ',
     >'(mu_s) ',' (rad) ','       ',' (MeV) ',9*'       ',
     > 4*'  (T)  ', 4*'(eV/m) ' ,
     >'  (m)  ',3*' ',
     >'  (m)  ','      ','  (m)  ','      ','      ','      ','  (m)  ',
     >22*'   '/

      END
