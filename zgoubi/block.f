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
      BLOCK DATA BLOCK
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
C----- Average orbit
      PARAMETER (MXPUD=9,MXPU=1000)
      COMMON/CO/ FPU(MXPUD,MXPU),KCO,NPU,NFPU,IPU
C----- CONSTANTES
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      COMMON/CONST2/ ZERO, UN
C--------
      COMMON/DEPL/ XF(3),DXF(3),DBRO,DTAR
      PARAMETER (MDR=9)
      COMMON/DROITE/ AM(MDR),BM(MDR),CM(MDR),IDRT
      COMMON/EFBS/ AFB(2), BFB(2), CFB(2), IFB
      CHARACTER  KAR(41)
      COMMON/KAR/ KAR
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      LOGICAL ZSYM
      COMMON/OPTION/ KFLD,MG,LC,ML,ZSYM
      COMMON/ORDRES/ KORD,IRD,IDS,IDB,IDE,IDZ
      COMMON/PTICUL/ AAM,Q,G,TO
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      COMMON/RIGID/ BORO,DPREF,DP,BR
      INCLUDE 'MXFS.H'
      COMMON/SCAL/SCL(MXF,MXS),TIM(MXF,MXS),NTIM(MXF),KSCL
      COMMON/STEP/ TPAS(3), KPAS
C      COMMON/STEP/ KPAS, TPAS(3) 
C----- CONVERSION DES COORD. (CM,MRD) -> (M,RD)
      COMMON/UNITS/ UNIT(MXJ)
      PARAMETER (MXV=40) 
      COMMON/VARY/NV,IR(MXV),NC,I1(MXV),I2(MXV),V(MXV),IS(MXV),W(MXV),
     >IC(MXV),IC2(MXV),I3(MXV),XCOU(MXV),CPAR(MXV,7)
 
      PARAMETER (MDR3= 3*MDR)
    
      DATA LF,LST / 2 * 0 /
      DATA KCO / 0 /
      DATA CL, QE, AMPROT / 2.99792458D8, 1.60217733D-19, 938.2723D0 /
      DATA XF,DXF,DBRO / 3*0.D0, 3*0.D0, 0.D0 /
      DATA ZERO, UN / 0.D0, 1.D0 /
      DATA IDRT, AM, BM, CM / 0, MDR3*0.D0 /
      DATA IFB / 0 /
      DATA (KAR(I),I=1,41) /
     > 'O','A','B','C','D','E','F','G','H','I','J','K','L','M','N'
     >,'P','Q','R','U','V','W','X','Y','Z','2','3','4','5','6'
     >,'7','8','9','0','(',')','+','-','/','=','"','*'/
      DATA KOBJ /0/
      DATA MG,LC,ML,ZSYM/ 1,2,3,.TRUE./
      DATA IDS, KORD / 4, 2 /
      DATA AAM, Q / 2 * 0.D0 /
      DATA NRBLT,IPASS/ 0, 1/
      DATA DPREF / 1.D0 /
      DATA KPAS / 0 /
      DATA CPAR / 280* 0.D0 /
C----- To yield MKSA units : 
C                                1      2     3     4     5    6     7
C                                Y      T     Z     P     S    D    time
C                                m     rad    m    rad    m    1     s
      DATA (UNIT(I),I=1,MXJ) / 1.D-2,1.D-3,1.D-2,1.D-3,1.D-2,1.D0,1.D-6/
      END
