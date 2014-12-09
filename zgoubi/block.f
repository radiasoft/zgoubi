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
C  Brookhaven National Laboratory  
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  -------
      BLOCK DATA BLOCK
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
C----- Pick-ups
      PARAMETER (MXPUD=9,MXPU=5000)
      COMMON/CO/ FPU(MXPUD,MXPU),KCO,NPU,NFPU,IPU
C----- CONSTANTES
      COMMON/CONST/ CL9,CL,PI,RAD,DEG,QEL,AMPROT,CM2M
      COMMON/CONST2/ ZRO, ONE
C--------
      COMMON/DEPL/ XF(3),DXF(3),DQBRO,DTAR
      PARAMETER (MDR=9)
      COMMON/DROITE/ AM(MDR),BM(MDR),CM(MDR),IDRT
      COMMON/EFBS/ AFB(2), BFB(2), CFB(2), IFB
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      CHARACTER(1) KAR(41)
      COMMON/KAR/ KAR
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      LOGICAL ZSYM
      COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      COMMON/ORDRES/ KORD,IRD,IDS,IDB,IDE,IDZ
      COMMON/PTICUL/ AAM,Q,G,TO
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      INCLUDE 'MXFS.H'
      INCLUDE 'MXSCL.H'
      COMMON/SCAL/ SCL(MXF,MXS,MXSCL),TIM(MXF,MXS),NTIM(MXF),KSCL
C      COMMON/SCAL/SCL(MXF,MXS),TIM(MXF,MXS),NTIM(MXF),JPA(MXF,MXP),KSCL
      COMMON/STEP/ TPAS(3), KPAS
C----- CONVERSION DES COORD. (CM,MRD) -> (M,RD)
      COMMON/UNITS/ UNIT(MXJ)
      PARAMETER (MXV=60) 
      COMMON/VARY/NV,IR(MXV),NC,I1(MXV),I2(MXV),V(MXV),IS(MXV),W(MXV),
     >IC(MXV),IC2(MXV),I3(MXV),XCOU(MXV),CPAR(MXV,7)
 
      PARAMETER (MDR3= 3*MDR)
    
      DOUBLE PRECISION QE0
      PARAMETER (QE0=1.602176487D-19)

      DATA LF,LST / 2 * 0 /
      DATA KCO / 0 /
      DATA CL, QEL, AMPROT / 2.99792458D8, QE0, 938.27203D0 /
      DATA XF,DXF,DQBRO / 3*0.D0, 3*0.D0, 0.D0 /
      DATA ZRO, ONE / 0.D0, 1.D0 /
      DATA IDRT, AM, BM, CM / 0, MDR3*0.D0 /
      DATA IFB / 0 /
      DATA AMQLU,PABSLU/6*.FALSE./
      DATA (KAR(I),I=1,41) /
     > 'O','A','B','C','D','E','F','G','H','I','J','K','L','M','N'
     >,'P','Q','R','U','V','W','X','Y','Z','2','3','4','5','6'
     >,'7','8','9','0','(',')','+','-','/','=','"','*'/
      DATA KOBJ /0/
      DATA MG,LC,ML,ZSYM/ 1,2,3,.TRUE./
      DATA IDS, KORD / 4, 2 /
      DATA AAM, Q / 0.D0, 1D0 /
      DATA NRBLT,IPASS/ 0, 1/
      DATA DPREF / 1.D0 /
      DATA KPAS / 0 /
      DATA CPAR / 280* 0.D0 /
      DATA IDMAX / 1 /

C----- To yield MKSA units : 
C                                1      2     3     4     5    6     7
C                                Y      T     Z     P     S    D    time
C                                m     rad    m    rad    m    1     s
      DATA (UNIT(I),I=1,MXJ) / 1.D-2,1.D-3,1.D-2,1.D-3,1.D-2,1.D0,1.D-6/
      DATA (NTIM(I),I=1,MXF) / MXF * 0 /
      END
