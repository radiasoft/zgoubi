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
C  Upton, NY, 11973, USA
C  -------
      BLOCK DATA BLOCK
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
C----- Pick-ups
      PARAMETER (MXPUD=9,MXPU=5000)
       INCLUDE "C.CO.H"     ! COMMON/CO/ FPU(MXPUD,MXPU),KCO,NPU,NFPU,IPU
C----- CONSTANTES
      INCLUDE "C.CONST_2.H"     ! COMMON/CONST/ CL9,CL,PI,RAD,DEG,QEL,AMPROT,CM2M
      INCLUDE "C.CONST2.H"     ! COMMON/CONST2/ ZERO, UN
C--------
      INCLUDE "C.DEPL.H"     ! COMMON/DEPL/ XF(3),DXF(3),DQBRO,DTAR
C      PARAMETER (MDR=9)
      INCLUDE "C.DROITE_2.H"     ! COMMON/DROITE/ AM(MDR),BM(MDR),CM(MDR),IDRT
      INCLUDE "C.EFBS.H"     ! COMMON/EFBS/ AFB(2), BFB(2), CFB(2), IFB
      INCLUDE "MAXTRA.H"
      INCLUDE "C.DESIN.H"     ! COMMON/DESIN/ FDES(7,MXT),IFDES,KINFO,IRSAR,IRTET,IRPHI,NDES
C     >,AMS,AMP,AM3,TDVM,TETPHI(2,MXT)
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      INCLUDE "C.GASC.H"     ! COMMON/GASC/ AI, DEN, KGA
      CHARACTER(1) KAR(41)
      INCLUDE "C.KAR.H"     ! COMMON/KAR/ KAR
      INCLUDE 'MXLD.H'
      PARAMETER (LBLSIZ=20)
      CHARACTER(LBLSIZ) LABEL
      INCLUDE "C.LABEL.H"     ! COMMON/LABEL/ LABEL(MXL,2)
      INCLUDE "C.OBJET.H"     ! COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT,KZOB
      INCLUDE "C.ORDRES.H"     ! COMMON/ORDRES/ KORD,IRD,IDS,IDB,IDE,IDZ
      INCLUDE 'C.PDATA.H'  ! COMMON /PDATA/ QELM,AMLEC,GLEC,AMPRO,GPRO,AMMU,GMU,TAUMU,Q3HE,AM3HE,G3HE
      INCLUDE "C.PTICUL_2.H"     ! COMMON/PTICUL/ AAM,Q,G,TO
      INCLUDE "C.REBELO.H"   ! COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,HDPRF,DP,QBR,BRI
      INCLUDE 'MXFS.H'
      INCLUDE 'MXSCL.H'
      INCLUDE "C.SCAL.H"     ! COMMON/SCAL/ SCL(MXF,MXS,MXSCL),TIM(MXF,MXS),NTIM(MXF),KSCL
      LOGICAL TSPCH
      COMMON/SPACECHA/ TLAMBDA,RBEAM(2),XAVE(2),EMITT(2),TAVE,BUNCH_LEN,
     >                EMITTZ, BTAG, SCKX, SCKY, TSPCH
      INCLUDE "C.STEP.H"     ! COMMON/STEP/ TPAS(3), KPAS
      INCLUDE "C.SYNRA.H"     ! COMMON/SYNRA/ KSYN
      INCLUDE "C.TYPFLD.H"     ! COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
C----- CONVERSION DES COORD. (CM,MRD) -> (M,RD)
      INCLUDE "C.UNITS.H"     ! COMMON/UNITS/ UNIT(MXJ)
      PARAMETER (MXV=60)
      INCLUDE "C.VARY.H"  ! COMMON/VARY/ NV,IR(MXV),NC,I1(MXV),I2(MXV),V(MXV),IS(MXV),W(MXV),
                          !     >IC(MXV),IC2(MXV),I3(MXV),XCOU(MXV),CPAR(MXV,27)

      PARAMETER (MDR3= 3*MDR)

      DOUBLE PRECISION QE0
      PARAMETER (QE0=1.602176487D-19)
      PARAMETER (PROTM=938.27203D0)

      DATA LF,LST / 2 * 0 /
      DATA KCO / 0 /
      DATA CL, QEL, AMPROT / 2.99792458D8, QE0, PROTM /
      DATA XF,DXF,DQBRO / 3*0.D0, 3*0.D0, 0.D0 /
      DATA ZERO, UN / 0.D0, 1.D0 /
      DATA IDRT, AM, BM, CM / 0, MDR3*0.D0 /
      DATA IFB / 0 /
      DATA IFDES / 0 /
      DATA KGA / 0 /
      DATA AMQLU,PABSLU/6*.FALSE./
      DATA (KAR(I),I=1,41) /
     > 'O','A','B','C','D','E','F','G','H','I','J','K','L','M','N'
     >,'P','Q','R','U','V','W','X','Y','Z','2','3','4','5','6'
     >,'7','8','9','0','(',')','+','-','/','=','"','*'/
      DATA KZOB, KOBJ / 0, 0 /
      DATA KFLD,MG,LC,ML,ZSYM/ 1,1,2,3,.TRUE./
      DATA IDS, KORD / 4, 2 /
      DATA AAM, Q / 0.D0, 1D0 /
      DATA NRBLT,IPASS/ 0, 1/
Changed dpref to dp/dp_ref - integer part. FM Mar 2018
C     DATA DPREF / 1.D0 /
      DATA DPREF, HDPRF / 0.D0, 1.D0 /
      DATA KPAS / 0 /
      DATA CPAR / 1620*0.D0 /
      DATA IDMAX / 1 /
      DATA KSYN / 0 /
C Particle data
      DATA QELM / QE0 /
      DATA AMLEC / 0.510998928D0 /
      DATA GLEC / 1.15965218076D-3 /
      DATA AMPRO / PROTM /
      DATA GPRO / 1.79284735D0 /
      DATA AMMU / 105.6583715D0 /
      DATA GMU / 11659208.9D-10 /
      DATA TAUMU / 2.1969811D-6 /
      PARAMETER (DQE0= 2.D0*QE0 )
      DATA Q3HE / DQE0 /
      DATA AM3HE / 2808.39148D0 /
      DATA G3HE / -4.1841538D0 /

      DATA TSPCH / .FALSE. /

C----- To yield MKSA units :
C                                1      2     3     4     5    6     7
C                                Y      T     Z     P     S    D    time
C                                m     rad    m    rad    m    1     s
      DATA (UNIT(I),I=1,MXJ) / 1.D-2,1.D-3,1.D-2,1.D-3,1.D-2,1.D0,1.D-6/
      DATA (NTIM(I),I=1,MXF) / MXF * 0 /
      PARAMETER (MXL2=MXL*2)
      DATA LABEL / MXL2*' ' /
      END
