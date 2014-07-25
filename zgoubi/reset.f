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
      SUBROUTINE RESET
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "MAXTRA.H"
      COMMON/CHAMBR/ LIMIT,IFORM,YLIM2,ZLIM2,SORT(MXT),FMAG,BMAX
     > ,YCH,ZCH
      PARAMETER (MXPUD=9,MXPU=5000)
      COMMON/CO/ FPU(MXPUD,MXPU),KCO,NPU,NFPU,IPU
      COMMON/DESIN/ FDES(7,MXT),IFDES,KINFO,IRSAR,IRTET,IRPHI,NDES
     >,AMS,AMP,AM3,TDVM,TETPHI(2,MXT)
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      PARAMETER (JH=24,KH=5)
      COMMON/HISTO/ ICTOT(JH,KH),MOYC(JH,KH) ,CMOY(JH,KH),JMAX(JH,KH)
      COMMON/HISTOG/ NC(JH,120,KH),J,NH,XMI(JH,KH),XMO(JH,KH),XMA(JH,KH)
      COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      PARAMETER (LBLSIZ=10)
      INCLUDE 'MXLD.H'
      CHARACTER(LBLSIZ) LABEL
      COMMON /LABEL/ LABEL(MXL,2)
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      LOGICAL ZSYM
      COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      COMMON/ORDRES/ KORD,IRD,IDS,IDB,IDE,IDZ
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      INCLUDE 'MXFS.H'
      INCLUDE 'MXSCL.H'
      COMMON/SCAL/ SCL(MXF,MXS,MXSCL),TIM(MXF,MXS),NTIM(MXF),KSCL
      PARAMETER (KSIZ=10)
      CHARACTER(KSIZ) FAM
      CHARACTER(LBLSIZ) LBF
      COMMON/SCALT/ FAM(MXF),LBF(MXF,MLF)
      COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)
 
      DO I=1,MXF
        FAM(I) = ' ' 
        DO J=1,MLF
          LBF(I,J) = ' ' 
        ENDDO
      ENDDO

C----- Pick-up signal calculation switched off
      KCO = 0

C     ....  Defaults: ORDRE=2 when KALC=3; 
C                     fields have z-symmetry when KALC=1
      KORD=2
      ZSYM=.TRUE.
 
C     ....  REBELOTE
      KOBJ = 0
C      IPASS = 1
C     .... COMPTEUR DE PASSAGE PAR 'REBELOTE'
      NRBLT = 0
 
C----- Decay simulation switched off
      IFDES = 0
C     ... COMPTEUR DE DESINTEGRATION
      NDES = 0

C     ... COORD AU POINT DE DESINTEGRATION
      CALL RAZ(FDES,6*MXT)
 
C     ... TEMPS DE VOL AU POINT DE SORTIE
      CALL RAZ(SORT,MXT)
 
C       ... RAZ LES HISTOS
        CALL IRAZ(NC,JH*120*KH)
        CALL IRAZ(JMAX,JH*KH)
C       ... RESET LES SOMMES
        CALL RAZ(XMO,JH*KH)
        DO 3 K=1,KH
          DO 3 J=1,JH
            XMI(J,K)=1.D10
            XMA(J,K)=-1.D10
 3      CONTINUE
 
C     ....  SPIN TRACKING
      KSPN = 0
 
C     .... SCALING
      KSCL=0
      CALL IRAZ(NTIM,MXF)

      CALL RAZ(F ,MXJ*MXT)
      CALL RAZ(FO,MXJ*MXT)
 
      ENTRY RESET2
C     ... # part. rejected in SBR INTEG
C     ... # part. out of acceptance (collimator, chamber)
C      NRJ = 0
C      NOUT = 0
      CALL CNTRST

C     .... OPTION CHAMBRE => LIMIT=1
      LIMIT=0

C     .... FIELD MAP FILE NAMES
      CALL TOSCA1

C Field map counter reset to 0
      CALL KSMAP0

C-----  structure length ------------
      CALL SCUMS(ZERO)
      CALL SCUMT(ZERO)
C------------------------------------------------

      RETURN
      END
