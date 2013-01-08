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
      SUBROUTINE COFIN(KART,NL,LST,DS,KEX,IT,AMT,QT,EVNT,
     >                                           Y,T,Z,P,X,SAR,TAR,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ----------------------------------
C     COORDONNEES ABSOLUES EN FIN DE DS
C      APPELE PAR INTEGR
C     ----------------------------------
      LOGICAL EVNT
      COMMON/AIM/ AE,AT,AS,RM,XI,XFF,EN,EB1,EB2,EG1,EG2
      INCLUDE "MAXTRA.H"
      COMMON/CHAMBR/ LIMIT,IFORM,YL2,ZL2,SORT(MXT),FMAG,BMAX
     > ,YC,ZC
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      COMMON/CONST2/ ZERO, UN
      INCLUDE 'MXSTEP.H'
      INCLUDE 'CSR.H'
C      COMMON/CSR/ KTRA,KCSR,YZXB(MXSTEP,41,36),DWC(MXT)
      COMMON/DEPL/ XF(3),DXF(3),DQBRO,DTAR
      COMMON/DESIN/ FDES(7,MXT),IFDES,KINFO,IRSAR,IRTET,IRPHI,NDES
     >,AMS,AMP,AM3,TDVM,TETPHI(2,MXT)
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      CHARACTER LET
      COMMON/FAISCT/ LET(MXT)
      COMMON/GASC/ AI, DEN, KGA
C      COMMON/MARK/ KART,KALC,KERK,KUASEX
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      LOGICAL ZSYM
      COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      COMMON/PTICUL/ AM,Q,G,TO
      COMMON/REBELO/ NPASS,IPASS,KWRT,NNDES,STDVM
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)
      COMMON/SYNRA/ KSYN

      COMMON/CHAVE/ B(5,3),V(5,3),E(5,3)
 
C----- Write position M0, field at M0,  etc.
      IF(LST .EQ. 2) CALL IMPPLT(NL,Y,T,Z,P,X,SAR,TAR,DS,AMT,QT,KEX,IT)

      SAR=SAR+DS
      Y=Y+XF(2)
      Z=Z+XF(3)

      T=ATAN2(DXF(2),DXF(1))
      P=ATAN( DXF(3)/SQRT(DXF(1)*DXF(1)+DXF(2)*DXF(2)) )
C      P=ATAN2( DXF(3),SQRT(DXF(1)*DXF(1)+DXF(2)*DXF(2)) )
      TAR = TAR + DTAR

      IF(KFLD.GE.LC) THEN
         QBR=QBR + DQBRO
         BRI = QT/QBR
      ENDIF

      IF(KART .EQ. 2) THEN
C-------- Polar coordinates
         DX=ATAN(XF(1)/Y)
         Y=SQRT(XF(1)*XF(1)+Y*Y)
         X=X+DX
         T=T+DX

      ELSE
         X=X+XF(1)
      ENDIF

      IF(EVNT) THEN 
        IF(KSPN .EQ. 1) THEN
C--------- Spin tracking
          IF(KART .EQ. 2) CALL SPNROT(IT,ZERO,ZERO,-DX)
          CALL SPNTRK(DS)
C----- Coherent synchrotron radiation
        ELSEIF(KCSR .EQ. 1) THEN
          IF    (IPASS.EQ.1) THEN
C----------- Store reference trajectory
C            IF(IT .EQ. KTRA) THEN
              CALL IMPCSR(DS,Y,T,Z,P,X,SAR,TAR,IT)
C            ENDIF
          ELSEIF(IPASS.EQ.2) THEN
C----------- Compute interaction undergone by current particle, due to all earlier CSR emittors
            CALL CSRINT(DS,IMAX)
          ENDIF
        ELSE        
          CALL EVENT(
     >    DS,Y,T,Z,P,X,ZERO,QBR,SAR,TAR,KEX,IT,
     >    AMT,Q,BORO,KART,IFDES,KGA,KSYN,IMAX,*99)
        ENDIF
      ENDIF
      RETURN
 99   RETURN 1
      END
