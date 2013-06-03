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
      SUBROUTINE PCKUP(NOEL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     -----------------------------------------------------
C     Pick-up signal (multiturn AND multiparticle summmation)
C     at labeled elements.
C     MPULAB = max number of LABEL's. MXPU = max number of 
C     pick-ups (virtual pick-ups, positionned at indicated 
C     labeled elements!) for CO measurments.
C     -----------------------------------------------------
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      PARAMETER (MXPUD=9,MXPU=1000)
      COMMON/CO/ FPU(MXPUD,MXPU),KCO,NPU,NFPU,IPU
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM

      DIMENSION FPUL(MXPUD,MXPU), FPUL2(MXPUD-2,MXPU), NOELPU(MXPU)
      SAVE FPUL, FPUL2, NOELPU

      DIMENSION TBT(MXPUD,2)
      SAVE TBT

      DATA NOELPU / MXPU*0 /

C----- Pick-up number. Reset to 0 in SBR PICKUP
      IPU = IPU + 1

      IF(IPASS.EQ.1) NOELPU(IPU) = NOEL

      IF(IPU .GT. MXPU) 
     >CALL ENDJOB('SBR PCKUP.  Too many c.o. pick-ups,  max is ',MXPU)

      NT = 0
      DO 2 I = 1, IMAX
        TBT(IPU,1) = TBT(IPU,1) + F(2,I)
        TBT(IPU,2) = TBT(IPU,2) + F(4,I)
        IF( IEX(I) .GT. 0) THEN
          NT = NT+1
          DO 1 J = 1, 7
            FJI =  F(J,I)
            IF(J.EQ.1)  FJI = FJI -1.D0
            FPUL(J,IPU) = FPUL(J,IPU) + FJI
            FPUL2(J,IPU) = FPUL2(J,IPU) + FJI*FJI
 1        CONTINUE
        ENDIF
 2    CONTINUE
      FPUL(8,IPU) = NT

          DO 3 J = 1, 5
            FPU(J,IPU) = FPU(J,IPU) + FPUL(J,IPU)
 3        CONTINUE
C          CUMULATED PATH LENGTH
          FPU(6,IPU) = FPUL(6,IPU)
C          cumulated time
          FPU(7,IPU) = FPUL(7,IPU)
C          cumulate turn #
          FPU(8,IPU) = FPU(8,IPU) + NT

C------- Record pick-up position (cm)
      IF(IPASS .EQ. 1) FPU(9,IPU) = F(6,1)

      RETURN

      ENTRY PCKUP1
      DO 4 I = 1, IPU
        NT = NINT(FPUL(8,I))
        WRITE(NFPU,FMT= 
     >    '(1P,2X,I4,1X,E15.7,7(1X,E12.4),I9,I7,1X,7(1X,E12.4))')
     >  I,FPU(9,I),(FPUL(J,I),J=2,6),FPUL(1,I),FPUL(7,I),NT,IPASS
     >  , (FPUL2(J,I),J=2,7), FPUL2(1,I)
 4    CONTINUE
      RETURN

      ENTRY PCKUP2
      IPU = 0
      CALL RAZ(TBT,MXPUD*2)
      CALL RAZ(FPUL,MXPUD*MXPU)
      CALL RAZ(FPUL2,(MXPUD-2)*MXPU)
      RETURN

      ENTRY PCKUP3(NOELI)
      IPU1 = 1
      DO IL = 1, NOELI
        IF(IPU1.LE.MXPU) THEN
          IF(NOELPU(IPU1).EQ.0) THEN
            GOTO 10
          ELSE 
            IF(NOELPU(IPU1).LT.NOELI) IPU1= IPU1+1
          ENDIF
        ELSE
          CALL ENDJOB('SBR PCKUP. Problem : IPU should be <= ',MXPU)
        ENDIF
      ENDDO
 10   CONTINUE
      IPU = IPU1-1            
      CALL RAZ(FPUL,MXPUD*MXPU)
      CALL RAZ(FPUL2,(MXPUD-2)*MXPU)
      IF(NRES.GT.0) THEN
        WRITE(NRES,*) 
     >   ' Found ',IPU,' pick-ups prior to element #',NOELI,'.'
        IF(IPU.GT.0) WRITE(NRES,*) ' Of which last PU is # ',
     >  IPU,' (located at  element NOEL = ',NOELPU(IPU),')'
      ENDIF
      RETURN

      ENTRY PCKUP5(MPU
     >                ,YCM,ZCM)
      if(   FPUL(8,MPU) .le.0 )
     >  call endjob(' SBR pckup. All particles fucked up !',-99)
      YCM = TBT(MPU,1) /    FPUL(8,MPU) 
      ZCM = TBT(MPU,2) /    FPUL(8,MPU) 
      RETURN

      END
