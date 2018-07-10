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
      SUBROUTINE SVDPCK
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     -----------------------------------------------------
C     Pick-up signal (multiturn AND multiparticle summmation)
C     at labeled elements.
C     MPULAB = max number of LABEL's. MXPU = max number of 
C     pick-ups (virtual pick-ups, positionned at indicated 
C     labeled elements!) for CO measurments.
C     -----------------------------------------------------
      INCLUDE 'C.CDF.H'     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      PARAMETER (MXPUD=9,MXPU=5000)
      INCLUDE 'C.CO.H'     ! COMMON/CO/ FPU(MXPUD,MXPU),KCO,NPU,NFPU,IPU
      INCLUDE 'MXLD.H'
      INCLUDE 'C.DON.H'     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE 'MAXTRA.H'
      INCLUDE 'MAXCOO.H'
      LOGICAL AMQLU(5),PABSLU
      INCLUDE 'C.FAISC.H'     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      PARAMETER (LBLSIZ=20)
      CHARACTER(LBLSIZ) LABEL
      INCLUDE 'C.LABEL.H'     ! COMMON/LABEL/ LABEL(MXL,2)
      INCLUDE 'C.REBELO.H'   ! COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM

      DIMENSION NOELPU(MXPU), IQPU(MXL)
      SAVE NOELPU, IQPU

      PARAMETER (KSIZ=10)
C      CHARACTER(KSIZ)  KLOBJ

      CHARACTER(KSIZ) KLEO
      LOGICAL CHKPU
      SAVE IPSS
      
      DATA NOELPU / MXPU*0 /
      DATA CHKPU / .TRUE. /
      DATA IPSS / 1 /
      
C-----Pick-up number. Reset to 0 via ENTRY PICKP2 below, by SBR PICKUP,
C or by SVD procedure in REBELOTE      
      IPU = IPU + 1

      IF(IPASS.EQ.2) THEN  ! Will never be 1, as REBELOTE is at end of pass
        NOELPU(IPU) = NOEL
        IQPU(NOEL) = IPU
      ENDIF

      IF(IPU .GT. MXPU) 
     >CALL ENDJOB('SBR PCKUP.  Too many c.o. pick-ups,  max is ',MXPU)

      NT = 0
      DO I = 1, IMAX
        IF( IEX(I) .GT. 0) THEN
          NT = NT+1
          DO J = 1, 7
            FJI =  F(J,I)
            IF(J.EQ.1)  FJI = FJI -1.D0
            FPU(J,IPU) = FJI
          ENDDO
        ENDIF
      ENDDO
      FPU(8,IPU) = NT
      
C------- Record pick-up position (cm)
      IF(IPASS .EQ. 2) FPU(9,IPU) = F(6,1)
      
      RETURN

      ENTRY SVDPC0
      IPU = 0
      RETURN

      END
