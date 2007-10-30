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
      SUBROUTINE PCKUP
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     -----------------------------------------------------
C     Average orbit (multiturn AND multiparticle summmation)
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
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),IMAX,IEX(MXT),IREP(MXT)
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM

      DIMENSION FPUL(MXPUD,MXPU), FPUL2(MXPUD-2,MXPU)
      SAVE FPUL, FPUL2

C----- Pick-up number. Reset to 0 in SBR PICKUP
      IPU = IPU + 1
      IF(IPU .GT. MXPU) 
     >CALL ENDJOB('Too  many  c.o.  pick-ups,  max  is ',MXPU)

      NT = 0
      DO 2 I = 1, IMAX
        IF( IEX(I) .GT. 0) THEN
          NT = NT+1
          DO 1 J = 1, 7
            FJI =  F(J,I)
            FPUL(J,IPU) = FPUL(J,IPU) + FJI
            FPUL2(J,IPU) = FPUL2(J,IPU) + FJI*FJI
 1        CONTINUE
        ENDIF
 2    CONTINUE
      FPUL(8,IPU) = NT

          DO 3 J = 1, 5
            FPU(J,IPU) = FPU(J,IPU) + FPUL(J,IPU)
 3        CONTINUE
C          cumulated path length
          FPU(6,IPU) = FPUL(6,IPU)
C          cumulated time
          FPU(7,IPU) = FPUL(7,IPU)
          FPU(8,IPU) = FPU(7,IPU) + NT

C------- Record pick-up position (cm)
      IF(IPASS .EQ. 1) FPU(9,IPU) = F(6,1)

      RETURN

      ENTRY PCKUP1
      DO 4 I = 1, IPU
        NT = NINT(FPUL(8,I))
        WRITE(NFPU,FMT= '(1P,2X,I4,1X,G15.7,7G12.4,I9,I7)' )
     >  I,FPU(9,I),(FPUL(J,I),J=2,7),FPUL(1,I),NT,IPASS
        WRITE(NFPU,FMT= '(22X,7G12.4)') (FPUL2(J,I),J=2,7),FPUL2(1,I)
 4    CONTINUE

      RETURN

      ENTRY PCKUP2
      CALL RAZ(FPUL,MXPUD*MXPU)
      CALL RAZ(FPUL2,(MXPUD-2)*MXPU)

      RETURN

      END
