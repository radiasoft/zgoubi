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
      SUBROUTINE PICKUP
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     -----------------------------------------------------
C     Average orbit (multiturn AND multiparticle summmation)
C     at labeled elements.
C     MPULAB = max number of LABEL's. MPX = max number of 
C     pick-ups (virtual pick-ups, positionned at indicated 
C     labeled elements!) for CO measurments.
C     -----------------------------------------------------
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      PARAMETER (MXPUD=9,MXPU=1000)
      COMMON/CO/ FPU(MXPUD,MXPU),KCO,NPU,NFPU,IPU
      PARAMETER (MPULAB=5)
      CHARACTER*10 PULAB
      COMMON/COT/ PULAB(MPULAB)
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      CHARACTER*80 TITRE
      COMMON/TITR/ TITRE 
       
      LOGICAL IDLUNI

      IF(KCO .EQ. 0) GOTO 98

      IF(NRES .GT. 0) THEN
        IF(KCO .EQ. 1) THEN
          WRITE(NRES,110) 
 110      FORMAT(/,10X,' Average  orbit  calculation  requested')

          WRITE(NRES,111) 
 111      FORMAT(/,15X,' Particle coordinates will be averaged',
     >      1X,'at elements labeled:')
          WRITE(NRES,112) (PULAB(I), I=1,NPU) 
 112      FORMAT(20X,A)
        ELSE
          WRITE(NRES,FMT='(/,10X,'' Average orbit calculation *OFF*'')')
        ENDIF
      ENDIF

      IF(IPASS .EQ. 1) THEN
        CALL RAZ(FPU,MXPUD*MXPU)
        IF(IDLUNI(
     >            NFPU)) THEN
          OPEN(UNIT=NFPU,FILE='zgoubi.co',ERR=99)
        ELSE
          GOTO 99
        ENDIF 
      ENDIF

      IF(IPASS .EQ. 1) THEN
        WRITE(NFPU,FMT='(A15,A40,5(1X,A10))') ' Average  orbit',
     >  ' - storage file.  PU are positionned at ',(PULAB(I), I=1,NPU)
        WRITE(NFPU,FMT='(A80)') TITRE
      ENDIF

C----- some reset actions at start of each new pass
C      Total pick-up number
      IPU = 0
      CALL PCKUP2

 98   RETURN
 99   CALL ENDJOB('*** Error, SBR PICKUP -> can`t open strage file',-99)
      RETURN
      END
