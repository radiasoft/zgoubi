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
      SUBROUTINE PICKUP
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     -----------------------------------------------------
C     Pick-up signal (multiturn AND multiparticle summmation)
C     at labeled elements.
C     MPULAB = max number of LABEL's. MPX = max number of 
C     pick-ups (virtual pick-ups, positionned at indicated 
C     labeled elements!) for CO measurments.
C     -----------------------------------------------------
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      PARAMETER (MXPUD=9,MXPU=5000)
      COMMON/CO/ FPU(MXPUD,MXPU),KCO,NPU,NFPU,IPU
      PARAMETER (LBLSIZ=10)
      CHARACTER(LBLSIZ) PULAB
      PARAMETER (MPULAB=5)
      COMMON/COT/ PULAB(MPULAB)
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
       
      LOGICAL IDLUNI, OPN
      SAVE OPN, KPCKUP
      LOGICAL FITING
C      INTEGER DEBSTR, FINSTR
 
      DATA OPN / .FALSE. /

      NPU = NINT(A(NOEL,1))
      IF(NPU .GE. 1) THEN
        KCO = 1
      ELSE
        KCO = 0
      ENDIF
      OPN = NINT(A(NOEL,2)) .EQ. 1

      KPCKUP = 1

      IF(KCO .EQ. 0) THEN
        KPCKUP = 0
        IF(NRES .GT. 0) WRITE(NRES,100)
 100      FORMAT(/,20X,' ++++  PICKUPS command is inactive  ++++',/)
        GOTO 98
      ENDIF

      IF(NRES .GT. 0) THEN
        WRITE(NRES,110) 
 110    FORMAT(/,10X,' Pick-up  signal  calculation  requested')
        WRITE(NRES,111) 
 111    FORMAT(/,15X,' Particle coordinates will be averaged',
     >    1X,'at elements labeled:')
        WRITE(NRES,112) (PULAB(I), I=1,NPU) 
 112    FORMAT(20X,A)
      ENDIF

      CALL FITSTA(5,
     >              FITING)

      IF(.NOT. FITING) THEN

        IF(IPASS .EQ. 1) THEN

          CALL RAZ(FPU,MXPUD*MXPU)

 10       CONTINUE
          IF(OPN) THEN
            INQUIRE(FILE='zgoubi.PICKUP.out',ERR=11,NUMBER=LN)
            IF(NRES.GT.0) WRITE(NRES,*) 
     >          ' Pick-up storage file zgoubi.PICKUP.out '
     >         ,' already open under logical unit number ', LN
          ELSE
            INQUIRE(FILE='zgoubi.PICKUP.out',
     >                      ERR=11,OPENED=OPN,NUMBER=LN)
            IF(OPN) GOTO 10
            IF(IDLUNI(
     >                NFPU)) THEN
              OPEN(UNIT=NFPU,FILE='zgoubi.PICKUP.out',ERR=99)
            ELSE
              GOTO 99
            ENDIF 
          ENDIF 
 11       CONTINUE

          WRITE(NFPU,FMT='(A,A,5(1X,A10))') '# PICKUPS ',
     >    '- storage file.  PU are positionned at ',(PULAB(I), I=1,NPU)
        ENDIF
      ENDIF

C----- some reset actions at start of each new pass
C      Total pick-up number
C      IPU = 0
      CALL PCKUP2
      GOTO 98

 99   CONTINUE
      KPCKUP = 0
      CALL ENDJOB('*** Error, SBR PICKUP -> can`t open strage file',-99)
 98   CONTINUE
      CALL REBEL2(KPCKUP)
      RETURN

      ENTRY PICKU1(
     >             KPU)
      KPU = KPCKUP
      RETURN

      END
