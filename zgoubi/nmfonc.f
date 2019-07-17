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
      DOUBLE PRECISION FUNCTION NMFONC(XX)

      IMPLICIT NONE
      DOUBLE PRECISION XX(*)
      DOUBLE PRECISION NMFINI, NMFON2
      INTEGER NN
      INTEGER MXV
      PARAMETER (MXV=60)
      INCLUDE "C.VAR.H"     ! COMMON /VAR/ X(3*MXV), P(MXV)
      DOUBLE PRECISION FF
      EXTERNAL FF
      INTEGER I
      DOUBLE PRECISION FR
      INTEGER N,NI
      SAVE N,NI
      LOGICAL SYSOUT
      SAVE SYSOUT

      DO 1000 I=1,N
         IF (X(I).LE.X(I+N).OR.X(I+2*N).LE.X(I)) THEN
            WRITE(6,1100) I
 1100       FORMAT(/,10('*'), ' ALARM ! VARIABLE N0',I3,' AT LIMIT ')
            NMFONC=3.40282347D+38
            RETURN
         ENDIF
         X(I) = XX(I)
 1000 CONTINUE
      NI = NI+1
      CALL CPTFCT(FF,FR)
      IF(NI.EQ.1 .OR. INT(NI/10)*10.EQ.NI) THEN
        CALL IMPVAR(6,NI)
        CALL IMPCTR(6,FR)
        CALL CPTWRT(I)
        IF(SYSOUT) WRITE (6,2000) FR
 2000   FORMAT(/,' Xi2 =',1P,E24.16,'   Busy...',/)
      ENDIF

      NMFONC = FR

      RETURN

      ENTRY NMFINI(NN)

      NI=0
      N=NN
      NMFINI = 0D0

      RETURN

      ENTRY NMFON2
      SYSOUT = .FALSE.
      RETURN
      END
