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
      SUBROUTINE MINO1(FONC,N,X,P,V,Y,XI,F0,FINI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(*),P(*),V(*),Y(*),XI(*)

      LOGICAL ITSENS

      PARAMETER (NITER=30)

      PARAMETER (VALMAX=1.D38)
      PARAMETER (D=3.0D0)
      EXTERNAL FONC

      SAVE PNLTY, ICPTMA
      DATA PNLTY, ICPTMA / 1D-10, 1000 /

      CALL CPTINI
      NI=0
      FINI=VALMAX
1     CONTINUE
      NI=NI+1
      CALL CPTFCT(FONC,F0)
      CALL IMPVAR(6,NI)
      CALL IMPCTR(6,F0)
      CALL CPTWRT(
     >            I)
      CALL BUTEE(N,X,X(N+1),X(2*N+1),XI)
      WRITE(6,1040) FINI
1040  FORMAT(/,' Xi2 =',1P,G20.5,'   Busy...',/)
C1040  FORMAT(/,' TAPER Ctrl C POUR REPRENDRE LE CONTROLE ',/)

C Stop test
      IF(FINI .LT. PNLTY
     >.OR. NI .GT. NITER
     >.OR. I .GT. ICPTMA) THEN
         IVAR=0
         DO I=1,N
            IF(XI(I) .EQ. X(I)) THEN
               P(I)=P(I)/D
               IVAR=IVAR+1
            ENDIF
         ENDDO
         IF(IVAR .NE. 0) THEN
            IF(NI.LE.3) THEN
               DO I=1,N
                  X(I)=XI(I)
               ENDDO
               FINI=VALMAX
            ENDIF
         ELSE
            DO I=1,N
               P(I)=P(I)/D
            ENDDO
         ENDIF
         CALL CPTFCT(FONC,F0)
         CALL IMPVAR(6,NI)
         CALL IMPCTR(6,F0)
         CALL CPTWRT(
     >               I)

        RETURN

      ENDIF

      IF(ITSENS(FONC,N,X,X(N+1),X(2*N+1),P,V,F0,D)) THEN
         CALL ITAVAN(FONC,N,X,X(N+1),X(2*N+1),Y,V,F0)
         GOTO 1
      ELSE
         IVAR=0
         DO 4 I=1,N
            IF(XI(I) .EQ. X(I)) THEN
               P(I)=P(I)/D
               IVAR=IVAR+1
            ENDIF
4        CONTINUE
         IF(IVAR .NE. 0) THEN
            IF(NI.LE.3) THEN
               DO 25 I=1,N
                  X(I)=XI(I)
25             CONTINUE
               FINI=VALMAX
            ENDIF
         ELSE
            DO 129 I=1,N
               P(I)=P(I)/D
129         CONTINUE
         ENDIF
         CALL CPTFCT(FONC,F0)
         CALL IMPVAR(6,NI)
         CALL IMPCTR(6,F0)
         CALL CPTWRT(
     >               I)

         IF(NI.LE.NITER) GOTO 1
C         WRITE(ABS(NRES),*) 
         WRITE(6,*) 
     >     ' SBR mino1 : Out of mino1 upon NITER = ',NI

      ENDIF

      RETURN

      ENTRY MINO12(PNLTI,ICPTM)
      PNLTY = PNLTI
      ICPTMA = ICPTM
      RETURN

      ENTRY MINO13(
     >             NITERO)
      NITERO = NITER
      RETURN

      END
