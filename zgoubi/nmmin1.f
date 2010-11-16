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
C  Brookhaven National Laboratory               és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      DOUBLE PRECISION FUNCTION NMMIN1(F,X,DX,N)

      IMPLICIT NONE
      DOUBLE PRECISION F
      EXTERNAL F
      INTEGER N
      DOUBLE PRECISION X(N),DX(N)

      INTEGER NMAX
      PARAMETER (NMAX=100)
      DOUBLE PRECISION FF(NMAX+3),XX(NMAX,NMAX+3)
      DOUBLE PRECISION XA(NMAX),XMIN(NMAX),XMAX(NMAX)
      INTEGER P(NMAX+3)
      DOUBLE PRECISION FMIN,FMAX,KN,XN
      INTEGER I,J,K
      LOGICAL SHRINK

      DOUBLE PRECISION DR,DE,DOC,DIC
      PARAMETER (DR=1D0,DE=2D0,DOC=0.5D0,DIC=-0.5D0)

      P(1) = 1
      DO 1000 I=1,N
         XN = -1D0/I
         KN = -SQRT(I*(I+2D0))*XN
         DO 1100 J=1,I
            DO 1110 K=1,I-1
               XX(K,J) = XX(K,J)*KN
 1110       CONTINUE
            XX(I,J) = XN
 1100    CONTINUE
         DO 1200 K=1,I-1
            XX(K,I+1) = 0D0
 1200    CONTINUE
         XX(I,I+1) = 1D0
         P(I+1) = I+1
 1000 CONTINUE

      DO 2000 I=1,N+1
         DO 2100 J=1,N
            XX(J,P(I)) = X(J)+DX(J)*XX(J,P(I))
 2100    CONTINUE
         FF(P(I)) = F(XX(1,P(I)))
         CALL LSTINS(FF,I,P,I-1)
 2000 CONTINUE

      P(N+2) = N+2
      P(N+3) = N+3

      FMIN = FF(P(1))
      FMAX = FF(P(N+1))

      K=0
 3000 IF (K.GT.N.OR.FF(P(N+1)).LE.FF(P(1))) GOTO 4000
      DO 3100 I=1,N
         XA(I) = 0D0
         DO 3110 J=1,N
            XA(I) = XA(I)+XX(I,P(J))
 3110    CONTINUE
         XA(I) = XA(I)/N
         XX(I,P(N+2)) = XA(I) + DR*(XA(I)-XX(I,P(N+1)))
 3100 CONTINUE
      FF(P(N+2)) = F(XX(1,P(N+2)))
      IF (FF(P(1)).LE.FF(P(N+2)).AND.FF(P(N+2)).LT.FF(P(N))) THEN
         CALL LSTINS(FF,N+2,P,N+1)
      ELSE IF (FF(P(N+2)).LT.FF(P(1))) THEN
         DO 3200 I=1,N
            XX(I,P(N+3)) = XA(I) + DE*(XA(I)-XX(I,P(N+1)))
 3200    CONTINUE
         FF(P(N+3)) = F(XX(1,P(N+3)))
         IF (FF(P(N+3)).LT.FF(P(N+2))) THEN
            CALL LSTINS(FF,N+3,P,N+1)
         ELSE
            CALL LSTINS(FF,N+2,P,N+1)
         ENDIF
      ELSE
         SHRINK= .TRUE.
         IF (FF(P(N+2)).LT.FF(P(N+1))) THEN
            DO 3300 I=1,N
               XX(I,P(N+3)) = XA(I) + DOC*(XA(I)-XX(I,P(N+1)))
 3300       CONTINUE
            FF(P(N+3)) = F(XX(1,P(N+3)))
            IF (FF(P(N+3)).LE.FF(P(N+2))) THEN
               CALL LSTINS(FF,N+3,P,N+1)
               SHRINK = .FALSE.
            ENDIF
         ELSE
            DO 3400 I=1,N
               XX(I,P(N+3)) = XA(I) + DIC*(XA(I)-XX(I,P(N+1)))
 3400       CONTINUE
            FF(P(N+3)) = F(XX(1,P(N+3)))
            IF (FF(P(N+3)).LE.FF(P(N+1))) THEN
               CALL LSTINS(FF,N+3,P,N+1)
               SHRINK = .FALSE.
            ENDIF
         ENDIF
         IF (SHRINK) THEN
            DO 3500 I=2,N+1
               DO 3510 J=1,N
                  XX(J,P(I)) = 0.5D0*(XX(J,P(I))+XX(J,P(1)))
 3510          CONTINUE
               FF(P(I)) = F(XX(1,P(I)))
               CALL LSTINS(FF,I,P,I-1)
 3500       CONTINUE
         ENDIF
      ENDIF
      IF (FF(P(1)).GE.FMIN.AND.FF(P(N+1)).GE.FMAX) THEN
         K = K+1
      ELSE
         K = 0
      ENDIF
      IF (FF(P(1)).LT.FMIN) THEN
         FMIN = FF(P(1))
         DO 3600 I=1,N
            XMIN(I) = XX(I,P(1))
            XMAX(I) = XX(I,P(1))
            DO 3610 J=2,N+1
               IF (XX(I,P(J)).LT.XMIN(I)) XMIN(I)=XX(I,P(J))
               IF (XX(I,P(J)).GT.XMAX(I)) XMAX(I)=XX(I,P(J))
 3610       CONTINUE
            XN = XMAX(I)-XMIN(I)
            IF (XN.GT.0D0) DX(I) = XN
 3600    CONTINUE
      ENDIF
      IF (FF(P(N+1)).LT.FMAX) FMAX = FF(P(N+1))
      GOTO 3000

 4000 DO 4100 I=1,N
         X(I) = XX(I,P(1))
 4100 CONTINUE
      NMMIN1 = FF(P(1))

      END
