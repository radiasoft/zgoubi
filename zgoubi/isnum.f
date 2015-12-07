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
      FUNCTION ISNUM(S)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL ISNUM
      CHARACTER(*) S
      INTEGER I, K, N
      INTEGER DEBSTR, FINSTR
      LOGICAL ALRDY

      ALRDY = .FALSE. 

      ISNUM = .FALSE.
      S = S(DEBSTR(S):FINSTR(S))
      N = LEN(S(DEBSTR(S):FINSTR(S)))

      I = 1
      K = ICHAR(S(I:I))
      IF (K .LT. 48) THEN
        IF (.NOT.(K.EQ.43 .OR. K.EQ.45 .OR. K.EQ.46)) THEN   ! + OR - OR .
          ISNUM = .FALSE.
          GOTO 88
        ELSEIF(K.EQ.46) THEN
          ALRDY = .TRUE.
        ENDIF
      ELSEIF (K .GT. 57) THEN     ! 49-57 : 1-9
        ISNUM = .FALSE.
        GOTO 88
      ENDIF

      DO I = 2, N
        K = ICHAR(S(I:I))
        IF (K .LT. 48) THEN
          IF (K.NE.46) THEN    !  46 : DOT
            ISNUM = .FALSE.
            GOTO 88
          ELSE
            IF(ALRDY) THEN
              ISNUM = .FALSE.
              GOTO 88
            ELSE
              ALRDY = .TRUE.
            ENDIF 
          ENDIF
        ELSEIF (K .GT. 57) THEN
          IF(  S(I:I).EQ.'e' .OR. S(I:I).EQ.'E' 
     >    .OR. S(I:I).EQ.'d' .OR. S(I:I).EQ.'D') THEN
            GOTO 77
          ELSE
            ISNUM = .FALSE.
            GOTO 88
          ENDIF
        ENDIF
C      write(*,*) ' n, i, char, isnum ',n, i, ICHAR(S(I:I)),s(i:i)
      ENDDO
      ISNUM = .TRUE.
      goto 88

 77   CONTINUE
      I = I+1
      K = ICHAR(S(I:I))
C      write(*,*) ' n, i, char, isnum ',n, i, ICHAR(S(I:I)),s(i:i)
      IF (K .LT. 48) THEN
        IF (.NOT.(K.EQ.43 .OR. K.EQ.45)) THEN   ! + OR - 
          ISNUM = .FALSE.
          GOTO 88
        ENDIF
      ELSEIF (K .GT. 57) THEN
        ISNUM = .FALSE.
        GOTO 88
      ENDIF
      J = I + 1
      DO I = J, N
        K = ICHAR(S(I:I))
C      write(*,*) ' n, i, char, isnum ',n, i, ICHAR(S(I:I)),s(i:i)
        IF (K .LT. 48) THEN
            ISNUM = .FALSE.
            GOTO 88
        ELSEIF (K .GT. 57) THEN
            ISNUM = .FALSE.
            GOTO 88
        ENDIF
C      write(*,*) ' n, i, char, isnum ',n, i, ICHAR(S(I:I)),s(i:i)
      ENDDO
      ISNUM = .TRUE.


 88   CONTINUE

      RETURN
      END 
