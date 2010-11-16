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
      SUBROUTINE RPARTI(NDAT,NOEL,
     >                            A)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'MXLD.H'
      DIMENSION A(MXL,*)
      CHARACTER*132 TXT
      CHARACTER STRA(6)*20
      INTEGER DEBSTR, FINSTR

      READ(NDAT,FMT='(A)') TXT
      IFT = FINSTR(TXT)
      TXT = TXT(DEBSTR(TXT):IFT)
      CALL STRGET(TXT,6,
     >                      IDUM,STRA)
      IFIN = FINSTR(STRA(1))
      IDT = IFIN+1
      IF(STRA(1)(1:1) .EQ. '{') THEN
C------- Read 2 masses
        CALL PARTII('MASS_CODE_1')
        STRA(1)=STRA(1)(2:IFIN) 
c        IFIN = FINSTR(STRA(1))
c        IF(STRA(1)(IFIN:IFIN) .EQ. '}') THEN
c          STRA(1)=STRA(1)(1:IFIN-1) 
c        ELSEIF(STRA(1)(IFIN-1:IFIN) .EQ. '},') THEN
c          STRA(1)=STRA(1)(1:IFIN-2) 
c        ELSE
c          GOTO 99
c        ENDIF
        STRA(2)=STRA(2)(1:IFIN-1) 

        READ(STRA(1),*) A(NOEL,1)
        READ(STRA(2),*) A(NOEL,2)

        TXT = TXT(IDT:IFT)
        IFT = FINSTR(TXT)
        IF(TXT(1:1) .EQ. ',') TXT = TXT(2:IFT)
        
C------- Read Q, G, tau, dum
        READ(stra(3),*) A(NOEL,3)
c          write(*,*) stra(1), stra(2), a(noel,1),a(noel,2),a(noel,3)
c             stop
C        READ(TXT,*) A(NOEL,3),A(NOEL,4),A(NOEL,5),A(NOEL,6)
      ELSE
C------- Old method 
        CALL PARTII('NONE')
        READ(TXT,*) (A(NOEL,I),I=1,5)
      ENDIF
      RETURN
C 99   STOP ' *** DATA ERROR : in PARTICUL, while reading M1, M2 ***'
      END
