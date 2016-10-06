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
      SUBROUTINE RPARTI(NDAT,NOEL,
     >                            A)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'MXLD.H'
      DIMENSION A(MXL,*)
      INCLUDE "C.DONT.H"     ! COMMON/DONT/ TA(MXL,MXTA)
      INCLUDE 'C.PDATA.H'      ! COMMON /PDATA/ QELM,AMLEC,GLEC,AMPRO,GPRO,AMMU,GMU,TAUMU
      CHARACTER(132) TXT
      CHARACTER(80) STRA(6)
      INTEGER DEBSTR, FINSTR
      LOGICAL STRCON

      READ(NDAT,FMT='(A)') TXT
      IF(STRCON(TXT,'!',
     >                  II)) TXT = TXT(DEBSTR(TXT):II-1)
      IFT = FINSTR(TXT)
      TXT = TXT(DEBSTR(TXT):IFT)
      CALL STRGET(TXT,6,
     >                  IDUM,STRA)
      IFIN = FINSTR(STRA(1))
      IDT = IFIN+1
      IF(STRA(1)(1:1) .EQ. '{') THEN
C------- Read 2 masses
        CALL PARTII('MASS_CODE_1')
        STRA(1)=STRA(1)(2:IFIN) 
        STRA(2)=STRA(2)(1:IFIN-1) 

        READ(STRA(1),*) A(NOEL,1)
        READ(STRA(2),*) A(NOEL,2)

        TXT = TXT(IDT:IFT)
        IFT = FINSTR(TXT)
        IF(TXT(1:1) .EQ. ',') TXT = TXT(2:IFT)
        
C------- Read Q, G, tau, dum
        READ(stra(3),*) A(NOEL,3)
c        write(*,*) stra(1), stra(2), a(noel,1),a(noel,2),a(noel,3)
c           stop
C        READ(TXT,*) A(NOEL,3),A(NOEL,4),A(NOEL,5),A(NOEL,6)
      ELSE
C------- Default method
        CALL PARTII('NONE')
        IF    (STRA(1) .EQ. 'ELECTRON') THEN
          A(NOEL,1) = AMLEC
          A(NOEL,2) = -QE0
          A(NOEL,3) = GLEC
          A(NOEL,4) = 1D99
        ELSEIF(STRA(1) .EQ. 'MUON+') THEN
          A(NOEL,1) = AMMU
          A(NOEL,2) = QE0
          A(NOEL,3) = GMU
          A(NOEL,4) = TAUMU
        ELSEIF(STRA(1) .EQ. 'MUON-') THEN
          A(NOEL,1) = AMMU
          A(NOEL,2) = -QE0
          A(NOEL,3) = GMU
          A(NOEL,4) = TAUMU / 1.000024
        ELSEIF(STRA(1) .EQ. 'PROTON') THEN
          A(NOEL,1) = AMPRO
          A(NOEL,2) = QE0
          A(NOEL,3) = GPRO
          A(NOEL,4) = 1D99
        ELSE
          READ(TXT,*) (A(NOEL,I),I=1,5)
          STRA(1) = ' ' 
        ENDIF
        TA(NOEL,1) = STRA(1)
        A(NOEL,5) = 0.d0        
      ENDIF
      RETURN
C 99   STOP ' *** DATA ERROR : in PARTICUL, while reading M1, M2 ***'
      END
