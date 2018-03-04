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
      SUBROUTINE MCOBJ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     **************************
C      CONSTITUTION DE L'OBJET PAR TIRAGE MONTE-CARLO
C      DOUBLE TIRAGE AVEC MELANGE DES CARTES
C     **************************
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
      INCLUDE "C.OBJET.H"     ! COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      INCLUDE "C.PTICUL_2.H"     ! COMMON/PTICUL/ AAM,Q,G,TO
      INCLUDE "C.REBELO.H"   ! COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,HDPRF,DP,QBR,BRI
 
      DIMENSION CENTRE(MXJ)
      PARAMETER (MXJ1=MXJ-1)

      CHARACTER(15)  KTIR(MXJ)
      CHARACTER(14) KOUV
      LOGICAL CINE
      CHARACTER(80) TXT
      LOGICAL FITING
 
      DATA KTIR, KOUV / MXJ*' ', ' ' / 
      DATA CENTRE / MXJ*0.D0 / 

      AMQLU(1) = .FALSE.
      AMQLU(2) = .FALSE.
      AMQLU(3) = .FALSE.
      AMQLU(4) = .FALSE.
      AMQLU(5) = .FALSE.
      PABSLU = .FALSE.

C----- Magnetic rigidity  (Kg*cm)
      BORO = A(NOEL,1)

      CALL REBELR(
     >            KREB3,KREB31,KDUM)
      IF(KREB3 .EQ. 99) THEN
C------- SET TO 99 IN SBR REBELOTE
        IF(NRES.GT.0) WRITE(NRES,133)
 133    FORMAT(//,15X,'Final  coordinates  of  previous  run',1X
     >  ,' taken  as  initial  coordinates')
        IF(KREB31 .NE. 0) THEN
C--------- add new beamlet next to the previous one(s), e.g. for multiturn injection
          IF(IPASS .LE. 1+KREB31) THEN
            IF(NRES.GT.0) WRITE(NRES,FMT='(
     >           15X,'' Injection run ; new beamlet launched'',/)') 
          ELSE
            GOTO 99
          ENDIF
        ELSE
          GOTO 99
        ENDIF
      ENDIF
 
      CALL FITSTA(5,
     >              FITING)

C FM, Mar 2018      
C      IF(FITING) DPREF = 1.D0
      IF(FITING) THEN
        DPREF = 0.D0
        HDPRF = 1.D0
      ENDIF

C------ Type of support : 
C      KOBJ=1 : window
C      KOBJ=2 : grid
C      KOBJ=3 : ellipses
      KOBJ = NINT(A(NOEL,10))
      KOBJ2 = NINT(A(NOEL,11))

      IF(KOBJ2.EQ.1) THEN
         KNRM = 1
      ELSE
         KNRM = 0
      ENDIF

      IF(KOBJ .EQ. 1) THEN
        KOUV='Window'
      ELSEIF(KOBJ .EQ. 2) THEN
        KOUV='Grid'
      ELSEIF(KOBJ .EQ. 3) THEN
        KOUV='Ellipse'
        IF(KOBJ2.EQ.1) KNRM = 1
      ENDIF

      IMAX  = NINT(A(NOEL,20))
      IF(KREB31 .EQ. 0) THEN
        IMI  = 1 
        IMA = IMAX
        IF(IMAX .GT. MXT) GOTO 98
      ELSE
C------- Multiturn injection
        IF(IMAX*(KREB31+1) .GT. MXT) GOTO 98
        IMI  = 1 + IMAX*(IPASS-1)
        IMA = IMAX*IPASS
        IMAX=IMA
      ENDIF

      IDMAX=1
      IMAXT=IMAX

      IF(KOUV.EQ.'Ellipse') THEN
        DO 3 I= 2, 6, 2
          KPDF = NINT(A(NOEL,30+I-1))
          IF(KPDF .EQ. 1) THEN
            KTIR(I) = 'Uniform'
          ELSEIF(KPDF .EQ. 2) THEN
            KTIR(I) = 'Gaussian'
          ELSEIF(KPDF .EQ. 3) THEN
            KTIR(I) = 'Parabolic'
          ENDIF
 3      CONTINUE
      ELSE
        DO 2 I= 1, MXJ1-1
          KPDF = NINT(A(NOEL,30+I-1))
          IF(KPDF .EQ. 1) THEN
            KTIR(I+1) = 'Uniform'
          ELSEIF(KPDF .EQ. 2) THEN
            KTIR(I+1) = 'Gaussian'
          ELSEIF(KPDF .EQ. 3) THEN
            KTIR(I+1) = 'Parabolic'
          ENDIF
 2      CONTINUE
        CINE=.FALSE.
        KPDF = NINT(A(NOEL,30+MXJ1-1))
        IF(KPDF .EQ. 1) THEN
          KTIR(1) = 'Uniform'
        ELSEIF(KPDF .EQ. 2) THEN
          KTIR(1) = 'Exponentl'
        ELSEIF(KPDF .EQ. 3) THEN
          KTIR(1) = 'Kinematic'
          CINE=.TRUE.
        ENDIF
      ENDIF

      IF(NRES.GT.0) WRITE(NRES,103) BORO,IMAX
 103  FORMAT(25X,' Reference  magnetic rigidity ='
     >,F15.3,' KG*CM'
     >   ,//,25X,' Object  built  up  of ',I7,'  particles')
 
C     .... Central values of the sorting
C       2-Y, 3-T, 4-Z, 5-P, 6-X, 1-D
      CENTRE(2) = A(NOEL,40)
      CENTRE(3) = A(NOEL,41)
      CENTRE(4) = A(NOEL,42)
      CENTRE(5) = A(NOEL,43)
      CENTRE(6) = A(NOEL,44)
      CENTRE(1) = A(NOEL,45)

      IF(KOBJ.LE.2) THEN
        CALL MCOBJ1(KTIR,KOUV,CINE,CENTRE)
      ELSEIF(KOBJ.EQ.3) THEN
        CALL MCOBJ3(KTIR,CENTRE,KNRM,IMI,IMA)
      ENDIF

 99   CONTINUE
      DO 991 I=1,IMAX
        IF(IEX(I) .LT. -1) GOTO 991 
        IF(F(1,I) .EQ. 0.D0) THEN
          TXT = '       zero momentum found.'
          CALL OBJERR(ABS(NRES),1,MXT,TXT)
          CALL ENDJOB('*** Error, SBR MCOBJ -> momentum  cannot  be ',0)
        ENDIF
        AMQ(1,I) = AAM
        AMQ(2,I) = Q
 991  CONTINUE 
      IF(IPASS .EQ. 1) CALL CNTMXW(IMAX)
      RETURN

 98   TXT = '    Too many particles'
      CALL OBJERR(ABS(NRES),2,MXT,TXT)
      CALL ENDJOB('Too many particles, max is ',MXT/(KREB31+1))
      RETURN

      END
