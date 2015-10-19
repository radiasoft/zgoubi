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
      SUBROUTINE RMCOBJ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     *******************************************
C     READS DATA FOR OBJECT DEFINITION BY 'OBJET'
C     *******************************************
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      INCLUDE "C.OBJET.H"     ! COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      INTEGER FINSTR
      LOGICAL STRCON
      CHARACTER(132) TXT132
      PARAMETER (MXJ2=MXJ-2)
      PARAMETER (I0 = 0)

C----- BORO
      LINE = 1
      READ(NDAT,*,ERR=90,END=90) A(NOEL,1)
C----- KOBJ - may be of the form "KOBJ.K2"
      LINE = 2
      READ(NDAT,*,ERR=90,END=90) TXT132
      IF(STRCON(TXT132,'!',
     >                     IS)) TXT132 = TXT132(1:IS-1)
      A(NOEL,11) = I0
      IF(STRCON(TXT132,'.',
     >                      IS)) THEN
        READ(TXT132(1:IS-1),*,ERR=90,END=90) KOBJ
        READ(TXT132(IS+1:FINSTR(TXT132)),*,ERR=90,END=90) K2
        A(NOEL,11) = K2
      ELSE
        READ(TXT132,*,ERR=90,END=90) KOBJ
      ENDIF
      A(NOEL,10) = KOBJ
C----- IMAX
      LINE = 3
      READ(NDAT,*,ERR=90,END=90) A(NOEL,20)
C----- Type of p.d.f, KY, KT, KZ, KP, KX, KD
      LINE = 4
      READ(NDAT,*,ERR=90,END=90) (A(NOEL,I),I=30,30 + MXJ2)
C----- Central values Yo, To, Zo, Po, Xo, Do
      LINE = 5
      READ(NDAT,*,ERR=90,END=90) (A(NOEL,I),I=40,40 + MXJ2)

      GOTO (1,2,3) KOBJ
      CALL ENDJOB('KOBJ should be 1, 2 or 3. No such object KOBJ=',KOBJ)
 
 1    CONTINUE
C----- dY, dT, dZ,...
      LINE = 6
      READ(NDAT,*,ERR=90,END=90) (A(NOEL,I),I=50,50 + MXJ2)
C----- Cut-offs Nsy, Nst, Nsz,...
      LINE = 7
      READ(NDAT,*,ERR=90,END=90) (A(NOEL,I),I=60,60 + MXJ2)
C----- No, Co, C1, C2 C3
      LINE = 8
      READ(NDAT,*,ERR=90,END=90) (A(NOEL,I),I=70,74)
C----- IR1, IR2, IR3
      LINE = 9
      READ(NDAT,*,ERR=90,END=90) (A(NOEL,I),I=80,82)
      RETURN
 
 2    CONTINUE
C----- # of bars IY, IT , IZ,...
      LINE = 6
      READ(NDAT,*,ERR=90,END=90) (A(NOEL,I),I=50,50 + MXJ2)
C----- Distances between bars PY, PT, PZ,...
      LINE = 7
      READ(NDAT,*,ERR=90,END=90) (A(NOEL,I),I=60,60 + MXJ2)
C----- Widths of the bars dY, dT, dZ,...
      LINE = 8
      READ(NDAT,*,ERR=90,END=90) (A(NOEL,I),I=70,70 + MXJ2)
C----- Cut-offs Nsy, Nst, Nsz,...
      LINE = 9
      READ(NDAT,*,ERR=90,END=90) (A(NOEL,I),I=80,80 + MXJ2)
C----- No, Co, C1, C2 C3
      LINE = 10
      READ(NDAT,*,ERR=90,END=90) (A(NOEL,I),I=90,94)
C----- IR1, IR2, IR3
      LINE = 11
      READ(NDAT,*,ERR=90,END=90) (A(NOEL,I),I=95,97)
      RETURN
 
 3    CONTINUE
      LINE = 6
      READ(NDAT,*,ERR=90,END=90) (A(NOEL,I),I=50,53)
      IF(A(NOEL,53) .LT. 0.D0) THEN
        BACKSPACE(NDAT,ERR=99)
        READ(NDAT,*,ERR=90,END=90) (A(NOEL,I),I=50,54)
        A(NOEL,53) = -A(NOEL,53) 
      ELSE
        A(NOEL,54) = 0.D0
      ENDIF 
      LINE = 7
      READ(NDAT,*,ERR=90,END=90) (A(NOEL,I),I=60,63)
      IF(A(NOEL,63) .LT. 0.D0) THEN
        BACKSPACE(NDAT,ERR=99)
        READ(NDAT,*,ERR=90,END=90) (A(NOEL,I),I=60,64)
        A(NOEL,63) = -A(NOEL,63) 
      ELSE
        A(NOEL,64) = 0.D0
      ENDIF 
      LINE = 8
      READ(NDAT,*,ERR=90,END=90) (A(NOEL,I),I=70,73)
      IF(A(NOEL,73) .LT. 0.D0) THEN
        BACKSPACE(NDAT,ERR=99)
        READ(NDAT,*,ERR=90,END=90) (A(NOEL,I),I=70,74)
        A(NOEL,73) = -A(NOEL,73) 
      ELSE
        A(NOEL,74) = 0.D0
      ENDIF 
      LINE = 9
      READ(NDAT,*,ERR=90,END=90) (A(NOEL,I),I=80,82)
      RETURN

 99   WRITE(NRES,FMT='(5X,A)') 
     >'*** Pgm rmcobj. error upon BACKSPACE'

 90   CALL ENDJOB('*** Pgm rmcobj. Input data error, at line ',line)

      RETURN
      END
