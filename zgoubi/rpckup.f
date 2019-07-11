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
C  USA
C  -------
      SUBROUTINE RPCKUP
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'C.CDF.H'     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      PARAMETER (MPULAB=5)
      INCLUDE 'C.COT.H'     ! COMMON/COT/ PULAB(MPULAB)
      INCLUDE 'MXLD.H'
      INCLUDE 'C.DON.H'     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      PARAMETER (LBLSIZ=20)
      CHARACTER(LBLSIZ) PULAB

      CHARACTER(132) TXT
      INTEGER DEBSTR
      LOGICAL STRCON
      PARAMETER (KSIZ=10)
      CHARACTER(KSIZ) KLE

C----- NPU = 0 (OFF) or NPU > 0 (# of distinct lmnt LABEL's)
      LINE = 1
      READ(NDAT,FMT='(A)',ERR=99,END=99) TXT
      TXT = TXT(DEBSTR(TXT):LEN(TXT))
      IF    (STRCON(TXT,'!',
     >                      IS)) TXT = TXT(1:IS-1)

      READ(TXT,*,ERR=99,END=99) NPU

      IF(NPU .GT. MPULAB) CALL ENDJOB('Pgm rpckup. '
     >//'Too  many  PU families, maximum allowed is  ',MPULAB )

      A(NOEL,1) = NPU

      IF    (STRCON(TXT,'PRINT',
     >                          IS)) A(NOEL,2)=1.D0

C----- LABEL's
      LINE = LINE + 1
      READ(NDAT,FMT='(A)',ERR=99,END=99) TXT
      TXT = TXT(DEBSTR(TXT):LEN(TXT))
      IF(NPU.GE.1) READ(TXT,*,ERR=98) (PULAB(I),I=1,NPU)

      RETURN

 98   WRITE(6,FMT='(A,1X,I0)')
     >'Data error met while reading pickup list. Check # of PUs. '
     >//' Expected # of itmes in list is ',NPU
      WRITE(NRES ,*)
     >'Data error met while reading pickup list. Check # of PUs. '
     >//' Expected # of itmes in list is ',NPU
      GOTO 90

 99   WRITE(6,*)
     >  ' *** Execution stopped upon READ : invalid input in PICKUPS'
      WRITE(NRES ,*)
     >  ' *** Execution stopped upon READ : invalid input in PICKUPS'
      GOTO 90

 90   CONTINUE
      CALL ZGKLEY(
     >            KLE)
      CALL ENDJOB('*** Pgm rpckup, keyword '//KLE//' : '//
     >'input data error, at line #',LINE)
      RETURN
      END
