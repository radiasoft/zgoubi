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
      SUBROUTINE RSVDOC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ***************************************
C     READS DATA FOR FIT PROCEDURE WITH 'FIT'
C     ***************************************
      INCLUDE 'MXLD.H'
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL

      PARAMETER(I300=300)
      CHARACTER(I300) TXT300, TXTMP

      INTEGER DEBSTR
      LOGICAL STRCON
      PARAMETER (I2=2)
      CHARACTER(I300) STRA(I2)

      PARAMETER (LBLSIZ=20)
      PARAMETER (MPULAB=5)
      parameter (IMON=MPULAB/3)
      parameter(mxpuh =IMON, mxpuv =IMON, mxpuhv =IMON)
      CHARACTER(LBLSIZ) HPNA(mxpuh), VPNA(mxpuv), HVPNA(mxpuhv)
      parameter (mxcoh=5, mxcov=5)
      CHARACTER(LBLSIZ) HCNA(mxcoh), VCNA(mxcov)

      LINE = LINE + 1
      READ(NDAT,FMT='(A)',ERR=98,END=98) TXT300
      IF(STRCON(TXT300,'!',
     >                   II))
     >TXT300 = TXT300(DEBSTR(TXT300):II-1)
      IF(STRCON(TXT300,'PUH',
     >                       II)) THEN
        TXTMP = TXT300(II:)
         IF(STRCON(TXTMP,'{',
     >                       II)) THEN
          TXTMP = TXTMP(II+1:)
          IF(STRCON(TXTMP,'}',
     >                        JJ)) THEN

            CALL STRGET(TXTMP(1:jj-1),MXPUH,
     >                                      NSR,STRA)
            IF(NSR.GT.mxpuH) THEN
              IF(NRES .NE. 0) WRITE(ABS(NRES),*)
     >        'SBR RREBEL. Too many H-PUs.'
     >        //' Maximum allowed is ',mxpuh
              GOTO 98
            ENDIF
            DO I = 1, NSR
              READ(STRA(I),*,ERR=98) HPNA(I) !  H PU name list
            ENDDO
          ENDIF
        ENDIF
      ENDIF
      IF(STRCON(TXT300,'PUV',
     >                       II)) THEN
        TXTMP = TXT300(II:)
         IF(STRCON(TXTMP,'{',
     >                       II)) THEN
          TXTMP = TXTMP(II+1:)
          IF(STRCON(TXTMP,'}',
     >                         JJ)) THEN

            CALL STRGET(TXTMP(1:jj-1),mxpuv,
     >                                 NSR,STRA)
            IF(NSR.GT.mxpuV) THEN
              IF(NRES .NE. 0) WRITE(ABS(NRES),*)
     >        'SBR RREBEL. Too many V-PUs.'
     >        //' Maximum allowed is ',mxpuv
              GOTO 98
            ENDIF
            DO I = 1, NSR
              READ(STRA(I),*,ERR=98) VPNA(I) !  V PU name list
            ENDDO
          ENDIF
        ENDIF
      ENDIF
      IF(STRCON(TXT300,'PUHV',
     >                       II)) THEN
        TXTMP = TXT300(II:)
         IF(STRCON(TXTMP,'{',
     >                       II)) THEN
          TXTMP = TXTMP(II+1:)
          IF(STRCON(TXTMP,'}',
     >                        JJ)) THEN
            CALL STRGET(TXTMP(1:jj-1),mxpuhv,
     >                                   NSR,STRA)
            IF(NSR.GT.mxpuHV) THEN
              IF(NRES .NE. 0) WRITE(ABS(NRES),*)
     >        'SBR RREBEL. Too many HV-PUs.'
     >        //' Maximum allowed is ',mxpuhv
              GOTO 98
            ENDIF
            DO I = 1, NSR
              READ(STRA(I),*,ERR=98) HVPNA(I) !  HV PU name list
            ENDDO
          ENDIF
        ENDIF
      ENDIF

      LINE = LINE + 1
      READ(NDAT,FMT='(A)',ERR=98,END=98) TXT300
      IF(STRCON(TXT300,'!',
     >                   II))
     >TXT300 = TXT300(DEBSTR(TXT300):II-1)
      IF(STRCON(TXT300,'CRH',
     >                       II)) THEN
        TXTMP = TXT300(II:)
        IF(STRCON(TXTMP,'{',
     >                       II)) THEN
          TXTMP = TXTMP(II+1:)
          IF(STRCON(TXTMP,'}',
     >                        JJ)) THEN
            CALL STRGET(TXTMP(1:jj-1),mxcoh,
     >                           NSRH,STRA)
            IF(NSRH.GT.mxcoh) THEN
              IF(NRES .NE. 0) WRITE(ABS(NRES),*)
     >        'SBR RREBEL. Too many H-correctors.'
              GOTO 98
            ENDIF
            DO I = 1, NSRH
              READ(STRA(I),*,ERR=98) HCNA(I)    !  H-corr name list 
            ENDDO
          ENDIF
        ENDIF
      ENDIF
      IF(STRCON(TXT300,'CRV',
     >                       II)) THEN
        TXTMP = TXT300(II:)
        IF(STRCON(TXTMP,'{',
     >                       II)) THEN
          TXTMP = TXTMP(II+1:)
          IF(STRCON(TXTMP,'}',
     >                         JJ)) THEN
            CALL STRGET(TXTMP(1:jj-1),mxcov,
     >                           NSRV,STRA)
            IF(NSRV.GT.mxcov) THEN
              IF(NRES .NE. 0) WRITE(ABS(NRES),*)
     >        'SBR RREBEL. Too many V-correctors.'
              GOTO 98
            ENDIF
            DO I = 1, NSRV
              READ(STRA(I),*,ERR=98) VCNA(I)    !  V-corr name list 
            ENDDO
          ENDIF
        ENDIF
      ENDIF

      LINE = LINE + 1
      READ(NDAT,*,ERR=98) (HKI,i=1,NSRH), (VKI,i=1,NSRv)
      A(NOEL,10) = HKI                 !  H-corrector kick value
      A(NOEL,20) = VKI                 !  V-corrector kick value

      CALL SVDOC2(HPNA,VPNA,HVPNA,HCNA,VCNA)
      NOELA = 1
      NOELB = NOEL
      CALL REBLT6(NOELA, NOELB)

      RETURN

 98   CONTINUE
      CALL ENDJOB('*** Pgm robjet, keyword ''SVDOC'' : '// 
     >'input data error, at line #',LINE)
      RETURN
 
      END
