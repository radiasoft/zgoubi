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
      SUBROUTINE RTOSCA(IDIM,
     >                       ND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     -----------------------------------
C     Read input data for 'TOSCA' keyword
C     -----------------------------------
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
C      PARAMETER (LNTA=132) ; CHARACTER(LNTA) TA
C      PARAMETER (MXTA=45)
      INCLUDE "C.DONT.H"     ! COMMON/DONT/ TA(MXL,MXTA)

      CHARACTER(132) TXT
      PARAMETER (MXSTR=20)
      CHARACTER(LNTA) STRA(MXSTR)
      INTEGER DEBSTR
      LOGICAL STRCON
      LOGICAL ISNUM
      PARAMETER (KSIZ=10)
      CHARACTER(KSIZ) KLE
      PARAMETER (MXC = 4)

C     ... IC, IL
      LINE = 1
      READ(NDAT,*,ERR=99,END=98) A(NOEL,1),A(NOEL,2)

C----- BNORM & X-,Y-,Z-NORM
      DO I = 1, IDIM    ! For compatibility with earlier versions
        A(NOEL,10+I) = 1.D0
      ENDDO
      LINE = LINE + 1
      READ(NDAT,FMT='(A)',ERR=99,END=98) TXT
      IF( STRCON(TXT,'!',
     >                   IS)) TXT = TXT(DEBSTR(TXT):IS-1)
      CALL STRGET(TXT,4,
     >                  NSTR,STRA)
      IF(NSTR .GT. MXSTR) CALL ENDJOB
     >('Pgm rtosca. Prblm with normalization coefficients list : '
     >//' expected is a list of  1 to 4 data at line ',LINE)

      I = 1
      DOWHILE (I .LE. NSTR)
        IF(ISNUM(STRA(I)))   ! In case unexpectedly a ! is missing at that line...
     >  READ(STRA(I),*,ERR=99,END=98) A(NOEL,9+I)
        I = I + 1
      ENDDO

C----- TITLE - Start TITLE with FLIP to get map flipped (implemented with TOSCA... to
C                       be completed for others)
      LINE = LINE + 1
      READ(NDAT,200,ERR=99,END=98) TA(NOEL,1)
 200  FORMAT(A)
      IF(STRCON(TA(NOEL,1),'!',
     >      IS)) TA(NOEL,1) = TA(NOEL,1)(DEBSTR(TA(NOEL,1)):IS-1)

      ND = 50
      IF    (IDIM .EQ. 1) THEN
        LINE = LINE + 1
        READ(NDAT,*,ERR=99,END=98) A(NOEL,20)
        NFIC = 1
      ELSEIF(IDIM .EQ. 2) THEN
        LINE = LINE + 1
        READ(NDAT,FMT='(A)',ERR=99,END=98) TXT
        IF(STRCON(TXT,'!',
     >                    IS)) TXT = TXT(DEBSTR(TXT):IS-1)
        CALL STRGET(TXT,4,
     >                    NSTR,STRA)
        READ(STRA(1),*,ERR=99,END=98) NX
        READ(STRA(2),*,ERR=99,END=98) NY
        A(NOEL,20)=NX
        A(NOEL,21)=NY

C Old style only has NX, NY for 2D maps. New style gets KZMA here for further use.
        IF(NSTR.GE.3) THEN
          IF    (NSTR.EQ.3) THEN
            READ(STRA(3),*,ERR=99,END=98) AMOD
            NFIC = 1
          ELSEIF(NSTR.GE.4) THEN
            READ(STRA(3),*,ERR=99,END=98) KZMA
            READ(STRA(4),*,ERR=99,END=98) AMOD
            MOD=INT(AMOD)
            MOD2 = NINT (10 *  (AMOD -  DBLE(MOD)) )
            NFIC = MOD2
          ENDIF
        ELSE
          AMOD = 0.D0
          KZMA = 1
          NFIC = 1
        ENDIF
        A(NOEL,22)=AMOD

      ELSEIF(IDIM .EQ. 3) THEN
C------- TOSCA, 2-D or 3-D maps, either cartesian mesh (MOD.le.19), or polar mesh (MOD.ge.20)
        LINE = LINE + 1
        READ(NDAT,FMT='(A)',ERR=99,END=98) TXT
        IF(STRCON(TXT,'!',
     >                    IS)) TXT = TXT(DEBSTR(TXT):IS-1)
        CALL STRGET(TXT,4,
     >                    NSTR,STRA)
        READ(STRA(1),*,ERR=99,END=98) NX
        READ(STRA(2),*,ERR=99,END=98) NY
        IF(NSTR.GE.3) THEN
          READ(STRA(3),*,ERR=99,END=98) KZMA
          IF(NSTR.GE.4) THEN
            READ(STRA(4),*,ERR=99,END=98) AMOD
          ELSE
            AMOD = 0
          ENDIF
        ELSE
          KZMA = 1
C------- To ensure compatibility with version 3 of Zgoubi
          AMOD = 0.D0
        ENDIF

        A(NOEL,20)=NX
        A(NOEL,21)=NY
        A(NOEL,22)=KZMA
        A(NOEL,23)=AMOD
        MOD=INT(AMOD)
        MOD2 = NINT (10 *  (AMOD -  DBLE(MOD)) )

        IF(KZMA .EQ. 1) THEN
          IDIM = 2
        ELSE
          IDIM = 3
        ENDIF

C--------------------------------------------------------
        IF    (MOD .LE. 19) THEN
C--------- Cartesian mesh
          KART = 1
          IF    (MOD .EQ. 0) THEN
C----------- 3-D map generated by  symmetrysing wrt horizontal plane
            NFIC = (KZMA/2) + 1
          ELSEIF(MOD .EQ. 3) THEN
C--------- AGS mid-plane 2-D maps
            NFIC = 1
            IF(MOD2.EQ.1) READ(TXT,*,ERR=99,END=98)
     >      IDUM1,IDUM2,IDUM3,DUM4,DB1,DB2
            A(NOEL,24)=DB1
            A(NOEL,25)=DB2
          ELSEIF(MOD .EQ.10) THEN
C----------- 3-D map generated by  symmetrysing wrt horizontal plane
            NFIC = (KZMA/2) + 1
          ELSEIF(MOD .EQ. 1) THEN
C----------- no symmetrization, 3-D map generated by taking 2-D maps as is
            NFIC = KZMA
          ELSEIF(MOD .EQ. 12) THEN
C--------- The all 3D map is contained in a single file
            NFIC = 1
          ELSEIF(MOD .EQ. 15) THEN
C--------- There are 'MOD2' files, they will be combined linearly
C          Each single file contains the all 3D volume
            NFIC = MOD2
            IF(4+NFIC .GT. MXSTR) THEN
              WRITE(ABS(NRES),*) 'At input data line ',LINE
              CALL ENDJOB('Pgm rtosca. Too many'
     >        //'  field map scaling factors. Max is ',MXC)
            ENDIF

            DO J = 1, MXC
              A(NOEL,23+J) = 1.D0
            ENDDO

            CALL STRGET(TXT,4+NFIC,
     >                           IDUM,STRA)
            DO J = 1, NFIC
              READ(STRA(4+j),*,ERR=99,END=98) A(NOEL,23+J)
            ENDDO

          ELSEIF(MOD .EQ. 16) THEN
C--------- There are 'MOD2' files.
C          Each single file contains the all 3D volume
C          Their field contributions at particle location will add
            NFIC = MOD2
            IF(4+NFIC .GT. MXSTR) THEN
              WRITE(ABS(NRES),*) 'At input data line ',LINE
              CALL ENDJOB('Pgm rtosca. Too many'
     >        //'  field map scaling factors. Max is ',MXC)
            ENDIF

            DO J = 1, MXC
              A(NOEL,23+J) = 1.D0
            ENDDO

            CALL STRGET(TXT,4+NFIC,
     >                           IDUM,STRA)
            DO J = 1, NFIC
              READ(STRA(4+j),*,ERR=99,END=98) A(NOEL,23+J)
            ENDDO

            ND = ND + 10*MOD2

          ELSE
            NFIC = 1
          ENDIF

        ELSEIF(MOD .GE. 20) THEN
C--------- Cylindrical mesh. Axis is Z. 3D maps are contained in a single file
C--------- Some MOD values allow superimposing maps.
          KART = 2
          IF(MOD2 .LE. 1) THEN
            NFIC = 1
          ELSE
            NFIC = MOD2
          ENDIF

          IF(4+NFIC .GT. MXSTR) THEN
            WRITE(ABS(NRES),*) 'At input data line ',LINE
            CALL ENDJOB('Pgm rtosca. Too many'
     >      //'  field map scaling factors. Max is ',MXC)
          ENDIF
          CALL STRGET(TXT,4+NFIC,
     >                           ITMP,STRA)

          DO J = 1, ITMP-4
            READ(STRA(4+J),*,ERR=99,END=98) A(NOEL,23+J)
          ENDDO
          DO J = ITMP-4+1, NFIC
            A(NOEL,23+J)  = 1.D0
          ENDDO

        ENDIF
      ENDIF

C----- MAP FILE NAME(S)
        DO IFIC=1,NFIC
          LINE = LINE + 1
          READ(NDAT,FMT='(A)',ERR=99,END=98) TXT
          IF(STRCON(TXT,'!',
     >                    IS)) TXT = TXT(DEBSTR(TXT):IS-1)
          CALL STRGET(TXT,6,
     >                      IDUM,STRA)
          TA(NOEL,1+IFIC) = STRA(1)
          IF(MOD .EQ. 16) THEN
            I5 = 5
            IF(IDUM .LT. 6) CALL ENDJOB('Pgm rtosca. A list of 5 '
     >      //'data expected for positioning, line ',LINE)
            DO I = 1, MIN(I5,IDUM)
              READ(STRA(I+1),*,ERR=90,END=90) A(NOEL,20+IFIC*10+I-1)
            ENDDO
          ENDIF
        ENDDO

C----- DROITE(S) DE COUPURE (IA=-1, 1, 2 or 3)
      LINE = LINE + 1
      READ(NDAT,FMT='(A)',ERR=99,END=98) TXT
      IF(STRCON(TXT,'!',
     >                    IS)) TXT = TXT(DEBSTR(TXT):IS-1)
      READ(TXT,*,ERR=99,END=98) IA
      JD = ND - 20
      IF(IA.GE.1) READ(TXT,*,ERR=99,END=98)
     >IA,(A(NOEL,I),I=JD+1,JD+3*ABS(IA))

      IF(IA .GT. 3) CALL ENDJOB(
     >  '*** Error, Pgm rtosca -> input data IDRT should be < ',4)

      A(NOEL,JD) = IA
C     ... IRD
      LINE = LINE + 1
      JD = ND - 10
      READ(NDAT,*,ERR=99,END=98) A(NOEL,JD)
C     ... XPAS
      LINE = LINE + 1
      READ(NDAT,*,ERR=99,END=98) A(NOEL,ND)

      IF(KART .EQ. 1) THEN
C       ... Cartesian map frame
C           KP, XCE, YCE, ALE
        LINE = LINE + 1
        READ(NDAT,*,ERR=99,END=98) IA,(A(NOEL,I),I=ND+10+1,ND+10+3)
        A(NOEL,ND+10) = IA
      ELSEIF(KART .EQ. 2) THEN
C       ... Polar map frame
        LINE = LINE + 1
        READ(NDAT,*,ERR=99,END=98) KP
        A(NOEL,ND+10) = KP
        IF( KP .EQ. 2 ) THEN
          LINE = LINE + 1
          READ(NDAT,*,ERR=99,END=98) (A(NOEL,I),I=ND+20,ND+20+3)
        ELSE
          LINE = LINE + 1
          READ(NDAT,*,ERR=99,END=98) A(NOEL,ND+20)
        ENDIF
      ENDIF

      RETURN

 99   WRITE(6,*)
     >  ' *** Execution stopped upon READ ERR : invalid input in rtosca'
      WRITE(NRES ,*)
     >  ' *** Execution stopped upon READ ERR : invalid input in rtosca'
      GOTO 90

 98   WRITE(6,*)
     >  ' *** Execution stopped upon READ END : invalid input in rtosca'
      WRITE(NRES ,*)
     >  ' *** Execution stopped upon READ END : invalid input in rtosca'
      IF(LINE.EQ.2) WRITE(NRES,FMT='(A,I0,A)')
     >'Pgm rtosca. Prblm with normalization coefficients list : '
     >//' expected is a list of  1 to ',MXC,' data at line ',LINE
      IF(LINE.EQ.6) WRITE(NRES,*)
     >'Expecting more data after MOD=3, at line ',LINE

 90   CONTINUE
      CALL ZGKLEY(
     >            KLE)
      CALL ENDJOB('*** Pgm rtosca, keyword '//KLE//' : '//
     >'input data error, at line #',LINE)
      RETURN
      END
