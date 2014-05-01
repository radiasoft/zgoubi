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
      SUBROUTINE MATIMP(R,F0,YNU,ZNU,CMUY,CMUZ,NMAIL,PRDIC,IREF) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL PRDIC
      DIMENSION F0(6,6)
      DIMENSION R(6,*) , T(6,6,*)
      DIMENSION  T3(5,*)

      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU

      LOGICAL KWRI, KWRMAT, IDLUNI
      CHARACTER FNAME*17
      LOGICAL EXS, OPN, OK
      SAVE KWRMAT
      CHARACTER*15 TXTYNU, TXTZNU

      DATA FNAME / 'zgoubi.MATRIX.out' /
  
      DETY=R(1,1)*R(2,2)-R(1,2)*R(2,1)
      DETZ=R(3,3)*R(4,4)-R(3,4)*R(4,3)
      RIJ = R(2,2)
      IF(RIJ .EQ. 0.D0) RIJ = 1.D-10
      SFH = - R(1,2)/RIJ
      RIJ = R(4,4)
      IF(RIJ .EQ. 0.D0) RIJ = 1.D-10
      SFZ = - R(3,4)/RIJ
 
      I=1
      IF(NRES.GT.0) THEN
        WRITE(NRES,103) I
 103    FORMAT(//,18X,'TRANSFER  MATRIX  ORDRE',I3,'  (MKSA units)',/)
        WRITE(NRES,104) (( R(IA,IB) , IB=1,6) , IA=1,6)
 104    FORMAT(6X,1P,6G16.6)
        WRITE(NRES,112) DETY-1.D0,DETZ-1.D0
112     FORMAT(/,10X,'DetY-1 = ',F18.10,',',4X,'DetZ-1 = ',F18.10)
        WRITE(NRES,FMT='(/,10X,''R12=0 at '',G12.4,'' m, '',7X, 
     >                       ''R34=0 at '',G12.4,'' m'')') SFH,SFZ
      ENDIF

      CALL SYMPL(R)

      IF(NMAIL.LE.0) THEN
c        IF(NRES.GT.0) 
c     >  WRITE(NRES,*) ' Pgm matimp. NUMBER OF PERIODS IS ERRONEOUS !'
c        RETURN
      ELSE
        IF(NRES.GT.0) WRITE(NRES,106) NMAIL
 106    FORMAT(//,15X,' TWISS  parameters,  periodicity  of',
     >         I4,'  is  assumed ',/
     >           ,35X,' - UNCOUPLED -')
      ENDIF

      TXTYNU = '-99999.'
      TXTZNU = '-99999.'

      IF(PRDIC) THEN
        IF(NRES.GT.0) THEN
          WRITE(NRES,113)
 113      FORMAT(/,6X,
     >    ' Beam  matrix  (beta/-alpha/-alpha/gamma)',
     >    ' and  periodic  dispersion  (MKSA units)',/)
          WRITE(NRES,114) (( F0(IA,IB) , IB=1,6) , IA=1,6)
 114      FORMAT(6X,6F13.6)
          WRITE(NRES,FMT='(/,35X,''Betatron  tunes'',/)') 

c         IF    (ABS(CMUY).LT.1.D0 .AND. ABS(CMUZ).LT.1.D0) THEN
c           WRITE(TXTYNU,FMT='(G15.8)') YNU
c           WRITE(TXTZNU,FMT='(G15.8)') ZNU
c         ELSEIF(ABS(CMUY).LT.1.D0 .OR. ABS(CMUZ).LT.1.D0) THEN
           IF(CMUY*CMUY .LT. 1.D0) THEN
             WRITE(TXTYNU,FMT='(G15.8)') YNU
           ELSE
             WRITE(TXTYNU,FMT='(A)') 'undefined'
           ENDIF
           IF(CMUZ*CMUZ .LT. 1.D0) THEN
             WRITE(TXTZNU,FMT='(G15.8)') ZNU
           ELSE
             WRITE(TXTZNU,FMT='(A)') 'undefined'
           ENDIF
c         ENDIF

         WRITE(NRES,FMT='(15X,2(5X,A,A))') 
     >            'NU_Y = ', TXTYNU, 'NU_Z = ', TXTZNU

        ENDIF
      ENDIF

      IF(KWRMAT) THEN
        INQUIRE(FILE=FNAME,EXIST=EXS,OPENED=OPN,IOSTAT=I)
        IF(EXS) THEN
          IF(.NOT. OPN) THEN
            OK = IDLUNI(
     >                  LNWRT)
            OPEN(UNIT=LNWRT, FILE=FNAME, status='OLD',ERR=96)
            WRITE(LNWRT,FMT='(A)') 
     >      '% R11 R12 R13 R14 R15 R16, R21 R22 R23 R24 ... R65 R66, '
     >      //' ALFY, BETY, ALFZ, BETZ, DY, DYP, DZ, DZP, PHIY, PHIZ,'
     >      //' F(1,IREF),F(2,IREF),F(3,IREF),F(4,IREF),F(5,IREF)'
     >      //' F(6,IREF),F(7,IREF), CMUY, CMUZ, QY, QZ, XCE, YCE, ALE'
c     >      '% R11 R12 R13 R14,  R21 R22 R23 R24, R31 ... R43 R44, '
c     >      //' R16, R26, R36, R46,  R51 R52 R53 R54 R55,  R66, '
c     >      //' ALFY, BETY, ALFZ, BETZ, DY, DYP, DZ, DZP, PHIY, PHIZ,'
c     >      //' F(1,IREF),F(2,IREF),F(3,IREF),F(4,IREF),F(5,IREF)'
            WRITE(LNWRT,FMT='(A)') 
     >     '% 1   2   3   4   5   6    7   8   9   10  ... 35  36   '
     >     //' 37    38    39    40    41  42   43  44   45    46   '
     >     //' 47        48        49        50        51       '
     >     //' 52        53         54    55    56  57  58   59   60'
          ENDIF
        ELSE
          OK = IDLUNI(
     >                LNWRT)
          OPEN(UNIT=LNWRT, FILE=FNAME, status='NEW',ERR=96)
          WRITE(LNWRT,FMT='(A)') 
     >     '% R11 R12 R13 R14 R15 R16, R21 R22 R23 R24 ... R65 R66, '
     >     //' ALFY, BETY, ALFZ, BETZ, DY, DYP, DZ, DZP, PHIY, PHIZ,'
     >     //' F(1,IREF),F(2,IREF),F(3,IREF),F(4,IREF),F(5,IREF)'
     >      //' F(6,IREF),F(7,IREF), CMUY, CMUZ, QY, QZ, XCE, YCE, ALE'
c     >     '% R11 R12 R13 R14 R21 R22 R23 ... R43 R44'
c     >     //' R16 R26 R36 R46 R51 R52 R53 R54 R55 R66'
c     >     //' ALFY, BETY, ALFZ, BETZ, DY, DYP, DZ, DZP, PHIY, PHIZ,'
c     >     //' F(1,IREF),F(2,IREF),F(3,IREF),F(4,IREF),F(5,IREF)'
            WRITE(LNWRT,FMT='(A)') 
     >     '% 1   2   3   4   5   6    7   8   9   10  ... 35  36   '
     >     //' 37    38    39    40    41  42   43  44   45    46   '
     >     //' 47        48        49        50        51       '
     >     //' 52        53         54    55    56  57  58   59   60'
        ENDIF

        CALL REFER3(
     >              XCE,YCE,ALE)
        WRITE(LNWRT,FMT='(1P,55(1X,E16.8),2(2X,A),3(1X,E16.8))') 
     >  ((R(IA,IB),IB=1,6),IA=1,6), 
     >  F0(1,1), -F0(1,2), F0(3,3), -F0(3,4), 
     >  F0(1,6), F0(2,6), F0(3,6), F0(4,6), YNU, ZNU,
     >  F(1,IREF),F(2,IREF),F(3,IREF),F(4,IREF),F(5,IREF),
     >  F(6,IREF),F(7,IREF),CMUY,CMUZ,TXTYNU,TXTZNU,xce,yce,ale
c        WRITE(LNWRT,FMT='(1P,41(1X,E16.8))') ((R(IA,IB),IB=1,4),IA=1,4), 
c     >  R(1,6), R(2,6), R(3,6), R(4,6), 
c     >  R(5,1), R(5,2), R(5,3), R(5,4), R(5,5), R(6,6), 
c     >  F0(1,1), -F0(1,2), F0(3,3), -F0(3,4), 
c     >  F0(1,6), F0(2,6), F0(3,6), F0(4,6), YNU, ZNU,
c     >  F(1,IREF),F(2,IREF),F(3,IREF),F(4,IREF),F(5,IREF)
C        CLOSE(LNWRT)
C        KWRMAT = .FALSE.

        CALL FLUSH2(LNWRT,.FALSE.)

      ENDIF


      RETURN

      ENTRY MATIM2(R,T,T3)

C MODIFIED, FM, 04/97
C       ** CHANGE MATRICE TRIANGULAIRE EN CARREE SYMMETRIQUE/DIAG
        DO 11 IA=1,6
          DO 11 IB=1,6
            IC1=IB+1
            DO 11 IC=IC1,6
              T(IA,IC,IB)=T(IA,IB,IC)
 11     CONTINUE

      I=2
      IF(NRES.GT.0) THEN
        WRITE(NRES,103) I
        DO 16 IA=1,6
          IF(IA.GT.1) WRITE(NRES,107)
 107      FORMAT(/)
          DO 16 IB=1,6
C          WRITE(NRES,108) ( IA,IC,IB, T(IA,IC,IB)  , IC=1,IB )
C MODIFIED, FM, 04/97
            WRITE(NRES,108) ( IA,IC,IB, T(IA,IC,IB)  , IC=1,6 )
 108        FORMAT( 6(I4,I2,I1,1P,G11.3) )
 16     CONTINUE
      ENDIF

      CALL SYMPL2(R,T)
 
      IF(NRES.GT.0) THEN
        WRITE(NRES,123) T3(1,1),T3(1,2),T3(1,3),T3(1,4)
 123    FORMAT(//,15X,'COEFFICIENTS  D''ORDRE  SUPERIEUR  ( MKSA ):'
     >  ,//,10X,' Y/Y3   ',5X,1P,G14.5
     >  , /,10X,' Y/T3   ',5X,   G14.5
     >  , /,10X,' Y/Z3   ',5X,   G14.5
     >  , /,10X,' Y/P3   ',5X,   G14.5,/)
        WRITE(NRES,124) T3(2,1),T3(2,2),T3(2,3),T3(2,4)
 124    FORMAT(
     >     10X,' T/Y3   ',5X,1P,G14.5
     >  ,/,10X,' T/T3   ',5X,   G14.5
     >  ,/,10X,' T/Z3   ',5X,   G14.5
     >  ,/,10X,' T/P3   ',5X,   G14.5,/)
        WRITE(NRES,125) T3(3,1),T3(3,2),T3(3,3),T3(3,4)
 125    FORMAT(
     >     10X,' Z/Y3   ',5X,1P,G14.5
     >  ,/,10X,' Z/T3   ',5X,   G14.5
     >  ,/,10X,' Z/Z3   ',5X,   G14.5
     >  ,/,10X,' Z/P3   ',5X,   G14.5,/)
        WRITE(NRES,126) T3(4,1),T3(4,2),T3(4,3),T3(4,4)
 126    FORMAT(
     >     10X,' P/Y3   ',5X,1P,G14.5
     >  ,/,10X,' P/T3   ',5X,   G14.5
     >  ,/,10X,' P/Z3   ',5X,   G14.5
     >  ,/,10X,' P/P3   ',5X,   G14.5)
 
        WRITE(NRES,101) IEX(1),(F(J,1),J=1,7)
  101   FORMAT(' TRAJ 1 IEX,D,Y,T,Z,P,S,time :',I3,1P,5G12.4,2G17.5)
      ENDIF

      RETURN

      ENTRY MATIM6(KWRI)
      KWRMAT = KWRI
      RETURN

 96   KWRMAT = .FALSE.
      CALL ENDJOB('ERROR upon  open  old  file '//FNAME,-99)
      RETURN

      END
