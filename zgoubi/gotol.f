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
      SUBROUTINE GOTOL(IPASS,MXKLE,KLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER(*) KLE(*)
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
C      PARAMETER (LNTA=132) ; CHARACTER(LNTA) TA
C      PARAMETER (MXTA=45)
      INCLUDE "C.DONT.H"     ! COMMON/DONT/ TA(MXL,MXTA)

      CHARACTER(132) TXT132
      PARAMETER (MX2=2)
      CHARACTER(20) STRA(MX2)
      INTEGER DEBSTR, FINSTR
      LOGICAL STRCON
      PARAMETER (MST=35)
      PARAMETER (LBLSIZ=10)
      CHARACTER(LBLSIZ) LBLST(MST)
      DIMENSION NLTO(MST), NLBCK(MST)
      CHARACTER(LEN=30) FRMT
      LOGICAL FIRST, OK, IDLUNI

      PARAMETER (KSIZ=10)
      CHARACTER(KSIZ) KEY, LBL1, LBL2

      SAVE NLTO, MLST, NLBCK, LUN, IQCNT

      DATA FIRST / .TRUE. /
      DATA IQCNT / 1 /

      IF(IPASS.GT.MST) CALL ENDJOB('Pgm gotol. Return address is too '//
     >'large. Increase MST - now MST = ',MST)

      IF(TA(NOEL,1) .NE. 'GOBACK') THEN
C GOTO 
        IF(TA(NOEL,1) .EQ. 'PASS#') THEN

C          IF(FIRST) THEN
          CALL STRGET(TA(NOEL,2),MST,
     >                              MLST,LBLST)
          IF(MLST.GT.MST) WRITE(NRES,FMT='(10X,A,I0)')
     >    'Pgm gotol. Warning : list of addresses shortened to '//
     >    'maximum allowed : ',MST

          IF(NRES.GT.0) THEN
            WRITE(NRES,FMT='(/,10X,
     >      ''Will switch to next element according to pass #. '',
     >      ''List of '',I0,'' switches is as follows : '',/)') MLST
c                 write(*,*) ' gotol frmt ta1 : ',
c     >             ta(noel,1)(DEBSTR(TA(NOEL,1)):FINSTR(TA(NOEL,1)))
c                 write(*,*) ' gotol frmt ta2 : ',
c     >             ta(noel,2)(DEBSTR(TA(NOEL,2)):FINSTR(TA(NOEL,2)))
c                 write(*,*) ' gotol frmt ',frmt
c                  read(*,*)
            WRITE(FRMT,FMT='(A,I0,A)') '(10X,A,', MLST, 'I10)' 
            WRITE(NRES,FMT=FRMT)'Pass #        : ',(I,I=1,MLST)
            WRITE(FRMT,FMT='(A,I0,A)') '(10X,A,4X,', MLST, 'A)' 
            WRITE(NRES,FMT=FRMT)'Element label : ',(LBLST(I),I=1,MLST)
            WRITE(NRES,FMT='(/,15X,''Present pass is  # '',I0)') IPASS
          ENDIF

          CALL SYSTEM('cp zgoubi.dat zg_temp_gotol.dat')
          OK = IDLUNI(
     >                  LUN)
          OPEN(UNIT=LUN,FILE='zg_temp_gotol.dat')

c               write(*,*) ' GOTOL : ok open '
C          ENDIF

C          REWIND(LUN)
          CALL ZGNBLM( 
     >                NBLMN)

          NUEL = 0
 11       CONTINUE
          TXT132 = ' ' 
          DOWHILE(TXT132(1:1).NE.'''' .AND. NUEL.LT.NBLMN)
c            IF(NUEL.EQ.NBLMN) CALL ENDJOB('Pgm gotol. Could not goto.'
c     >      //' Check GOTO list (too short ?). Now number of elements'
c     >      //' in GOTO list is ',-99)
C Read keyword [/ label1 [/ label2]]
            READ(LUN,FMT='(A)',ERR=10,END=10) TXT132
            TXT132 = TXT132(DEBSTR(TXT132):FINSTR(TXT132))
          ENDDO

          IF(NUEL.GE.NBLMN) then
            WRITE(ABS(NRES),*) 'Branching tag at this GOTO is ''',
     >      LBLST(IPASS)(DEBSTR(LBLST(IPASS)):FINSTR(LBLST(IPASS))),''''
            CALL ENDJOB('Pgm gotol, keyword GOTO : '//
     >      ' could not fing branching tag. Reached lmnt # ',NUEL-1)
          ENDIF

          NUEL = NUEL + 1
          CALL STRGET(TXT132,MX2,
     >                           IDUM,STRA)
          IF(STRA(2) .EQ. LBLST(IPASS)) THEN
            NLBCK(IPASS) = NOEL+1
            NOEL = NUEL
            CALL GO2KEY(NOEL,IQCNT,MXKLE,KLE,
     >                                       KEY, LBL1, LBL2) 
          ELSE
            GOTO 11
          ENDIF

          IF(NRES.GT.0) THEN
            WRITE(NRES,FMT='(/,15X,''This GOTO switches to'',
     >      '' element # '',I0,'' with label #1 : '',A)') NOEL,STRA(2)
            WRITE(NRES,FMT='(15X,''Return address will be'',
     >      '' element # '',I0)') NLBCK(IPASS)
          ENDIF 
          CLOSE(LUN,status='delete')
          NOEL = NOEL - 1

        ELSE

          WRITE(NRES,FMT='(10X,A,A)') 'Pgm goto. No such option ',
     >    TA(NOEL,1)(DEBSTR(TA(NOEL,1)):FINSTR(TA(NOEL,1)))
          CALL ENDJOB('Pgm goto. Check input data list ',-99)

        ENDIF

      ELSE
C GOBACK
        NOEL = NLBCK(IPASS) 
        CALL GO2KEY(NOEL,iqcnt,mxkle,kle,
     >                                   KEY, LBL1, LBL2) 

c            WRITE(*,*) ' GOBACK NLBCK(IPASS)  ',NLBCK(IPASS) 
c                  read(*,*)

          IF(NRES.GT.0) THEN
            WRITE(NRES,FMT='(/,15X,''Present pass is  # '',I0)') IPASS
            WRITE(NRES,FMT='(/,15X,''This GOTO/GOBACK switches to '',
     >      ''return address, namely :   element # '',I0)') NLBCK(IPASS)
          ENDIF 

          NOEL = NOEL - 1
      ENDIF

      RETURN

 10   CONTINUE
      CALL ENDJOB('Pgm gotol. Error : could not find goto label.'
     >//' Current NOEL = ',NOEL)
      RETURN
      END
