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
      SUBROUTINE RFIT(KLEY,
     >                     PNLTY,ITRMA,ICPTMA,FITFNL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ***************************************
C     READS DATA FOR FIT PROCEDURE WITH 'FIT'
C     ***************************************
      CHARACTER(*) KLEY
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      PARAMETER (MXV=40) 
      COMMON/MIMA/ DX(MXV),XMI(MXV),XMA(MXV)
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER(80) TA
      PARAMETER (MXTA=45)
      COMMON/DONT/ TA(MXL,MXTA)
      COMMON/VARY/NV,IR(MXV),NC,I1(MXV),I2(MXV),V(MXV),IS(MXV),W(MXV),
     >IC(MXV),IC2(MXV),I3(MXV),XCOU(MXV),CPAR(MXV,7)

      CHARACTER(132) TXT132
      LOGICAL STRCON, CMMNT, FITFNL
      CHARACTER(40) STRA(10)

      PARAMETER (ICPTM1=1000, ICPTM2=1000)
      LOGICAL FITSAV
      INTEGER DEBSTR, FINSTR
      LOGICAL OK
      CHARACTER(80) FNAME
      LOGICAL EMPTY
      LOGICAL FIRST 

      DATA FIRST / .TRUE. /
c      DATA PNLTY, ITRMA, ICPTMA / 1.D-10, 90, 1000 /
c      DATA FITFNL / .TRUE. /

C  READ NV [,'nofinal','save' [FileName]]
      READ(NDAT,FMT='(A)') TXT132
      IF(STRCON(TXT132,'!',
     >                     IIS)) TXT132 = TXT132(1:IIS-1)
      READ(TXT132,*) NV
      IF(NV.LT.1) RETURN
      FITFNL = .NOT. STRCON(TXT132,'nofinal',
     >                                       IIS) 
      FITSAV = STRCON(TXT132,'save',
     >                             JJS) 
      IF(FITSAV) THEN
        TXT132 = TXT132(JJS+4:FINSTR(TXT132))
        TXT132 = TXT132(DEBSTR(TXT132):FINSTR(TXT132))

        IF(.NOT. EMPTY(TXT132)) THEN
          IF(TXT132(DEBSTR(TXT132):DEBSTR(TXT132)+6).NE.'nofinal') THEN
            READ(TXT132(DEBSTR(TXT132):FINSTR(TXT132)),*) FNAME
          ELSE
            FNAME = 'zgoubi.FITVALS.out'
          ENDIF
        ELSE
          FNAME = 'zgoubi.FITVALS.out'
        ENDIF
        IF(FIRST) CALL FITNU6(FNAME)
        FIRST = .FALSE.  
      ENDIF

      DO I=1,NV
        READ(NDAT,FMT='(A)') TXT132
        CMMNT = STRCON(TXT132,'!',
     >                            III) 
        IF(CMMNT) THEN
          III = III - 1
        ELSE
          III = 132
        ENDIF
        IF(STRCON(TXT132(1:III),'[',
     >                              II)) THEN
C--------- New method
          READ(TXT132(1:II-1),*) IR(I),IS(I),XCOU(I)
          IF(STRCON(TXT132,']',
     >                      II2)) THEN
            READ(TXT132(II+1:II2-1),*) XMI(I),XMA(I)
          ELSE
            CALL ENDJOB(' SBR RFIT, wrong input data / variables',-99)
          ENDIF
        ELSE
C--------- Old method
 
           READ(TXT132,*) IR(I),IS(I),XCOU(I),DX(I)
           XI = A(IR(I),IS(I))
           XMI(I)=XI-ABS(XI)*DX(I)
           XMA(I)=XI+ABS(XI)*DX(I)

        ENDIF
      ENDDO

C  READ NC [,PNLTY [,ITRMA [,ICPTMA]]]
      READ(NDAT,FMT='(A)') TXT132
      IF(STRCON(TXT132,'!',
     >                     III)) TXT132 = TXT132(1:III-1) 
      CALL STRGET(TXT132,3,
     >                     NSTR,STRA)
      IF(NSTR .GT. 3) NSTR = 3
      IF(NSTR.GE.1) THEN
        READ(STRA(1),*,ERR=97,END=97) NC
        IF(NSTR.GE.2) THEN
          READ(STRA(2),*,ERR=43,END=43) PNLTY
          IF(NSTR.EQ.3) THEN
            READ(STRA(3),*,ERR=44,END=44) ITRMA
            IF(NSTR.EQ.4) THEN
              READ(STRA(4),*,ERR=44,END=44) ICPTMA
            ELSE
              GOTO 44
            ENDIF
          ELSE
            GOTO 42
          ENDIF
        ELSE
          GOTO 43
        ENDIF
      ENDIF
      GOTO 45
 43   CONTINUE
      PNLTY = 1.D-10
 42   CONTINUE
      ITRMA = 90
 44   CONTINUE
      IF(KLEY .EQ. 'FIT') THEN
        ICPTMA = ICPTM1
      ELSEIF(KLEY .EQ. 'FIT2') THEN
        ICPTMA = ICPTM2
      ENDIF
 45   CONTINUE

      IF(NC.LT.1) RETURN
      DO 4 I=1,NC
        READ(NDAT,*,ERR=41,END=41) XC,I1(I),I2(I),I3(I),V(I),W(I),
     >  CPAR(I,1),(CPAR(I,JJ),JJ=2,NINT(CPAR(I,1))+1)
 41     CONTINUE
        IC(I) = INT(XC)
        IC2(I) = NINT(10.D0*XC - 10*IC(I))
 4    CONTINUE  

C----- Looks for possible parameters, with values in CPAR, and action
      DO 5 I=1,NC
        IF    (IC(I).EQ.3) THEN 
C--------- Traj coord
          IF(I1(I).EQ.-3)  THEN
            CALL DIST2W(CPAR(I,2), CPAR(I,3), CPAR(I,4))
          ENDIF
        ELSEIF(IC(I).EQ.5) THEN 
C--------- Numb. particls
          IF(I1(I).GE.1) THEN
            IF(I1(I).LE.3) THEN
              CALL ACCENW(CPAR(I,2))
            ELSEIF(I1(I).LE.6) THEN
              CALL ACCEPW(CPAR(I,2))
            ENDIF
          ENDIF
        ENDIF
 5    CONTINUE  

      RETURN

 97   CALL ENDJOB('SBR rfit, error input data at NC',-99)
      RETURN
      END
