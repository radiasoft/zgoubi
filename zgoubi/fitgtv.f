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
      SUBROUTINE FITGTV(NOMFIC,
     >                         FITGET)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER NOMFIC*(*)
      LOGICAL FITGET

      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
 
      CHARACTER*60 NAMFIC
      INTEGER DEBSTR,FINSTR

      LOGICAL OPN, IDLUNI

C      PARAMETER (NNEWV=MXL*MXD)
      CHARACTER TXT132*132, TXT1*1

      INCLUDE 'MXFS.H'
      PARAMETER (LBLSIZ=8)
      PARAMETER (KSIZ=10)
      CHARACTER FAM*(KSIZ),LBF*(LBLSIZ),KLEY*(KSIZ),LABEL*(LBLSIZ)
      COMMON/SCALT/ FAM(MXF),LBF(MXF,2),KLEY,LABEL(MXL,2)

      DIMENSION AFIT(MXL), IPRM(MXL)
C      DIMENSION AFIT(MXL,MXD), NEWVAL(MXL,MXD)
C      SAVE AFIT, NMI, NMA, NEWVAL
      SAVE AFIT
C      SAVE AFIT, NEWVAL
      CHARACTER KLE*(KSIZ),LBL1*(LBLSIZ),LBL2*(LBLSIZ)
      CHARACTER KLEFIT(MXL)*(KSIZ)
      CHARACTER LB1FIT(MXL)*(LBLSIZ),LB2FIT(MXL)*(LBLSIZ)
      SAVE KLEFIT, LB1FIT, LB2FIT, KREAD

C      DATA NEWVAL / NNEWV*0/

      FITGET = .FALSE.

      NAMFIC=NOMFIC(DEBSTR(NOMFIC):FINSTR(NOMFIC))

      IF(NRES .GT. 0) WRITE(NRES,100) NAMFIC
 100  FORMAT(/,10X, 
     >'Parameter values to be refreshed from FIT storage file : ',A,/)

      IF(NAMFIC.EQ.'NONE' .OR. NAMFIC.EQ.'none') GOTO 97

      INQUIRE(FILE=NAMFIC,ERR=98,OPENED=OPN,NUMBER=LUN)
      IF(OPN) CLOSE(UNIT=LUN)
        
      IF(IDLUNI(
     >          LUN)) OPEN(UNIT=LUN,FILE=NAMFIC,ERR=99)

      WRITE(NRES,101) NAMFIC
 101  FORMAT(/,10X,' Ok, opened storgae file ',A,/)

C      CALL HEADER(LUN,NRES,4,.FALSE.,
      CALL HEADER(LUN,NRES,2,.FALSE.,
     >                              *999)
      WRITE(NRES,*) ' ' 
      
C      NMI = 999999
C      NMA = -999999
      KREAD = 0
 1    CONTINUE
        READ(LUN,FMT='(A)',ERR=10,END=10) TXT132
        TXT1 = TXT132(DEBSTR(TXT132):DEBSTR(TXT132))
        IF(TXT1 .EQ. '%' .OR. TXT1 .EQ. '#' .OR. TXT1 .EQ. '!') GOTO 1
        READ(TXT132,*,end=10) NUML,I,ISI,XK,XII,XI,XJ,PI,KLE,LBL1,LBL2
C        write(*,*) 'fitgtv ',NUML,I,ISI,XK,XII,XI,XJ,PI,KLE,LBL1,LBL2
        IPRM(I) = ISI
        AFIT(I) = XI
        KLEFIT(I) = KLE
        LB1FIT(I) = LBL1
        LB2FIT(I) = LBL2
C        NEWVAL(NUML,ISI) = 1
        KREAD = KREAD+1

        K=I+NV
        J=K+NV
        IF(NRES.GT.0) 
     >  WRITE(NRES,400) NUML,I,ISI,XK,XII,XI,XJ,PI,KLE,LBL1,LBL2
400     FORMAT(1P, 
     >  2X,I3,3X,I2,4X,I3,2(2X,G10.3),2X,G17.10,2(1X,G10.3),3(1X,A))

        GOTO 1

 10   CONTINUE
      IF(KREAD.GE.1) FITGET = .TRUE.
      IF(NRES .GT. 0) WRITE(NRES,103) KREAD,NAMFIC
 103  FORMAT(/,10X, I3, 
     >  ' variables have been read from FIT storage file ',A,/)

C----- OBJECT coordinates are redefined according to FIT variables, if any
C To be completed

      RETURN

 97   CONTINUE
      RETURN 
 98   CONTINUE
      CALL ENDJOB('*** Error, SBR FITGTV -> at INQUIRE ',-99)
      RETURN 
 99   CONTINUE
      CALL ENDJOB('*** Error, SBR FITGTV -> can`t open strage file',-99)
      RETURN 
 999  CONTINUE
      CALL ENDJOB('*** Error, SBR FITGTV ->  Read error at header',-99)
      RETURN 

      ENTRY FITGT1
        DO 2 IV = 1, KREAD
C          IF(NEWVAL(NOEL,IV) .EQ. 1) THEN
          IF(KLEFIT(IV) .EQ. KLEY) THEN
            IF(
     >        (LB1FIT(IV) .EQ. '*' 
     >        .OR. LB1FIT(IV) .EQ. LABEL(NOEL,1)) 
     >        .AND.
     >        (LB2FIT(IV) .EQ. '*' 
     >        .OR. LB2FIT(IV) .EQ. LABEL(NOEL,2)) ) THEN 
                 TEMP = A(NOEL,IPRM(IV)) 
                 A(NOEL,IPRM(IV)) = AFIT(IV)
                 IF(NRES .GT. 0) WRITE(NRES,
     >           FMT='('' GETFITVAL procedure.  Former  A('',I3,
     >           '','',I3,'') = '',1P,G12.4,'' changed to new value '', 
     >           G12.4)') NOEL,IPRM(IV),TEMP,A(NOEL,IPRM(IV))
            ENDIF
          ENDIF
 2      CONTINUE
      RETURN
      END
