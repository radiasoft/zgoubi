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
      SUBROUTINE FITGTV(FITING,NOMFIC,
     >                                FITGET)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL FITING
      CHARACTER(*) NOMFIC
      LOGICAL FITGET

      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     >     IREP(MXT),AMQLU,PABSLU
      PARAMETER (LBLSIZ=10)
      CHARACTER(LBLSIZ) LABEL
      COMMON /LABEL/ LABEL(MXL,2)
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      INCLUDE 'MXFS.H'
      PARAMETER (KSIZ=10)
      CHARACTER FAM*(KSIZ),LBF*(LBLSIZ)
      COMMON/SCALT/ FAM(MXF),LBF(MXF,MLF)
 
      CHARACTER(60) NAMFIC
      INTEGER DEBSTR,FINSTR

      LOGICAL OPN, IDLUNI, EMPTY, EXS

      CHARACTER(200) TXT200
      CHARACTER(1) TXT1

      DIMENSION AFIT(MXL), IPRM(MXL)
      SAVE AFIT
      CHARACTER KLE*(KSIZ),LBL1*(LBLSIZ),LBL2*(LBLSIZ)
      CHARACTER(KSIZ) KLEFIT(MXL)
      CHARACTER LB1FIT(MXL)*(LBLSIZ),LB2FIT(MXL)*(LBLSIZ)
      SAVE KLEFIT, LB1FIT, LB2FIT, KREAD

      CHARACTER(KSIZ) KLEY
      LOGICAL GTTEXT, OK, STRCON
      LOGICAL FITFNL
      LOGICAL FIRST
      LOGICAL SKPLIN

      DATA SKPLIN / .FALSE. /

      IF(FITING) THEN
        IF(NRES.GT.0) WRITE(NRES,FMT='(
     >  /,10X,''Now fitting, hence GETFITVAL is inhibited'')')
        RETURN
      ENDIF

      FITGET = .FALSE.

      NAMFIC=NOMFIC(DEBSTR(NOMFIC):FINSTR(NOMFIC))

      IF(NRES .GT. 0) WRITE(NRES,100) 
     >  NOMFIC(DEBSTR(NOMFIC):FINSTR(NOMFIC))
 100  FORMAT(/,10X, 
     >'Parameter values to be refreshed from FIT storage file : ',A,/)

      IF(NAMFIC.EQ.'NONE' .OR. NAMFIC.EQ.'none') GOTO 97

      INQUIRE(FILE=NAMFIC,ERR=2,OPENED=OPN,NUMBER=LUN,EXIST=EXS)
      IF(.NOT. EXS) GOTO 2
      IF(OPN) CLOSE(UNIT=LUN)
      GOTO 3
 2    IF(NRES .GT. 0) WRITE(NRES,FMT='(3A,/)')     
     >'    Could not find file ',NOMFIC(DEBSTR(NOMFIC):FINSTR(NOMFIC)),
     >' ;  proceeding without it ... '
      RETURN

 3    CONTINUE        
      IF(IDLUNI(
     >    LUN)) OPEN(UNIT=LUN,FILE=NAMFIC,status='OLD',ERR=99)

      IF(NRES.GT.0) WRITE(NRES,101) NAMFIC
 101  FORMAT(/,10X,' Ok, opened storage file ',A,/)

      OK = GTTEXT(NRES,LUN,'STATUS OF VARIABLES'
     >                                         ,TXT200)
      IF(.NOT. OK) GOTO 98

      KREAD = 0
 1    CONTINUE
        READ(LUN,FMT='(A)',ERR=1,END=10) TXT200
        IF(EMPTY(TXT200)) GOTO 1
        IF(STRCON(TXT200,'STATUS OF CONSTRAINTS',
     >                                           IS)) GOTO 10
        TXT1 = TXT200(DEBSTR(TXT200):DEBSTR(TXT200))
        IF(TXT1 .EQ. '%' .OR. TXT1 .EQ. '#' .OR. TXT1 .EQ. '!') THEN
          IF(NRES.GT.0) WRITE(NRES,FMT='(A)') TXT200(DEBSTR(TXT200):132)
          GOTO 1
        ENDIF
        READ(TXT200,*,end=10,err=1) 
     >      NUML,I,ISI,XK,XII,XI,XJ,PI,KLE,LBL1,LBL2

        KREAD = KREAD + 1
        IPRM(KREAD) = ISI
        AFIT(KREAD) = XI
        KLEFIT(KREAD) = KLE
        LB1FIT(KREAD) = LBL1
        LB2FIT(KREAD) = LBL2
C        NEWVAL(NUML,ISI) = 1

C        K=KREAD+NV
C        J=K+NV
        IF(NRES.GT.0) 
     >  WRITE(NRES,400) NUML,I,ISI,XK,XII,XI,XJ,PI,KLE,LBL1,LBL2
400     FORMAT(1P, 
     >  2X,I3,3X,I2,4X,I3,2(2X,G10.3),2X,G17.10,2(1X,G10.3),3(1X,A))

        GOTO 1

 10   CONTINUE
      IF(KREAD.GE.1) FITGET = .TRUE.
      IF(NRES .GT. 0) THEN
        WRITE(NRES,103) KREAD,NAMFIC
 103    FORMAT(/,10X, I3, 
     >  ' variables have been read from FIT storage file ',A,' :',/)
        WRITE(NRES,FMT='(A)')
     >  ' Var#  Param       Val     Kley      Lbl1   Lbl2 '
        DO I = 1, KREAD
          WRITE(NRES,409)I,IPRM(I),AFIT(I),KLEFIT(I),LB1FIT(I),LB2FIT(I)
 409      FORMAT(1P,2(2X,I3),3X,E13.4,2X,3(1X,A))
        enddo
      ENDIF
C----- OBJECT coordinates are redefined according to FIT variables, if any
C To be completed

C FM 14-08-01
      CLOSE(LUN)
      RETURN

 97   CONTINUE
      IF(NRES.GT.0) THEN
        WRITE(NRES,FMT='(20X,''GETFITVAL :  file name is '''''',A
     >  ,'''''''')')
     >  NOMFIC(DEBSTR(NOMFIC):FINSTR(NOMFIC))
        WRITE(NRES,FMT='(20X,''Thus no updating requested, 
     >  carrying on without that...'',/)') 
      ENDIF
      RETURN 

 98   CONTINUE
      IF(NRES.GT.0) THEN
        WRITE(NRES,FMT='(20X,''GETFITVAL :  sorry, could not'',
     >  '' find string ''''STATUS OF VARIABLES'''' in file '',A)') 
     >  NOMFIC(DEBSTR(NOMFIC):FINSTR(NOMFIC))
        WRITE(NRES,FMT='(20X,''Thus, no initialization of'',
     >  '' variables performed, carrying on without that...'',/)') 
      ENDIF
      RETURN 

 99   CONTINUE
      IF(NRES.GT.0) WRITE(NRES,FMT='(3A)') ' Could not open '
     >  ,NOMFIC(DEBSTR(NOMFIC):FINSTR(NOMFIC))
     >  ,', give ''none'' as name so to allow proceeding.  '
      CALL ENDJOB('*** Error, SBR FITGTV -> can`t open strage file',-99)
      RETURN 

      ENTRY FITGT1
      CALL ZGKLEY(
     >            KLEY)
      CALL FITST5(
     >            FITFNL)     

        FIRST = .TRUE. 
        DO IV = 1, KREAD
C          IF(NEWVAL(NOEL,IV) .EQ. 1) THEN
          IF(KLEFIT(IV) .EQ. KLEY) THEN
            IF(.NOT. EMPTY(LABEL(NOEL,1))) THEN
              IF(
     >        (LB1FIT(IV) .EQ. '*' 
     >        .OR. LB1FIT(IV) .EQ. LABEL(NOEL,1)) 
     >        .AND.
     >        (EMPTY(LABEL(NOEL,2)) .OR. LB2FIT(IV) .EQ. '*' 
     >        .OR. LB2FIT(IV) .EQ. LABEL(NOEL,2)) ) THEN 
                SKPLIN = .TRUE.
                IF(.NOT. FITFNL) THEN
                  TEMP = A(NOEL,IPRM(IV)) 
                  A(NOEL,IPRM(IV)) = AFIT(IV)
                  IF(NRES .GT. 0) WRITE(NRES,
     >            FMT='('' GETFITVAL procedure. ''
     >            ,'' Former  value  A(NOEL=''
     >            ,I3,'','',I3,'') = '',1P,G15.6,'',   changed to '', 
     >            G15.6)') NOEL,IPRM(IV),TEMP,A(NOEL,IPRM(IV))
                ELSE
                  IF(FIRST) THEN
                    IF(NRES .GT. 0) WRITE(NRES,FMT='('' Pgm fitgtv. ''
     >              ,''GETFITVAL procedure will not be applied, since ''
     >              ,'' this is a last run with ''
     >              ,/,'' variable values''
     >              ,'' following from just completed FIT.'')')
                  ENDIF
                  FIRST = .FALSE.
                  IF(NRES .GT. 0) WRITE(NRES,FMT='(
     >            '' Variable concerned, and its present value :  ''
     >            ,'' A(NOEL='',I3,'','',I3,'') = '',1P,G15.6)') 
     >            NOEL,IPRM(IV),A(NOEL,IPRM(IV))
                ENDIF          
              ENDIF
            ELSE
              IF(NRES .GT. 0) WRITE(NRES,
     >        FMT='('' GETFITVAL procedure, cannot be applied : ''
     >        ,'' Elements need be labeled. Provide label #1.'')')
            ENDIF
          ENDIF
        ENDDO
        IF(SKPLIN) THEN
          IF(NRES .GT. 0) WRITE(NRES,FMT='(A)') ' '
          SKPLIN = .FALSE.
        ENDIF

      RETURN
      END
