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
CDECK GRAPH1 
      SUBROUTINE GRAPH7(NLOG,LM,NOMFIC,wrkDir)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER(*) NOMFIC
      CHARACTER(*) WRKDIR
      
      LOGICAL OKECH, OKVAR, OKBIN
      COMMON/ECHL/ OKECH, OKVAR, OKBIN
      INCLUDE 'MXVAR.H'
      CHARACTER(7) KVAR(MXVAR), KDIM(MXVAR)
      CHARACTER(9) KPOL(2)
      COMMON/INPVR/ KVAR, KPOL, KDIM
      COMMON/LUN/NDAT,NRES,NPLT,NFAI,NMAP,NSPN
      COMMON/REBELO/ NRBLT,IPASS,KWRI,NNDES,STDVM
      INCLUDE 'MAXNTR.H'
      COMMON/TRACKM/COOR(NTRMAX,9),NPTS,NPTR
      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB

      LOGICAL OKOPN, INPECH, CHANGE, KLIPS, KHIST, OK24, OK23
      LOGICAL TYLAB,TYLABR
      CHARACTER(16) TXTL
      LOGICAL FIRST

      CHARACTER(1) KLET

      CHARACTER(1) RLIPS,RLIPS0,RLIPS2,RSNPT,RSNPT0,RSNPT2
      CHARACTER(1) RHIST,RHIST0,RHIST2

      SAVE FIRST
      SAVE NL
      SAVE NPASS

      LOGICAL OKREW
     
      LOGICAL IDLUNI
      CHARACTER(11) FNAM
      LOGICAL FIRST2
      SAVE FNAM, NLIPS, FIRST2

      PARAMETER (ITWO=2)
     
      INTEGER DEBSTR, FINSTR

      DATA  OKOPN / .FALSE. /
      DATA CHANGE / .TRUE./
      DATA FIRST / .TRUE. /
      DATA RSNPT,RSNPT2 / 'N','Y' /
      DATA RLIPS,RLIPS2 / 'Y','N' /
      DATA RHIST,RHIST2 / 'Y','N' /
      DATA OKREW / .TRUE. /

      DATA FNAM /'zpop.lips'/
      DATA FIRST2 /.TRUE./

C       nlips = 77
C         goto 921
C         write(*,*) 'Pgm graph7 ' 
C          read(*,*)

      IF(FIRST2) THEN
        IF (IDLUNI(NLIPS)) THEN
          OPEN(UNIT=NLIPS,FILE=FNAM,ERR=997)
          CLOSE(UNIT=NLIPS,STATUS='DELETE')
          OPEN(UNIT=NLIPS,FILE=FNAM,STATUS='NEW',ERR=997)
          WRITE(6,*)
     >      ' Pgm graph7. Opened ',FNAM(debstr(FNAM):finstr(FNAM))
        ELSE
          WRITE(6,*) ' *** Problem in graph7 : No idle unit number ! '
          WRITE(6,*) '     No idle unit number !  Forced to 77 '
          WRITE(6,*) '                            file will be fort.77 '
          NLIPS=77
        ENDIF
        FIRST2 = .FALSE.
      ENDIF

      GOTO 921

 920  CONTINUE      
      CALL FBGTXT
      I=IDLG('('' Press RETURN for more :'')','    ',1)

 921  CONTINUE
      CALL FBGTXT
      CALL HOMCLR

      WRITE(6,100) '[...]'//wrkDir(finstr(wrkDir)-60:finstr(wrkDir))
     >,NOMFIC,KVAR(KX),KX,KVAR(KY),KY
 100  FORMAT(//,3X,60('*'),//,20X,' MENU - Graphic processing :' ,//
     > ,9X,'        ' ,/
     > ,4x,' Working directory now : ',/,5x,A,  //
     > ,9X,'        ' ,/
     1 ,9X,'  1    OPEN  FILE - current is ',A ,/
     2 ,9X,'  2    VARIABLES  TO  PLOT ' ,/
     2 ,9X,'            PRESENT:  H-AXIS :',A,'(KX=',I2,')',/
     2 ,9X,'                      V-AXIS :',A,'(KY=',I2,')',/  
     3 ,9X,'  3    PLOT  OPTIONS',/
     4 ,9X,'  4    MANUAL  SCALES ',/
     5 ,9X,'  5    AUTOMATIC  SCALES ',/
     7 ,9X,'  7    **  PLOT  Y V.S. X  **         ',/
     8 ,9X,'  8    Save  graphic  screen ',/
     9 ,9X,'  9    EXIT  THIS  MENU     ',/
     > ,9X,' 12    ERASE  DISPLAY    ',/
     > ,9X,' 15    KILL/CREATE  GRAPHIC  SUB-PROCESS      ',/
     > ,9X,' 20    Superpose a curve',/
     > ,9X,' 77    **  PLOT  Y V.S. X, turn by turn  **    ',/
     > ,9X,' 88    HELP',/
     > ,3X,60('*'),/)

      IF(FIRST) THEN
        CALL PLINIT('////MENU07////',
     >                               NL,OKOPN,CHANGE,NOMFIC)
        FIRST = .FALSE.
      ENDIF
      IF(.NOT. OKOPN) CALL OPNWRN(1)

      WRITE(6,101) ' * Option  number : '
 101  FORMAT(A20)
      READ(5,201,ERR=921) IOPT
 201  FORMAT(I2)
      IF(IOPT.EQ.88) GOTO 88
      IF(IOPT.EQ.77) GOTO 77
      GOTO ( 1, 2, 3, 4, 5,921, 7, 8,99,921,921,12,
     > 921,921,15,921,921,921,921,20 ) IOPT  
      GOTO 921

 1    CONTINUE  
        IOP = 0
        CALL OPNMN(IOP,
     >                 NL,OKOPN,CHANGE,NOMFIC)
        IF(CHANGE) THEN
          IF(KY .EQ. 28) OKBIN=.FALSE.
        ENDIF
      GOTO 921

 2    CONTINUE
        KX0=KX
        KY0=KY
        CALL INPVAR
        OKVAR = .TRUE.
      GOTO 921 

 3    CONTINUE                      
        CALL PLTOPT(
     >              LM,OKECH,OKBIN,KOPT)
      GOTO 921

 4    CONTINUE
        OKECH=INPECH()
      GOTO 921
        
 5    CONTINUE
        IF(OKOPN) THEN
          OKECH = .FALSE.
          WRITE(6,*)
          WRITE(6,*) ' Busy: calculating scales'
          CALL CALECH(NL,
     >                   NOC,NPASS)
          CALL READC3(KL1,KL2)
          IF(KL1.EQ.-1) THEN
            WRITE(TXTL,FMT='(A5)') '* all'
          ELSE
            WRITE(TXTL,FMT='(I5,A,I6)') KL1,' to ',KL2
          ENDIF
C          IF(LM.EQ.-1) THEN
C             WRITE(TXTL,FMT='(A5)') '* all'
C          ELSE
C             WRITE(TXTL,FMT='(I5)') LM
C          ENDIF
          CALL READC5(KT1,KT2)
          CALL READC1(
     >                KP1,KP2,KP3)
          CALL READC9(KKEX,KLET)
          WRITE(6,*)
          WRITE(6,*) ' Particle  # :',KT1,' to ',KT2,
     >                              ' ;  element  # : ',TXTL,' ;'
          WRITE(6,*) ' pass # ',KP1,' to ',KP2,', ipass-modulo ',KP3
     >    ,' ;  case  KEX/KTag =', KKEX,'/',KLET
          WRITE(6,*) ' calculation  of  scales  from ',NOC,'  points;'
          WRITE(6,*) ' Xmin-max :',XMI,XMA,KDIM(KX)
          WRITE(6,*) ' Ymin-max :',YMI,YMA,KDIM(KY)
          IF(NOC .EQ. 0) THEN
            WRITE(6,*) 
            WRITE(6,*) '   ***  No  data  could  be  read  ***'
            WRITE(6,*) '     Check particle or element number ' 
            WRITE(6,*) '        -  Main menu, option 2  -'
          ENDIF
          IF(XMI.LT. XMA .AND. YMI.LT. YMA) THEN
C            CALL TXTFBG
            CALL LINTYP(-1)
          DDX=ABS(XMA-XMI)*0.05D0
          DDY=ABS(YMA-YMI)*0.05D0
          CALL TRAXES(XMI,XMA,YMI,YMA,ITWO)
C          CALL TRAXES(XMI-DDX,XMA+DDX,YMI-DDY,YMA+DDY,ITWO) 
            OKECH = .TRUE.
          ENDIF
        ELSE
          CALL OPNWRN(1)
        ENDIF
      GOTO 920
     
 7    CONTINUE 
        IF(NOMFIC.EQ.'none') THEN
          IOP=1
        ELSE
          IOP=5
        ENDIF
        IF(TYLAB(KX,KY)) THEN
          RSNPT0=RSNPT
 75       CONTINUE      
          WRITE(6,FMT='(/,
     >    ''  Plot synoptic ? ('',A1,''/'',A1,'') : ''
     >    )') RSNPT,RSNPT2
          READ(5,FMT='(A1)',ERR=75) RSNPT
          IF(RSNPT.EQ.'n') RSNPT='N'
          IF(RSNPT.EQ.'y') RSNPT='Y'
          IF(RSNPT.NE.RSNPT2) THEN
            RSNPT=RSNPT0
          ELSE
            IF(RSNPT2.EQ.'Y') THEN
              RSNPT2='N'
            ELSE
              RSNPT2='Y'
            ENDIF
          ENDIF
          IF(RSNPT.EQ.'Y') THEN
            CALL OPNMN(IOP,
     >                     NL,OKOPN,CHANGE,NOMFIC)
                 IF(TYLABR()) CALL TRACSY(2,*79)
 79            CONTINUE
          ENDIF
        ENDIF

        IF(OKOPN) THEN

          IF( .NOT. OKECH) THEN
            WRITE(6,*) 
            WRITE(6,*) ' Busy: calculating scales'
            CALL CALECH(NL,
     >                     NOC,NPASS)
            CALL READC3(KL1,KL2)
            IF(KL1.EQ.-1) THEN
              WRITE(TXTL,FMT='(A5)') '* all'
            ELSE
              WRITE(TXTL,FMT='(I6,A,I6)') KL1,' to ',KL2
            ENDIF
C            IF(LM.EQ.-1) THEN
C              WRITE(TXTL,FMT='(A5)') '* all'
C            ELSE
C              WRITE(TXTL,FMT='(I5)') LM
C            ENDIF
            CALL READC5(KT1,KT2)
            CALL READC1(
     >                  KP1,KP2,KP3)
            CALL READC9(KKEX,KLET)
            WRITE(6,*)
            WRITE(6,*) ' Particle  # :',KT1,' to ',KT2,
     >                              ' ;  element  # : ',TXTL,' ;'
            WRITE(6,*) ' pass # ',KP1,' to ',KP2,', ipass-modulo ',KP3
     >      ,' ;  case  KEX/KTag =', KKEX,'/',KLET
            WRITE(6,*) ' calculation  of  scales  from ',NOC,'  points;'
            WRITE(6,*) ' Xmin-max :',XMI,XMA,KDIM(KX)
            WRITE(6,*) ' Ymin-max :',YMI,YMA,KDIM(KY)
            IF(NOC .EQ. 0) THEN
              WRITE(6,*) 
              WRITE(6,*) '   ***  No  data  could  be  read  ***'
              WRITE(6,*) '     Check particle or element number ' 
              WRITE(6,*) '        -  Main menu, option 2  -'
            ENDIF
            IF(XMI.LT. XMA .AND. YMI.LT. YMA) THEN
C              CALL TXTFBG
              CALL LINTYP(-1)
              DDX=ABS(XMA-XMI)*0.05D0
              DDY=ABS(YMA-YMI)*0.05D0
              CALL TRAXES(XMI,XMA,YMI,YMA,ITWO)
C              CALL TRAXES(XMI-DDX,XMA+DDX,YMI-DDY,YMA+DDY,ITWO)
              OKECH = .TRUE.
            ENDIF
          ENDIF

          IF(OKECH) THEN
            IF(OKVAR) THEN
              IF(KY.EQ. 28) THEN
                CALL PLTHIS(NLOG,NL,NOC,OKBIN,OKECH)
              ELSE
                OK23=((KX.EQ.2  .AND. KY.EQ.3)  .OR.
C                    xxp phase-space
     >               (KX.EQ.12 .AND. KY.EQ.13) .OR.
C                    xxp_o initial phase-space
     >               (KX.EQ.4  .AND. KY.EQ.5)  .OR.
C                    zzp initial phase-space
     >               (KX.EQ.14 .AND. KY.EQ.15) .OR.    
C                    zzp_0 initial phase-space
     >              ((KX.EQ.6  .OR.  KX.EQ.7  .OR. KX.EQ.18)  .AND. 
     >                     (KY.EQ.1 .OR. KY.EQ.19 .OR. KY.EQ.20)) .OR.     
C                   dp/Ekin vs. s/time/phase 
     >              ((KX.EQ.16 .OR.  KX.EQ.17 .OR. KX.EQ.18
     >                         .OR. KX.EQ.6 .OR. KX.EQ.7) .AND. 
     >                     (KY.EQ.11 .OR. KY.EQ.19 .OR. KY.EQ.20 )))
C                   initial dp/Ekin vs. s/time/phase
                OK24=((KX.EQ.2  .AND. KY.EQ.4)  .OR.
C                      Y-Z cross section
     >               (KX.EQ.12 .AND. KY.EQ.14))
C                      Y_o-Z_o cross section
                IF(KX.LT.10 .OR. KX.EQ.18) THEN
                  KPS=1
                ELSE
                  KPS=0
                ENDIF
 
                IF(OK23 .OR. OK24) THEN
                  RLIPS0=RLIPS
 710              CONTINUE      
                  WRITE(6,FMT='(/,
     >              ''  Plot matched ellipse ? ('',A1,''/'',A1,'') : ''
     >              )') RLIPS,RLIPS2
                  READ(5,FMT='(A1)',ERR=710) RLIPS
                  IF(RLIPS.EQ.'n') RLIPS='N'
                  IF(RLIPS.EQ.'y') RLIPS='Y'
                  IF(RLIPS.NE.RLIPS2) THEN
                    RLIPS=RLIPS0
                  ELSE
                    IF(RLIPS2.EQ.'Y') THEN
                      RLIPS2='N'
                    ELSE
                      RLIPS2='Y'
                    ENDIF
                  ENDIF
                  KLIPS = ( (OK23   .AND.  (RLIPS.EQ.'Y')) 
     >                .OR.  (OK24   .AND.  (RLIPS.EQ.'Y'))  
     >                )
                  RHIST0=RHIST
 711              CONTINUE      
                  WRITE(6,FMT='(/,
     >              ''  Plot histograms ? ('',A1,''/'',A1,'') : ''
     >              )') RHIST,RHIST2
                  READ(5,FMT='(A1)',ERR=711) RHIST
                  IF(RHIST.EQ.'n') RHIST='N'
                  IF(RHIST.EQ.'y') RHIST='Y'
                  IF(RHIST.NE.RHIST2) THEN
                    RHIST=RHIST0
                  ELSE
                    IF(RHIST2.EQ.'Y') THEN
                      RHIST2='N'
                    ELSE
                      RHIST2='Y'
                    ENDIF
                  ENDIF
                  KHIST = (  (OK23    .AND. (RHIST.EQ.'Y'))
     >                 .OR.  (OK24      .AND. (RHIST.EQ.'Y'))  )
                ELSE
                  KLIPS = .FALSE.
                  KHIST = .FALSE.
                ENDIF

                CALL PLOT6(KLIPS)
                
                OKREW = .TRUE. 
                CALL PLOTER(NLOG,NL,KPS,NPTS,NPTR,OKREW)
                IF(KLIPS) CALL PLLIPS(NLOG,NLIPS,LM,IPASS,KX,KY,*921)  
                IF(KHIST) CALL PLHIST(NLOG,NL)

              ENDIF
            ELSE
              WRITE(6,*) ' DEFINE  VARIABLES !'
     >        ,' (MENU , OPTION 2)'
            ENDIF
          ENDIF

        ELSE

          CALL OPNWRN(1)

        ENDIF

      GOTO 920
    
 8    CONTINUE
C        CALL MENVCF
        CALL SAVPLT(*920)
      GOTO 921

 12   CONTINUE
        CALL CLSCR
        OKECH=.FALSE.
      GOTO 921

 15   CONTINUE
      WRITE(6,*) 
      WRITE(6,*) ' KILL/CREATE GRAPHIC SUB-PROCESS (0/1):'
      READ(5,*,ERR=15) KC
      IF(KC.EQ. 0) THEN
        CALL ENDVCF
      ELSEIF(KC.EQ. 1) THEN
        CALL BEGVCF
      ENDIF   
      GOTO 921

 20   CONTINUE 
C Superpose a curve
        CALL SUPERP(OKECH,NLOG,*920)
      GOTO 920

 77   CONTINUE 

        IF(OKOPN) THEN

          IF(.NOT.OKVAR) THEN          
            WRITE(6,*) ' DEFINE  VARIABLES !'
     >        ,' (MENU , OPTION 2)'
            GOTO 920
          ENDIF

          IF( .NOT. OKECH) THEN
            WRITE(6,*) 
            WRITE(6,*) ' Busy: calculating scales'
            CALL CALECH(NL,
     >                     NOC,NPASS)
            CALL READC3(KL1,KL2)
            IF(KL1.EQ.-1) THEN
              WRITE(TXTL,FMT='(A5)') '* all'
            ELSE
              WRITE(TXTL,FMT='(I6,A,I6)') KL1,' to ',KL2
            ENDIF
            CALL READC5(KT1,KT2)
            CALL READC1(
     >                  KP1,KP2,KP3)
            CALL READC9(KKEX,KLET)
            WRITE(6,*)
            WRITE(6,*) ' Particle  # :',KT1,' to ',KT2,
     >                              ' ;  element  # : ',TXTL,' ;'
            WRITE(6,*) ' pass # ',KP1,' to ',KP2,', ipass-modulo ',KP3
     >      ,' ;  case  KEX/KTag =', KKEX,'/',KLET
            WRITE(6,*) ' calculation  of  scales  from ',NOC,'  points;'
            WRITE(6,*) ' Xmin-max :',XMI,XMA,KDIM(KX)
            WRITE(6,*) ' Ymin-max :',YMI,YMA,KDIM(KY)
            IF(NOC .EQ. 0) THEN
              WRITE(6,*) 
              WRITE(6,*) '   ***  No  data  could  be  read  ***'
              WRITE(6,*) '     Check particle or element number ' 
              WRITE(6,*) '        -  Main menu, option 2  -'
            ENDIF
            IF(XMI.LT. XMA .AND. YMI.LT. YMA) THEN
              CALL LINTYP(-1)
              DDX=ABS(XMA-XMI)*0.05D0
              DDY=ABS(YMA-YMI)*0.05D0
              CALL TRAXES(XMI,XMA,YMI,YMA,ITWO)
              OKECH = .TRUE.
            ENDIF
          ENDIF

          IF(OKECH) THEN
            IF(KY.EQ. 28) THEN
              CALL PLTHIS(NLOG,NL,NOC,OKBIN,OKECH)
            ELSE
              OK23=((KX.EQ.2  .AND. KY.EQ.3)  .OR.
C                  xxp phase-space
     >             (KX.EQ.12 .AND. KY.EQ.13) .OR.
C                  xxp initial phase-space
     >             (KX.EQ.4  .AND. KY.EQ.5)  .OR.
C                  zzp initial phase-space
     >             (KX.EQ.14 .AND. KY.EQ.15) .OR.    
C                  xxp initial phase-space
     >            ((KX.EQ.6  .OR.  KX.EQ.7  .OR. KX.EQ.18)  .AND. 
     >                   (KY.EQ.1 .OR. KY.EQ.19 .OR. KY.EQ.20)) .OR.     
C                 s/time/phase-dp/Ekin  phase-space
     >            ((KX.EQ.16 .OR.  KX.EQ.17 .OR. KX.EQ.18) .AND. 
     >                   (KY.EQ.11 .OR. KY.EQ.19 .OR. KY.EQ.20 )))
C                 s/time/phase-dp/Ekin initial phase-space
              OK24=((KX.EQ.2  .AND. KY.EQ.4)  .OR.
C                    YZ cross section
     >             (KX.EQ.12 .AND. KY.EQ.14))
C                    Y_oZ_o cross section
              IF(KX.LT.10 .OR. KX.EQ.18) THEN
                KPS=1
              ELSE
                KPS=0
              ENDIF
 
              IF(OK23 .OR. OK24) THEN
                RLIPS0=RLIPS
C 770            CONTINUE      
C                WRITE(6,FMT='(/,
C     >            ''  Plot matched ellipse ? ('',A1,''/'',A1,'') : ''
C     >            )') RLIPS,RLIPS2
C                READ(5,FMT='(A1)',ERR=770) RLIPS
                RLIPS='Y'
                IF(RLIPS.EQ.'n') RLIPS='N'
                IF(RLIPS.EQ.'y') RLIPS='Y'
                IF(RLIPS.NE.RLIPS2) THEN
                  RLIPS=RLIPS0
                ELSE
                  IF(RLIPS2.EQ.'Y') THEN
                    RLIPS2='N'
                  ELSE
                    RLIPS2='Y'
                  ENDIF
                ENDIF
                KLIPS = ( (OK23 .AND. (RLIPS.EQ.'Y')) 
     >              .OR.  (OK24   .AND. (RLIPS.EQ.'Y'))  )
                RHIST0=RHIST
C 771            CONTINUE      
C                WRITE(6,FMT='(/,
C     >            ''  Plot histograms ? ('',A1,''/'',A1,'') : ''
C     >            )') RHIST,RHIST2
C                READ(5,FMT='(A1)',ERR=771) RHIST
                RHIST='Y'
                IF(RHIST.EQ.'n') RHIST='N'
                IF(RHIST.EQ.'y') RHIST='Y'
                IF(RHIST.NE.RHIST2) THEN
                  RHIST=RHIST0
                ELSE
                  IF(RHIST2.EQ.'Y') THEN
                    RHIST2='N'
                  ELSE
                    RHIST2='Y'
                  ENDIF
                ENDIF
                KHIST = (  (OK23    .AND. (RHIST.EQ.'Y'))
     >               .OR.  (OK24    .AND. (RHIST.EQ.'Y'))  )
              ELSE
                KLIPS = .FALSE.
                KHIST = .FALSE.
              ENDIF

              CALL PLOT6(KLIPS)

              KPU1 = KP1
              KPU2 = KP2
              IF(KPU2 .GT.NPASS) KPU2 = NPASS
              KPU3 = KP3                
              KPU = KPU1
              CALL READC2B(KPU,KPU,1)
              CALL CLSCR
              CALL TRAXES(XMI,XMA,YMI,YMA,ITWO)
              OKREW = .TRUE. 
              CALL PLOTER(NLOG,NL,KPS,NPTS,NPTR,OKREW)
              IF(KLIPS) CALL PLLIPS(NLOG,NLIPS,LM,IPASS,KX,KY,*921)  
C              IF(KHIST) CALL PLHIST(NLOG,NL)
 772          CONTINUE
               KPU = KPU + 1
               IF(KPU.LT.KPU2) THEN
                 IF((KPU/KPU3)*KPU3 .EQ. KPU) THEN
                   CALL READC2B(KPU,KPU,1)
                   CALL CLSCR
                   CALL TRAXES(XMI,XMA,YMI,YMA,ITWO)
C                   OKREW = .FALSE. 
                   CALL PLOTER(NLOG,NL,KPS,NPTS,NPTR,OKREW)
                   IF(KLIPS) CALL PLLIPS(NLOG,NLIPS,LM,IPASS,KX,KY,*921)
C                   IF(KHIST) CALL PLHIST(NLOG,NL)
                 ENDIF
                 GOTO 772
               ENDIF
               CALL READC2B(KPU,KPU,1)
               CALL CLSCR
               CALL TRAXES(XMI,XMA,YMI,YMA,ITWO)
               CALL PLOTER(NLOG,NL,KPS,NPTS,NPTR,OKREW)
               IF(KLIPS) CALL PLLIPS(NLOG,NLIPS,LM,IPASS,KX,KY,*921)  
C               IF(KHIST) CALL PLHIST(NLOG,NL)

C Reset KP1-3
               CALL READC2B(KP1,KP2,KP3)

            ENDIF
          ENDIF

        ELSE

          CALL OPNWRN(1)
        ENDIF

      GOTO 920
    
 88   CONTINUE 
        CALL HELP('1/7')
      GOTO 920

 99   CONTINUE
C      CALL CLMITV       
      IF(OKOPN) THEN
        CLOSE(NL) 
        OKOPN = .FALSE.
      ENDIF
      RETURN 

 997  WRITE(6,*) ' Error upon OPEN ',FNAM
      RETURN 

      END
