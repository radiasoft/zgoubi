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
C  François Méot <meot@lpsc.in2p3.fr>
C  Service Accélerateurs
C  LPSC Grenoble
C  53 Avenue des Martyrs
C  38026 Grenoble Cedex
C  France
      SUBROUTINE SCALIN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER*80 TA
      COMMON/DONT/ TA(MXL,20)
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      INCLUDE 'MXFS.H'
      COMMON/SCAL/SCL(MXF,MXS),TIM(MXF,MXS),NTIM(MXF),KSCL
      CHARACTER FAM*8,LBF*8,KLEY*10,LABEL*8
      COMMON/SCALT/ FAM(MXF),LBF(MXF,2),KLEY,LABEL(MXL,2)
 
      PARAMETER(MST=MXLF+1)
      LOGICAL EMPTY, IDLUNI

      PARAMETER (ND=200000)
      DIMENSION XM(ND), YM(ND), freq(nd), ekin(nd), turn(nd)

      CHARACTER*9 OPT(2)

      DATA OPT/ '++ OFF ++','         ' /
 
      IOPT = A(NOEL,1)
      NFAM = A(NOEL,2)

      KSCL=0
      IF(IOPT*NFAM.NE.0)  KSCL=1
 
      IF(NRES .GT. 0) WRITE(NRES,100) NFAM, OPT(IOPT+1)
 100  FORMAT(/,15X,' Scaling  applied  on',I3
     >,'  families  of  optical  elements',/,T40,A9)

      IF(IOPT .EQ. 0) RETURN
 
      DO 1 IF=1,NFAM

        IF(NRES .GT. 0) THEN
          WRITE(NRES,101) FAM(IF),(LBF(IF,KL),KL=1,MXLF)
 101      FORMAT(/,10X,'Element [/label] to be scaled :',/,
     >                                   15X,A,5X,2('/',A))

          IF( EMPTY(LBF(IF,1)) ) THEN
             WRITE(NRES,FMT='(15X,
     >       ''Family not labeled ;  this scaling will be apply'',
     >       '' to all elements  "'',A,''"'')') FAM(IF)
          ENDIF
        ENDIF

        IF(NTIM(IF) .GE. 0) THEN

          DO IT = 1, NTIM(IF) 
            SCL(IF,IT) = A(NOEL,10*IF+IT-1)
            TIM(IF,IT) = A(NOEL,10*IF+NTIM(IF)+IT-1)
          ENDDO

          IF(NRES .GT. 0) THEN
            WRITE(NRES,102) NTIM(IF), (SCL(IF,IT) ,IT=1,NTIM(IF))
 102        FORMAT(20X,' # TIMINGS :', I2
     >          ,/,20X,' SCALING   : ',10(F11.4,1X))
            WRITE(NRES,103) ( INT(TIM(IF,IT)), IT=1,NTIM(IF))
 103        FORMAT(20X,' TURN #    : ',10(I8,1X))
          ENDIF

        ELSEIF(NTIM(IF) .EQ. -1) THEN
C--------- Scaling is taken from CAVITE (ENTRY CAVIT1)
C          Starting value is SCL(IF,1)
          SCL(IF,1) = A(NOEL,10*IF)

          IF(NRES .GT. 0) THEN
            WRITE(NRES,FMT='(15X,''Scaling of fields follows ''
     >      ,''increase of rigidity taken from CAVITE'')')
            WRITE(NRES,FMT='(15X,''Starting scaling value is ''
     >      ,1P,E14.6)') SCL(IF,1)
          ENDIF

        ELSEIF(NTIM(IF) .EQ. -2) THEN
C--------- Field law for scaling FFAG, LPSC, Sept. 2007

          IF(IPASS.EQ.1) THEN
            IF(IDLUNI(
     >                LUN)) THEN
              OPEN(UNIT=LUN,
     >            FILE='zgoubi.freqLaw.In',STATUS='OLD',ERR=597)
            ELSE
              WRITE(NRES,*) ' Tried  to  open  zgoubi.freqLaw.In...'
              GOTO 596
            ENDIF

            I = 1
 11         CONTINUE
C Array XM must be  increasing function of index
C              READ(LUN,*,ERR=599,END=598) XM(I),YM(I)
C Read                 turn#,  freq,   phi,   oclock, Ekin 
              READ(LUN,*,ERR=599,END=598) 
     >                 turn(i),freq(i),YM(I), XM(I),  ekin(i)
              IF(I.GT.ND) STOP ' SBR SCALIN, too many data'
c              WRITE(88,fmt='(1p,4e14.6,2x,a)') turn,freq, yM(I), xM(I), 
c     >                   '  sbr scalin turn freq phase oclock '
              I = I+1
            GOTO 11

 598        WRITE(*,*) '  READ ended upon EOF, # data read is ',I-1
            GOTO 21
 599        WRITE(*,*) '  READ ended upon ERROR, # data read is ',I-1
            GOTO 21

 21         CONTINUE
         
              N = I-1
              WRITE(*,*) ' # of data to be interpolated :  N= ',N
              IF(N.EQ.0) STOP ' SBR SCALIN, no data to be interpolated'

              IF(XM(2).LT.XM(1)) STOP 
     >        'SBR SCALIN, array X must be increasing function of index'

             DUM = SCALE6(XM,YM,turn,freq,ekin,N)           

          ENDIF

        ELSEIF(NTIM(IF) .EQ. -77) THEN
C--------- Field law proton driver, FNAL, Nov.2000
          SCL(IF,1) = A(NOEL,10*IF)
          SCL(IF,2) = A(NOEL,10*IF+1)
          SCL(IF,3) = A(NOEL,10*IF+2)
          SCL(IF,4) = A(NOEL,10*IF+3)
          TIM(IF,1) = A(NOEL,10*IF+4)
          TIM(IF,2) = A(NOEL,10*IF+5)
          IF(NRES .GT. 0) THEN
            WRITE(NRES,112) (SCL(IF,IC) ,IC=1,4)
 112        FORMAT(20X,' Brho-min, -max, -ref,  Frep  : ', 4(F17.8,1X))
            WRITE(NRES,113) ( INT(TIM(IF,IT)), IT=1,2)
 113        FORMAT(20X,' TURN #    : ',I8,'  TO  ',I8)
          ENDIF

        ELSEIF(NTIM(IF) .EQ. -88) THEN
C--------- AC dipole at BNL
          SCL(IF,1) = A(NOEL,10*IF)       ! C
          SCL(IF,2) = A(NOEL,10*IF+1)     ! Q1
          SCL(IF,3) = A(NOEL,10*IF+2)     ! Q2
          SCL(IF,4) = A(NOEL,10*IF+3)     ! P
          TIM(IF,1) = A(NOEL,10*IF+4)     ! Nramp
          TIM(IF,2) = A(NOEL,10*IF+5)     ! Nflat
          TIM(IF,3) = A(NOEL,10*IF+6)     ! Ndown
          IF(NRES .GT. 0) THEN
            WRITE(NRES,FMT='(5X,1P,''C, Q1, Q2, P :'',4(1X,E14.6))') 
     >          (SCL(IF,IC) ,IC=1,4)
            WRITE(NRES,FMT='(5X,'' N-ramp, -flat, -down : '',3I8)') 
     >          (NINT(TIM(IF,IT)), IT=1,3)
          ENDIF

        ELSEIF(NTIM(IF) .EQ. -87) THEN
C--------- Scaling is taken from CAVITE (ENTRY CAVIT1), as for NTIM=-1. 
C          Starting value is SCL(IF,1). 
C          Switch-on is triggered by G.gamma value : trigger is at G.gamma=Integer+dN (dN=0.75-frac(Qx), 
C                with normally frac(Qx)~0.73 hence dN~0.25)
          SCL(IF,1) = A(NOEL,10*IF)
          TIM(IF,1) = A(NOEL,10*IF+1)     ! G.gamma value at start of jump 
          TIM(IF,2) = A(NOEL,10*IF+2)     ! dN
          TIM(IF,3) = A(NOEL,10*IF+3)     ! # of turns on up and on down ramps (~40)

          IF(NRES .GT. 0) THEN
            WRITE(NRES,FMT='(15X,''Scaling of field in Q-jump quads ''
     >      ,''follows increase of rigidity taken from CAVITE'')')
            WRITE(NRES,FMT='(15X,''Starting scaling value is ''
     >      ,1P,E14.6)') SCL(IF,1)
            WRITE(NRES,FMT='(5X,
     >      ''Quad jump series is started at G.gamma = N + dN ='',
     >      F10.4,'' + '',F10.4)') TIM(IF,1),TIM(IF,2)
            WRITE(NRES,FMT=
     >      '(5X,''Number of turns of up- and down-ramps is : '',I4)') 
     >      NINT(TIM(IF,3))

          ENDIF

        ENDIF

 1    CONTINUE
 
      RETURN

 597  CONTINUE
      WRITE(6,FMT='(3A)') ' SBR SCALIN : COULD NOT OPEN FILE ',
     > 'zgoubi.freqLaw.In','.  Check existence...'
      STOP 
 596  CONTINUE
      STOP ' SBR SCALIN :  IDLE  UNIT  PROBLEM  AT  OPEN  FILE '
      RETURN
      END
