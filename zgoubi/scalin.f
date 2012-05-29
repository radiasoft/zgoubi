C  ZGOUBI, a program for computing the trajectories of charged particles
C  in electric and magnetic fields
C  Copyright (C) 1988-2007  Fran�ois M�ot
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
C  Fran�ois M�ot <fmeot@bnl.gov>
C  Brookhaven National Laboratory                    �s
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE SCALIN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER*80 TA
      COMMON/DONT/ TA(MXL,40)
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      INCLUDE 'MXFS.H'
      COMMON/SCAL/SCL(MXF,MXS),TIM(MXF,MXS),NTIM(MXF),KSCL
      PARAMETER (LBLSIZ=8)
      PARAMETER (KSIZ=10)
      CHARACTER FAM*(KSIZ),LBF*(LBLSIZ),KLEY*(KSIZ),LABEL*(LBLSIZ)
      COMMON/SCALT/ FAM(MXF),LBF(MXF,2),KLEY,LABEL(MXL,2)
 
      PARAMETER(MST=MXLF+1)
      LOGICAL EMPTY, IDLUNI

      PARAMETER (ND=200000)
      DIMENSION XM(ND), YM(ND), freq(nd), ekin(nd), turn(nd)

      CHARACTER*9 OPT(2)
      DIMENSION MODSCI(MXF),MODSCL(MXF)
      SAVE MODSCL

      DIMENSION SCL2(MXF,MXD), TIM2(MXF,MXD), NTIM2(MXF)
      DIMENSION SCL2I(MXF,MXD), TIM2I(MXF,MXD), NTIM2I(MXF)
      SAVE NTIM2, SCL2, TIM2

      DATA OPT/ '++ OFF ++','         ' /
 
      IOPT = NINT(A(NOEL,1))
      NFAM = NINT(A(NOEL,2))

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
     >       ''Family not labeled ;  this scaling will apply'',
     >       '' to all elements  "'',A,''"'')') FAM(IF)
          ENDIF
        ENDIF

        IF(NTIM(IF) .GE. 0) THEN

          IF    (MODSCL(IF) .LT. 10) THEN

            DO IT = 1, NTIM(IF) 
               SCL(IF,IT) = A(NOEL,10*IF+IT-1)
               TIM(IF,IT) = A(NOEL,10*IF+NTIM(IF)+IT-1)
            ENDDO

          ELSEIF(MODSCL(IF) .EQ. 10) THEN

            IF(IPASS.EQ.1) THEN

              DO IT = 1, NTIM(IF) 
                SCL(IF,IT) = SCL(IF,IT) * SCL(IF,MXS) 
              ENDDO

            ENDIF

          ELSEIF(MODSCL(IF) .EQ. 11) THEN

            IF(IPASS.EQ.1) THEN

              IIT = NTIM2(IF) 
              DO IT = 1, 10
                 SCL2(IF,IT) = A(NOEL,10*(2*IF-1)+IT-1)
                 TIM2(IF,IT) = A(NOEL,10*(2*IF  )+IT-1)
C        write(66,*) ' scalin ',scl2(iF,it),if,it
              ENDDO

              if(iit.eq.1 .and. TIM2(IF,1).eq.0.d0) 
     >                         TIM2(IF,1) =  BORO 

cC               if(if.eq.1 .or. if.eq.2)   then
c               if(if.eq.1 .or. if.eq.2)   then
c                  write(*,*) ' scalin '
c                  write(*,*) ' A(10*if / scl2) /  IF :  ',if,noel
c                  write(*,*)   10*(2*IF-1),'-',10*(2*IF-1)+IIT-1
c                  write(*,*)  (A(NOEL,10*(2*IF-1)+IT-1),IT=1,IIT)
c                  write(*,*) ' A(10*if... / tim2) : '
c                  write(*,*)   10*(2*IF  ),'-',10*(2*IF  )+IIT-1
c                  write(*,*)  (A(NOEL,10*(2*IF  )+IT-1),IT=1,IIT)
c                   write(*,*) ' '
c                endif

              call scale2(
     >                    SCL2,TIM2,NTIM2,IF)


            ENDIF

          ENDIF

          IF(NRES .GT. 0) THEN

            IF(MODSCL(IF) .LT. 10) THEN

              WRITE(NRES,102) NTIM(IF), (SCL(IF,IT) ,IT=1,NTIM(IF))
 102          FORMAT(20X,' # TIMINGS :', I2
     >          ,/,20X,' SCALING   : ',10(F11.4,1X))
              WRITE(NRES,103) ( INT(TIM(IF,IT)), IT=1,NTIM(IF))
 103          FORMAT(20X,' TURN #    : ',10(I12,1X))

            ELSE

              WRITE(NRES,104) 
     >        TA(NOEL,IF)(1:LEN(TA(NOEL,IF))),
     >        NTIM(IF), TIM(IF,1), TIM(IF,NTIM(IF)),
     >        SCL(IF,1), SCL(IF,NTIM(IF))
 104          FORMAT(15X,'Scaling of field follows the law'
     >        ,' function of Rigidity taken from file : ', A
     >        ,/,20X,' # TIMINGS : ', I5
     >        ,/,20X,' From :      ', ES9.2, ' To ', ES9.2, ' kG.cm'
     >        ,/,20X,' Scal from : ', ES9.2, ' To ', ES9.2)

              IF    (MODSCL(IF) .EQ. 10) THEN
                WRITE(NRES,FMT='(10X,
     >          ''Scaling law from file is further multiplied by ''
     >          ,''constant factor = '',1P,E15.6)') SCL(IF,MXS) 

              ELSEIF(MODSCL(IF) .EQ. 11) THEN
                WRITE(NRES,FMT='(10X,
     >          ''Scaling law from file is further multiplied by ''
     >          ,''Brho-dependent scaling SCL2'')')
              ENDIF

            ENDIF

          ENDIF

        ELSEIF(NTIM(IF) .EQ. -1) THEN
C--------- Scaling is taken from CAVITE (ENTRY CAVIT1)
C          Starting value is either SCL(IF,1) or BORO
          SCL(IF,1) = A(NOEL,10*IF) 
CCCCCCCCCCCCCCC          IF(SCL(IF,1) .EQ. 0.D0)  SCL(IF,1) =  BORO * 1.D-3

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

        ELSEIF(NTIM(IF) .EQ. -60) THEN
C--------- Scaling is taken from CAVITE (ENTRY CAVIT1)
C          Starting value is SCL(IF,1)
          SCL(IF,1) = A(NOEL,10*IF)

          IF(NRES .GT. 0) THEN
            WRITE(NRES,FMT='(15X,''Scaling of fields follows ''
     >      ,''increase of rigidity taken from CAVITE'')')
            WRITE(NRES,FMT='(15X,''Starting scaling value is ''
     >      ,1P,E14.6)') SCL(IF,1)
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
 113        FORMAT(20X,' TURN #    : ',I12,'  TO  ',I12)
          ENDIF

        ELSEIF(NTIM(IF) .EQ. -88) THEN
C--------- AC dipole at  BNL
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

C      ENTRY SCALI4(
C     >             NTIM2I,JF)
C      NTIM2(JF) = NTIM2I(JF)
C      RETURN

      ENTRY SCALI4(
     >             SCL2I,TIM2I,NTIM2I,JF)
      NTIM2(JF) = NTIM2I(JF)
      DO IT = 1, NTIM2(JF)
        SCL2(JF,IT) = SCL2I(JF,IT)
        TIM2(JF,IT) = TIM2I(JF,IT)
      ENDDO
      RETURN

      ENTRY SCALI6(MODSCI)
      DO I = 1, MXF 
        MODSCL(i) = modsci(i)
      ENDDO
      RETURN

      ENTRY SCALI5(
     >              MODSCI)
      DO I = 1, MXF 
        MODSCI(I) = MODSCL(I)
      ENDDO
      RETURN

      END
