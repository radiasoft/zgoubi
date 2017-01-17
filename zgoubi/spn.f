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
      SUBROUTINE SPN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "C.DONT.H"     ! COMMON/DONT/ TA(MXL,MXTA)
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      INCLUDE "C.OBJET.H"     ! COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      INCLUDE "C.PTICUL.H"     ! COMMON/PTICUL/ AM,Q,G,TO
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      INCLUDE "C.SPIN.H"     ! COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)
 
      CHARACTER(1) KAX(3)
 
      SAVE SXM,SYM,SZM
 
      LOGICAL ONCE
      SAVE ONCE

      INTEGER DEBSTR, FINSTR
      LOGICAL IDLUNI, OK
      CHARACTER(300) TXT,TXT2
      LOGICAL GTTEXT

      DATA KAX / 'X' , 'Y' , 'Z' /
      DATA ONCE / .FALSE. /
 
      KSO = NINT(A(NOEL,1))
      KSO2 = NINT( 10.D0*A(NOEL,1) - 10.D0*DBLE(KSO) )
 
      ONCE = KSO .GE. 1  .OR. ONCE
 
      CALL REBELR(KREB3,KDUM,KDUM)
      IF(  KREB3 .EQ. 99
     >.AND.  KSO .NE .0
     >.AND.  KSO .NE. -1
     >                   ) THEN
C       ... SET TO 99 IN SBR REBELOTE - FOR PERIODIC MACHINES
        IF(NRES.GT.0) WRITE(NRES,103)
 103    FORMAT(//,15X,
     >  'Final  spins  of  last  run  taken  as  initial  spins')
        RETURN
      ENDIF
 
 
      IF(NRES.GT.0) THEN
        IF    (KSO .EQ. 0) THEN
          WRITE(NRES,107)
 107      FORMAT(/,15X,' KSO=0 :  Spin  tracking  off. ',/)
        ELSEIF(KSO .EQ. -1) THEN
          IF(ONCE) THEN
            WRITE(NRES,FMT=
     >     '(/,15X,''  KSO=-1 :  Spin  tracking  resumes. '',/)')
          ELSE
            WRITE(NRES,FMT=
     >     '(/,15X,''  KSO=-1  is inoperant in this context.'',
     >     '' It needs be preceded at least once by'',/,15X,
     >     '' KSO .ge. 1. That is necessary for launching the spin '',
     >     '' tracking machinery. '',/)')
          ENDIF
        ELSE
          WRITE(NRES,110) AM, G
 110      FORMAT(/,15X,' SPIN  TRACKING  REQUESTED  ',1P
     >    ,//,25X,' PARTICLE  MASS          = ',G15.7,' MeV/c2'
     >    , /,25x,' GYROMAGNETIC  FACTOR  G = ',G15.7)
 
          WRITE(NRES,111) KSO
 111      FORMAT(/,25X,' INITIAL SPIN CONDITIONS TYPE ',I2,' :')
 
          IF    (KSO .LE. 3) THEN
            WRITE(NRES,101) KAX(KSO)
 101        FORMAT(30X,' ALL PARTICLES HAVE SPIN PARALLEL TO  ',A1
     >      ,'  AXIS')
          ELSEIF(KSO .EQ. 4) THEN
            IF    (KSO2 .EQ.0) THEN
              WRITE(NRES,104) NINT(A(NOEL,9))
 104          FORMAT(
     >           30X,'All spins entered particle by particle'
     >        ,/,30X,'Particles # 1 to ',I7,' may be subjected to spin '
     >        ,      'matching using FIT procedure')
            ELSEIF(KSO2 .EQ.1) THEN
              WRITE(NRES,113) NINT(A(NOEL,9))
 113          FORMAT(
     >           30X,'Same spin for all particles'
     >        ,/,30X,'Particles # 1 to ',I7,' may be subjected to spin '
     >        ,      'matching using FIT procedure')
            ELSEIF(KSO2 .EQ.2) THEN
              WRITE(NRES,112) 
     >        TA(NOEL,1)(DEBSTR(TA(NOEL,1)):FINSTR(TA(NOEL,1)))
 112          FORMAT(
     >           30X,'Same spin for all particles, read from file ',A)
            ENDIF
          ELSEIF(KSO .EQ. 5) THEN
            WRITE(NRES,108)
 108        FORMAT(15X,' OPTION 5 UNAVAILABLE IN THIS VERSION',/)
            KSO=0
            RETURN
         ENDIF
 
          P = BORO*CL9*Q
          BE = P/SQRT(P*P + AM*AM)
          GG = G/SQRT(1.D0-BE*BE)
          WRITE(NRES,102) BORO,BE,GG
 102      FORMAT(
     >    //,25X,' PARAMETRES  DYNAMIQUES  DE  REFERENCE :'
     >    ,/,30X,' BORO   =  ',F12.3,' KG*CM'
     >    ,/,30X,' BETA   =  ',F10.6,/,30X,' GAMMA*G = ',F10.6)
 
        ENDIF
      ENDIF
 
      IF    (KSO.EQ.-1) THEN
        IF(ONCE) THEN
          KSPN = 1
        ELSE
          KSPN = 0
        ENDIF
      ELSEIF(KSO.EQ.0) THEN
        KSPN = 0
      ELSEIF(KSO.GE.1) THEN
        IF(AM .EQ. 0.D0) THEN
          IF(NRES.GT.0) WRITE(NRES,106)
 106      FORMAT(//,15X,' SVP  INDIQUER  LA  MASSE  DES  PROJECTILES !'
     >         ,/,15X,' - UTILISER  LE  MOT-CLE  ''PARTICUL''',/)
          CALL ENDJOB('SBR SPN. Need to provide particle mass.',-99)
        ENDIF
        KSPN = 1
      ENDIF
 
      GOTO(1,1,1,4,5) KSO
      RETURN
 
 1    CONTINUE
      SX = 0D0
      SY = 0D0
      SZ = 0D0
      IF(KSO .EQ. 1) SX = 1.D0
      IF(KSO .EQ. 2) SY = 1.D0
      IF(KSO .EQ. 3) SZ = 1.D0
      DO 11 I=1,IMAX
        SI(1,I) = SX
        SI(2,I) = SY
        SI(3,I) = SZ
        SI(4,I) = 1.D0
        SF(1,I) = SX
        SF(2,I) = SY
        SF(3,I) = SZ
        SF(4,I) = 1.D0
 11   CONTINUE
      GOTO 98
 
 4    CONTINUE
      IF    (KSO2.EQ.0) THEN
        IM = IMAX
        IF(IM.GT.MXD/10) IM=MXD/10
        IA = 0
        DO I=1,IM
          IA = IA+10
          SX = A(NOEL,IA)
          SY = A(NOEL,IA+1)
          SZ = A(NOEL,IA+2)
          SI(1,I) = SX
          SI(2,I) = SY
          SI(3,I) = SZ
          SI(4,I) = SQRT(SX*SX+SY*SY+SZ*SZ)
          SF(1,I) = SX
          SF(2,I) = SY
          SF(3,I) = SZ
          SF(4,I) = SI(4,I)
c           write(*,*) ' spn ',i,sx,sy,sz 
       ENDDO
c        DO I=IM+1,IMAX
c        DO I=1,IMAX
c          SX = SI(1,I)
c          SY = SI(2,I)
c          SZ = SI(3,I)
c          SI(4,I) = SQRT(SX*SX+SY*SY+SZ*SZ)
c          SF(1,I) = SX
c          SF(2,I) = SY
c          SF(3,I) = SZ
c          SF(4,I) = SI(4,I)
c        ENDDO
      ELSEIF(KSO2.EQ.1) THEN
        SX = A(NOEL,10)
        SY = A(NOEL,11)
        SZ = A(NOEL,12)
        DO I=1,IMAX
          SI(1,I) = SX
          SI(2,I) = SY
          SI(3,I) = SZ
          SI(4,I) = SQRT(SX*SX+SY*SY+SZ*SZ)
          SF(1,I) = SX
          SF(2,I) = SY
          SF(3,I) = SZ
          SF(4,I) = SI(4,I)
        ENDDO
      ELSEIF(KSO2.EQ.2) THEN
        OK = IDLUNI(
     >              LR)
        OPEN(UNIT=LR,FILE=TA(NOEL,1))
        ok = gttext(NRES,LR,'(deg)       (deg)',
     >              txt2)
c        do i = 1, 11
c            write(*,*) ' spn ',txt2
c                  read(*,*)
          READ(LR,FMT='(A)') txt
c        enddo
c            write(*,*) ' spn ',txt
c                  read(*,*)
        READ(txt,*) txt2,SXi, SYi, SZi,smi,SX, SY, SZ,sm
c            write(*,*) ' spn ',txt
c            write(*,*) ' SX, SY, SZ ',SX, SY, SZ
c                  read(*,*)
        CLOSE(LR)
        DO I=1,IMAX
          SI(1,I) = SX
          SI(2,I) = SY
          SI(3,I) = SZ
          SI(4,I) = SQRT(SX*SX+SY*SY+SZ*SZ)
          SF(1,I) = SX
          SF(2,I) = SY
          SF(3,I) = SZ
          SF(4,I) = SI(4,I)
        ENDDO
      ELSE
        CALL ENDJOB('*** Pgm spn, keyword SPNTRK : '// 
     >' No such option KSO2 =',KSO2)
      ENDIF
      GOTO 98
 
 5    CONTINUE
 
      TO = A(NOEL,10)
      PO = A(NOEL,11)
      AL = A(NOEL,20)
      DA = A(NOEL,21)
 
      SX = COS(PO)*COS(TO)
      SY = COS(PO)*SIN(TO)
      SZ = SIN(PO)
      DO 15 I=1,IMAX
        SI(1,I) = SX
        SI(2,I) = SY
        SI(3,I) = SZ
        SI(4,I) = 1.D0
        SF(1,I) = SX
        SF(2,I) = SY
        SF(3,I) = SZ
        SF(4,I) = 1.D0
 15   CONTINUE
 
 
 98   CONTINUE
 
      SXM = 0.D0
      SYM = 0.D0
      SZM = 0.D0
      DO 20 I=1,IMAX
        SXM = SXM + SI(1,I)
        SYM = SYM + SI(2,I)
        SZM = SZM + SI(3,I)
 20   CONTINUE
 
      IF(NRES.GT.0) THEN
        RIMX = DBLE(IMAX)
        SM = SQRT(SXM*SXM+SYM*SYM+SZM*SZM)/RIMX
        WRITE(NRES,120) IMAX, SXM/RIMX, SYM/RIMX, SZM/RIMX, SM
 120    FORMAT(//,25X,' POLARISATION  INITIALE  MOYENNE  DU'
     >  ,2X,'FAISCEAU  DE  ',I7,'  PARTICULES :'
     >  ,/,30X,' <SX> = ',F10.4
     >  ,/,30X,' <SY> = ',F10.4
     >  ,/,30X,' <SZ> = ',F10.4
     >  ,/,30X,' <S>  = ',F10.4)
      ENDIF
 
      RETURN
      END
