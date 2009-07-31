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
      SUBROUTINE MCDESI(*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ---------------------------------------------
C     Initialise parametres necessairy for managing
C      decays in flight. 
C     --------------------------------------------
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST/ CL9,CEL,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE "MAXTRA.H"
      COMMON/DESIN/ FDES(7,MXT),IFDES,KINFO,IRSAR,IRTET,IRPHI,NDES
     >,AMS,AMP,AM3,TDVM,TETPHI(2,MXT)
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),IMAX,IEX(MXT),IREP(MXT),AMQLU
      CHARACTER LET
      COMMON/FAISCT/ LET(MXT)
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IMAXD,IMAXT
      COMMON/PTICUL/ AM,Q,G,TO
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      COMMON/RIGID/ BORO,DPREF,DP,BR
      COMMON/UNITS/ UNIT(MXJ)
 
      SAVE CTOLAB,E3STAR,B3STAR
 
C----- Mass of secondary particle to be tracked
      AMS= A(NOEL,1)
C----- Mass of secondary particle to be abandonned
      AM3= A(NOEL,2)
C----- Life time (s) of secondary particle. Unused if zero
      TOS= A(NOEL,3)

      IF(IPASS .EQ. 1) THEN
C       *** Initialisation at 1st pass only
        IFDES = 1
        IRSAR = A(NOEL,10)
        IRTET = A(NOEL,11)
        IRPHI = A(NOEL,12)
        IRSAR =(IRSAR /2)*2+1
        IRTET =(IRTET /2)*2+1
        IRPHI =(IRPHI /2)*2+1
      ENDIF 

      IF(AMS.LT.0.D0) THEN
C       ** To WRITE decay related info
        KINFO = 1
        AMS = -AMS
      ELSE
        KINFO=0
      ENDIF
 
      AMS2 = AMS*AMS
      AM32 = AM3*AM3
 
C------- Sort life time for each parent particle, 
C      in a law   N/N0 = EXP(-S/S0) , S=fligth distance in cm
      TDVM=0.D0
      SE = 0.D0
      SB  = 0.D0
      STO = 0.D0
      SCTO = 0.D0
      IT = 0
      DO 1 I=1,IMAX
        IF(IEX(I).LT.-1) GOTO 1
C------- Life time (s) of secondary particle
        FDES(7,I) = TOS
        IF(LET(I).EQ.'S') GOTO 1
        AMP = AMQ(1,I)
        TO = AMQ(4,I)
        IF(AMP*TO .EQ. 0.D0) THEN
          WRITE(NRES,106)
 106      FORMAT(//,15X,' Give  mass  &  TO  of  projectiles !'
     >           ,/,15X,' - use  keyword  ''PARTICUL''',/)
          RETURN 1
        ENDIF
        IT = IT + 1
        U=RNDM()
        IF(U .EQ. 0.D0) U=1.D0
        AMP2 = AMP * AMP
        ENSTAR = .5D0*(AMP2+AMS2-AM32)/AMP
        BSTAR  = SQRT(1.D0-AMS2/ENSTAR/ENSTAR)
        SE = SE + ENSTAR
        SB  = SB + BSTAR
        STO = STO + TO
        CTOLAB = 1.D-7*CEL*CEL * BORO * TO / AMP
        SCTO = SCTO + CTOLAB
        FDES(6,I) = -CTOLAB*F(1,I)* LOG(U)
        TETPHI(1,I) = ACOS(1.D0 - 2.D0*RNDM())
        TETPHI(2,I) = 2.D0*PI* RNDM()
        TDVM=TDVM+FDES(6,I)
CC---------- Test densities ("ligne verte")//
C        ppi = F(1,I)*BORO*CL9 *Q
C        ppi2 = ppi * ppi
C        Epi = SQRT(ppi2+AMP2)
C        btapi = ppi / Epi
C        qapa = (amp2+ams2)/(2.D0 * amp2)
C        ctmu1 =-(qapa * Epi + ams/sqrt(1.d0-btapi**2))/((1.d0-qapa)*ppi)
C        ctmu2 =-(qapa * Epi - ams/sqrt(1.d0-btapi**2))/((1.d0-qapa)*ppi)
C        TETPHI(1,I) = ACOS(ctmu2)
C         write(77,*) i, TETPHI(1,I), ppi, ctmu1,ctmu2, 'SBR MCDESI'
CC-----------------------------------------//
 1    CONTINUE
      TDVM=TDVM/IT
      ENSTAR = SE/IT
      BSTAR  = SB/IT
      TO = STO/IT
      CTOLAB =SCTO/IT
      NDES=0
 
      IF(NRES.GT.0) THEN
        WRITE(NRES,100)
 100    FORMAT(20X,'            IN  FLIGHT  DECAY  SET ',
     >  //,20X,' Mean life path, Theta & phi sorted at random',/)
 
        IF(KOBJ .EQ. 1) WRITE(NRES,104)
 104    FORMAT(20X,' !! ATTENTION,  OBJET/KOBJ=1 :  Po<0',2X,
     >  '&  Zo<0  are  obtained by symmetry, not  computed  !!',/)
 
        E3STAR = .5D0*(AMP2+AM32-AMS2)/AMP
        B3STAR = SQRT(1.D0-AM32/E3STAR/E3STAR)
        WRITE(NRES,101) AMP, TO*CEL, CTOLAB*UNIT(5), AMS, ENSTAR ,BSTAR
     >  ,AM3, E3STAR, B3STAR
 101    FORMAT(/,10X,
     >  '  Parent  particle,  theoretical  data :  ',/,
     >  15X,1P,'(#1)   M1 =',
     >  G12.4,' MeV/c2 ;  Average life dist. in com  = ',G12.4,' m',
     >  25X,'i.e.,  average life dist. in lab  =  ',G12.4,' m  ',/,
     >  10X,'  Secondary  particles,  theoretical  data :',
     >  /,15X,'(#2)   M2 =',G12.4,' MeV/c2',6X,'E* =',G12.4,' MEV',
     >  6X,'BETA* =',G12.5,
     >  /,15X,'(#3)   M3 =',G12.4,' MeV/c2',6X,'E* =',G12.4,' MEV',
     >  6X,'BETA* =',G12.5)
C
        WRITE(NRES,102) IT,TDVM*UNIT(5)
102     FORMAT(/,10X,'  Monte Carlo  data :  ',/,
     >  15X,1P,'Average  mean  life  time  of  ',I5,
     >  '  parent  particles : ',G12.4,'  m',/)
C
        IF(KINFO .EQ. 1) THEN
          WRITE(NRES,103) IT
103       FORMAT(/,' Mean  life  time  of  the  ',I5,
     >                           '  particles  (m) :',/)
          WRITE(NRES,105) (FDES(6,I)*UNIT(5),I=1,IMAX)
105       FORMAT(5X,1P,10G12.4)
        ENDIF

        IF(TOS.NE.0.D0) WRITE(NRES, FMT='(/,10X,
     >  ''  Secondary  particles have limited life-time,  '',1P,G12.4,
     >  '' s (com) ;   they wiil be subject to decay process'')') TOS

      ENDIF
C
      RETURN
      END
