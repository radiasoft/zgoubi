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
      SUBROUTINE MCDESI(*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ---------------------------------------------
C     Initialise parametres necessairy for managing
C      decays in flight.
C     --------------------------------------------
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.CONST_2.H"     ! COMMON/CONST/ CL9,CL,PI,RAD,DEG,QEL,AMPROT,CM2M
      INCLUDE "MAXTRA.H"
      INCLUDE "C.DESIN.H"     ! COMMON/DESIN/ FDES(7,MXT),IFDES,KINFO,IRSAR,IRTET,IRPHI,NDES
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
      CHARACTER(1) LET
      INCLUDE "C.FAISCT.H"     ! COMMON/FAISCT/ LET(MXT)
      INCLUDE "C.OBJET.H"     ! COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT,KZOB
      INCLUDE "C.REBELO.H"   ! COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,HDPRF,DP,QBR,BRI
      INCLUDE "C.UNITS.H"     ! COMMON/UNITS/ UNIT(MXJ)

      SAVE VTOLAB,E3STAR,B3STAR

C----- Mass of secondary particle to be tracked
      AMS= A(NOEL,1)
      IF(AMS.LE.0.D0) CALL ENDJOB('Wrong M2 value :  should be > 0',-99)
C----- Mass of secondary particle to be abandonned
      AM3= A(NOEL,2)
C----- Life time (s) of secondary particle. Zero if to be abandonned.
      TOS= A(NOEL,3)
      KINFO = NINT(A(NOEL,4))

      IF(IPASS .EQ. 1) THEN
C       *** Initialisation at 1st pass only
        IFDES = 1
        IRSAR = NINT(A(NOEL,10))
        IRTET = NINT(A(NOEL,11))
        IRPHI = NINT(A(NOEL,12))
        IRSAR =(IRSAR /2)*2+1
        IRTET =(IRTET /2)*2+1
        IRPHI =(IRPHI /2)*2+1
      ELSE
        CALL REBELR(KREB3,KREB31,KDUM)
        IF(KREB3.EQ.99) RETURN
      ENDIF

      AMS2 = AMS*AMS
      AM32 = AM3*AM3

C------- Sort life time for each parent particle,
C      in a law   N/N0 = EXP(-S/S0) , S=fligth distance in cm
      TDVM=0.D0
      SE = 0.D0
      SB  = 0.D0
      STO = 0.D0
      SVTO = 0.D0
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
        VTOLAB = 1.D-7*CL*CL * BORO * TO / AMP
        SVTO = SVTO + VTOLAB
        FDES(6,I) = -VTOLAB*F(1,I)* LOG(U)
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
      VTOLAB =SVTO/IT
      NDES=0

      IF(NRES.GT.0) THEN
        WRITE(NRES,100)
 100    FORMAT(20X,'            IN-FLIGHT  DECAY  SWITCHED ON ! ',/)

        IF(KOBJ .EQ. 1) WRITE(NRES,104)
 104    FORMAT(20X,'!! ATTENTION, if using OBJET/KOBJ=1 :'
     >  ,/,20X,'trajectories with initial  Po<0'
     >  ,'&  Zo<0  are derived by symmetry, not from  ray-tracing  !!')

        E3STAR = .5D0*(AMP2+AM32-AMS2)/AMP
        B3STAR = SQRT(1.D0-AM32/E3STAR/E3STAR)
        WRITE(NRES,101) AMP, TO, VTOLAB*UNIT(5)
     >  , AMS, ENSTAR ,BSTAR, TOS, AM3, E3STAR, B3STAR
 101    FORMAT(1P,//,10X,' Parent  particle,  theoretical  data :  '
     >  ,/,15X,' M1 =',G12.4,' MeV/c2 '
     >  ,/,15X,' theoretical life time at rest, tau = ',G12.4,' s'
     >  ,/,15X,' i.e.,  life  distance  in  lab,  v*tau  = ',G12.4,' m'
     >  ,/,10X,' Daughter  particle,  theoretical  data :'
     >  ,/,15X,' M2 =',G12.4,' MeV/c2, ',6X,'E* =',G12.4,' MeV'
     >  ,            6X,'beta* =',G12.5
     >  ,/,15X,' theoretical life time at rest, tau_2 = ',G12.4,' s'
     >  ,/,15X,' M3 =',G12.4,' MeV/c2, ',6X,'E* =',G12.4,' MeV'
     >  ,'      beta* =',G12.5)

        WRITE(NRES,102) IT,TDVM*UNIT(5), TDVM/VTOLAB
102     FORMAT(/,10X,'  Monte Carlo  outcome  :  '
     >  ,/,15X,1P,'Average  life  distance  of  the ',I5,
     >  '  parent  particles : ',G12.4,' (m)     / theor. life dist. ->'
     >  ,G12.4,/)

        IF(KINFO .EQ. 1) THEN
          SFD6 = 0.D0
          II = 0
          DO I = 1, IMAX
            IF(IEX(I).GE.1) THEN
              II = II+1
              SFD6 = SFD6 + FDES(6,I)
            ENDIF
          ENDDO
          SFD6 = SFD6 * UNIT(5) / II
          WRITE(NRES,103) IT
103       FORMAT(/,'  Individual  life  distances  of  the  ',I5,
     >                           '  particles  (m) :',/)
          WRITE(NRES,105) (I,FDES(6,I)*UNIT(5),I=1,IMAX)
105       FORMAT(5(2X,I6,'/',2X,F13.4))
          WRITE(NRES,FMT='(1P,''   Average : '',G12.4,'' m'')') SFD6
        ENDIF

        WRITE(NRES, FMT='(/
     >  ,/,10X,'' Daughter particle M2 will be tracked and subject to''
     >  ,         '' decay process, ''
     >  ,/,10X,'' Daughter particle M3 will be abandonned. '')')

      ENDIF
C
      RETURN
      END
