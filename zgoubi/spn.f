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
      SUBROUTINE SPN(*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU,PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      COMMON/PTICUL/ AM,Q,G,TO
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)
 
      CHARACTER KAX(3)
 
      SAVE SXM,SYM,SZM
 
      DATA KAX / 'X' , 'Y' , 'Z' /
 
      CALL REBELR(KREB3,KREB31)
      IF(KREB3 .EQ. 99) THEN
C       ... SET TO 99 IN SBR REBELOTE - FOR PERIODIC MACHINES
        IF(NRES.GT.0) WRITE(NRES,103)
 103    FORMAT(//,15X,
     >  'FINAL  SPINS  OF  LAST  RUN  TAKEN  AS  INITIAL  SPINS')
        RETURN
      ENDIF
 
      KSPN = 1
      KSO= A(NOEL,1)
 
      IF(NRES.GT.0) THEN
 
        IF(KSO .EQ. 0) THEN
          WRITE(NRES,107)
 107      FORMAT(/,15X,' SPIN  TRACKING  OFF ',/)
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
            WRITE(NRES,104)
 104        FORMAT(30X,' ALL SPINS ENTERED PARTICLE BY PARTICLE')
          ELSEIF(KSO .EQ. 5) THEN
            WRITE(NRES,108)
 108        FORMAT(15X,' OPTION 5 UNAVAILABLE IN THIS VERSION',/)
            KSO=0
            RETURN
         ENDIF
 
C          P = BORO*CL*1D-9*(Q/QE)
C          P = BORO*CL*1D-9*Q
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
 
      IF(AM .EQ. 0.D0) THEN
        WRITE(NRES,106)
 106    FORMAT(//,15X,' SVP  INDIQUER  LA  MASSE  DES  PROJECTILES !'
     >         ,/,15X,' - UTILISER  LE  MOT-CLE  ''PARTICUL''',/)
        RETURN 1
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
      DO 14 I=1,IMAX
        SX = SI(1,I)
        SY = SI(2,I)
        SZ = SI(3,I)
        SI(4,I) = SQRT(SX*SX+SY*SY+SZ*SZ)
        SF(1,I) = SX
        SF(2,I) = SY
        SF(3,I) = SZ
        SF(4,I) = SI(4,I)
 14   CONTINUE
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
        SM = SQRT(SXM*SXM+SYM*SYM+SZM*SZM)/IMAX
        WRITE(NRES,120) IMAX,SXM/IMAX,SYM/IMAX,SZM/IMAX,SM
 120    FORMAT(//,25X,' POLARISATION  INITIALE  MOYENNE  DU'
     >  ,2X,'FAISCEAU  DE  ',I3,'  PARTICULES :'
     >  ,/,30X,' <SX> = ',F10.4
     >  ,/,30X,' <SY> = ',F10.4
     >  ,/,30X,' <SZ> = ',F10.4
     >  ,/,30X,' <S>  = ',F10.4)
      ENDIF
 
      RETURN
      END
