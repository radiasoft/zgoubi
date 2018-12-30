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
C  USA
C  -------
      SUBROUTINE ESL(IOPT,XL,I1,I2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ------------------
C     SECTION SANS Champ
C     ------------------
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      INCLUDE "C.CHAMBR.H"     ! COMMON/CHAMBR/ LIMIT,IFORM,YLIM2,ZLIM2,SORT(MXT),FMAG,YCH,ZCH

      INCLUDE "C.CONST_2.H"     ! COMMON/CONST/ CL9,CL,PI,RAD,DEG,QEL,AMPROT,CM2M
      INCLUDE "C.DESIN.H"     ! COMMON/DESIN/ FDES(7,MXT),IFDES,KINFO,IRSAR,IRTET,IRPHI,NDES
C     >,AMS,AMP,AM3,TDVM,TETPHI(2,MXT)
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      INCLUDE "C.GASC.H"     ! COMMON/GASC/ AI, DEN, KGA
      INCLUDE "C.OBJET.H"     ! COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT,KZOB
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,HDPRF,DP,QBR,BRI
      INCLUDE "C.SPIN.H"     ! COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)
C----- Conversion  coord. (cm,mrd) -> (m,rd)
      INCLUDE "C.UNITS.H"     ! COMMON/UNITS/ UNIT(MXJ)

      DL=XL
      CALL SCUMW(XL)
      CALL SCUMR(
     >           XL,SCUM,TCUM)

      IF(IOPT .EQ. 1)THEN
        IF(NRES.GT.0) THEN
          WRITE(NRES,109) DL
 109      FORMAT(/,30X,'Drift,  length =',F12.5,'  cm',/)
C 109      FORMAT(/,30X,'ESPACE  LIBRE =',F12.5,'  CM',/)
        ENDIF
      ENDIF

C----- COMPTAGE / LIMITES CHAMBRE
C ?? Do not test limit here due to possible decay etc. Only at end of ESL. ??
      IF(LIMIT .EQ. 1) THEN
        IF(DL.GE.0.D0) CALL CHMBR(I1,I2)
      ENDIF

      IF( DL .NE. 0D0 ) THEN

C        DO 1 I=1,IMAX
        DO 1 I=I1, I2
C---------- IEX < -1 <=> PARTICLE STOPPEE
           IF(IEX(I) .LT. -1) GOTO 1

           F(2,I)=F(2,I)+DL*TAN(F(3,I)*1.D-3)
           F(4,I)=F(4,I)+DL*TAN(F(5,I)*1.D-3)*(1.D0/COS(F(3,I)*1.D-3))
           DS = DL/(COS(F(3,I)*1.D-3)*COS(F(5,I)*1.D-3))
           F(6,I)= F(6,I) + DS
           P = BORO*F(1,I)*CL9*AMQ(2,I)
           AMI = AMQ(1,I)
           IF(AMI*P.NE.0.D0) F(7,I) = F(7,I) +
     >         (DS*1.D4 / (P/SQRT(P*P+AMI*AMI)*CL))
 1      CONTINUE

      ENDIF

C-------- DESINTEGRATION EN COURS DE VOL
C         DOIT PRECEDER LE TEST CHMBR !!
      IF(IFDES .EQ. 1) THEN
        IF(DL.GT.0.D0) THEN
          CALL MCDESL(DL,I1,I2)
          IF(KGA .EQ. 1) THEN
            IF(NRES.GT.0) THEN
              WRITE(NRES,*)
              WRITE(NRES,*) ' *** WARNING, subroutine ESL'
              WRITE(NRES,*) '   mixing of gas-scattering and decay in '
              WRITE(NRES,*) '   flight not correctly implemented'
              WRITE(NRES,*)
            ENDIF
          ENDIF
        ELSE
          IF(NRES.GT.0) THEN
            WRITE(NRES,*)
            WRITE(NRES,*) ' *** WARNING, subroutine ESL'
            WRITE(NRES,*) '    No correction to decay in negative drift'
            WRITE(NRES,*)
          ENDIF
        ENDIF
      ENDIF

C-------- Gas-scattering
      IF(KGA .EQ. 1) THEN
        IF(DL.GT.0D0) CALL GASESL(DL,I1,I2)
      ENDIF

C-------- COMPTAGE / LIMITES CHAMBRE
      IF(LIMIT .EQ. 1) THEN
        IF(DL.GE.0.D0) CALL CHMBR(I1,I2)
      ENDIF

      IF(IOPT .EQ. 1)THEN
       IF(NRES.GT.0) THEN
         WRITE(NRES,101) IEX(1),-1.D0+F(1,1),(F(J,1),J=2,7)
 101     FORMAT('TRAJ #1 IEX,D,Y,T,Z,P,S,time :',
     >   I3,1P,5E14.6,1X,E15.7,1X,E13.5)
         IF(KSPN.EQ.1) WRITE(NRES,102) IEX(1),(SF(I,1),I=1,4)
 102     FORMAT('TRAJ #1 SX, SY, SZ, |S| :',1X,I2,2X,1P,4(E14.6,1X))
        ENDIF
      ENDIF

      IF(NRES .GT. 0)
     >  WRITE(NRES,FMT='(/,'' Cumulative length of optical axis = '',
     >  1P,G17.9,'' m  '',
     >  '' ;  Time  (for reference rigidity & particle) = '',
     >  1P,G14.6,'' s '')')  SCUM*UNIT(5), TCUM

      RETURN
      END
