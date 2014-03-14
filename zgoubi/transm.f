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
      SUBROUTINE TRANSM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ***********************************
C     TRANSFERT MATRICIEL AU SECOND ORDRE
C     APPELE PAR KLEY 'TRANSMAT'
C     ***********************************
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      COMMON/PTICUL/ AM,Q,G,TO
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      COMMON/SYNCH/ PH(MXT), DPR(MXT), PS
      COMMON/TRNSM/ R(6,6),T(6,6,6)
      COMMON/UNITS/ UNIT(MXJ)

      PARAMETER (MXJ1=MXJ-1)
      DIMENSION FF(MXJ1),FC(MXJ1)
      SAVE P0
 

      INCLUDE 'MXFS.H'
      COMMON/SCALP/ VPA(MXF,MXP),JPA(MXF,MXP)
      INCLUDE 'MXSCL.H'
      dimension kfm(MXSCL)

      DUM = SCALER(1, NOEL, 
     >     DUM)
      CALL SCALE9(
     >            KFM )
      do ifm = 1, MXSCL
            IF(KFM(ifm) .le. 0) goto 121
        DO I= 1 , JPA(KFM(ifm),MXP)
           A(NOEL,JPA(KFM(ifm),I)) = VPA(KFM(IFM) ,I)
        ENDDO
      enddo
 121  continue

C     *** R EST AU FORMAT 6 x 6
CModified, FM, 04/97
C         T EST AU FORMAT 6 x 6 x 6
C         UNITES  MKSA
 
      IORD = A(NOEL,1)
      XL = A(NOEL,10)
      IIB = 20
      DO 5  IA=1, 6
        DO 5 IB=1, 6
          R(IA,IB) = A(NOEL,IIB)
          IIB = IIB + 1
 5     CONTINUE
 
C Modified, FM, 04/97
C      IF(IORD .EQ. 2) THEN
CC       ** CHANGE MATRICE TRIANGULAIRE EN CARREE
C        DO 11 IA=1,6
C          DO 11 IB=1,5
C            IB1=IB+1
C            DO 11 IC=IB1,6
C              T(IA,IC,IB)=T(IA,IB,IC)
C 11     CONTINUE
C      ENDIF
 
      IF(NRES.GT.0) THEN
        WRITE(NRES,105) IORD,XL
 105    FORMAT(/,20X,8('+'),3X,'TRANSFERT  MATRICIEL  A  L''ORDRE  ',
     >  I1,3X,8('+'), 1P,
     >  /,20X,'Length  of  element  :',G12.4,' m',/)
        I=1
        WRITE(NRES,100) I
 100    FORMAT(//,15X,'TRANSFER  MATRIX  ORDRE',I3,'  (MKSA units)',/)
        WRITE(NRES,101) (( R(IA,IB) , IB=1,6) , IA=1,6)
 101    FORMAT(6X,6F13.6)
        IF(IORD .EQ. 2) THEN
          WRITE(NRES,100) IORD
          DO 16 IA=1,6
            IF(IA.GT.1) WRITE(NRES,106)
 106        FORMAT(/)
            DO 16 IB=1,6
              WRITE(NRES,102) ( IA,IC,IB, T(IA,IB,IC) , IC=1,6)
 102        FORMAT( 6(I4,I2,I1,1PE11.3) )
 16       CONTINUE
        ENDIF
      ENDIF
 
C     ... ACCELERATION WITH 'CAVITE': PS=SYNCHRONOUS MOMENTUM
      IF(IPASS .EQ. 1) THEN
C        P0 = BORO*CL*1.D-9 *Q/QE
C        P0 = BORO*CL*1.D-9 *Q
        P0 = BORO*CL9 *Q
        PS=P0
C       ... ELSE: PS IS UP DATED BY 'CAVITE'
      ENDIF
 
      DO 1 I=1,IMAX
        FF(6)=F(1,I)*P0/PS-1.D0
        FF(1)=F(2,I)
        FF(2)=F(3,I)
        FF(3)=F(4,I)
        FF(4)=F(5,I)
        DO 2 IA=1,6
          FC(IA)=0.D0
          DO 2 IB=1,6
            FC(IA)=FC(IA) + FF(IB) * R(IA,IB)/UNIT(IA)*UNIT(IB)
            IF(IORD .EQ. 1) GOTO 2
C            DO 3 IC=1,IB
            DO 3 IC=1,6
              FC(IA)=FC(IA) +
     >          FF(IB)*FF(IC)*T(IA,IB,IC)/UNIT(IA)*UNIT(IB)*UNIT(IC)
 3          CONTINUE
 2      CONTINUE
        F(1,I)=(1.D0+FC(6))*PS/P0
        F(2,I)=FC(1)
        F(3,I)=FC(2)
        F(4,I)=FC(3)
        F(5,I)=FC(4)
C       ... Path length :
        F(6,I)=F(6,I) + FC(5) + XL/UNIT(5)
C       ... Time of flight update ????
C        F(7,I)= ??

 1    CONTINUE
 
      RETURN
      END
