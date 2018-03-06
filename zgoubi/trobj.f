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
C  Brookhaven National Laboratory                    és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE TROBJ (IOPT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ------------------
C     TRANSFORMATION
C     ------------------
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      INCLUDE "C.CHAMBR.H"     ! COMMON/CHAMBR/ LIMIT,IFORM,YLIM2,ZLIM2,SORT(MXT),FMAG,YCH,ZCH
 
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE "C.DESIN.H"     ! COMMON/DESIN/ FDES(7,MXT),IFDES,KINFO,IRSAR,IRTET,IRPHI,NDES
C     >,AMS,AMP,AM3,TDVM,TETPHI(2,MXT)
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      INCLUDE "C.OBJET.H"     ! COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT,KZOB
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,HDPRF,DP,QBR,BRI
 
      I1=INT(A(NOEL,1))
      I2=INT(A(NOEL,2))
      ALP=A(NOEL,3)
      BET=A(NOEL,4)
      GAM=A(NOEL,5)
      IF(IOPT .EQ. 1) THEN
        IF(NRES.GT.0) WRITE(NRES,109) I1,I2,ALP,BET,GAM
 109    FORMAT(/,10X,'TRANSFORMATION ',I2,1X,I2,' COEF ',
     S         F12.5,1X,F12.5,1X,F12.5/)
      ENDIF
 
      DO 1 I=1,IMAX
 
         FO(I1,I)=ALP*FO(I1,I)+BET*FO(I2,I)+GAM
         F(I1,I)=FO(I1,I)
 
 1    CONTINUE
 
 
      IF(NRES.GT.0) WRITE(NRES,101) IEX(1),(F(J,1),J=1,7)
  101 FORMAT(' TRAJ 1 IEX,D,Y,T,Z,P,S,time :',I3,1P,5G12.4,2G17.5)

      RETURN
      END
