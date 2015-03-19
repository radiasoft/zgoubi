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
      SUBROUTINE PATH(SL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ------------------------------------------------------
C     Each particle flies a drift distance  to 'PATH' length
C     ------------------------------------------------------
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
      CHARACTER(1) LET
      COMMON/FAISCT/ LET(MXT)
      COMMON/GASC/ AI, DEN, KGA
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
C----- Conversion  coord. (cm,mrd) -> (m,rd)
      INCLUDE "C.UNITS.H"     ! COMMON/UNITS/ UNIT(MXJ)
  
      IF( SL .NE. 0D0 ) THEN

        DO 1 I=1,IMAX
C---------- IEX < -1 <=> PARTICLE STOPPEE
           IF(IEX(I) .LT. -1) GOTO 1
            
           IF(LET(I).NE.'S') THEN
C------------ Particle is not yet decayed
             UL = FDES(6,I)-F(6,I)
             IF(UL .LE. SL) THEN
C-------------- Particle will decay...
               DL = UL + 1.D-10
               XL = DL*(COS(F(3,I)*UNIT(2))*COS(F(5,I)*UNIT(4)))
               CALL ESL(0,XL,I,I)    
C-------------- Particle has decayed !
               XL = (SL-DL)*(COS(F(3,I)*UNIT(2))*COS(F(5,I)*UNIT(4)))
             ELSE
               XL =  SL*(COS(F(3,I)*UNIT(2))*COS(F(5,I)*UNIT(4)))
             ENDIF
           ELSE
             XL =  SL*(COS(F(3,I)*UNIT(2))*COS(F(5,I)*UNIT(4)))
           ENDIF
           CALL ESL(0,XL,I,I)
 
 1      CONTINUE

      ENDIF

      RETURN
      END
