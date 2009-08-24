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
      SUBROUTINE SPNPRN(KPR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER*80 TA
      COMMON/DONT/ TA(MXL,20)
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      CHARACTER LET
      COMMON/FAISCT/ LET(MXT)
      COMMON/PTICUL/ AM,Q,G,TO
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)
 
      IF(KPR .GT. 1) THEN
        IF(IPASS.GT.1) THEN
          IF(KPR*(IPASS/KPR) .NE. IPASS) RETURN
        ENDIF
      ENDIF
 
      IF(IPASS .EQ. 1) CALL OPEN2('SPNPRN',NSPN,TA(NOEL,1))
 
      DO 10 I=1,IMAX
        P = BORO*CL9 *F(1,I) * AMQ(2,I)
        GA = SQRT(P*P/(AMQ(1,I)*AMQ(1,I)) + 1.D0)
        WRITE(NSPN,101) LET(I),IEX(I),(SI(J,I),J=1,4),(SF(J,I),J=1,4)
     >  ,(GA-1.D0)*AMQ(1,I),I,IMAX,IPASS,NOEL            
 101    FORMAT(1X,A1,1X,I2,1X,1P,8(1X,E15.7)
     >  ,/,E15.7,3(1X,I7),1X,I5)
 10   CONTINUE
 
      RETURN
      END
