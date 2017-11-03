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
      SUBROUTINE DPKICK(DPKCK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ------------------
C     SECTION SANS Champ
C     ------------------
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      INCLUDE "C.OBJET.H"     ! COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      INCLUDE "C.REBELO.H"   ! COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      INCLUDE 'MXFS.H'
      INCLUDE 'MXSCL.H'
      INCLUDE "C.SCAL.H"     ! COMMON/SCAL/ SCL(MXF,MXS,MXSCL),TIM(MXF,MXS),NTIM(MXF),KSCL
C      COMMON/SCAL/SCL(MXF,MXS),TIM(MXF,MXS),NTIM(MXF),JPA(MXF,MXP),KSCL
 
      DATA DTA1 / 0.D0 /

C     ... FACTEUR D'ECHELLE DES ChampS. UTILISE PAR 'SCALING'
      SCAL = SCAL0()
      IF(KSCL .EQ. 1) SCAL = SCAL0()*SCALER(IPASS,NOEL,
     >                                                 DTA1)

        IF(NRES.GT.0) THEN
          WRITE(NRES,109) DPKCK
 109      FORMAT(/,30X,'dP kick =',F12.5,'  (relative)',/)
        ENDIF
        DO 1 I=1,IMAX
C---------- IEX <-1 <=> PARTICULE STOPPEE
           IF(IEX(I) .LT. -1) GOTO 1
 
           F(1,I)=F(1,I)+DPKCK*SCAL
 1      CONTINUE
      RETURN
      END
