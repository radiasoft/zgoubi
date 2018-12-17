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

      SUBROUTINE SPACH
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.CDF.H"      ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.CONST.H"    ! COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      LOGICAL AMQLU(5),PABSLU
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"      ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      INCLUDE "C.FAISC.H"    ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT)
      INCLUDE "C.PTICUL.H"   ! COMMON/PTICUL/ AM,Q,G,TO
      INCLUDE "C.REBELO.H"   ! COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      INCLUDE "C.RIGID.H"    ! COMMON/RIGID/ BORO,DPREF,DP,QBR,BR
      logical TSPCH
      INCLUDE "C.SPACECHA.H" ! COMMON/SPACECHA/ TLAMBDA,RBEAM(2),XAVE(2),EMITT(2),TAVE,BUNCH_LEN,
C                             >                EMITTZ, BTAG, SCKX, SCKY, TSPCH
          dimension XPave(2)


      pi = 4.d0 * atan (1.d0)
      c = 299 792 458 !m / s

c      CALL SPA1(
c     >             TLAMBDA1,TSPCH1)


      IF(NRES.GT.0) THEN
         IF  (TSPCH)  THEN
          WRITE(NRES,FMT='(/,A,ES15.8,A,/)')
     >    'The linear charge density is: ', tlambda,'*10**(10) e.m-1'

         ELSE
          WRITE(NRES,107)
 107      FORMAT(/,15X,' Space charge is off. ',/)
         ENDIF
      ENDIF


C----- P0, AM  are  in  MEV/c, /c^2
       P0 = BORO*CL9*Q

C----- The momentum of the reference particle is used to give the space charge kick:
C      YOU SHOULD TRACK ONLY ONE DISTRIBUTION WHEN YOU WANT TO COMPUTE SPACE CHARGE EFFECTS.
C----- Be careful if you have RF you have to check the particle momentum

      P = P0*F(1,1)
      AM2 = AM*AM
      ENRG = SQRT(P*P+AM2)
      BTA = P/ENRG
      BTA2 = BTA**2
      GAMMA=1/SQRT(1-BTA2)
      BTAG = BTA*GAMMA
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccc    compute the perveance term    cccccccccccccccccccccc

!      epsilon0 = 8.854187817*10**(-12)   in Faraday
      GAMMA3 = GAMMA**3
      fact1=5526.347
!      Q_perv = Q*tlambda*(10**6)/(2*pi*epsilon0*AM*GAMMA3*BTA2)
      Q_perv = Q*tlambda/(2*pi*fact1*AM*GAMMA3*BTA2)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c here I need to apply the space charge kick to each slice: the assumption is that the beam size does not change a lot, in order to apply a uniform space charge kick to each slice -- to be updated


      DO 10 JJ = 1, 2
         CALL RMS(JJ,
     >                           EMIT,ALP,BET,XM,XPM,rx)
         Xave(JJ)=XM
         XPave(JJ)=XPM
         Rbeam(JJ)=rx*100.0   ! zgoubi unity is in cm
         Emitt(JJ)=EMIT
 10   CONTINUE

cccccccccccccccccc  the longitudinal rms beamsize is computed here -- JJ=3 activates this calculation
         JJ=3
         CALL RMS(JJ,
     >                           Emittz,ALP2,BET2,Tave,XPM2,Bunch_len)


cccccccccccccccccccccccccccccc   compute the linear space charge kick here   cccccccccccccccccccccccccccc
      if ((Rbeam(1) .GT. 0.) .AND. (Rbeam(2) .GT. 0.) .AND.
     >      (TSPCH)) THEN
       SCkx=2*Q_perv/((Rbeam(1)+Rbeam(2))*Rbeam(1))!*1/(1+Bunch_len*1000)
       SCky=2*Q_perv/((Rbeam(1)+Rbeam(2))*Rbeam(2))!*1/(1+Bunch_len*1000)
       open(unit=lunWW,file='perv.out')    ! malek
             write(lunWW,*) Q_perv
       close(lunWW)
c         WRITE(NRES,FMT='(/,A,/)')
c     > 'Summary: Space charge calculation ended up normally'
      else
         SCkx=0.0
         SCky=0.0
c         WRITE(NRES,FMT='(/,A,/)')
c     > 'Space charge is not applied: the beam does not have any radius'

      endif

      RETURN

      END
