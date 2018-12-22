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
      SUBROUTINE PRSR(LNSR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      INCLUDE "C.CDF.H"       ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL,PI,RAD,DEG,QE,AMPROT,CM2M
      INCLUDE "C.DESIN.H"     ! COMMON/DESIN/ FDES(7,MXT),IFDES,KINFO,IRSAR,IRTET,IRPHI,NDES
C     >,AMS,AMP,AM3,TDVM,TETPHI(2,MXT)
C     1,AMS ,AMP,ENSTAR,BSTAR,TDVM,TETPHI(2,MXT)
      PARAMETER (LBLSIZ=20)
      CHARACTER(LBLSIZ) LABEL
      INCLUDE 'MXLD.H'
      INCLUDE "C.LABEL.H"     ! COMMON/LABEL/ LABEL(MXL,2)
      INCLUDE "MAXCOO.H"
      INCLUDE "C.OBJET.H"     ! COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT,KZOB
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      CHARACTER(1) LET
      INCLUDE "C.FAISCT.H"    ! COMMON/FAISCT/ LET(MXT)
C      LOGICAL ZSYM
      INCLUDE "C.TYPFLD.H"    ! COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      INCLUDE "C.PTICUL.H"    ! COMMON/PTICUL/ AMASS,Q,G,TO
      INCLUDE "C.REBELO.H"    ! COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,HDPRF,DP,QBR,BRI
      INCLUDE "C.UNITS.H"     ! COMMON/UNITS/ UNIT(MXJ)

C      DIMENSION SIG(4,4)
      CHARACTER(5) TXT(3)
      CHARACTER(10) UU(3)

      DIMENSION EMIT(3),ALP(3),BET(3),XM(3),XPM(3)
      DIMENSION NLIV(3),NINL(3),RATIN(3)
      DIMENSION SIGX(3), SIGXP(3)
      PARAMETER (KSIZ=10)
      CHARACTER(KSIZ) KLE
      SAVE DEAV,PHEAV,PHERMS

      CHARACTER(LBLSIZ) TLBL1, TLBL2
      LOGICAL EMPTY

      DATA TXT / '(Y,T)', '(Z,P)', '(t,K)' /
      DATA UU / '(cm,rd)', '(cm,rd)', '(mu_s,MeV)' /
      DATA DE, SIGE2 / 0.d0, 0.d0 /

Compute rms ellipse
      PI = 4.D0 * ATAN(1.D0)
      DO JJ = 1, 3
        CALL LPSFIT(JJ,
     >                 EMIT(jj),ALP(jj),BET(jj),XM(jj),XPM(jj))
Compute number of particles alive and numberinside ellipse
        CALL CNTINL(JJ,PI*EMIT(jj),ALP(jj),BET(jj),XM(jj),XPM(jj),
     >                                        NLIV(jj),NINL(jj))
        RATIN = DBLE(NINL(jj))/DBLE(IMAX)
        SIGX(JJ) = SQRT(EMIT(JJ) * BET(JJ))
        SIGXP(JJ) = SQRT(EMIT(JJ) * (1.D0+ALP(JJ)**2)/BET(JJ))
      ENDDO

      CALL ZGNOEL(
     >            NOEL)
      CALL ZGKLEY(
     >             KLE)
      TEMP = EMIT(3)/PI * (1.D0+ALP(3)**2)/BET(3)

      TLBL1=LABEL(NOEL,1)
      IF(EMPTY(TLBL1)) TLBL1 = 'none'
      TLBL2=LABEL(NOEL,2)
      IF(EMPTY(TLBL2)) TLBL2 = 'none'

      WRITE(LNSR,FMT='(1P,
     >3(1X,A),2(1X,I6),5(1X,E14.6),1X,I6,
     >5(1X,E14.6),2(1X,I6),1X,E14.6,
     >5(1X,E14.6),2(1X,I6),1X,E14.6,
     >5(1X,E14.6),2(1X,I6),1X,E14.6,
     >1X,E14.6,
     >5(1X,E14.6)
     >)')
     >KLE,TLBL1,TLBL2,NOEL,IPASS,BORO,DPREF+HDPRF,AMASS,Q,G,IMAX,
C Beam properties
     >PI*EMIT(1),ALP(1),BET(1),XM(1),XPM(1),NLIV(1),NINL(1),RATIN(1),
     >PI*EMIT(2),ALP(2),BET(2),XM(2),XPM(2),NLIV(2),NINL(2),RATIN(2),
     >PI*EMIT(3),ALP(3),BET(3),XM(3),XPM(3),NLIV(3),NINL(3),RATIN(3),
C local contribution to dE
     >XPM(3) - DE,
C local contribution to sig_E
     >SQRT(TEMP -SIGE2),
C Theoretical data, from field specs
     >DEAV,PHEAV,PHERMS

      DE = XPM(3)
      SIGE2 = TEMP

      RETURN

      ENTRY PRSR2(DEAVI,PHEAVI,PHERMI)
      DEAV = DEAVI
      PHEAV = PHEAVI
      PHERMS = PHERMI
      RETURN

      END
