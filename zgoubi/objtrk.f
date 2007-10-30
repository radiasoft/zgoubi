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
      SUBROUTINE OBJTRK(VECT) 
C     ------------------------------------------------
C     An object for developping the TRACKING procedure
C     ------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VECT(6)
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      COMMON/CHAMBR/ LIMIT,IFORM,YLIM2,ZLIM2,SORT(MXT),FMAG,BMAX
     > ,YCH,ZCH
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IIP(MXL),NB,NOEL
      CHARACTER*80 TA
      COMMON/DONT/ TA(MXL,20)
      INCLUDE "MAXCOO.H"
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),IMAX,IEX(MXT),IREP(MXT)
      CHARACTER LET
      COMMON/FAISCT/ LET(MXT)
      CHARACTER  KAR(41)
      COMMON/KAR/ KAR
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      COMMON/SYNCH/ RET(MXT), DPR(MXT),PS

      CHARACTER LETI

      ntr = 1

      Q = 1.D0
      P0 = BORO*CL9*Q

      YFAC = 1D0
      TFAC = 1D0
      ZFAC = 1D0  
      PFAC = 1D0
      SFAC = 1D0
      DPFAC= 1D0
      TIFAC= 1D0
      YREF = 0D0
      TREF = 0D0
      ZREF = 0D0
      PREF = 0D0
      SREF = 0D0
      DPREF= 0D0
      TIREF= 0D0
      DPREF=1.D0

C----- Traj. counter
      I = 0
      IT1 = 0
      IKAR = 0
      IPASS1 = 0
      IEND = 0
  17  CONTINUE
        I = I + 1
        IT1 = IT1 + 1
            IKAR = IKAR+1
            IF(IKAR.GT.41)  IKAR=1
c            READ(NL,*,ERR=97,END=95) Y,T,Z,P,S, DP
            Y = VECT(IT1)
            T = VECT(IT1)
            Z = VECT(IT1)
            P = VECT(IT1)
            S = VECT(IT1)
            DP = VECT(IT1)
            TIM = 0.D0
            LETI=KAR(IKAR)
            IEXI=1
            IT = IT1
            IREPI = IT
            BRO = BORO
            YO= Y
            TO= T
            ZO= Z
            PO= P
            SO= S
            DPO= DP
            TIMO= TIM

        LET(IT1)=LETI
        IEX(IT1)=IEXI
C        FO(1,IT1)=1.D0 + DPO
C        FO(1,IT1)=(1.D0 + DPO) * BRO/BORO
        FO(1,IT1)= DPO * BRO/BORO
        FO(2,IT1)=YO
        FO(3,IT1)=TO
        FO(4,IT1)=ZO
        FO(5,IT1)=PO
        FO(6,IT1)=SO
        FO(7,IT1)=TIMO
        F(1,IT1)=(DP*DPFAC + DPREF) * BRO/BORO
        F(2,IT1)=  Y*YFAC  + YREF
        F(3,IT1)=  T*TFAC  + TREF
        F(4,IT1)=  Z*ZFAC  + ZREF
        F(5,IT1)=  P*PFAC  + PREF
        F(6,IT1)=  S*SFAC  + SREF
        F(7,IT1)=TIM*TIFAC + TIREF 
        RET(IT1)=RETI
        DPR(IT1)=DPRI
C        IREP(IT1) = IREPI
        IREP(IT1) = IT1
        SORT(IT1) = SORTI
        AMQ(1,IT1) = AMQ1
        AMQ(2,IT1) = AMQ2
        AMQ(3,IT1) = AMQ3
        AMQ(4,IT1) = AMQ4
        AMQ(5,IT1) = AMQ5

        IF(IT1 .EQ. NTR) GOTO  169
      GOTO  17
 
 169  CONTINUE
C Changed for linear FFAG sudies
C      IMAX=IT
      IMAX=IT1
C-----
      IF(NRES .GT. 0) 
     >   WRITE(NRES,*) ntr,' PARTICLE COORDINATES RECORDED'
      IDMAX = 1
      IMAXT=IMAX

      RETURN
      END
