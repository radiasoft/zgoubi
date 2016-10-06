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
C  François Meot <fmeot@bnl.gov>
C  Brookhaven National Laboratory  
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  -------
      SUBROUTINE AGSBLW(MOD2,NOEL,DEV,NBLW,NB,WN,WA,
     >                                              BM)
C     ----------------------------------------------
C     NBLW : # of blwg ; WN(i=1,nblw) : # of turns ; 
C     WA(i=1,nblw) : amperes.   
C     ----------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION WN(*),WA(*),BM(NB)
      INCLUDE 'MXFS.H'
      INCLUDE 'MXLD.H'
      PARAMETER (LBLSIZ=20)
      CHARACTER(LBLSIZ) LABEL
      INCLUDE "C.LABEL.H"     ! COMMON/LABEL/ LABEL(MXL,2)
      PARAMETER (KSIZ=10)
      CHARACTER(KSIZ) FAM
      CHARACTER(LBLSIZ) LBF
      INCLUDE "C.SCALT.H"     ! COMMON/SCALT/ FAM(MXF),LBF(MXF,MLF)

      CHARACTER(LBLSIZ) MMNM 
      PARAMETER(TURNMM=16.D0)
      INTEGER DEBSTR, FINSTR

      MMNM = LABEL(NOEL,1)
      MMNM = MMNM(DEBSTR(MMNM):FINSTR(MMNM))

      IF    (MOD2 .EQ.0) THEN
C User defined values. 

      ELSEIF(MOD2 .EQ.1) THEN
C Actual AGS values

c/========================================================================
c/ Cold Snake a16-19_blw bump.  Backlegs on A16, A17, -A18, -A19 in series.
c/========================================================================

C        WN(2) = 0.D0

        IF    (MMNM .EQ. 'MM_A16AD') THEN
          NBLW = 1
          SIGN = 1.D0
          NWL = 10
          WN(1) = DBLE(NWL) * SIGN
        ELSEIF(MMNM .EQ. 'MM_A17CF') THEN
          NBLW = 1
          SIGN = 1.D0
          NWL = 10
          WN(1) = DBLE(NWL) * SIGN
        ELSEIF(MMNM .EQ. 'MM_A18CF') THEN
          NBLW = 1
          SIGN = -1.D0
          NWL = 10
          WN(1) = DBLE(NWL) * SIGN
        ELSEIF(MMNM .EQ. 'MM_A19BD') THEN
          NBLW = 1
          SIGN = -1.D0
          NWL = 12
          WN(1) = DBLE(NWL) * SIGN


C/================================================================
C/ COLD SNAKE A20 BLW. BACKLEG WIRED ONLY on A20, 12 turns.
c/	Note: There are three power supplies wired on A20, but the controls
c/	      system (currently) combines them into one total current which
c/	      we use here.
c/================================================================
        ELSEIF(MMNM .EQ. 'MM_A20BD') THEN
          NBLW = 1
          SIGN = 1.D0
          NWL = 12
          WN(1) = DBLE(NWL) * SIGN

c/====================================================================
c/ Cold Snake b2-5_blw bump. Backlegs wired on  B2, B3, -B4, -B5 in series
c/====================================================================

c B02 magnet has 2 windings
        ELSEIF(MMNM .EQ. ' MM_B02BF') THEN
          NBLW = 2
          SIGN = 1.D0
          NWL = 12
          WN(1) = DBLE(NWL) * SIGN
          SIGN = 1.D0
          NWL = 6
          WN(2) = DBLE(NWL) * SIGN
        ELSEIF(MMNM .EQ. 'MM_B03CD') THEN
          NBLW = 1
          SIGN = 1.D0
          NWL = 10
        ELSEIF(MMNM .EQ. 'MM_B04CD ') THEN
          NBLW = 1
          SIGN = -1.D0
          NWL = 10
        ELSEIF(MMNM .EQ. 'MM_B05AF') THEN
          NBLW = 1
          SIGN = -1.D0
          NWL = 10

c/====================================================================
c/  L20 Position Bump
c/   BLWs on K19, K20, L13, L14, A7, A8, B1, B2
c/   with      6,   6,  -5,  -5, -5, -5,  6,  6 turns respectively.
c/====================================================================
        ELSEIF(MMNM .EQ. 'MM_K19BD ') THEN
          NBLW = 1
          SIGN = 1.D0
          NWL = 6
          WN(1) = DBLE(NWL) * SIGN
        ELSEIF(MMNM .EQ. 'MM_K20BD') THEN
          NBLW = 1
          SIGN = 1.D0
          NWL = 6
          WN(1) = DBLE(NWL) * SIGN
        ELSEIF(MMNM .EQ. 'MM_L13CF ') THEN
          NBLW = 1
          SIGN = -1.D0
          NWL = 5
          WN(1) = DBLE(NWL) * SIGN
        ELSEIF(MMNM .EQ. 'MM_L14CF') THEN
          NBLW = 1
          SIGN = -1.D0
          NWL = 5
          WN(1) = DBLE(NWL) * SIGN
        ELSEIF(MMNM .EQ. 'MM_A07CD ') THEN
          NBLW = 1
          SIGN = -1.D0
          NWL = 5
          WN(1) = DBLE(NWL) * SIGN
        ELSEIF(MMNM .EQ. 'MM_A08CD') THEN
          NBLW = 1
          SIGN = -1.D0
          NWL = 5
          WN(1) = DBLE(NWL) * SIGN
        ELSEIF(MMNM .EQ. 'MM_B01BF') THEN
          NBLW = 1
          SIGN = 1.D0
          NWL = 6
          WN(1) = DBLE(NWL) * SIGN
C          BLIL20=BLI

c/===============================================================
c/ L20 Angle Bump (1 lambda)
c/  BLW windings around L6, L7, A14, A15
c/               with    5,  5,  -5,  -5 turns, respectively.
c/===============================================================
        ELSEIF(MMNM .EQ. 'MM_L06AF') THEN
          NBLW = 1
          SIGN = 1.D0
          NWL = 5
          WN(1) = DBLE(NWL) * SIGN
        ELSEIF(MMNM .EQ. 'MM_L07CD') THEN
          NBLW = 1
          SIGN = 1.D0
          NWL = 5
          WN(1) = DBLE(NWL) * SIGN
        ELSEIF(MMNM .EQ. 'MM_A14CF') THEN
          NBLW = 1
          SIGN = -1.D0
          NWL = 5
          WN(1) = DBLE(NWL) * SIGN
        ELSEIF(MMNM .EQ. 'MM_A15AD') THEN
          NBLW = 1
          SIGN = -1.D0
          NWL = 5
          WN(1) = DBLE(NWL) * SIGN

c/================================================================
c/ F10 BLW Bump
c/  Windings around E6, E7, E20, F1, F14, F15, G8, G9
c/           with   -5, -5,  6,   6,   5,   5, -6, -6 turns, respectively.
c/================================================================
        ELSEIF(MMNM .EQ. 'MM_E06AF') THEN
          NBLW = 1
          SIGN = -1.D0
          NWL = 5
          WN(1) = DBLE(NWL) * SIGN
        ELSEIF(MMNM .EQ. 'MM_E07CD') THEN
          NBLW = 1
          SIGN = -1.D0
          NWL = 5
          WN(1) = DBLE(NWL) * SIGN
        ELSEIF(MMNM .EQ. 'MM_E20BD') THEN
          NBLW = 1
          SIGN = 1.D0
          NWL = 6
          WN(1) = DBLE(NWL) * SIGN
        ELSEIF(MMNM .EQ. 'MM_F01BF') THEN
          NBLW = 1
          SIGN = 1.D0
          NWL = 6
          WN(1) = DBLE(NWL) * SIGN
        ELSEIF(MMNM .EQ. 'MM_F14CF') THEN
          NBLW = 1
          SIGN = 1.D0
          NWL = 5
          WN(1) = DBLE(NWL) * SIGN
        ELSEIF(MMNM .EQ. 'MM_F15AD') THEN
          NBLW = 1
          SIGN = 1.D0
          NWL = 5
          WN(1) = DBLE(NWL) * SIGN
        ELSEIF(MMNM .EQ. 'MM_G08CD') THEN
          NBLW = 1
          SIGN = -1.D0
          NWL = 5
          WN(1) = DBLE(NWL) * SIGN
        ELSEIF(MMNM .EQ. 'MM_G09BF') THEN
          NBLW = 1
          SIGN = -1.D0
          NWL = 6
          WN(1) = DBLE(NWL) * SIGN

c/============================================================
c/ J10 Dump Bump
c/  BLWs on I10, I11, J04, J05, J18, J19, K12, K13
c/  with      6,   6,   5,   5,   5,   6,   6,   5
c/============================================================
c$$$      write(lunout,*) 'MULTIPOL MM_I10BF'       
c$$$      write(lunout,*) factor(scalfact,-1.d0,06.d0,blI,AMM)
c$$$      write(lunout,*) 'MULTIPOL MM_I11BD'
c$$$      write(lunout,*) factor(scalfact,-1.d0,06.d0,blI,AMM)
c$$$      write(lunout,*) 'MULTIPOL MM_J04CD'
c$$$      write(lunout,*) factor(scalfact, 1.d0,05.d0,blI,AMM)
c$$$      write(lunout,*) 'MULTIPOL MM_J05AF'
c$$$      write(lunout,*) factor(scalfact, 1.d0,05.d0,blI,AMM)
c$$$      write(lunout,*) 'MULTIPOL MM_J18CF'
c$$$      write(lunout,*) factor(scalfact, 1.d0,05.d0,blI,AMM)
c$$$      write(lunout,*) 'MULTIPOL MM_J19BD'
c$$$      write(lunout,*) factor(scalfact, 1.d0,06.d0,blI,AMM)
c$$$      write(lunout,*) 'MULTIPOL MM_K12BD'
c$$$      write(lunout,*) factor(scalfact,-1.d0,06.d0,blI,AMM)
c$$$      write(lunout,*) 'MULTIPOL MM_K13CF'
c$$$      write(lunout,*) factor(scalfact,-1.d0,05.d0,blI,AMM)

c/=================================================
c/ "G9" Backleg winding xtractn bump for G10 kicker
c/   Windings  F8  F9  -G2  -G3  -G16  -G17  
c/=================================================
        ELSEIF(MMNM .EQ. 'MM_F08CD') THEN
c$$$      write(lunout,*) factor(scalfact, 1.d0,05.d0,blI,AMM)
          NBLW = 1
          SIGN = 1.D0
          NWL = 5
          WN(1) = DBLE(NWL) * SIGN
        ELSEIF(MMNM .EQ. 'MM_F09BF') THEN
c$$$      write(lunout,*) factor(scalfact, 1.d0,06.d0,blI,AMM)
          NBLW = 1
          SIGN = 1.D0
          NWL = 19 ! ?????????????????????????
          WN(1) = DBLE(NWL) * SIGN
        ELSEIF(MMNM .EQ. 'MM_G02BF') THEN
c$$$      write(lunout,*) factor(scalfact,-1.d0,06.d0,blI,AMM)
          NBLW = 1
          SIGN = -1.D0
          NWL = 6
          WN(1) = DBLE(NWL) * SIGN
        ELSEIF(MMNM .EQ. 'MM_G03CD') THEN
c$$$      write(lunout,*) factor(scalfact,-1.d0,05.d0,blI,AMM)
          NBLW = 1
          SIGN = -1.D0
          NWL = 5
          WN(1) = DBLE(NWL) * SIGN
        ELSEIF(MMNM .EQ. 'MM_G16AD') THEN
c$$$      write(lunout,*) factor(scalfact,-1.d0,10.d0,blI,AMM)
          NBLW = 1
          SIGN = -1.D0
          NWL = 10
          WN(1) = DBLE(NWL) * SIGN
        ELSEIF(MMNM .EQ. 'MM_G17CF') THEN
c$$$      write(lunout,*) factor(scalfact,-1.d0,10.d0,blI,AMM)
          NBLW = 1
          SIGN = -1.D0
          NWL = 10
          WN(1) = DBLE(NWL) * SIGN

c/==================================================
c/ "H11" Backleg winding xtractn bump for H10 septum
c/   Windings  -H4 -H5 -H18 -H19 I12 I13
c/==================================================
        ELSEIF(MMNM .EQ. 'MM_H04CD') THEN
          NBLW = 1
          SIGN = -1.D0
          NWL = 10
          WN(1) = DBLE(NWL) * SIGN
        ELSEIF(MMNM .EQ. 'MM_H05AF') THEN
          NBLW = 1
          SIGN = -1.D0
          NWL = 10
          WN(1) = DBLE(NWL) * SIGN
        ELSEIF(MMNM .EQ. 'MM_H18CF') THEN
          NBLW = 1
          SIGN = -1.D0
          NWL = 5
          WN(1) = DBLE(NWL) * SIGN
        ELSEIF(MMNM .EQ. 'MM_H19BD') THEN
          NBLW = 1
          SIGN = -1.D0
          NWL = 6
          WN(1) = DBLE(NWL) * SIGN
        ELSEIF(MMNM .EQ. 'MM_I12BD') THEN
          NBLW = 1
          SIGN = 1.D0
          NWL = 6
          WN(1) = DBLE(NWL) * SIGN
        ELSEIF(MMNM .EQ. 'MM_I13CF') THEN
          NBLW = 1
          SIGN = 1.D0
          NWL = 5
          WN(1) = DBLE(NWL) * SIGN

        ELSE

          NBLW = 0
          SIGN = 0.d0
          NWL = 0
          WN(1) = 0.D0
          WN(2) = 0.D0

        ENDIF

      ELSE

        CALL ENDJOB('SBR AGSBLW. No such possibility MOD2',-99)
      ENDIF

      BLWI = 0.D0
      DO I = 1, NBLW
        BLWI = BLWI + WA(I)*WN(I)
      ENDDO

c         IF(MMNM .EQ. 'MM_F08CD') write(*,*) ' MM_F08CD'
c         IF(MMNM .EQ. 'MM_F08CD') write(*,*) ' blwi  ',blwi

c      if(blwi .gt. 1e-6) then
      DO I = 1, NB
        BM(I) = BM(I) * ( 1.D0 + BLWI/ AGSMMA(DEV) )
c        IF(I.EQ.1) THEN
c         IF(MMNM .EQ. 'MM_F08CD') write(*,*) ' MM_F08CD'
c        IF(MMNM.EQ.'MM_F08CD')write(*,*)bm(i),blwi,AGSMMA(DEV),dev,wa(1) 
c         IF(MMNM .EQ. 'MM_F08CD') write(*,*) ' '
cc         IF(MMNM .EQ. 'MM_F08CD') read(*,*)
c         IF(MMNM .EQ. 'MM_F09BF') write(*,*) ' MM_F09BF'
c         IF(MMNM .EQ. 'MM_F09BF') write(*,*) bm(i),blwi,AGSMMA(DEV),dev 
cc         IF(MMNM .EQ. 'MM_F09BF') read(*,*)
c        ENDIF
      ENDDO
c      endif
               
      RETURN
      END
