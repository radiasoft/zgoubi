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
C  -------
      SUBROUTINE AGSBLW(MOD2,NOEL,DEV,NBLW,WN,WA,
     >                                           BM,NB)
C     ----------------------------------------------
C     NBLW : # of blwg ; WN(i=1,nblw) : # of turns ; 
C     WA(i=1,nblw) : amperes.   
C     ----------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION WN(*),WA(*),BM(NB)
      INCLUDE 'MXFS.H'
      INCLUDE 'MXLD.H'
      PARAMETER (LBLSIZ=10)
      CHARACTER(LBLSIZ) LABEL
      COMMON /LABEL/ LABEL(MXL,2)
      PARAMETER (KSIZ=10)
      CHARACTER(KSIZ) FAM
      CHARACTER(LBLSIZ) LBF
      COMMON/SCALT/ FAM(MXF),LBF(MXF,MLF),JPA(MXF,MXP)

      CHARACTER MMNM*(LBLSIZ)
      parameter(turnMM=16.d0)
      integer debstr, finstr

      MMNM = LABEL(NOEL,1)
      MMNM = MMNM(DEBSTR(MMNM):FINSTR(MMNM))

      IF(MOD2 .EQ.1) THEN
C Default values are taken

c/========================================================================
c/ Cold Snake a16-19_blw bump.  Backlegs on A16, A17, -A18, -A19 in series.
c/========================================================================
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
          BLIL20=BLI

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

c/=============================
c/ G09 Backleg winding bump
c/   Windings around F8  F9  -G2  -G3  -G16  -G17  
c/==================================
c$$$      write(lunout,*) 'MULTIPOL MM_F08CD'
c$$$      write(lunout,*) factor(scalfact, 1.d0,05.d0,blI,AMM)
c$$$      write(lunout,*) 'MULTIPOL MM_F09BF'
c$$$      write(lunout,*) factor(scalfact, 1.d0,06.d0,blI,AMM)
c$$$      write(lunout,*) 'MULTIPOL MM_G02BF'
c$$$      write(lunout,*) factor(scalfact,-1.d0,06.d0,blI,AMM)
c$$$      write(lunout,*) 'MULTIPOL MM_G03CD'
c$$$      write(lunout,*) factor(scalfact,-1.d0,05.d0,blI,AMM)
c$$$      write(lunout,*) 'MULTIPOL MM_G16AD'
c$$$      write(lunout,*) factor(scalfact,-1.d0,10.d0,blI,AMM)
c$$$      write(lunout,*) 'MULTIPOL MM_G17CF'
c$$$      write(lunout,*) factor(scalfact,-1.d0,10.d0,blI,AMM)

        ENDIF
      ENDIF

      BLWI = 0.D0
      DO I = 1, NBLW
        BLWI = BLWI + WA(I)*WN(I)
      ENDDO
      DO I = 1, NB
        BM(I) = BM(I) * ( 1.D0 + BLWI/ AGSMMA(DEV) )
      ENDDO

      RETURN
      END
