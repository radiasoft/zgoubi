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
      SUBROUTINE DRVG(IG,C,S,G,DG,D2G,D3G,D4G,D5G,D6G)
C     >                                   D7G,D8G,D9G,D10G)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION C(*)
      COMMON/TRAJ/ Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT

      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG

      D2G = 0.D0
      D3G = 0.D0
      D4G = 0.D0
      D5G = 0.D0
      D6G = 0.D0
      D7G = 0.D0
      D8G = 0.D0
      D9G = 0.D0
      D10G = 0.D0

      POL = C(1)+(C(2)+(C(3)+(C(4)+(C(5)+C(6)*S)*S)*S)*S)*S
C      WRITE(NRES,*) ' SBR DRVG C, S, POL :',(C(I),I=1,6),S,POL
      EP = EXP(POL)
      G  = 1.D0/(1.D0+EP)
C----- U=1+EP, GU=1

      IF(IG .LE. 0) RETURN

      DP  = C(2)+(2.D0*C(3)+(3.D0*C(4) +(4.D0*C(5) +5.D0*C(6)*S)*S)*S)*S
      DU  = DP*EP
C----- (GU)' = 0
      DG =-DU*G*G
 
      IF(IG .LE. 1) RETURN

      D2P = 2.D0*C(3) +(6.D0*C(4) +(12.D0*C(5) +20.D0*C(6)*S)*S)*S
      D2U = D2P*EP + DP*DU
C----- (GU)'' = 0
      D2G=-G*(2.D0*DG*DU+G*D2U)
 
      IF(IG .LE. 2) RETURN
 
      D3P = 6.D0*C(4) +(24.D0*C(5) + 60.D0*C(6)*S)*S
      D3U = D3P*EP + 2.D0*D2P*DU + DP*D2U
C----- (GU)''' = 0
      D3G=-G*(3.D0*(D2G*DU+DG*D2U)+G*D3U)
 
      IF(IG .LE. 3) RETURN
 
      D4P = 24.D0*C(5) +120.D0*C(6)*S
      D4U = D4P*EP + 3.D0*D3P*DU + 3.D0*D2P*D2U + DP*D3U
C----- (GU)'''' = 0
      D4G=-G*(4.D0*(D3G*DU+DG*D3U)+6.D0*D2G*D2U+G*D4U)
 
      IF(IG .LE. 4) RETURN
 
      D5P = 120.D0*C(6)
      D5U = D5P*EP + 4.D0*D4P*DU + 6.D0*D3P*D2U + 4.D0*D2P*D3U + DP*D4U
C----- (GU)''''' = 0
      D5G=-G*(5.D0*(D4G*DU+DG*D4U)+10.D0*(D3G*D2U+D2G*D3U)+G*D5U)
 
      IF(IG .LE. 5) RETURN
 
      D6U = 5.D0*D5P*DU+ 10.D0*D4P*D2U+
     >  10.D0*D3P*D3U+ 5.D0*D2P*D4U+ DP*D5U
C----- (GU)""" = 0
      D6G=-G*( 6.D0*(D5G*DU+ DG*D5U)+ 15.D0*( D4G*D2U+ D2G*D4U)+
     >  20.D0*D3G*D3U+ G*D6U)
 
      IF(IG .LE. 6) RETURN
 
      D7U =15.D0*D5P*D2U+20.D0*D4P*D3U+15.D0*D3P*D4U+6.D0*D2P*D5U+DP*D6U
C----- (GU)"""' = 0
      D7G=-G*( 7.D0*(D6G*DU+ DG*D6U)+ 21.D0*( D5G*D2U+ D2G*D5U)+
     >  35.D0*(D4G*D3U + D3G*D4U)+ G*D7U)
 
      IF(IG .LE. 7) RETURN
 
      D8U =70.D0*D5P*D3U+56.D0*D4P*D4U+28.D0*D3P*D5U+8.D0*D2P*D6U+DP*D7U
C----- (GU)"""" = 0
      D8G=-G*( 8.D0*(D7G*DU+ DG*D7U)+ 28.D0*( D6G*D2U+ D2G*D6U)+
     >  56.D0*(D5G*D3U + D3G*D5U)+ 70.D0*D4G*D4U + G*D8U)
 
      IF(IG .LE. 8) RETURN
 
      D9U =126.D0*D5P*D4U+84.D0*D4P*D5U+36.D0*D3P*D6U
     >                                              +9.D0*D2P*D7U+DP*D8U
C----- (GU)""""' = 0
      D9G=-G*( 9.D0*(D8G*DU+ DG*D8U)+ 36.D0*( D7G*D2U+ D2G*D7U)+
     >  84.D0*(D6G*D3U + D3G*D6U)+ 126.D0*(D5G*D4U + D4G*D5U) + G*D9U)
 
      RETURN
      END
