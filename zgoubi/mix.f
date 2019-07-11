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
      SUBROUTINE MIX(IMAX,FO,FAISC1,FAISC2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "MAXCOO.H"
      DIMENSION FO(MXJ,*),FAISC1(MXJ,*),FAISC2(MXJ,*)
      PARAMETER (MXJ1=MXJ-1)

C     ***  CONSTITUTION DU FAISCEAU A PARTIR DE FAISC1 ET
C       FAISC2 PAR MELANGE DES 2 TIRAGES ( ON BAT LES CARTES )
C       ON TIRE DANS L INTERVALLE (0,1) :
C     DE (0 ,.5)  FAISCEAU PRIS DANS 1
C     DE (.5,1.)  FAISCEAU PRIS DANS 2
C      GENERATEUR IR3
C       TRAJECTOIRES

C      JRND11=2
C      JRND21=2
C      JRND31=2
C      JRND41=2
C      JRND51=2
C      JRND12=2
C      JRND22=2
C      JRND32=2
C      JRND42=2
C      JRND52=2
C      JRND1=2
C      JRND2=2

      JRND11=1
      JRND21=1
      JRND31=1
      JRND41=1
      JRND51=1
      JRND12=1
      JRND22=1
      JRND32=1
      JRND42=1
      JRND52=1
      JRND1=1
      JRND2=1

      DO 3 I=1,IMAX
        XARPHA=RNDM()
        IF(XARPHA .LE. .5D0) THEN
          FO(2,I)=FAISC1(2,JRND21)
          JRND21=JRND21+1
        ELSE
          FO(2,I)=FAISC2(2,JRND22)
          JRND22=JRND22+1
        ENDIF
        XARPHA=RNDM()
        IF(XARPHA.LE.5D0) THEN
          FO(3,I)=FAISC1(3,JRND31)
          JRND31=JRND31+1
        ELSE
          FO(3,I)=FAISC2(3,JRND32)
          JRND32=JRND32+1
        ENDIF
        XARPHA=RNDM()
        IF(XARPHA.LE.5D0) THEN
          FO(4,I)=FAISC1(4,JRND41)
          JRND41=JRND41+1
        ELSE
          FO(4,I)=FAISC2(4,JRND42)
          JRND42=JRND42+1
        ENDIF
        XARPHA=RNDM()
        IF(XARPHA.LE.5D0) THEN
          FO(5,I)=FAISC1(5,JRND51)
          JRND51=JRND51+1
        ELSE
          FO(5,I)=FAISC2(5,JRND52)
          JRND52=JRND52+1
        ENDIF
        XARPHA=RNDM()
        IF(XARPHA.LE.5D0) THEN
          FO(1,I)=FAISC1(1,JRND11)
          JRND11=JRND11+1
        ELSE
          FO(1,I)=FAISC2(1,JRND12)
          JRND12=JRND12+1
        ENDIF
        XARPHA=RNDM()
        IF(XARPHA.LE.5D0) THEN
          JRND1=JRND1+1
        ELSE
          JRND2=JRND2+1
        ENDIF
 3    CONTINUE
      RETURN
      END
