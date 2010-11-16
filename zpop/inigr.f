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
C  Brookhaven National Laboratory                                                               és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE INIGR(
     >                 LM, NOMFIC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*(*) NOMFIC

      LOGICAL OKECH, OKVAR, OKBIN
      COMMON/ECHL/OKECH, OKVAR, OKBIN

      INCLUDE 'MXVAR.H'
      CHARACTER KVAR(MXVAR)*7, KPOL(2)*9, KDIM(MXVAR)*7
      COMMON/INPVR/ KVAR, KPOL, KDIM

      PARAMETER (NCANAL=2500)
      COMMON/SPEDF/BORNE(6),SPEC(NCANAL,3),PMAX(3),NC0(3)

      INCLUDE 'MAXNTR.H'
      COMMON/TRACKM/COOR(NTRMAX,9),NPTS,NPTR

      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB

      CHARACTER * 9   DMY
      CHARACTER*80 TXT
      CHARACTER LOGOT*18, TEMP*80

      SAVE KARSIZ, NLOG

      DATA KARSIZ / 3 /

C----- Histograms 
      IF(KX .EQ. 28) THEN
C        Counts are always on Y axis
        KX=KY
        KY=28
      ENDIF
      IF(KX*KY .LT. 0) THEN
C        Will plot mean value of Y vs X
        MBIN = 1
        IF(KX .LT. 0) THEN
          KYB = -KX 
          KX = KY
          KY = 28 
        ELSEIF(KY .LT. 0) THEN
          KYB = -KY 
          KY = 28 
        ENDIF
      ELSE
        KYB = 0    
        MBIN = 0
      ENDIF
C      Defines vertical axis of histograms. Default (MBIN=0) is COUNTS, 
C      otherwise (MBIN=1) can be mean value of any KY type coordinate
      CALL BIN3W(MBIN,KYB)

C----- Number of the lmnt concerned by the plot (-1 for all)
      LM = -1

C----- zpop log unit
      NLOG = 30
C----- Input data file name
C      Normally, default is zgoubi.fai, .plt, .spn, .map...
      NOMFIC = 'none'

C----- Line type for plot, e.g., 9=dots, 1=solid line...
C      CALL LINTYW(1)
      CALL LINTYW(9)

      ENTRY DFKSIZ
C----- Graphic character size
      CALL DEFCAR(KARSIZ,0,0)  

      RETURN

      ENTRY LOGO
      CALL DEFCAR(1,0,0)
      WRITE(TXT,FMT='(A12)') 'Zgoubi|Zpop '
      CALL TRTXT(38.D0,251.D0,TXT,0)
      CALL DATE2(DMY)
      WRITE(TXT,FMT='(A)') DMY
      CALL TRTXT(38.D0,244.D0,TXT,0)
      CALL DEFCAR(KARSIZ,0,0)  
      RETURN

      ENTRY LOGO2(TEMP)
      LOGOT = TEMP(51:69)
      LOGOT = 'Zgoubi'            
      RETURN

      ENTRY INIGR1(
     >             NLOGO)
      NLOGO = NLOG
      RETURN

      END
