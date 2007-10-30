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
      SUBROUTINE INIGR(
     >                 NLOG, LM, NOMFIC)
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

      DATA  OKECH ,OKVAR, OKBIN / .FALSE., .TRUE., .FALSE.  /

      DATA KVAR/
     >' dp/p ','   Y  ','  Y'' ','   Z  ','  Z'' ','   s  ',' Time ',
     >'   X  ',' Step ','   r  ',
     >'dp/p|o',' Y_o  ',' Y''_o ',' Z_o  ',' Z''_o  ',' s_o  ',' Time ',
     >'Phase ',' dp/p ','KinEnr',
     >'  Sx  ','  Sy  ','  Sz  ',' <S>  ',
     >' <Sx> ',' <Sy> ',' <Sz> ','COUNT ','      ',
     >'  Bx  ','  By  ','  Bz  ','  Br  ',
     >'  Ex  ','  Ey  ','  Ez  ','  Er  ',
     >' S_out',' Pass#'  ,2*'      ',
     >' Y_Lab','      ',' Z_Lab','      ','      ','      ',' X_Lab',
     >8*' ',
     >' lmnt#' ,
     >13*' '
     >/
      DATA KPOL/ 'cartesian' , 'cylindr.' /

C      DATA BORNE/ .01D0, .99D0, .01D0, .99D0, .001D0, .999D0 /
      DATA BORNE/ .0D0, .5D0, .0D0, .5D0, .001D0, .999D0 /
      DATA NC0/ 2000, 2000, 2000 /

      DATA NPTS / NTRMAX /

      DATA KDIM/
     >'       ','  (m)  ',' (rad) ','  (m)  ',' (rad) ','  (m)  ',
     >'(mu_s) ','  (m)  ','  (m)  ','  (m)  '                  ,
     >'       ','  (m)  ',' (rad) ','  (m)  ',' (rad) ','  (m)  ',
     >'       ',' (rad) ','       ',' (MeV) ',9*'       ',
     > 4*'  (T)  ', 4*'(eV/m) ' ,
     >'  (m)  ',3*' ',
     >'  (m)  ','      ','  (m)  ','      ','      ','      ','  (m)  ',
     >22*'   '/


      DATA KARSIZ / 3 /
      SAVE KARSIZ 

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
      CALL TRTXT(2.D0,251.D0,TXT,32,0)
      CALL DATE2(DMY)
      WRITE(TXT,FMT='(A)') DMY
      CALL TRTXT(2.D0,245.D0,TXT,9,0)
      CALL DEFCAR(KARSIZ,0,0)  
      RETURN

      ENTRY LOGO2(TEMP)
      LOGOT = TEMP(51:69)
      LOGOT = 'Zgoubi'
            
      RETURN
      END
