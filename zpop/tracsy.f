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
C  USA
C  -------
      SUBROUTINE TRACSY(NN,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C----------------------------------------------------------------------------
C
C     TRACE DU SYNOPTIQUE
C     
C----------------------------------------------------------------------------
      LOGICAL OKECH, OKVAR, OKBIN
      COMMON/ECHL/OKECH, OKVAR, OKBIN
      COMMON/VEN/Y1,S1,D,DQ,DS,DM,DSS,DPU,TETA1,
     >           SMIN,SMAX,YMIN,YMAX,SB0,YB0
      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB

      CALL INSY(IPLAN,SB0,YB0,TETA0,D,DQ,DS,DM,DSS,DPU,*99)
      CALL CALTRA(0,0)
C      CALL DEFCAR(1,0,0) 
         YMIN=0.D0
         YMAX=0.D0
         SMIN=0.D0
         SMAX=0.D0

      N=1
100   CONTINUE
        IF(N.EQ.2) THEN 
          WRITE(6,*)'------ LATTICE SYNOPTIC ---------'
          WRITE(6,*)'             plotting... '
C          CALL TRAXPI                     ! Set proportional axis

          IF(.NOT.OKECH) THEN
C            CORR has been introduced so to make sure plot window does encompass xy-min/max
            XMI=SMIN
            XMA=SMAX
            CORR = (XMA-XMI)* 2.D-2
            XMI=XMI - CORR
            XMA=XMA + CORR
            YMI=YMIN
            YMA=YMAX
            CORR = (YMA-YMI)* 1.D-1
            IF(ABS(CORR) .LT. 1E-8) CORR = .1D0
            YMI=YMI - CORR
            YMA=YMA + CORR
            CALL TRAXES(XMI,XMA,YMI,YMA,1)
C            CALL TRAXES(XMI,XMA,YMI,YMA,3)
            OKECH = .TRUE.
          ENDIF
          CALL VECTPL(SB0,YB0,4)
        ENDIF

        TETA1=TETA0
        CALL CALTRI(N,0.D0,0.D0)
        CALL CALTRA(N,IPLAN)

        IF(N.EQ.2) THEN
C           CALL TRAXPJ                                
        ELSE
           N=N+1
           IF(N.LE.NN) GOTO 100
        ENDIF

      RETURN

 99   CONTINUE
      WRITE(6,*) 
      WRITE(6,*) ' Sorry, synoptic will not be plotted,'
      WRITE(6,*) '           Lab frame is undefined, '
      WRITE(6,*) '           Lab coordinates probably meaningless...'
C      CALL DEFCAR(2,0,0) 
      RETURN 1

      END
