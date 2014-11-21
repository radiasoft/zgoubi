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
      SUBROUTINE SPEANA(YM,BORNE,NC0,
     >                               YNU,SPEC,PMAX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NCANAL=2500)
      DIMENSION YM(3), BORNE(6), NC0(3), YNU(3), SPEC(NCANAL,3), PMAX(3)

      INCLUDE 'MAXNTR.H'
      COMMON/TRACKM/COOR(NTRMAX,9),NPTS,NPTR
 
      PARAMETER ( PI=3.1415926536 , DEUXPI=2.0*PI )

      DO 1 INU = 1, 5, 2
        JNU = 1 + INU/2
        ANUI = BORNE(INU)
        ANUF = BORNE(INU+1)
        DELNU=(ANUF - ANUI) / NC0(JNU)
        PAS=DEUXPI * DELNU
        VAL=DEUXPI *(ANUI - 0.5d0 * DELNU)
        PMAX(JNU)=0.D0
        PMIN=1.D12
        DO 20 NC=1,NC0(JNU)
          VAL=VAL+PAS
          SR=0.D0
          SI=0.D0
          DO 10 NT=1,NPTS
C ** ERR, FM Jan/04            FF = COOR(NT,INU) - YM(INU)
              FF = COOR(NT,INU)  - YM(JNU)
              SR=SR + FF * COS(NT*VAL)
              SI=SI + FF * SIN(NT*VAL)
 10       CONTINUE
          PP=SR*SR+SI*SI
          IF(PP.GT. PMAX(JNU)) THEN
            PMAX(JNU)=PP
            KMAX=NC
          ELSEIF(PP.LT. PMIN) THEN
            PMIN=PP
          ENDIF
          SPEC(NC,JNU)=PP
 20     CONTINUE

c           write(88,*) ' speana spec(nc,jnu) -  kmax : ',kmax

        IF (PMAX(JNU) .GT. PMIN) THEN
           IF (KMAX .LT. NCANAL) THEN
             DEC=0.5D0 * (SPEC(KMAX-1,JNU)-SPEC(KMAX+1,JNU))
     >       /(SPEC(KMAX-1,JNU) - 2.D0 *SPEC(KMAX,JNU)+SPEC(KMAX+1,JNU))
           ELSE
             DEC=0.5D0 
           ENDIF
           YNU(JNU)= ANUI + (DBLE(KMAX)+DEC-0.5D0) * DELNU
        ELSE
           YNU(JNU) = 0.D0
        ENDIF
 1    CONTINUE
      RETURN
      END
