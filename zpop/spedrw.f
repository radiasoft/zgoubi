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
C  Brookhaven National Laboratory                                               és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE SPEDRW(NT,BORNE,YNU,PMAX,SPEC,NC0,OKECH)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL OKECH
      PARAMETER (NCANAL=2500)
      DIMENSION BORNE(6),YNU(3),PMAX(3),SPEC(NCANAL,3),NC0(3)

      COMMON/VXPLT/ DU1,DU2,DU3,DU4,KX,KY,IAX,LIS,NB
 
C      CALL TXTFBG
      CHARACTER REP
      INTEGER FINSTR

      DATA ISCA / 1 /

      KX0=KX
      KY0=KY
      KX=99
      KY=99

      GOTO 21

 20   CONTINUE      
      CALL FBGTXT
      WRITE(*,FMT='(/,''  Press RETURN for more'')') 
      READ(*,200,ERR=20) REP 
 200  FORMAT(A1)

 21   CONTINUE
      CALL FBGTXT
      CALL HOMCLR

      WRITE(*,104) 
 104  FORMAT(5X,' OPTIONS : ')
      WRITE(*,109) 
 109  FORMAT(/,5X,'  1,  2,  3  ** Plot spectrum  (Y, Z, D) **')
      WRITE(*,110)
 110  FORMAT(  
     4         5X,'          4  Normalized amplitude'
     8      ,/,5X,'          8  Print screen'
     9      ,/,5X,'          9  EXIT  THIS  MENU '
     2      ,/,5X,'         12  ERASE  DISPLAY    '
     >,/)

      WRITE(*,100)
 100  FORMAT('$  Option  number : ')
      READ(*,101,ERR=21) IOPT
101   FORMAT(I2)
      GOTO ( 1, 1, 1, 4,21,21,21, 8, 9,21,21,12) IOPT  
      GOTO 21

 1    CONTINUE
      JNU = IOPT
      IF    (ISCA.EQ. 1) THEN
C------- Vertical scale normalized to 1
        FAC = 1.D0 / PMAX(JNU)
        YMAX = 1.05D0       
      ELSEIF(ISCA.EQ. 2) THEN
        FAC = 1.D0
        YMAX = PMAX(JNU) * 1.05D0
      ENDIF

      ANUI = BORNE(2*JNU-1)
      ANUF = BORNE(2*JNU)

      IF( .NOT. OKECH ) THEN
        CALL TRAXES(ANUI,ANUF,0.D0,YMAX,1)
        OKECH = .TRUE.
      ENDIF

      DELNU=(ANUF - ANUI)/NC0(JNU)
      Z = ANUI - 0.5D0 * DELNU
 
      CALL LINTYP(1)
      DO 15 K=1,NC0(JNU)
        Z=Z+DELNU
        CALL VECTPL(Z,0.D0,4)
        CALL VECTPL(Z,SPEC(K,JNU) * FAC,2)
        call fbgtxt 
        write(77,*) ' spedrw, Z,SPEC(K,JNU) : ',Z,SPEC(K,JNU),JNU,K
     >      ,'  nu,SPEC(K,JNU),JNU,iocc. '
 15   CONTINUE

      CALL SPEGR(NT,JNU,YNU,PMAX)
      CALL FBGTXT
      GOTO 20

 4    CONTINUE
        WRITE(*,*) ' Vertical  amplitude  normalized  to  1 (y/n) :'
        READ(*,FMT='(A)',ERR=4) REP
        IF(REP .EQ. 'Y') REP = 'y'
        IF(FINSTR(REP) .EQ. 0 .OR. REP .EQ. 'y') THEN
          ISCA = 1
        ELSE
          ISCA = 2
        ENDIF
      GOTO 21

 8    CONTINUE
C        CALL MENVCF
        CALL SAVPLT
      GOTO 21

 9    CONTINUE
      CALL LINTYP(-1)
      KX=KX0
      KY=KY0
      RETURN

 12   CONTINUE
      CALL CLSCR
      OKECH=.FALSE.
      GOTO 21

      END
