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
      SUBROUTINE INPVAR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL OKECH, OKVAR, OKBIN
      COMMON/ECHL/OKECH, OKVAR, OKBIN
      INCLUDE 'MXVAR.H'
      CHARACTER KVAR(MXVAR)*7, KPOL(2)*9, KDIM(MXVAR)*7
      COMMON/INPVR/ KVAR, KPOL, KDIM
      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB

      LOGICAL TYLAB

      GOTO 21

 20   CONTINUE      
      WRITE(6,FMT='(/,''  Press RETURN for more'')') 
      READ(5,FMT='(A1)',ERR=20) REP

 21   CONTINUE

      WRITE(6,104) KVAR(KX),KX,KVAR(KY),KY
 104  FORMAT(/,' Present variables,   Horiz.  : ',A,'(KX=',I2,'), '
     > /,      '                      Vertic. : ',A,'(KY=',I2,')',/)

      WRITE(6,116)
 116  FORMAT(
     >/,'  * Available  variables -  KX (Horiz.), KY (Vert.) :'
     >,//,'   CURRENT COORDINATES :'
     >, /,8X,'1, 2, 3, 4, 5, 6, 7, 8: dp/p, Y,  T,  Z,  P,  S,  time, X'
     >, /,8X,'                              m  rd   m  rd   m   mu-s  m'
     >, /,8X,'* Note : Y, X are in optical element frame.'
     >, /,8X,'         In polar frames, Y=Radius and X=Angle'
     >, /,8X,' 9 :   Step size (cm)  -  with .plt type files      '                    
     >, /,8X,'10 :   r = sqrt(y^2+z^2)      '                    
     >, /,8X,'20 :   Kinetic E (MeV)  -  with .fai type files '                    
     >,//,'   INITIAL COORDINATES :'
     >, /,8X,'11,12,13,14,15,16,17  : dp/p, Y,  T,  Z,  P,  S,  time'  
     >, /,8X,'                              m  rd   m  rd   m   mu-s'
     >, /,'   SYNCHROTRON MOTION:'
     >, /,8X,'18, 19 :   phi-phi_s (rad),  dp/p  '                    
     >,//,'   SPIN:'
     >, /,8X,'21,22,23,24:    Sx, Sy, Sz, <S>'
     >, /,8X,'25,26,27,:      SumSx/#turn, SumSy/#turn SumSz/#turn'
     >,//,'   HISTOGRAM :'
     >, /,8X,'28: Counts'
     >,//,'   E and B FIELDS :'
     >, /,8X,'30,31,32,33 :    Bx, By, Bz, Br = sqrt(By^2+Bz^2)  (T)'
     >, /,8X,'34,35,36,37 :    Ex, Ey, Ez, Er = sqrt(Ey^2+Ez^2)  (eV/m)'
     >,//,'   COORDINATES IN LABORATORY:'
     >, /,8X,'42,44,48:  Y_Lab,  Z_Lab,  X_Lab,  '
     >, /,8X,'             m       m       m      '
     >, /,'   OTHER :'
     >, /,8X,'38 :  S_lost (m) = Distance where particle is lost'
     >, /,8X,'39 :  IPASS      = Pass (or Turn) #  '    
     >, /,8X,'57 :  NOEL       = Element #  '
     >,/)

      KX0 = KX
      KY0 = KY
      WRITE(6,100) ' * Give desired variables  (KX=KY to quit) : '
 100  FORMAT(A45,$)
      READ(5,*,ERR=20) KX, KY

      IF(KX .NE. KY) THEN

C------- Histogram
        IF(KX.EQ. 28) THEN
C--------- COUNTS
          KX=KY
          KY=28
        ENDIF

C------- Histogram type of plot
        IF(KX*KY .LT. 0) THEN
C--------- Will plot mean value of Y vs X
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
C        Defines vertical axis of histograms. Default (MBIN=0) is COUNTS, 
C        otherwise (MBIN=1) can be mean value of any KY type coordinate
        CALL BIN3W(MBIN,KYB)
       
        IF(KY.EQ. 28) THEN
C          To force re-binning
          IF(KX.NE.KX0) OKBIN=.FALSE.
        ENDIF
        

C------  COORDINATES IN LABORATORY
        IF(TYLAB(KX,KY) ) THEN
          IF(KX.NE.48 .AND. KY.NE.48) THEN
            WRITE(6,FMT='(/,A/)')
     >       '*** Sorry, X or Y has to be variable X_Lab (#48) '
            GOTO 20
          ELSE
            IF(KX.NE.48) THEN
              KY=KX
              KX=48
            ENDIF
          ENDIF
        ELSEIF(KX/10.EQ.4 .OR. KY/10.EQ.4) THEN
          KX = KX0
          KY = KY0            
          WRITE(6,FMT='(/,A/)')
     >     '*** Sorry, X and Y must both be lab variables  '
          GOTO 20
c        ELSEIF(KX.EQ.68 .AND. KY.EQ.62) THEN
c          MOD = 0
c          RM = 3.482590E2 
c          CM2M = 0.01
c          CALL READCC(MOD,RM*CM2M,RM*CM2M)
        ENDIF

      ELSE

        KX = KX0
        KY = KY0

      ENDIF

      RETURN
      END
