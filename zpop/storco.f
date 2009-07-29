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
      SUBROUTINE STORCO(MODSTO,NL,LM,KPS,
     >                                          NPASS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ---------------------------------------------------
C     Read coordinates from zgoubi output file, and store  
C     ---------------------------------------------------
      INCLUDE 'MAXNTR.H'          
      COMMON/TRACKM/COOR(NTRMAX,9),NPTS,NPTR

      CHARACTER LET 

      INCLUDE 'MXVAR.H'
      DIMENSION YZXB(MXVAR),NDX(5)

      CALL RAZ(COOR,NTRMAX*9)

            CALL REWIN2(NL,*96)

      IF(MODSTO.EQ.2) THEN
C----- Case call by ellipse fit
        CALL READC1(KP1,KP2)
        IF(KP1.EQ.-2) THEN
C------- Only last pass wanted, hence first search value of NPASS
          CALL REWIN2(NL,*96)
          WRITE(6,*) '  Reading coordinates, looking for NPASS, wait...'
          NOC=0
          NRBLT = -1 
C--------- BOUCLE SUR READ FICHIER NL 
 45       CONTINUE
            CALL READCO(NL,LM,
     >                        KART,LET,YZXB,NDX,*12,*78)

C--------- NDX: 1->KEX, 2->IT, 3->IREP, 4->IMAX

            NOC=NOC+1

            IF(NINT(YZXB(39)) .GE. NRBLT+1) NRBLT = NINT(YZXB(39)) -1
            IF(NOC.EQ. NPTR) GOTO 13

          GOTO 45
C         ----------------------------------
 78       CONTINUE
          WRITE(6,*) ' *** Coordinates reading stopped : error during',
     >     ' read of event # ',NOC+1
          GOTO 13

 12       CONTINUE
          WRITE(6,*) ' READ  OK; END  OF  FILE  ENCOUNTERED'

 13       CONTINUE
          NPASS = NRBLT + 1
          WRITE(6,*) ' Found NPASS = ',NPASS
        ENDIF
      ENDIF

      NOC=0
      NRBLT = -1 
C----- BOUCLE SUR READ FICHIER NL 
 44   CONTINUE
              
c             write(*,*) ' storco  nl,lm : ',nl,lm
        CALL READCO(NL,LM,
     >                    KART,LET,YZXB,NDX,*10,*79)

C----- NDX: 1->KEX, 2->IT, 3->IREP, 4->IMAX

        NOC=NOC+1

             nrbltav = nrblt
        IF(NINT(YZXB(39)) .GE. NRBLT+1) NRBLT = NINT(YZXB(39)) -1
c             write(*,*) ' nrbltav, nrbltap ',nrbltav, nrblt

        IF(MODSTO.EQ.2) THEN
C------- Case call by ellipse fit
          IF(KP1.EQ.-2) THEN
C----------- Means only last pass wanted
            IF(NRBLT+1 .NE. NPASS) GOTO 44
          ELSE
C----------- Only pass in [KP1,KP2] wanted, or all passes if KP=-1 
C            Nothing to do, this is handled by OKKP in READCO
          ENDIF
        ENDIF
        CALL FILCOO(KPS,NOC,YZXB)
        IF(NOC.EQ. NPTR) GOTO 11

      GOTO 44             
C     ----------------------------------

 79   CONTINUE
      WRITE(6,*) ' *** Coordinates  storage  stopped: error during',
     > ' read of event # ',NOC+1
      GOTO 11

 10   CONTINUE
      WRITE(6,*) ' READ  OK; END  OF  FILE  ENCOUNTERED'

 11   CONTINUE
      NPASS = NRBLT + 1
      NPTR=NOC
      CALL READC5(KT1,KT2)
C        write(*,*) 'RRRRRRRRRRRRRRR  ',kt1, kt2
      IF(KT1 .EQ. -1 .OR. KT2 .GT. KT1) THEN
        WRITE(6,*) '  Analysis of particles from a set '
        IF(KPS.EQ. 0) WRITE(6,*) '  Initial  phase-space'
        IF(KPS.EQ. 1) WRITE(6,*) '  Final  phase-space'
        WRITE(6,*) '  # of turns  in the structure :',NPASS
      ELSEIF(KT1 .EQ. KT2) THEN
        WRITE(6,*) '  A single  particle  analized          :',KT1
        WRITE(6,*) '  # of turns in the structure   :',NPASS
      ENDIF

      WRITE(6,*) ' ',NOC,' points have been stored'

 96   RETURN                  
      END
