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
      SUBROUTINE STORCO(MODSTO,NL,KPS,
     >                                NPASS)
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
        CALL READC1(
     >              KP1,KP2,KP3)
        CALL LASTP(NL,
     >                NPASS)
      ENDIF

      NOC=0
      NRBLT = -1 
C----- BOUCLE SUR READ FICHIER NL 
 44   CONTINUE
              
        CALL READCO(NL,
     >                 KART,LET,YZXB,NDX,*10,*79)

C----- NDX: 1->KEX, 2->IT, 3->IREP, 4->IMAX

        NOC=NOC+1

             nrbltav = nrblt
        IF(NINT(YZXB(39)) .GE. NRBLT+1) NRBLT = NINT(YZXB(39)) -1

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
      WRITE(6,*) ' END OF FILE encountered, read ',NOC,' points'

 11   CONTINUE
      NPASS = NRBLT + 1
      NPTR=NOC
      CALL READC5(
     >            KT1,KT2)
      write(*,*) ' Pgm storco, particles kt1:kt2 : ',kt1,':',kt2
      IF(KT1 .EQ. -1 .OR. KT2 .GT. KT1) THEN
        WRITE(6,*) '  Analysis of particles from a set '
        IF(KPS.EQ. 0) WRITE(6,*) '  Initial  phase-space'
        IF(KPS.EQ. 1) WRITE(6,*) '  Current  phase-space'
        WRITE(6,*) '  # of turns  in the structure :',NPASS
      ELSEIF(KT1 .EQ. KT2) THEN
        WRITE(6,*) '  A single  particle  analysed          :',KT1
        WRITE(6,*) '  # of turns in the structure   :',NPASS
      ENDIF

      WRITE(6,*) ' ',NOC,' points have been stored'

 96   RETURN                  
      END
