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
      SUBROUTINE XYZBR(NL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      CHARACTER LET
      INCLUDE 'MXVAR.H'
      DIMENSION YZXB(MXVAR), NDX(5)

      CHARACTER LETO
      DIMENSION YZXBO(MXVAR), NDXO(5)

      INCLUDE 'MXSTEP.H'
      PARAMETER (MXC=16)
      DIMENSION YZXBU(MXSTEP,MXC)
      SAVE YZXBU

C------ read the all data trajectory file (e.g., zgoubi.plt)

      WRITE(6,*) 
      WRITE(6,*) ' Busy, reading and storing the all data file... '

      CALL REWIN2(NL,*12)
      NSTP=0
C      CALL READC5(KT,KT)
     
C----- Loop on data reading 
 44   CONTINUE                 

        NSTP=NSTP+1
        CALL READCO(NL,
     >                    KART,LET,YZXB,NDX,*11,*12)

        YZXBU(NSTP,1) = YZXB(9)           ! Size of next step (m)
        YZXBU(NSTP,2) = YZXB(6)          ! Path length (m)

        YZXBU(NSTP,3) = YZXB(30)        ! Bx, Tesla
        YZXBU(NSTP,4) = YZXB(31)        ! By
        YZXBU(NSTP,5) = YZXB(32)        ! Bz

C        CT = COS( YZXB(3) )
C        CP = COS( YZXB(5) )
        YZXBU(NSTP,6) = YZXB(3) 
        YZXBU(NSTP,7) = YZXB(5)

C------- Particle position in magnet frame  ( r(t') )
C        YZXBU(NSTP,8) = YZXB(7)
        YZXBU(NSTP,8) = YZXB(8)            ! X (m)
        YZXBU(NSTP,9) = YZXB(2)            ! Y (m)
        YZXBU(NSTP,10) = YZXB(4)           ! Z (m)

C        BRO = YZXB(36) * (1.D0 + YZXB(1))        ! Rigidity, T.m
        YZXBU(NSTP,11) = YZXB(36) 
        YZXBU(NSTP,12) = YZXB(1) 

        YZXBU(NSTP,13) = KART
C------------ IEX
        YZXBU(NSTP,14) = NDX(1)
        YZXBU(NSTP,15) = ICHAR(LET)
C-------- NEWL = NOEL .NE. NDX(5)
        YZXBU(NSTP,16) = NDX(5)

        IF(NSTP .EQ. MXSTEP) GOTO 13
      GOTO 44             
C     ------------------------------------------

 11   CONTINUE
      WRITE(6,*) ' Information : stopped reading at step # ',NSTP
      WRITE(6,*) '         upon end of data file'
      CALL SREFW(NSTP-1)
      RETURN 

 12   CONTINUE
      WRITE(6,*) ' Information : stopped reading at step # ',NSTP
      WRITE(6,*) '         upon read error'
      CALL SREFW(NSTP-1)
      RETURN 

 13   CONTINUE
      WRITE(6,*) ' Information : stopped reading at step # ',NSTP
      WRITE(6,*) '         upon reaching storage limit '
      CALL SREFW(NSTP)
      RETURN 

      ENTRY XYZBW(NRD,
     >               KARTO,LETO,YZXBO,NDXO)

C----- Get coordinates, fields, etc. at current trajectory step
      YZXBO(9) = YZXBU(NRD,1)
      YZXBO(6) = YZXBU(NRD,2)
      YZXBO(30) = YZXBU(NRD,3)
      YZXBO(31) = YZXBU(NRD,4)
      YZXBO(32) = YZXBU(NRD,5)
      YZXBO(3) = YZXBU(NRD,6)
      YZXBO(5) = YZXBU(NRD,7)
C      YZXBO(7) = YZXBU(NRD,8)
      YZXBO(8) = YZXBU(NRD,8)
      YZXBO(2) = YZXBU(NRD,9)
      YZXBO(4) = YZXBU(NRD,10)
      YZXBO(36) = YZXBU(NRD,11)
      YZXBO(1) = YZXBU(NRD,12)
      KARTO = INT(YZXBU(NRD,13) + 1.D-6)
      NDXO(1) = YZXBU(NRD,14)
      LETO = CHAR(INT(YZXBU(NRD,15)+ 1.D-6))
      NDXO(5) = YZXBU(NRD,16)
      RETURN

      END
