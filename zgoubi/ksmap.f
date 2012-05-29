C  ZGOUBI, a program for computing the trajectories of charged particles
C  in electric and magnetic fields
C  Copyright (C) 1988-2007  François Meot
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
C  François Meot <fmeot@bnl.gov>
C  BNL
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE KSMAP(
     >                 KMAPO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C------------------------------
C     Monitors field map number
C------------------------------
      INCLUDE 'PARIZ.H'
      CHARACTER*(80) NAMSAV(MMAP,IZ)
      CHARACTER*(*) NOMFIC(*)
      LOGICAL NEWFIC, OLDFIC
      PARAMETER (MIZ = MMAP*IZ)

      SAVE NBMAPS, KMAP

      DATA NBMAPS, KMAP / 0, 0 /
      DATA NAMSAV / MIZ*' ' /

C Read
      KMAPO = KMAP
      RETURN

C Reset
      ENTRY KSMAP0
        NBMAPS = 0
      RETURN

C Stack, count 
      ENTRY KSMAP4(NOMFIC,NFIC,
     >                         NEWFIC,NBMAPO,KMAPO)
      IF(NFIC.GT.IZ) CALL ENDJOB(' SBR KSMAP. NFIC should be <',IZ)
      II = 0
      DO I = 1, MMAP
        II = II+1
        OLDFIC = .TRUE.
        DO IFIC = 1, NFIC
          OLDFIC = OLDFIC .AND. (NOMFIC(IFIC).EQ.NAMSAV(I,IFIC))
        ENDDO
        IF(OLDFIC) GOTO 1
      ENDDO
 1    CONTINUE
      NEWFIC = .NOT. OLDFIC
      IF(NEWFIC) THEN
        IF(NBMAPS.GE.MMAP) THEN 
          WRITE(NRES,*) ' SBR ksmap(b) : too many different field maps'
          CALL ENDJOB(' In PARIZ.H :  have  MMAP >',MMAP)
        ENDIF
        NBMAPS = NBMAPS+1
        KMAP = NBMAPS
        DO IFIC = 1, NFIC
          NAMSAV(NBMAPS,IFIC) = NOMFIC(IFIC)
        ENDDO
      ELSE
        KMAP = II
      ENDIF
      NBMAPO = NBMAPS
      KMAPO = KMAP
      RETURN

      END
