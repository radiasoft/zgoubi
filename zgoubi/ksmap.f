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
     >                 IMAPO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C---------------------------------------------------------------------------------------
C Monitors field map number. 
C IMAP is the number of a field map stored in HC(ID,IA,IR,IZ,IMAP) array, as a result of
C  - either individual reading if 1! file is read (e.g., MOD.MOD2=15.1)
C  - or a combination of IFIC=1,NFIC field maps from as many files (e.g., MOD.MOD2=15.4)
C---------------------------------------------------------------------------------------
      INCLUDE 'PARIZ.H'
      PARAMETER (MXC = 4)
      PARAMETER (NFM = IZ+MXC)
      CHARACTER*(80) NAMSAV(MMAP,NFM)
C      CHARACTER*(80) NAMSAV(MMAP,IZ)
      CHARACTER(*) NOMFIC(*)
      LOGICAL NEWFIC, OLDFIC
      PARAMETER (MIZ = MMAP*IZ)
 
      SAVE NBMAPS, IMAP
 
      PARAMETER (MMNF=MMAP*NFM)
      DATA NBMAPS, IMAP / 0, 0 /
      DATA NAMSAV / MMNF*' ' /
C      DATA NAMSAV / MIZ*' ' /
 
C     16/01/14 For KSMAP4 to remember the maps coeffcients used
C      PARAMETER (MXC = 4)
      DIMENSION COEFS(MMAP,NFM), AA(NFM)
C      DIMENSION COEFS(MMAP,IZ), AA(IZ)
 
C Read
      IMAPO = IMAP
 
      RETURN
 
C Reset
      ENTRY KSMAP0
 
      NBMAPS = 0
 
      RETURN 

C Stack, count.
      ENTRY KSMAP4(NOMFIC,NFIC,AA,
     >                         NEWFIC,NBMAPO,IMAPO)
C      IF(NFIC.GT.IZ) CALL ENDJOB(' SBR KSMAP. NFIC should be <',IZ)
      II = 0
      DO I = 1, MMAP
        II = II+1
        OLDFIC = .TRUE.
        DO IFIC = 1, NFIC
C OLDFIC stays true iff (the 1-NFIC series has already been met 
C AND the linear combination coefficients have not been changed)
          OLDFIC = OLDFIC .AND. (NOMFIC(IFIC).EQ.NAMSAV(I,IFIC))
     >          .AND. (AA(IFIC).EQ.COEFS(I,IFIC))
        ENDDO
        IF(OLDFIC) GOTO 1
      ENDDO

 1    CONTINUE

      NEWFIC = .NOT. OLDFIC
      IF(NEWFIC) THEN
        IF(NBMAPS.GE.MMAP) THEN
          WRITE(NRES,*) ' SBR ksmap : too many different field maps'
          CALL ENDJOB(' In PARIZ.H :  have  MMAP >',MMAP)
        ENDIF
        NBMAPS = NBMAPS+1
        IMAP = NBMAPS

            write(*,*) ' ksmap imap ',imap
        DO IFIC = 1, NFIC
          NAMSAV(IMAP,IFIC) = NOMFIC(IFIC)
C     16/01/14 to save the coefficients
          COEFS(IMAP,IFIC) = AA(IFIC)
       ENDDO
      ELSE
        IMAP = II
      ENDIF
C Total number of maps stored in HC(ID,IA,IR,IZ,IMP=1,NBMAPS)
      NBMAPO = NBMAPS
C Number of this map, HC(ID,IA,IR,IZ,IMAP)
      IMAPO = IMAP
      RETURN
 
      ENTRY KSMAP5(
     >             IMAPO)
C Stack, count, case self-made maps : dipole-m, aimant, cartemes, etc.
C INSTALLATION HAS TO BE COMPLETED,
C It is not in the present state compatible with KSMAP4, neither does it allow stacking maps.
      NBMAPS = NBMAPS+1
      IMAP = NBMAPS
      IMAPO = IMAP
      RETURN
 
      END
