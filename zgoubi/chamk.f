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
C  Upton, NY, 11973, USA
C  -------
      SUBROUTINE CHAMK(A1,R1,Z1,*)
      USE DYNHC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (I3=3)
      PARAMETER (MXC = 4)
      DIMENSION NMPT(MXC), NMPTI(MXC)

C      SAVE NFIC, NMPT
      SAVE IOP

      DATA NFIC / 1 /
      DATA NMPT / MXC * 1 /

      IF    (NFIC .EQ. 1) THEN 
        CALL KSMAP(
     >           IMAP) 
        CALL CHAMK1(A1,R1,Z1,IMAP,*999)

      ELSEIF(NFIC .GE. 2) THEN 
C Case MOD=16. IMAP has been incrized by NFIC units in tosca, each field map stored in HC
        I = 1

c               write(*,*) 'chamk NMPT(I) ',i,NMPT(I)

        CALL CHAMK1(A1,R1,Z1,NMPT(I),*999)
        IOP = 1
        CALL ADPOL(ID,IOP,B,DB,DDB,D3BX,D3BY,D4BX,D4BY,BT)
        DO I = 2, NFIC-1

c               write(*,*) 'chamk NMPT(I) ',i,NMPT(I)

          CALL CHAMK1(A1,R1,Z1,NMPT(I),*999)          
          CALL ADPOL(ID,IOP,B,DB,DDB,D3BX,D3BY,D4BX,D4BY,BT)
        ENDDO
        I = NFIC

c               write(*,*) 'chamk NMPT(I) ',i,NMPT(I)

        CALL CHAMK1(A1,R1,Z1,NMPT(I),*999)
        CALL ADPOL(ID,I3,B,DB,DDB,D3BX,D3BY,D4BX,D4BY,BT)

      ELSE
        CALL ENDJOB('Pgm chamk. No such possibility NFIC = ',NFIC)
      ENDIF

c               write(*,*) ' ' 

      RETURN

 999  RETURN 1

      ENTRY CHAMKW(NFICI,NMPTI)
      NFIC = NFICI
      NMPT = NMPTI
      RETURN

      END
