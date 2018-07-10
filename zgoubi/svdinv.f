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
      SUBROUTINE SVDINV(LSVD,MR,NC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(MR,NC),U(MR,MR),S(NC),V(NC,NC)
      character(800) txt800

      a=0.d0
      
      WRITE(*,*) ' svdinv  mr, nc  ', mr, nc  
      read(*,*)
      REWIND(LSVD)
      DO I = 1, 8
        READ(LSVD,*)
      ENDDO

      DO I=1,MR
C       DO J=1,NC
          READ(LSVD,fmt='(a)',end=10,err=10) txt800
C          write(*,*) ' svdinv ',trim(txt800)
          READ(txt800,*,end=10,err=10) (A(I,J),j=1,nc)
C        END DO  
      END DO

 10   continue
      WRITE(*,*) 'A ='
      WRITE(*,fmt='(24(1x,e12.4))') ((A(I,J),J=1,nc),I=1,mr)

       CALL SVD(A,U,S,V,MR,NC)

      WRITE(*,*) 'U ='
      WRITE(*,'(24F10.4)') ((U(I,J),J=1,24),I=1,24)

      WRITE(*,*) 'V='
      WRITE(*,'(24F10.4)') V
      WRITE(*,*) 'S =' 
      WRITE(*,'(22F10.4)')S
      RETURN
      END
