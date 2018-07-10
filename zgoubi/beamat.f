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
      SUBROUTINE BEAMAT(R,PRDIC,OKCPLD,
     >                                 F0,PHY,PHZ,CSTRN,RPRM)
C Transport the beam matrix. Initial beam matrix is in FI, set by OBJ5 or by TWISS.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(6,6),F0(6,6)
      LOGICAL PRDIC, OKCPLD
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M

      PARAMETER (N4 = 4)
      PARAMETER (N6 = 6)
      DIMENSION R44(N4,N4)

      DIMENSION C(2,2), P(4,4)
      DIMENSION FII(6,6)
      DIMENSION FI(6,6)
      SAVE FI

C Betatron phase advance
      PHY = ATAN2(R(1,2) , ( R(1,1)*FI(1,1) - R(1,2)*FI(1,2)))
      IF(PHY.LT.0.D0) PHY = 2.D0*PI + PHY
      PHZ = ATAN2(R(3,4) , ( R(3,3)*FI(3,3) - R(3,4)*FI(3,4)))
      IF(PHZ.LT.0.D0) PHZ = 2.D0*PI + PHZ

      IF(.NOT. OKCPLD) THEN

        F0(1,1) = R(1,1)*R(1,1)*FI(1,1)-2.D0*R(1,1)*R(1,2)*FI(2,1) +
     >  R(1,2)*R(1,2)*FI(2,2)
        F0(1,2) = -(  R(1,1)*R(2,1)*FI(1,1) -
     >  ( 1.D0 + 2.D0*R(1,2)*R(2,1) )*FI(2,1)+ R(1,2)*R(2,2)*FI(2,2) )
        F0(2,1) = F0(1,2) 
        F0(2,2) = R(2,1)*R(2,1)*FI(1,1) - 2.D0*R(2,2)*R(2,1)*FI(2,1) +
     >  R(2,2)*R(2,2)*FI(2,2)

c            call zgnoel(noelo)
c           noel = noelo
c          if(noel.eq.99) then
c                write(*,*) 'R(1,1),R(1,2)'
c                write(*,*) 'FI(1,1),FI(2,1),FI(2,2)'
c                write(*,*) 'F0(1,1),F0(1,2) ,F0(2,1),F0(2,2)'
c                write(*,*) R(1,1),R(1,2)
c                write(*,*) FI(1,1),FI(2,1),FI(2,2)
c                write(*,*) F0(1,1),F0(1,2) ,F0(2,1),F0(2,2)
c                write(*,*) 
c                  read(*,*)
c          endif


        F0(3,3) = R(3,3)*R(3,3)*FI(3,3)-2.D0*R(3,3)*R(3,4)*FI(4,3) +
     >  R(3,4)*R(3,4)*FI(4,4)
        F0(3,4) = -(  R(3,3)*R(4,3)*FI(3,3) -
     >  ( 1.D0 + 2.D0*R(3,4)*R(4,3) )*FI(4,3)+ R(3,4)*R(4,4)*FI(4,4) )
        F0(4,3) = F0(3,4) 
        F0(4,4) = R(4,3)*R(4,3)*FI(3,3) - 2.D0*R(4,4)*R(4,3)*FI(4,3) +
     >  R(4,4)*R(4,4)*FI(4,4)
  
        F0(1,6) = R(1,1)*FI(1,6) +R(1,2)*FI(2,6) +R(1,6)
        F0(2,6) = R(2,1)*FI(1,6) +R(2,2)*FI(2,6) +R(2,6)
        F0(3,6) = R(3,3)*FI(3,6) +R(3,4)*FI(4,6) +R(3,6)
        F0(4,6) = R(4,3)*FI(3,6) +R(4,4)*FI(4,6) +R(4,6)

        CSTRN = 0.D0
        RPRM = 0.D0

      ELSEIF(OKCPLD) then


        IF(.NOT. PRDIC) THEN

          CALL ENDJOB(' Case non-periodic coupled beam-matrix '
     >    //' transport to be installed.',-99)

        ELSE

C R is the 6x6 matrix from OBJET down to here. Make it N4xN4
          DO J=1,N4
            DO I=1,N4
              R44(I,J) = R(I,J)
            ENDDO
          ENDDO

C Get the 4x4 P matrix as stored by first TWISS pass
          CALL TUNES1(
     >                P,C)

          CALL PROPAG(R44,P,
     >                      F0,C,RPRM)

          F0(1,6) = R(1,1)*FI(1,6) +R(1,2)*FI(2,6) + 
     >             R(1,3)*FI(3,6) +R(1,4)*FI(4,6) + R(1,6)
          F0(2,6) = R(2,1)*FI(1,6) +R(2,2)*FI(2,6) + 
     >             R(2,3)*FI(3,6) +R(2,4)*FI(4,6) + R(2,6)
          F0(3,6) = R(3,3)*FI(3,6) +R(3,4)*FI(4,6) + 
     >             R(3,1)*FI(1,6) +R(3,2)*FI(2,6) + R(3,6)
          F0(4,6) = R(4,3)*FI(3,6) +R(4,4)*FI(4,6) + 
     >             R(4,1)*FI(1,6) +R(4,2)*FI(2,6) + R(4,6)

        ENDIF
      ENDIF

      RETURN

      ENTRY BEAMA1(
     >             FII)
      DO J=1,6
        DO I=1,6 
          FII(I,J) = FI(I,J) 
        ENDDO                
      ENDDO
      RETURN

      ENTRY BEAMA2(FII,SIGN)
      DO J=1,4
        DO I=1,4 
          FI(I,J) = SIGN**(I+J) * FII(I,J)  ! FI is in the form beta, +alpha, gamma
        ENDDO                               ! by contrast w F0 in the form beta, -alpha, gamma
      ENDDO
      DO I=1,4
        FI(I,6) = FII(I,6)
      ENDDO      
      RETURN

      END
