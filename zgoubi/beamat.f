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
     >                                 F0,PHY,PHZ,CSTRN,RPARAM)
C Transport the beam matrix. Initial beam matrix is in FI, set by OBJ5 or by TWISS.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(6,*),F0(6,*)
      LOGICAL PRDIC, OKCPLD
      COMMON/BEAM/ FI(6,6)
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M

      DIMENSION FII(6,6)
      LOGICAL OK, IDLUNI
      CHARACTER(50) CMMND
      PARAMETER (N4 = 4)
      PARAMETER (N6 = 6)
      DIMENSION RT66(N6,N6), RT(N4,N4)
      DIMENSION R44(N4,N4), RINV(N4,N4), RINT(N4,N4)
      DIMENSION N0(N4), B(N4)

C      DIMENSION R66(6,6)

      IF(.NOT. OKCPLD) THEN

        F0(1,1) = R(1,1)*R(1,1)*FI(1,1)-2.D0*R(1,1)*R(1,2)*FI(2,1) +
     >  R(1,2)*R(1,2)*FI(2,2)
        F0(1,2) = -(  R(1,1)*R(2,1)*FI(1,1) -
     >  ( 1.D0 + 2.D0*R(1,2)*R(2,1) )*FI(2,1)+ R(1,2)*R(2,2)*FI(2,2) )
        F0(2,1) = F0(1,2) 
        F0(2,2) = R(2,1)*R(2,1)*FI(1,1) - 2.D0*R(2,2)*R(2,1)*FI(2,1) +
     >  R(2,2)*R(2,2)*FI(2,2)

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

C Betatron phase advance
c        PHY = atan2(R(1,2) , ( R(1,1)*F0(1,1) - R(1,2)*F0(1,2)))
c        PHZ = atan2(R(3,4) , ( R(3,3)*F0(3,3) - R(3,4)*F0(3,4)))
        PHY = atan2(R(1,2) , ( R(1,1)*FI(1,1) - R(1,2)*FI(1,2)))
        IF(PHY.LT.0.D0) PHY = 2.D0*PI + PHY
        PHZ = atan2(R(3,4) , ( R(3,3)*FI(3,3) - R(3,4)*FI(3,4)))
        IF(PHZ.LT.0.D0) PHZ = 2.D0*PI + PHZ

        CSTRN = 0.D0
        RPARAM = 0.D0

      ELSEIF(OKCPLD) then


       if(.not. PRDIC) then

          call endjob(' Case non-periodic coupled beam-matrix '
     >    //' transport to be installed.',-99)

       else


              stop ' //////////// ici '

C R is the matrix from OBJET down to here. Make it N4xN4
        DO J=1,N4
          DO I=1,N4
            R44(I,J) = R(I,J)
            RINV(I,J) = R44(I,J)
          ENDDO
        ENDDO

Compute the inverse of R
        CALL DLAIN(N4,N4,RINV,N0,B,IER)

C Get the 6*6 1-turn map as stored by first TWISS pass
        CALL TWISS1(
     >              RT66)
C Get the 4*4 part of it
        do j = 1, 4
          do i = 1, 4
            rt(i,j) = rt66(i,j)
          enddo
        enddo
        CALL ZGNOEL(
     >               NOEL)
        CALL PMAT(RT,RINV,RINT,N4,N4,N4)
        CALL PMAT(R44,RINT,RT,N4,N4,N4)
C RT now contains the local 1-turn map 

C        ELSEIF(OKCPLD) THEN

c          OK = IDLUNI(
c     >              LUNW)
c          IF(.NOT. OK) CALL ENDJOB(
c     >    'SBR BEAMAT. Problem open idle unit for WRITE. ',-99)
c          OPEN(lunW,FILE='transfertM.dat',STATUS='UNKNOWN',IOSTAT=IOS2)
c          WRITE(lunW,FMT='(//)')
c          WRITE(lunW,FMT='(
c     >    ''TRANSPORT MATRIX (written by Zgoubi, used by ETparam ):'')')
c          DO I=1,N4
c             WRITE(lunW,FMT='(6(F15.8,1X))') (RT(I,J),J=1,N4)
c          ENDDO
c          WRITE(lunW,FMT='(/,a,i6)') ' Element # ',noel
c          CLOSE(lunW,IOSTAT=IOS2)
cCompute coupled optics
c          cmmnd = '~/zgoubi/current/coupling/ETparam'
c          CALL SYSTEM(cmmnd)
c          OK = IDLUNI(
c     >                 LUNW)
c          call et2res(lunw)
c          CLOSE(lunW)  
c          call et2re1(
c     >             F011,f012,f033,f034,phy,phz,Cstrn)

          CALL TUNESC(RT, 
     >                   F0,YNU,ZNU,CMUY,CMUZ,IERY,IERZ,RPARAM,CSTRN)

c          F0(1,1) = F011
c          F0(1,2) = f012
c          F0(2,1) = F0(1,2) 
c          F0(2,2) = (1.d0 + f012*f012 ) / F011
c          F0(3,3) = f033
c          F0(3,4) = f034
c          F0(4,3) = F0(3,4) 
c          F0(4,4) = (1.d0 + f034*f034 ) / F033

c          IF(PHY.LT.0.D0) PHY = 2.D0*PI + PHY
c          IF(PHZ.LT.0.D0) PHZ = 2.D0*PI + PHZ

          F0(1,6) = -((RT(1,4)*
     -        (-((R(1,6)*(RT(2,4)*(-1 + RT(3,3)) - RT(2,3)*RT(3,4)) + 
     -         RT(1,4)*(R(2,6) - R(2,6)*RT(3,3) + RT(2,3)*R(3,6)) + 
     -                RT(1,3)*(R(2,6)*RT(3,4) - RT(2,4)*R(3,6)))*
     -              ((RT(1,4)*RT(2,3) - RT(1,3)*RT(2,4))*
     -                 (RT(1,4)*RT(4,2) - RT(1,2)*(-1 + RT(4,4))) - 
     -                (RT(1,4)*(-1 + RT(2,2)) - RT(1,2)*RT(2,4))*
     -                 (RT(1,4)*RT(4,3) - RT(1,3)*(-1 + RT(4,4))))) + 
     -           ((RT(1,4)*RT(2,3) - RT(1,3)*RT(2,4))*
     -               (RT(1,4)*RT(3,2) - RT(1,2)*RT(3,4)) - 
     -              (RT(1,4)*(-1 + RT(2,2)) - RT(1,2)*RT(2,4))*
     -               (RT(1,4)*(-1 + RT(3,3)) - RT(1,3)*RT(3,4)))*
     -          (R(1,6)*(RT(2,4)*RT(4,3) - RT(2,3)*(-1 + RT(4,4))) + 
     -             RT(1,4)*(-(R(2,6)*RT(4,3)) + RT(2,3)*R(4,6)) + 
     -         RT(1,3)*(R(2,6)*(-1 + RT(4,4)) - RT(2,4)*R(4,6)))))/
     -       (-(((RT(1,4)*RT(2,3) - RT(1,3)*RT(2,4))*
     -               (RT(1,4)*RT(3,1) - (-1 + RT(1,1))*RT(3,4)) - 
     -              (RT(1,4)*RT(2,1) - (-1 + RT(1,1))*RT(2,4))*
     -               (RT(1,4)*(-1 + RT(3,3)) - RT(1,3)*RT(3,4)))*
     -            ((RT(1,4)*RT(2,3) - RT(1,3)*RT(2,4))*
     -               (RT(1,4)*RT(4,2) - RT(1,2)*(-1 + RT(4,4))) - 
     -              (RT(1,4)*(-1 + RT(2,2)) - RT(1,2)*RT(2,4))*
     -               (RT(1,4)*RT(4,3) - RT(1,3)*(-1 + RT(4,4))))) + 
     -         RT(1,4)*((RT(1,4)*RT(2,3) - RT(1,3)*RT(2,4))*
     -             (RT(1,4)*RT(3,2) - RT(1,2)*RT(3,4)) - 
     -            (RT(1,4)*(-1 + RT(2,2)) - RT(1,2)*RT(2,4))*
     -             (RT(1,4)*(-1 + RT(3,3)) - RT(1,3)*RT(3,4)))*
     -          (RT(2,3)*(-1 + RT(1,1) + RT(1,4)*RT(4,1)) - 
     -            RT(1,4)*RT(2,1)*RT(4,3) + 
     -       RT(1,3)*(-(RT(2,4)*RT(4,1)) + RT(2,1)*(-1 + RT(4,4))) + 
     -       (-1 + RT(1,1))*(RT(2,4)*RT(4,3) - RT(2,3)*RT(4,4)))))

          F0(2,6) = ((-(R(1,6)*RT(2,4)) + RT(1,4)*R(2,6))*
     -        (RT(1,4)*(-1 + RT(3,3)) - RT(1,3)*RT(3,4)) - 
     -       (RT(1,4)*RT(2,3) - RT(1,3)*RT(2,4))*
     -        (-(R(1,6)*RT(3,4)) + RT(1,4)*R(3,6)) + 
     -       (RT(1,4)*((RT(1,4)*RT(2,3) - RT(1,3)*RT(2,4))*
     -             (RT(1,4)*RT(3,1) - (-1 + RT(1,1))*RT(3,4)) - 
     -            (RT(1,4)*RT(2,1) - (-1 + RT(1,1))*RT(2,4))*
     -             (RT(1,4)*(-1 + RT(3,3)) - RT(1,3)*RT(3,4)))*
     -          (-((R(1,6)*(RT(2,4)*(-1 + RT(3,3)) - RT(2,3)*RT(3,4)) +
     -       RT(1,4)*(R(2,6) - R(2,6)*RT(3,3) + RT(2,3)*R(3,6)) + 
     -                 RT(1,3)*(R(2,6)*RT(3,4) - RT(2,4)*R(3,6)))*
     -               ((RT(1,4)*RT(2,3) - RT(1,3)*RT(2,4))*
     -                  (RT(1,4)*RT(4,2) - RT(1,2)*(-1 + RT(4,4))) - 
     -                 (RT(1,4)*(-1 + RT(2,2)) - RT(1,2)*RT(2,4))*
     -                  (RT(1,4)*RT(4,3) - RT(1,3)*(-1 + RT(4,4))))) + 
     -            ((RT(1,4)*RT(2,3) - RT(1,3)*RT(2,4))*
     -                (RT(1,4)*RT(3,2) - RT(1,2)*RT(3,4)) - 
     -               (RT(1,4)*(-1 + RT(2,2)) - RT(1,2)*RT(2,4))*
     -                (RT(1,4)*(-1 + RT(3,3)) - RT(1,3)*RT(3,4)))*
     -             (R(1,6)*(RT(2,4)*RT(4,3) - RT(2,3)*(-1 + RT(4,4))) + 
     -               RT(1,4)*(-(R(2,6)*RT(4,3)) + RT(2,3)*R(4,6)) + 
     -          RT(1,3)*(R(2,6)*(-1 + RT(4,4)) - RT(2,4)*R(4,6)))))/
     -        (-(((RT(1,4)*RT(2,3) - RT(1,3)*RT(2,4))*
     -                (RT(1,4)*RT(3,1) - (-1 + RT(1,1))*RT(3,4)) - 
     -               (RT(1,4)*RT(2,1) - (-1 + RT(1,1))*RT(2,4))*
     -                (RT(1,4)*(-1 + RT(3,3)) - RT(1,3)*RT(3,4)))*
     -             ((RT(1,4)*RT(2,3) - RT(1,3)*RT(2,4))*
     -                (RT(1,4)*RT(4,2) - RT(1,2)*(-1 + RT(4,4))) - 
     -               (RT(1,4)*(-1 + RT(2,2)) - RT(1,2)*RT(2,4))*
     -                (RT(1,4)*RT(4,3) - RT(1,3)*(-1 + RT(4,4))))) + 
     -          RT(1,4)*((RT(1,4)*RT(2,3) - RT(1,3)*RT(2,4))*
     -              (RT(1,4)*RT(3,2) - RT(1,2)*RT(3,4)) - 
     -             (RT(1,4)*(-1 + RT(2,2)) - RT(1,2)*RT(2,4))*
     -              (RT(1,4)*(-1 + RT(3,3)) - RT(1,3)*RT(3,4)))*
     -           (RT(2,3)*(-1 + RT(1,1) + RT(1,4)*RT(4,1)) - 
     -             RT(1,4)*RT(2,1)*RT(4,3) + 
     -         RT(1,3)*(-(RT(2,4)*RT(4,1)) + RT(2,1)*(-1 + RT(4,4))) + 
     -         (-1 + RT(1,1))*(RT(2,4)*RT(4,3) - RT(2,3)*RT(4,4)))))/
     -     ((RT(1,4)*RT(2,3) - RT(1,3)*RT(2,4))*
     -        (RT(1,4)*RT(3,2) - RT(1,2)*RT(3,4)) - 
     -       (RT(1,4)*(-1 + RT(2,2)) - RT(1,2)*RT(2,4))*
     -        (RT(1,4)*(-1 + RT(3,3)) - RT(1,3)*RT(3,4)))

          F0(3,6) = (R(1,6)*RT(2,4)*RT(3,2) - RT(1,4)*R(2,6)*RT(3,2) + 
     -    R(1,6)*RT(3,4) - R(1,6)*RT(2,2)*RT(3,4) 
     -                           + RT(1,2)*R(2,6)*RT(3,4) - 
     -    RT(1,4)*R(3,6) + RT(1,4)*RT(2,2)*R(3,6) 
     -                             - RT(1,2)*RT(2,4)*R(3,6) + 
     -    (RT(1,4)*(RT(1,4)*(-((-1+RT(2,2))*RT(3,1)) + RT(2,1)*RT(3,2))+ 
     -            RT(1,2)*(RT(2,4)*RT(3,1) - RT(2,1)*RT(3,4)) - 
     -    (-1 + RT(1,1))*(RT(2,4)*RT(3,2) + RT(3,4) - RT(2,2)*RT(3,4)))*
     -          (-((R(1,6)*(RT(2,4)*(-1 + RT(3,3)) - RT(2,3)*RT(3,4)) + 
     -          RT(1,4)*(R(2,6) - R(2,6)*RT(3,3) + RT(2,3)*R(3,6)) + 
     -                 RT(1,3)*(R(2,6)*RT(3,4) - RT(2,4)*R(3,6)))*
     -               ((RT(1,4)*RT(2,3) - RT(1,3)*RT(2,4))*
     -                  (RT(1,4)*RT(4,2) - RT(1,2)*(-1 + RT(4,4))) - 
     -                 (RT(1,4)*(-1 + RT(2,2)) - RT(1,2)*RT(2,4))*
     -                  (RT(1,4)*RT(4,3) - RT(1,3)*(-1 + RT(4,4))))) + 
     -            ((RT(1,4)*RT(2,3) - RT(1,3)*RT(2,4))*
     -                (RT(1,4)*RT(3,2) - RT(1,2)*RT(3,4)) - 
     -               (RT(1,4)*(-1 + RT(2,2)) - RT(1,2)*RT(2,4))*
     -                (RT(1,4)*(-1 + RT(3,3)) - RT(1,3)*RT(3,4)))*
     -             (R(1,6)*(RT(2,4)*RT(4,3) - RT(2,3)*(-1 + RT(4,4))) + 
     -               RT(1,4)*(-(R(2,6)*RT(4,3)) + RT(2,3)*R(4,6)) + 
     -       RT(1,3)*(R(2,6)*(-1 + RT(4,4)) - RT(2,4)*R(4,6)))))/
     -        (-(((RT(1,4)*RT(2,3) - RT(1,3)*RT(2,4))*
     -                (RT(1,4)*RT(3,1) - (-1 + RT(1,1))*RT(3,4)) - 
     -               (RT(1,4)*RT(2,1) - (-1 + RT(1,1))*RT(2,4))*
     -                (RT(1,4)*(-1 + RT(3,3)) - RT(1,3)*RT(3,4)))*
     -             ((RT(1,4)*RT(2,3) - RT(1,3)*RT(2,4))*
     -                (RT(1,4)*RT(4,2) - RT(1,2)*(-1 + RT(4,4))) - 
     -               (RT(1,4)*(-1 + RT(2,2)) - RT(1,2)*RT(2,4))*
     -                (RT(1,4)*RT(4,3) - RT(1,3)*(-1 + RT(4,4))))) + 
     -          RT(1,4)*((RT(1,4)*RT(2,3) - RT(1,3)*RT(2,4))*
     -              (RT(1,4)*RT(3,2) - RT(1,2)*RT(3,4)) - 
     -             (RT(1,4)*(-1 + RT(2,2)) - RT(1,2)*RT(2,4))*
     -              (RT(1,4)*(-1 + RT(3,3)) - RT(1,3)*RT(3,4)))*
     -           (RT(2,3)*(-1 + RT(1,1) + RT(1,4)*RT(4,1)) - 
     -             RT(1,4)*RT(2,1)*RT(4,3) + 
     -          RT(1,3)*(-(RT(2,4)*RT(4,1)) + RT(2,1)*(-1 + RT(4,4))) + 
     -           (-1 + RT(1,1))*(RT(2,4)*RT(4,3) - RT(2,3)*RT(4,4)))))/
     -    (RT(2,4)*(-(RT(1,3)*RT(3,2)) + RT(1,2)*(-1 + RT(3,3))) + 
     -    RT(1,4)*(-1+RT(2,3)*RT(3,2)-RT(2,2)*(-1+ RT(3,3)) + RT(3,3))+
     -       (RT(1,3)*(-1 + RT(2,2)) - RT(1,2)*RT(2,3))*RT(3,4))

          F0(4,6) = 
     -    (-(RT(1,3)*R(2,6)*RT(3,2)*RT(4,1)) - RT(1,3)*R(3,6)*RT(4,1) + 
     -       RT(1,3)*RT(2,2)*R(3,6)*RT(4,1) 
     -               - R(2,6)*RT(4,2) + RT(1,1)*R(2,6)*RT(4,2) + 
     -       RT(1,3)*R(2,6)*RT(3,1)*RT(4,2) + R(2,6)*RT(3,3)*RT(4,2) - 
     -       RT(1,1)*R(2,6)*RT(3,3)*RT(4,2) 
     -               - RT(1,3)*RT(2,1)*R(3,6)*RT(4,2) - 
     -       RT(2,3)*R(3,6)*RT(4,2) + RT(1,1)*RT(2,3)*R(3,6)*RT(4,2) - 
     -       R(2,6)*RT(3,2)*RT(4,3) 
     -             + RT(1,1)*R(2,6)*RT(3,2)*RT(4,3) - R(3,6)*RT(4,3) + 
     -       RT(1,1)*R(3,6)*RT(4,3) + RT(2,2)*R(3,6)*RT(4,3) - 
     -       RT(1,1)*RT(2,2)*R(3,6)*RT(4,3) - 
     -       R(1,6)*(-((-1 + RT(2,3)*RT(3,2) 
     -               - RT(2,2)*(-1 + RT(3,3)) + RT(3,3))*RT(4,1)) + 
     -          RT(3,1)*(RT(2,3)*RT(4,2) + RT(4,3) - RT(2,2)*RT(4,3)) + 
     -          RT(2,1)*(RT(4,2) - RT(3,3)*RT(4,2) + RT(3,2)*RT(4,3))) + 
     -       (-1 + RT(1,3)*RT(3,1) + (RT(1,3)*RT(2,1) 
     -               + RT(2,3))*RT(3,2) + RT(3,3) - 
     -          RT(2,2)*(-1 + RT(1,3)*RT(3,1) + RT(3,3)) - 
     -          RT(1,1)*(-1 + RT(2,3)*RT(3,2) - RT(2,2)*(-1 
     -               + RT(3,3)) + RT(3,3)))*R(4,6) + 
     -     RT(1,2)*(R(2,6)*((-1 + RT(3,3))*RT(4,1) - RT(3,1)*RT(4,3)) + 
     -          RT(2,3)*(-(R(3,6)*RT(4,1)) + RT(3,1)*R(4,6)) + 
     -          RT(2,1)*(R(3,6)*RT(4,3) + R(4,6) - RT(3,3)*R(4,6))))/
     -     (-1 + RT(1,3)*RT(3,1) + RT(2,2)*(1 - RT(1,3)*RT(3,1)) 
     -               + RT(1,3)*RT(2,1)*RT(3,2) + 
     -       RT(2,3)*RT(3,2) + RT(3,3) + RT(1,4)*RT(4,1) 
     -               - RT(1,4)*RT(2,3)*RT(3,2)*RT(4,1) + 
     -       RT(1,3)*RT(2,4)*RT(3,2)*RT(4,1) - RT(1,4)*RT(3,3)*RT(4,1) + 
     -       RT(1,4)*RT(2,2)*RT(3,3)*RT(4,1) + RT(1,3)*RT(3,4)*RT(4,1) - 
     -       RT(1,3)*RT(2,2)*RT(3,4)*RT(4,1) - RT(2,2)*(RT(3,3)
     -                + RT(1,4)*RT(4,1)) + 
     -       RT(1,4)*RT(2,1)*RT(4,2) + RT(2,4)*RT(4,2) 
     -               + RT(1,4)*RT(2,3)*RT(3,1)*RT(4,2) - 
     -       RT(1,3)*RT(2,4)*RT(3,1)*RT(4,2) 
     -               - RT(1,4)*RT(2,1)*RT(3,3)*RT(4,2) - 
     -       RT(2,4)*RT(3,3)*RT(4,2) + RT(1,3)*RT(2,1)*RT(3,4)*RT(4,2) + 
     -       RT(2,3)*RT(3,4)*RT(4,2) + RT(1,4)*RT(3,1)*RT(4,3) - 
     -       RT(1,4)*RT(2,2)*RT(3,1)*RT(4,3) 
     -               + RT(1,4)*RT(2,1)*RT(3,2)*RT(4,3) + 
     -       RT(2,4)*RT(3,2)*RT(4,3) + RT(3,4)*RT(4,3) 
     -               - RT(2,2)*RT(3,4)*RT(4,3) + RT(4,4) - 
     -       (RT(2,3)*RT(3,2) + RT(1,3)*(RT(3,1) 
     -               + RT(2,1)*RT(3,2)) + RT(3,3) - 
     -          RT(2,2)*(-1 + RT(1,3)*RT(3,1) + RT(3,3)))*RT(4,4) - 
     -       RT(1,1)*(-((-1 + RT(3,3))*(-1 + RT(2,4)*RT(4,2))) + 
     -          (RT(2,4)*RT(3,2) + RT(3,4))*RT(4,3) + 
     -          RT(2,3)*(RT(3,4)*RT(4,2) - RT(3,2)*(-1 
     -               + RT(4,4))) + RT(4,4) - 
     -          RT(3,3)*RT(4,4) - RT(2,2)*
     -           (-1 + RT(3,4)*RT(4,3) - RT(3,3)*(-1 
     -               + RT(4,4)) + RT(4,4))) + 
     -       RT(1,2)*(RT(2,4)*(RT(4,1) - RT(3,3)*RT(4,1) 
     -               + RT(3,1)*RT(4,3)) - 
     -          RT(2,1)*(-1 + RT(3,4)*RT(4,3) - RT(3,3)*(-1 
     -               + RT(4,4)) + RT(4,4)) + 
     -      RT(2,3)*(RT(3,1) + RT(3,4)*RT(4,1) - RT(3,1)*RT(4,4))))

        ENDIF
      ENDIF

      RETURN

      ENTRY BEAMA2(FII,sign)
      DO J=1,4
        DO I=1,4 
          FI(I,J) = sign**(i+j) * FII(I,J)  ! FI is in the form beta, +alpha, gamma
        ENDDO                               ! by contrast w F0 in the form beta, -alpha, gamma
      ENDDO
      DO I=1,4
        FI(I,6) = FII(I,6)
      ENDDO      
      RETURN

      END
