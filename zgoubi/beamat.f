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
C  Upton, NY, 11973
C  -------
      SUBROUTINE BEAMAT(R,PRDIC,OKCPLD,
     >                                 F0,PHY,PHZ,Cstrn)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(6,*),F0(6,*)
      LOGICAL PRDIC, OKCPLD
      COMMON/BEAM/ FI(6,6)
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M

      DIMENSION FII(6,6)
      LOGICAL OK, IDLUNI
      INTEGER DEBSTR, FINSTR
      CHARACTER(50) CMMND
      PARAMETER (N4 = 4)
      DIMENSION RT(N4,N4)
      DIMENSION R44(N4,N4), RINV(N4,N4), RINT(N4,N4)
      DIMENSION N0(N4), B(N4)

      IF(.NOT. PRDIC) THEN
        F0(1,1) = R(1,1)*R(1,1)*FI(1,1)-2.D0*R(1,1)*R(1,2)*FI(2,1) +
     >  R(1,2)*R(1,2)*FI(2,2)
        F0(1,2) = -(  R(1,1)*R(2,1)*FI(1,1) -
     >  ( 1.D0 + 2.D0*R(1,2)*R(2,1) )*FI(2,1) +
     >  R(1,2)*R(2,2)*FI(2,2)  )
        F0(2,1) = F0(1,2) 
        F0(2,2) = R(2,1)*R(2,1)*FI(1,1) - 2.D0*R(2,2)*R(2,1)*FI(2,1) +
     >  R(2,2)*R(2,2)*FI(2,2)
        F0(3,3) = R(3,3)*R(3,3)*FI(3,3)-2.D0*R(3,3)*R(3,4)*FI(4,3) +
     >  R(3,4)*R(3,4)*FI(4,4)
        F0(3,4) = -(  R(3,3)*R(4,3)*FI(3,3) -
     >  ( 1.D0 + 2.D0*R(3,4)*R(4,3) )*FI(4,3) +
     >  R(3,4)*R(4,4)*FI(4,4)  )
        F0(4,3) = F0(3,4) 
        F0(4,4) = R(4,3)*R(4,3)*FI(3,3) - 2.D0*R(4,4)*R(4,3)*FI(4,3) +
     >  R(4,4)*R(4,4)*FI(4,4)
        F0(1,6) = R(1,1)*FI(1,6) +R(1,2)*FI(2,6) +R(1,3)*FI(3,6) +
     >  R(1,4)*FI(4,6) +R(1,5)*FI(5,6) +R(1,6)*FI(6,6) 
        F0(2,6) = R(2,1)*FI(1,6) +R(2,2)*FI(2,6) +R(2,3)*FI(3,6) + 
     >  R(2,4)*FI(4,6) +R(2,5)*FI(5,6) +R(2,6)*FI(6,6) 
        F0(3,6) = R(3,1)*FI(1,6) +R(3,2)*FI(2,6) +R(3,3)*FI(3,6) + 
     >  R(3,4)*FI(4,6) +R(3,5)*FI(5,6) +R(3,6)*FI(6,6) 
        F0(4,6) = R(4,1)*FI(1,6) +R(4,2)*FI(2,6) +R(4,3)*FI(3,6) + 
     >  R(4,4)*FI(4,6) +R(4,5)*FI(5,6) +R(4,6)*FI(6,6) 

C Betatron phase advance
c        PHY = atan2(R(1,2) , ( R(1,1)*F0(1,1) - R(1,2)*F0(1,2)))
c        PHZ = atan2(R(3,4) , ( R(3,3)*F0(3,3) - R(3,4)*F0(3,4)))
        PHY = atan2(R(1,2) , ( R(1,1)*FI(1,1) - R(1,2)*FI(1,2)))
        IF(PHY.LT.0.D0) PHY = 2.D0*PI + PHY
        PHZ = atan2(R(3,4) , ( R(3,3)*FI(3,3) - R(3,4)*FI(3,4)))
        IF(PHZ.LT.0.D0) PHZ = 2.D0*PI + PHZ

        Cstrn = 0.d0

      ELSE

C R is the matrix from OBJET down to here. Make it 4x4
        DO J=1,4
          DO I=1,4
            R44(I,J) = R(I,J)
            RINV(I,J) = R44(I,J)
          ENDDO
        ENDDO

Compute the inverse of R
        CALL DLAIN(N4,N4,RINV,N0,B,IER)

C Get the 4x4 1-turn map
        call twiss1(
     >                RT)
        call ZGNOEL(
     >               NOEL)
        call pmat(RT,RINV,RINT,4,4,4)
        call pmat(R44,RINT,RT,4,4,4)
C RT now contains the local 1-turn map 

        IF(OKCPLD) THEN

          OK = IDLUNI(
     >              LUNW)
          IF(.NOT. OK) CALL ENDJOB(
     >    'SBR BEAMAT. Problem open idle unit for WRITE. ',-99)

          OPEN(lunW,FILE='transfertM.dat',STATUS='UNKNOWN',IOSTAT=IOS2)
          WRITE(lunW,FMT='(//)')
          WRITE(lunW,FMT='(
     >    ''TRANSPORT MATRIX (written by Zgoubi, used by ETparam ):'')')
          DO I=1,4
             WRITE(lunW,FMT='(4(F15.8,1X))') (RT(I,J),J=1,4)
          ENDDO
          WRITE(lunW,FMT='(/,a,i6)') ' Element # ',noel

          CLOSE(lunW,IOSTAT=IOS2)

Compute coupled optics
          cmmnd = '~/zgoubi/current/coupling/ETparam'
c        write(6,*) ' Pgm beamat. Now doing ' 
c     >  // cmmnd(debstr(cmmnd):finstr(cmmnd))
          CALL SYSTEM(cmmnd)

          OK = IDLUNI(
     >                 LUNW)
          call et2res(lunw)
          CLOSE(lunW)
   
          call et2re1(
     >             F011,f012,f033,f034,phy,phz,Cstrn)

          F0(1,1) = F011
          F0(1,2) = f012
          F0(2,1) = F0(1,2) 
          F0(2,2) = (1.d0 + f012*f012 ) / F011
          F0(3,3) = f033
          F0(3,4) = f034
          F0(4,3) = F0(3,4) 
          F0(4,4) = (1.d0 + f034*f034 ) / F033

          IF(PHY.LT.0.D0) PHY = 2.D0*PI + PHY
          IF(PHZ.LT.0.D0) PHZ = 2.D0*PI + PHZ

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
     - R(1,6)*RT(3,4) - R(1,6)*RT(2,2)*RT(3,4) 
     -                           + RT(1,2)*R(2,6)*RT(3,4) - 
     - RT(1,4)*R(3,6) + RT(1,4)*RT(2,2)*R(3,6) 
     -                             - RT(1,2)*RT(2,4)*R(3,6) + 
     - (RT(1,4)*(RT(1,4)*(-((-1 + RT(2,2))*RT(3,1)) + RT(2,1)*RT(3,2)) + 
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
     -     (RT(2,4)*(-(RT(1,3)*RT(3,2)) + RT(1,2)*(-1 + RT(3,3))) + 
     - RT(1,4)*(-1+ RT(2,3)*RT(3,2)- RT(2,2)*(-1 + RT(3,3)) + RT(3,3))+ 
     -       (RT(1,3)*(-1 + RT(2,2)) - RT(1,2)*RT(2,3))*RT(3,4))

       F0(4,6) = 
     -  (-(RT(1,3)*R(2,6)*RT(3,2)*RT(4,1)) - RT(1,3)*R(3,6)*RT(4,1) + 
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

        ELSEIF(.NOT. OKCPLD) THEN

          CALL TUNES(R,F0,1,IERY,IERZ,.FALSE.,
     >                                       PHY,PHZ,CMUY,CMUZ)          

        ENDIF
      ENDIF

      RETURN

      ENTRY BEAMA1(FII)
      DO 2 J=1,6
        DO 2 I=1,6 
          FI(I,J) = FII(I,J)
 2    CONTINUE
      RETURN
      END
