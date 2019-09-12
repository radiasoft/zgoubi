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
      SUBROUTINE ELCYLM(MPOL,EM,QLEM,QLSM,QE,QS,A,R,Z,
     >                                          E,DE,DDE)
C Mid-plane field, cylindrical coordinates.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(MCOEF=6)
      DIMENSION QLEM(MPOL),QLSM(MPOL),QE(MPOL,MCOEF),QS(MPOL,MCOEF)
      DIMENSION EM(*)
      DIMENSION E(5,3),DE(3,3),DDE(3,3,3)

      INCLUDE "C.AIM.H"     ! COMMON/AIM/ AE,AT,AS,RM,XI,XF,EN,EB1,EB2,EG1,EG2
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.CHAFUI.H"     ! COMMON/CHAFUI/ XE,XS,CE(MCOEF),CS(MCOEF),QCE(MCOEF),QCS(MCOEF)
      INCLUDE "C.INTEG.H"     ! COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,HDPRF,DP,QBR,BRI

C------- Extent of Entrance/Exit fringe field region
C      EQUIVALENCE (EB1,XLE), (EB2,XLS), (EG1,V0), (EG2,u), (RM,rho)
      EQUIVALENCE (EB1,XLE), (EB2,XLS), (EG1,E0), (EG2,u), (RM,rho)

      LOGICAL CHFE, CHU, CHFS
      PARAMETER (I0=0, I1=1)

      DUM = EM(1)

      VN = V0*BRI

C----- LambdaE,LambdaS
      QLE = QLEM(1)
      QLS = QLSM(1)

      AA = A
      RR = R
      R0R = RM/R
      R0RA = R0R**ANDX
      ZR = Z / R
      ZR2 = ZR*ZR

      IF(QLE+QLS .EQ. 0D0) THEN
C------- Sharp edges approximation

        IF ( AA .GE. AE .AND. AA .LE. AE+AT ) THEN
C         Uniform region

CCCCCCCC          E(1,2) = - VN / RR
      ! in Zgoubi y is radial, x longitudinal, z vertical               !AUL:110214
c          print*,'E0,rho,E0/rho,y=',EM(1),rho,EM(1)/rho,y;stop

c          c = 0.d0   ! cylindrical
cc          c = 1.d0   ! spherical
c          u = (rr-rm)/rm
c          v = z/rm  !!!!! z is not known

C Er
c          E(1,2) =  V0 *( 1.d0 - u * (1.d0 + c))
C Ez
c          E(1,3) =  V0 * v * c
c
c          DE(2,2) =  -EM(1)*BRI * (1.d0 + c)
c          DE(3,3) =  EM(1)*BRI * c

c          E(1,2) = EM(1)*BRI * R0RA * (
c     >    1.d0 - 0.5d0*(andx*andx-1.d0) * zr2 +
c     >    1.d0/24.d0*(andx*andx-1.d0)*(andx+1.d0)*(andx+3.d0) * zr2*zr2
c     >    )
c          E(1,3) = EM(1)*BRI * R0RA * (
c     >    (andx-1.d0) * zr -
c     >    1.d0/6.d0*(andx*andx-1.d0)*(andx+1.d0)*(andx+3.d0) * zr2*zr
c     >    )

C          E(1,2) = EM(1)*BRI * R0R
C          E(1,2) = EM(1)*BRI
c          E(1,3) = 0.d0

c          e(2,2) = - e(1,2) / R
c          e(3,2) = - e(2,2) * 2.d0/R
c          e(4,2) = - e(3,2) * 3.d0/R
c          e(5,2) = - e(4,2) * 4.d0/R

c      write(*,fmt='(a,i2,1x,20e12.4)')
c     > '  ex, ey, ez : ',e(1,1),e(1,2),e(1,3)

      EI = E0 * BRI

      E(1,2) = (EI*(1.d0 - ((-1.d0 + u**2)*z**2)/(2.d0*r**2) +
     - ((1.d0 + u)*(3.d0 + u)*(-1.d0 + u**2)*z**4)/(24.d0*r**4)))/r**u

      E(1,3) = (EI*(((-1.d0 + u)*z)/r -
     -      ((1.d0 + u)*(-1.d0 + u**2)*z**3)/(6.d0*r**3)))/r**u


c      write(88,fmt='(a,1x,20e12.4)')
c     > ' elcylm  ex, ey, ez : ',a,r,e(1,1),e(1,2)
C              read(*,*)
C      dE2/dr
      DE(2,2)=  (EI*r**(-5.d0 - u)*(-24.d0*r**4*u +
     -    12.d0*r**2*(-1.d0 + u)*(1.d0 + u)*(2.d0 + u)*z**2 -
     -    (-1.d0 + u)*(1.d0 + u)**2*(3.d0 + u)*(4.d0 + u)*z**4))/24.d0
C      dEr/dz, dEz/dr
      DE(3,2)=(EI*r**(-4 - u)*(-1.d0 + u**2)
     -   *z*(-6.d0*r**2 + (1.d0 + u)*(3.d0 + u)*z**2))/6.d0

C      d2Erdr2
      DDE(2,2,2)=  (EI*r**(-6.d0 - u)*(1.d0+ u)*
     - (24.D0*r**4*u - 12.d0*r**2*(-1.d0+ u)*(2.d0 + u)*(3.d0+ u)*z**2 +
     - (-1.d0+ u)*(1.d0+ u)*(3.d0+ u)*(4.d0+ u)*(5.D0 + u)*z**4))/24.D0
C      d2Erdrdz
      DDE(3,2,2) = (EI*r**(-5.d0 - u)*(-1.d0 + u)*(1.d0 + u)*z*
     - (6.d0*r**2*(2.d0+ u) - (1.d0+ u)*(3.d0+ u)*(4.d0+ u)*z**2))/6.d0

C      d2Erdz2
      DDE(3,3,2)=(EI*r**(-4.d0 - u)*(-1.d0 + u**2)*
     -  (-2.d0 *r**2 + (1.d0 + u)*(3.d0 + u)*z**2))/2.d0


        ENDIF

      ELSE

C        Position w.r.t. Entrance/Exit EFB

        CHFE= AA .LE. 2D0*AE .AND. QLE.NE.0D0
        CHU=  AA .GT. 2D0*AE .AND. AA .LT. AT - 2D0*AS
        CHFS= AA .GE. AT - 2D0*AS .AND. QLS.NE.0D0

C        write(*,*) ' chfe, chu, chfs  :',chfe, chu, chfs

        IF(CHU) THEN
C---------- Central region

          E(1,2) = - VN / RR

        ELSE
C---------- Entrance or exit fringe field

          IF(.NOT. CHFE) THEN
C           *** Either beyond effect of entrance fringe
C               or     hard edge

            GE  =1D0
            DGE = 0D0

          ELSE
C           *** Entrance fringe

C------------ SE>0 outside
            SE=(AE-AA)*RR/QLE
            CALL DRVG(I1,QCE,SE,GE,DGE,D2GE,D3GE,D4GE,D5GE,D6GE)
C     >                                   D7GE,D8GE,D9GE,D10GE)

            DSEA = - 1.D0

          ENDIF

          IF(.NOT. CHFS) THEN
C           *** Either beyond effect of exit fringe
C               or     hard edge

            GS  =1D0
            DGS = 0D0

          ELSE
C           *** Exit fringe

C------------ SS>0 outside
            SS=(AA-(AT-AS))*RR/QLS
            CALL DRVG(I1,QCS,SS,GS,DGS,D2GS,D3GS,D4GS,D5GS,D6GS)
C     >                                   D7GS,D8GS,D9GS,D10GS)

            DSSA = 1.D0

          ENDIF

          G = ( GE+GS-1D0 )

          IF(G.LT.0D0) CALL ENDJOB(
     >      ' SBR CHAMC :  problem  with  Gradient  ->  GE+GS-1 < ',I0)
            stop
          DG  = (DGE  * QE(1,1)*DSEA + DGS*QS(1,1)*DSSA)  *BRI

C--------- Er
          ER = - G * VN / RR
          DRA = - DG * RM / RR

C--------- Ea, from Maxwell's equation
          RRM = (RR-RM) / RM
          DAR = - DG / RM * ( 1.D0 + ( - 3.D0 + ( 5.5D0
     >       - 50.D0/6.D0 * RRM) * RRM) * RRM)  ! Check that !!
          EA = DRA - RR * DAR  ! Spherical ?? Check that !!

          E(1,1) = EA*VN*0.d0                  !AUL:110214
          E(1,2) = ER

        ENDIF
      ENDIF

c      write(*,fmt='(a,i2,1x,20e12.4)')
c     >  ' elcyl qle,qls,AA, AE, at,  rr :',
c     >  qle,qls,AA, AE, at,  rr
c      write(*,fmt='(a,i2,1x,20e12.4)')
c     > '  ex, ey, ez : ',e(1,1),e(1,2),e(1,3)


      RETURN
      END
