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
C  USA
C  -------
      SUBROUTINE BNDTH(BKG,ALCM,WDGE,FINTE,GAPE,WDGS,FINTS,GAPS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      INCLUDE "C.CHAMBR.H"     ! COMMON/CHAMBR/ LIMIT,IFORM,YLIM2,ZLIM2,SORT(MXT),FMAG,YCH,ZCH

      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL,PI,RAD,DEG,QEL,AMPROT,CM2M
      INCLUDE "C.CONST2.H"     ! COMMON/CONST2/ ZERO, UN
      INCLUDE "C.DROITE.H"     ! COMMON/DROITE/ CA(9),SA(9),CM(9),IDRT
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      INCLUDE "C.MARK.H"     ! COMMON/MARK/ KART,KALC,KERK,KUASEX
C      LOGICAL ZSYM
      INCLUDE "C.TYPFLD.H"     ! COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      INCLUDE "C.PTICUL.H"     ! COMMON/PTICUL/ AMASS,Q,G,TOO
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,HDPRF,DP,QBR,BRI
      INCLUDE "C.SPIN.H"     ! COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)
      INCLUDE "C.TRAJ.H"     ! COMMON/TRAJ/ Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT
      DIMENSION VEC(3), OT(4)

      DUM=WDGE; DUM=FINTE; DUM=GAPE
      DUM=WDGS; DUM=FINTS; DUM=GAPS

C----- OPTION COMPTAGE DES TRAJECTOIRES
C      HORS DES LIMITES   Y - Z  DU  SEPARATEUR
      IF(LIMIT .EQ. 1) THEN
        IF(ALCM.GE.0.D0) CALL CHMBR(1,IMAX)
      ENDIF

            YO = Y*.01D0
            TO = T
            ZO = Z*.01D0
            PO = P
            XO = X*.01D0
            AL = ALCM*.01D0

c       write(nres,*) ' bndth   '
c       write(nres,*) ' bndth  y,t,z,p,x,s at start : ', y,t,z,p,x,sar

C need install ~ WEDGKI(1,T,Z,P,WDGE,FINTE,GAPE)

          C1 = SIN(TO)        !*COS(PO)
          C2 = COS(TO)        !*COS(PO)
          ZP= PO        !SIN(PO)
            R = QBR / BKG  *1.D-2
            IF(C1 .EQ. 0.D0) C1=1.D-10
            EPS = ATAN2( C2 , C1 )
c       write(nres,*) 'bndth atan(c1/c2), eps, y, zp :',
c     >                                    atan(c2/c1),eps,y,zp

Compute intersection of circle w exit face (exit face has eq. a/c*x+b/c*y+1=0)
      XC =  R * SIN(TO) + XO
      YC = -R * COS(TO) + YO
      AC = CA(2) / (CM(2) * 1.D-2)
      BC = SA(2) / (CM(2) * 1.D-2)
      AA = (BC*BC + AC*AC)
      BB = 2.D0*(BC*(-XC*BC + AC*YC) + AC)
      CC = BC*(BC*(XC*XC + YC*YC - R*R) + 2.D0*YC) + 1.D0

      DSC = BB*BB -4.D0*AA*CC
       IF(DSC.LT.0.D0) STOP '  DSC IS <0 '
      XF = (-BB + SQRT(DSC))/(2.D0*AA)
      XL = XF - XO

            ARG = ACOS( C1 - XL/R )
            ARCA = ARG - EPS
            ARCL =  ARCA * R

            COA = COS(ARG)
            SIA = SIN(ARG)
            Y =  (SIA - C2)*R + YO
c            write(nres,*) '    bndth  y : ',SIA,C2,(SIA-C2),YO,y
            X = (-COA + C1)*R + XO
            Z = tan(ZP)*ARCL + ZO     ! ZP*ARCL + ZO

            T = ATAN2(COA,SIA)
            P  = ZP      !ATAN2(ZP,SQRT(COA*COA+SIA*SIA))
            Y = Y *1.D2
            Z = Z *1.D2
            X = X *1.D2
            SAR = SAR + ARCL/cos(zp)*1.D+2

c       write(nres,*) ' bndth   '
c       write(nres,*) ' bndth  y,t,z,p,x,s at end of arc : ',
c     >    y,t,z,p,x,sar

C To take particle back or forward to plane normal to X at AL
       dx = alcm - x
       dy = dx * tan(t)
       du2 = dx*dx + dy*dy
       dz = sqrt(du2) * tan(P)
       ds = sqrt(du2 + dz*dz)
       if(dx.lt.0.d0) ds = -ds
       x = x + dx
       y = y + dy
       z = z + dz
       sar = sar + ds

c       write(nres,*) ' bndth   '
c       write(nres,*)'bndth  y,t,z,p,x,s on bend exit face : ',
c     >    y,t,z,p,x,sar

C need install ~   CALL WEDGKI(2,T,Z,P,WDGS,FINTS,GAPS)

            if(amt*qt.ne.0.d0) then
              PP = QBR*CL9
              energ = SQRT(PP*PP+AMT*AMT)
              BTA = PP/energ
              V = BTA*CL
              TAR = TAR + ( ARCL+ ds*1.d-2)/V * 1.D11
            endif

      IF( KSPN .EQ. 1 ) THEN
            gamma = energ/amt
            angle = G * gamma * arca
            S1 = SF(1,IT)
            S2 = SF(2,IT)
            S3 = SF(3,IT)
            VEC(1)=0.d0
            VEC(2)=0.d0
            VEC(3)=1.d0
C             write(*,*)' bndth ',arca*deg,angle*deg
            Call SPINRO(ANGLE*DEG,VEC,OT)
            SF(1,IT) = S1*(OT(1)**2+OT(2)**2-OT(3)**2-OT(4)**2)+
     +              S2*2.d0*(OT(2)*OT(3)-OT(1)*OT(4))+
     +              S3*2.d0*(OT(2)*OT(4)+OT(1)*OT(3))
            SF(2,IT) = S2*(OT(1)**2-OT(2)**2+OT(3)**2-OT(4)**2)+
     +              S1*2.d0*(OT(2)*OT(3)+OT(1)*OT(4))+
     +              S3*2.d0*(OT(3)*OT(4)-OT(1)*OT(2))
            SF(3,IT) = S3*(OT(1)**2-OT(2)**2-OT(3)**2+OT(4)**2)+
     +              S2*2.d0*(OT(3)*OT(4)+OT(1)*OT(2))+
     +              S1*2.d0*(OT(2)*OT(4)-OT(1)*OT(3))

      ENDIF

      RETURN
      END
