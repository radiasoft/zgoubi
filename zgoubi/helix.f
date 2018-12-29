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
      SUBROUTINE HELIX(SCAL,
     >                      XL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE "C.AIM.H"     ! COMMON/AIM/ BO,RO,FG,GF,XI,XF,EN,EB1,EB2,EG1,EG2
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.CONST_2.H"     ! COMMON/CONST/ CL9,CL,PI,RAD,DEG,QEL,AMPROT,CM2M
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "C.INTEG.H"     ! COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
C      LOGICAL ZSYM
      INCLUDE "C.TYPFLD.H"     ! COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      INCLUDE "C.ORDRES.H"     ! COMMON/ORDRES/ KORD,IRD,IDS,IDB,IDE,IDZ
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI

      SAVE NN, RESOL
      CHARACTER(24) TYPCAL(2)
      SAVE TYPCAL

      DIMENSION B(5,3),DB(3,3),DDB(3,3,3)
      DIMENSION FTAB3(3,3,3,3)

      LOGICAL OKX, LTXI, GTXF
      SAVE PITCHB, XROT, AK

      DATA OKX, LTXI, GTXF / .TRUE., .FALSE.   , .FALSE.   /

      DATA TYPCAL / ' analytic', ' numerical interpolation'/

          XL =A(NOEL,10)
          PITCHB =A(NOEL,11)
          BO =A(NOEL,12)*SCAL
          XROT = A(NOEL,13) * RAD

          AK = 2.D0 * PI / PITCHB
          IF(NRES.GT.0) THEN
            WRITE(NRES,100)'Helical magnet',XL,PITCHB,BO,XROT/RAD
 100        FORMAT(/,1P, 5X,' -----  ',A15,'  : '
     >      ,/,15X,' Length  of  element     : ',G12.4,'  cm'
     >      ,/,15X,' Twist  pitch            : ',G12.4,'  cm'
     >      ,/,15X,' Field                   : ',G12.4,'  kG'
     >      ,/,15X,' Initial  field  angle   : ',G12.4,'  deg')
          ENDIF

      IRDA = NINT(A(NOEL,20))
      IF    (IRDA.EQ.0) THEN
C    analytic. A(NOEL,21) = 2 or 4 = order of derivatives (not necessarily operational...)
        IDB = NINT(A(NOEL,21))
        IF(IDB.NE.4) IDB=2
        KAN = 0
        IRD = 3
      ELSEIF(IRDA.EQ.3) THEN
C 3D interpolation from 3D 3*3*3 points flying grid with mesh size integration step/Resol.
        RESOL=A(NOEL,21)
        KAN = 1
        NN=3
      ELSE
        CALL ENDJOB('ERROR - SBR HELIX, wrong value IRD',-99)
      ENDIF
      CALL CHAMC6(KAN)

          XI = 0.D0
          XLIM = XL
          XF = XLIM
          XS = XLIM
          IF(BO .EQ. 0.D0) KFLD=0

      IF(NRES.GT.0) THEN
        IF(IRDA .EQ. 0) THEN
          WRITE(NRES,FMT='(/,5X,'' Field & derivatives calculation''
     >    ,'' method :'',A)') TYPCAL(1)
          WRITE(NRES,FMT='(20X,
     >    ''Derivatives computed to order '',I1)') IDB
        ELSEIF(IRDA.EQ.3) THEN
          WRITE(NRES,FMT='(/,5X,'' Field & derivatives calculation''
     >    ,'' method :'',A)') TYPCAL(2)
          WRITE(NRES,FMT='(20X
     >        ,''3*3*3-point  interpolation, 2nd degree polynomial'')')
          WRITE(NRES,121) RESOL
 121      FORMAT(1P,20X
     >     ,'Size of flying mesh is :   integration step /',G12.4)
        ENDIF
      ENDIF

      RETURN

C-----------------------------------------------------------
C  Compute helix field derivatives from flying 3D field-mesh
      ENTRY HELIXF(X,Y,Z,
     >                   XX,YY,ZZ,DX,DY,DZ,FTAB3,XROTO)

      XROTO = XROT
      CALL INTEG5(
     >            STEP)

      CALL RAZ(FTAB3,3*3*3*3)
      DX = STEP/RESOL
      DY = DX
      DZ = DX

            XX = ZERO
            YY = ZERO
            ZZ = ZERO

      KMAG = 1

      NN2 = INT(NN/2)
      DO  1  KZ = 1,NN
        ZK = Z - DZ * DBLE(NN-KZ-NN2)

      DO  1  JY = 1,NN
        YJ = Y + DY * DBLE(NN-JY-NN2)

C        LTXI = X-DX .LT. XI
C        GTXF = X+DX .GT. XF
C        OKX = .NOT.(LTXI .OR. GTXF)
CC        write(*,*) ' dipi  okX, ltxi, gtxf ', okX, ltxi, gtxf
C        IF(OKX) THEN
C          XX = 0.D0
C        ELSEIF(LTXI) THEN
C          XX=-DX
C          DX = 2.D0 * DX
C        ELSEIF(GTXF) THEN
C          XX= DX
C          DX = 2.D0 * DX
C        ENDIF

        DO  1  IX = 1,NN
          XI = X - DX * DBLE(NN-IX-NN2)

            AK2 = AK*AK
            AK24 = AK2/4.D0
            AK28 = AK24/2.D0
            XK = XI * AK
            CKX = COS(XK)
            SKX = SIN(XK)
            Y2 = YJ*YJ
            Z2 = ZK*ZK
            YZ = YJ*ZK
            Y2Z2 = Y2+Z2
            Y23Z2 = Y2Z2 + 2.D0*Z2
            Z23Y2 = 2.D0*Y2 + Y2Z2

C--------- Components Bx, By, Bz of field
c        r = sqrt(Y2+Z2)
c        r2 = r*r
c        akr = ak*r
c        fp = BI1p(akr)
c        fx = x1BI1(akr)
c        fac =  fp - fx
c          if(r.eq.0.d0) r=1.d-6
c          byo = 2.d0/r2*( fac*y*z*CKX - (y2*fp + z2* fx)*SKX ) * BO
c          bzo = 2.d0/r2*( (z2*fp + y2* fx)*CKX - fac*y*z*SKX ) * BO
c          bxo = -2.d0*ak * fx * (y*CKX + z*SKX)                * BO

            BX= -BO*AK*(YJ*CKX+ZK*SKX)*(1.D0 + AK28*Y2Z2)
            BY= -BO*((1.D0 + AK28*Z23Y2)*SKX - AK24*YZ*CKX)
            BZ=  BO*((1.D0 + AK28*Y23Z2)*CXX - AK24*YZ*SKX)

            FTAB3(1,IX,JY,KZ) = BX
            FTAB3(2,IX,JY,KZ) = BY
            FTAB3(3,IX,JY,KZ) = BZ

 1    CONTINUE
C--------- end loop on   NN x DX,  NN x DY, NN x DZ

      RETURN

C---------------------------
C 3D analytic model of field
      ENTRY HELIXA(X,Y,Z,
     >                   B,DB,DDB,XROTO)

C      CALL ENDJOB
C     >    ('SBR HELIX : analytical  model not implemented. ',-99)

      XROTO = XROT

      BN = BO*BRI
      AK2 = AK*AK
      AK24 = AK2/4.D0
      AK28 = AK24/2.D0
      XK = X * AK
      CX = COS(XK)
      SX = SIN(XK)
      Y2 = Y*Y
      Z2 = Z*Z
      YZ = Y*Z
      Y2Z2 = Y2+Z2
      Y23Z2 = Y2Z2 + 2.D0*Z2
      Z23Y2 = 2.D0*Y2 + Y2Z2
      ay2z2 = AK28*Y2Z2
      AZ23Y2 = AK28*Z23Y2
      AY23Z2 = AK28*Y23Z2
C--------- Components Bx, By, Bz of field
          B(1,1)= -BN*AK*( Y*CX + Z*SX ) * ( 1.D0 + AY2Z2 )
          B(1,2)= -BN*((1.D0 + AZ23Y2)*SX - AK24*YZ*CX)
          B(1,3)=  BN*((1.D0 + AY23Z2)*CX - AK24*YZ*SX)

          cp = -ak * sx
          sp =  ak * cx
          cpp = -ak * sp
          spp =  ak * cp
          ay = AK28*Y
          az = AK28*Z
          uay2z2 = 1.D0 + ay2z2
          yczs =  (y*cx + z*sx)
          ycpzsp =  (y*cp + z*sp)
          day = 2.D0 * ay
          daz = 2.D0 * az
          u2ay = 1.D0 + day
          u2az = 1.D0 + daz

C         ... dBx/dX
c          uuu = -(ycpzsp * uay2z2                )* BN*AK
          DB(1,1) = ak * (y*cp+z*sp) * uay2z2  * (-bn)
c             write(*,*) ' uuu : ',uuu-DB(1,1),DB(1,1)
cC         ... dBx/dY = dBy/dX
c          uuu = -(cx * uay2z2 + yczs * day       )* BN*AK
          DB(2,1) = ak*(cx*(1.d0+ak28*Z23Y2) + ak24*yz*sx )  *(-bn)
c             write(*,*) ' uuu : ',uuu-DB(2,1),DB(2,1)
cC         ... dBy/dY
c          uuu = -(3.d0*day *sx - 2.d0*az * cx              )* BN  ********
          DB(2,2) = ak24 * (3.d0 *y *sx - z *cx )  *(-bn)
c             write(*,*) ' uuu : ',uuu-DB(2,2),DB(2,2)
cC         .. dBx/dZ = dBz/dX
c          uuu = -(sx * uay2z2 + yczs * daz       )* BN*AK
          DB(3,1) = ak * ( sx*(1.d0+ak28*Y23Z2) + ak24*yz*cx )  *(-bn)
c             write(*,*) ' uuu : ',uuu-DB(3,1),DB(3,1)
cC         .. dBy/dZ = dBZ/dY
c          uuu = -(2.D0 * az * sx - 2.*ay * cx       )* BN   ******
          DB(3,2) = ak24 * (z *sx - y *cx )  *(-bn)
c             write(*,*) ' uuu : ',uuu-DB(3,2),DB(3,2)
cC         ... d2Bx/dX2
c          uuu = -((y*cpp + z*spp) * uay2z2  )* BN*AK
          DDB(1,1,1) = ak * (y*cpp+z*spp) * uay2z2  *(-bn)
c             write(*,*) ' uuu : ',uuu-dDB(1,1,1),DdB(1,1,1)
cC         ... d2Bx/dXdY = d2By/dX2
c          uuu = -(cp * uay2z2 + ycpzsp*day   )* BN*AK
          DDB(2,1,1) = ak*(cp*(1.d0+ak28*Z23Y2) + ak24*yz*sp )  *(-bn)
c             write(*,*) ' uuu : ',uuu-dDB(2,1,1),DdB(2,1,1)
cC         ... d2Bx/dY2
c          uuu = - (2.D0 *cx* day  + yczs*2.D0 *AK28) * BN*AK  *******
          DDB(2,2,1) = ak * ak24 * (3.d0 *y *cx + z *sx)  *(-bn)
c             write(*,*) ' uuu : ',uuu-dDB(2,2,1),DdB(2,2,1)
cC         ... d2By/dY2
c          uuu = -(2.D0 * AK28 * sx     *3.d0       )* BN    ********
          DDB(2,2,2) = 3.d0 * ak24 * sx  *(-bn)
c             write(*,*) ' uuu : ',uuu-dDB(2,2,2),DdB(2,2,2)
cC         .. d2Bx/dXdZ = d2Bz/dX2
c          uuu = -(sp * uay2z2 + ycpzsp * daz )* BN*AK
          DDB(3,1,1) = ak * (sp*(1.d0+ak28*Y23Z2) + ak24*yz*cp )  *(-bn)
c             write(*,*) ' uuu : ',uuu-dDB(3,1,1),DdB(3,1,1)
cC         .. d2By/dXdZ = d2Bz/dXdY = d2Bx/dYdZ
c          uuu = -(2.D0* az * sp - 2.d0*ay * cp )* BN  *********
          DDB(3,2,1) = ak24 * (z *sp - y *cp )  *(-bn)
c             write(*,*) ' uuu : ',uuu-dDB(3,2,1),DdB(3,2,1)
cC         .. d2By/dYdZ = d2Bz/dY2
c          uuu = -(- AK28 * cx *2.D0           )* BN
          DDB(3,2,2) = -ak24 * cx  *(-bn)
c             write(*,*) ' uuu : ',uuu-dDB(3,2,2),DdB(3,2,2)

      RETURN
      END
