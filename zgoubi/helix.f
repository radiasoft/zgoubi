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
C  François Méot <meot@lpsc.in2p3.fr>
C  Service Accélerateurs
C  LPSC Grenoble
C  53 Avenue des Martyrs
C  38026 Grenoble Cedex
C  France
      SUBROUTINE HELIX(SCAL,
     >                      XL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON/AIM/ BO,RO,FG,GF,XI,XF,EN,EB1,EB2,EG1,EG2
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QEL,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      LOGICAL ZSYM
      COMMON/OPTION/ KFLD,MG,LC,ML,ZSYM
      COMMON/ORDRES/ KORD,IRD,IDS,IDB,IDE,IDZ
      COMMON/RIGID/ BORO,DPREF,DPPP,QBR,BRI
 
      SAVE NN, RESOL
      CHARACTER TYPCAL(2)*24
      SAVE TYPCAL

      DIMENSION B(5,3),DB(3,3),DDB(3,3,3)
      DIMENSION FTAB3(3,3,3,3)
 
      LOGICAL OKX, LTXI, GTXF
      SAVE PITCHB, ANG0, AK

      DATA OKX, LTXI, GTXF / .TRUE., .FALSE.   , .FALSE.   /

      DATA TYPCAL / ' analytic', ' numerical interpolation'/

          XL =A(NOEL,10)
          PITCHB =A(NOEL,11)
          BO =A(NOEL,12)*SCAL
          ANG0 = A(NOEL,13)
 
          AK = 2.D0 * PI / PITCHB
          IF(NRES.GT.0) THEN
            WRITE(NRES,100) 'Helical  magnet',XL,PITCHB,BO,ANG0
 100        FORMAT(/,1P, 5X,' -----  ',A15,'  : '
     >      ,/,15X,' Length  of  element     : ',G12.4,'  cm'
     >      ,/,15X,' Twist  pitch            : ',G12.4,'  cm'
     >      ,/,15X,' Field                   : ',G12.4,'  kG'
     >      ,/,15X,' Initial  field  angle   : ',G12.4,'  rad')

          ENDIF

      IRD2 = NINT(A(NOEL,20))
      IF    (IRD2.EQ.3) THEN
C 3D interpolation from 3D 3*3*3 points flying grid with mesh size integration step/Resol. 
        RESOL=A(NOEL,21)
        IRD = IRD2
        IRD2 = 1
        NN=3
      ELSEIF(IRD2.EQ.0) THEN
C    analytic. A(NOEL,ND) has the form 0.j, j=2 or 4. 
        IDB = NINT(10*A(NOEL,ND))
        IF(IDB.NE.4) IDB=2
      ELSE
        CALL ENDJOB('ERROR - SBR HELIX, wrong value IRD',-99)
      ENDIF
      CALL CHAMC6(IRD2)
    
          XI = 0.D0
          XLIM = XL
          XF = XLIM
          XS = XLIM
          IF(BO .EQ. 0.D0) KFLD=0

      IF(NRES.GT.0) THEN
        WRITE(NRES,FMT='(/,5X,'' Field & derivatives calculation''
     >  ,'' method :'',A)') TYPCAL(IRD2+1)
        IF(IRD .EQ. 3) THEN
          WRITE(NRES,FMT='(20X
     >        ,''3*3*3-point  interpolation, 2nd degree polynomial'')')
          WRITE(NRES,121) RESOL
 121      FORMAT(1P,20X
     >     ,'Size of flying mesh is :   integration step /',G12.4)
        ELSEIF(IRD2.EQ.0) THEN
          WRITE(NRES,FMT='(20X,
     >    ''Derivatives computed to order '',I1)') IDB
        ENDIF
      ENDIF

      RETURN

C----------------------------------------------------------------
C  Compute helix field  and derivatives from flying 3D field-mesh
      ENTRY HELIXF(X,Y,Z,
     >                   XX,YY,ZZ,DX,DY,DZ,FTAB3)

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
            BZ=  BO*((1.D0 + AK28*Y23Z2)*CKX - AK24*YZ*SKX)

c             write(88,*) x,r,bx-bxo,by-byo,bz-bzo,' x,r, ...'

            FTAB3(1,IX,JY,KZ) = BX
            FTAB3(2,IX,JY,KZ) = BY
            FTAB3(3,IX,JY,KZ) = BZ
C           write(90,*) ' helix b1L, L F :  1',ix,jy,kz,bx,xi,yj,zk
C           write(90,*) ' helix b1L, L F :  2',ix,jy,kz,by,xi,yj,zk
C           write(90,*) ' helix b1L, L F :  3',ix,jy,kz,bz,xi,yj,zk
          
 1    CONTINUE
C--------- end loop on   NN x DX,  NN x DY, NN x DZ

c      DO  lL = 1,3
c        DO  K=1,3
c          Kk=K-2
c          DO  J=1,3
c            Jj=2-J
c            DO  I=1,3
c              Ii=I-2
c                 write(89,*) FTAB3(lL,Ii,Jj,Kk), lL,Ii,Jj,Kk
c        enddo
c        enddo
c        enddo
c        enddo

      RETURN

C----------------------------
C 3D analytic model of field
      ENTRY HELIXA(X,Y,Z,
     >                   B,DB,DDB)

C      CALL ENDJOB
C     >    ('SBR HELIX : analytical  model not implemented. ',-99)

      BN = BO*BRI
      AK2 = AK*AK
      AK24 = AK2/4.D0
      AK28 = AK24/2.D0
      XK = X * AK
      CKX = COS(XK)
      SKX = SIN(XK)
      Y2 = Y*Y
      Z2 = Z*Z
      YZ = Y*Z
      Y2Z2 = Y2+Z2
      Y23Z2 = Y2Z2 + 2.D0*Z2
      Z23Y2 = 2.D0*Y2 + Y2Z2

C--------- Components Bx, By, Bz of field
          B(1,1)= -BN*AK*(Y*CKX+Z*SKX)*(1.D0 + AK28*Y2Z2)
          B(1,2)= -BN*((1.D0 + AK28*Z23Y2)*SKX - AK24*YZ*CKX)
          B(1,3)=  BN*((1.D0 + AK28*Y23Z2)*CKX - AK24*YZ*SKX)

C Just to calm down the compiler
          DB(1,1) = 0.D0
          DDB(1,1,1) = 0.D0

      RETURN
      END
