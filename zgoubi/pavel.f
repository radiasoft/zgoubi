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
      SUBROUTINE PAVEL(IORDE,A1,R1,Z1,
     >                                B,DB,DDB,KERK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION B(5,3),DB(3,3),DDB(3,3,3)
      INCLUDE 'PARIZ.H'
      INCLUDE "XYZHC.H"
C      COMMON//XH(MXX),YH(MXY),ZH(IZ),HC(ID,MXX,MXY,IZ,IMAP),IXMA,JYMA,KZMA
C      COMMON/CHAVE/ B(5,3),V(5,3),E(5,3)
C      COMMON/DDBXYZ/ DB(3,3),DDB(3,3,3)
 
      COMMON/SS1/X(MXX),Y(MXY),Z(IZ),JY(25),JX(25),JZ(25),NX,NY,NZ,NN
      COMMON/SSS/EE(15,25)
 
      SAVE INDEX, NT
C      DATA INDEX/0/
      DATA IMAP / 1 /

      CALL KSMAP(
     >           IMAP) 

      IF(IORDE .EQ. 2) GOTO 3001

C EXTRAPOLATION TO 4th ORDER IN COORDINATES
C EXTRAPOLATION IS DONE FROM THE NT PLANE 
C      NT=1

      IF(INDEX .EQ. 1) GOTO 102
      CALL SERV5
      CALL SERV2
      CALL SERV4
      INDEX=1

102   CONTINUE
      X0=A1
      Y0=R1
      Z0=Z1
      CALL ICENTR(X0,Y0,DX,DY,ICX,ICY,
     >                                KERK)

      DZ=Z0-Z(NT)

      DX2 = DX**2
      DX3 = DX2 * DX
      DX4 = DX2 * DX2
      DY2 = DY**2
      DY3 = DY2 * DY
      DY4 = DY2 * DY2

      A00=0.D0
      A10=0.D0
      A01=0.D0
      A20=0.D0
      A11=0.D0
      A02=0.D0
      A30=0.D0
      A21=0.D0
      A12=0.D0
      A03=0.D0
      A40=0.D0
      A31=0.D0
      A22=0.D0
      A13=0.D0
      A04=0.D0
C        I3=3 introduced to avoid compiler complainig when IZ=1...
      I3 = 3
      DO 500 M=1,NN
        IX=ICX+JX(M)
        IY=ICY+JY(M)
        TEMP = HC(I3,IX,IY,NT,IMAP)
        A00=A00+EE(1,M)*TEMP
        A10=A10+EE(2,M)*TEMP
        A01=A01+EE(3,M)*TEMP
        A20=A20+EE(4,M)*TEMP
        A11=A11+EE(5,M)*TEMP
        A02=A02+EE(6,M)*TEMP
        A30=A30+EE(7,M)*TEMP
        A21=A21+EE(8,M)*TEMP
        A12=A12+EE(9,M)*TEMP
        A03=A03+EE(10,M)*TEMP
        A40=A40+EE(11,M)*TEMP
        A31=A31+EE(12,M)*TEMP
        A22=A22+EE(13,M)*TEMP
        A13=A13+EE(14,M)*TEMP
500     A04=A04+EE(15,M)*TEMP

      ZF0=A00+A10*DX+A01*DY+A20*DX2+A11*DX*DY+A02*DY2+
     *A30*DX3+A21*DY*DX2+A12*DX*DY2+A03*DY3+
     *A40*DX4+A31*DY*DX3+A22*DX*DX*DY2+A13*DX*DY3+A04*DY4
      ZFY=A01+2.D0*A02*DY+A11*DX+3.D0*A03*DY2+A12*2.D0*DX*DY+A21*DX2
     *+4.D0*A04*DY3+3.D0*A13*DY*DY*DX+A22*2.D0*DY*DX2+A31*DX3
      ZFX=A10+2.D0*A20*DX+A11*DY+3.D0*A30*DX2+A21*2.D0*DX*DY+A12*DY2
     *+4.D0*A40*DX3+3.D0*A31*DX*DX*DY+A22*2.D0*DX*DY2+A13*DY3
      ZFYX=A11+A21*2.D0*DX+A12*DY*2.D0
     *+3.D0*A31*DX*DX+A22*4.D0*DX*DY+A13*DY2*3D0
      ZF2X=2.D0*A20+6.D0*A30*DX+2.D0*A21*DY+12.D0*A40*DX2+
     > 6.D0*A31*DX*DY+
     *2.D0*A22*DY2
      ZF2Y=2.D0*A02+6.D0*A03*DY+2.D0*A12*DX+12.D0*A04*DY2+
     > 6.D0*A13*DX*DY+
     *2.D0*A22*DX2
      ZF2XY=2.D0*A21+6.D0*A31*DX+
     *4.D0*A22*DY
      ZF2YX=2.D0*A12+6.D0*A13*DY+
     *4.D0*A22*DX
      ZF3X=6.D0*A30+24.D0*A40*DX+6.D0*A31*DY
      ZF3XY=6.D0*A31
      ZF3Y=6.D0*A03+24.D0*A04*DY+6.D0*A13*DX
      ZF3YX=6.D0*A13
      ZF4X=24.D0*A40
      ZF4Y=24.D0*A04
      ZF2X2Y=4.D0*A22
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      A00=0.D0
      A10=0.D0
      A01=0.D0
      A20=0.D0
      A11=0.D0
      A02=0.D0
      A30=0.D0
      A21=0.D0
      A12=0.D0
      A03=0.D0
      A40=0.D0
      A31=0.D0
      A22=0.D0
      A13=0.D0
      A04=0.D0
      DO 501 M=1,NN
        IX=ICX+JX(M)
        IY=ICY+JY(M)
        TEMP = HC(1,IX,IY,NT,IMAP)
        A00=A00+EE(1,M)*TEMP
        A10=A10+EE(2,M)*TEMP
        A01=A01+EE(3,M)*TEMP
        A20=A20+EE(4,M)*TEMP
        A11=A11+EE(5,M)*TEMP
        A02=A02+EE(6,M)*TEMP
        A30=A30+EE(7,M)*TEMP
        A21=A21+EE(8,M)*TEMP
        A12=A12+EE(9,M)*TEMP
        A03=A03+EE(10,M)*TEMP
        A40=A40+EE(11,M)*TEMP
        A31=A31+EE(12,M)*TEMP
        A22=A22+EE(13,M)*TEMP
        A13=A13+EE(14,M)*TEMP
501     A04=A04+EE(15,M)*TEMP

      XF0=A00+A10*DX+A01*DY+A20*DX2+A11*DX*DY+A02*DY2+
     *A30*DX3+A21*DY*DX2+A12*DX*DY2+A03*DY3+
     *A40*DX4+A31*DY*DX3+A22*DX*DX*DY2+A13*DX*DY3+A04*DY4
      XFY=A01+A11*DX+A02*DY*2.D0+
     *A21*DX2+A12*DX*DY*2.D0+3.D0*A03*DY2+
     *A31*DX3+A22*DX*DX*DY*2.D0+3.D0*A13*DX*DY2+4.D0*A04*DY3
      XFX=A10+2.D0*A20*DX+A11*DY+3.D0*A30*DX2+A21*2.D0*DX*DY+A12*DY2
     *+4.D0*A40*DX3+3.D0*A31*DX*DX*DY+A22*2.D0*DX*DY2+A13*DY3
      XF2X=2.D0*A20+6.D0*A30*DX+2.D0*A21*DY+12.D0*A40*DX2+
     > 6.D0*A31*DX*DY+
     *2.D0*A22*DY2
      XF2Y=2.D0*A02+6.D0*A03*DY+2.D0*A12*DX+12.D0*A04*DY2+
     > 6.D0*A13*DX*DY+
     *2.D0*A22*DX2
      XF3Y=6.D0*A03+24.D0*A04*DY+6.D0*A13*DX
      XF2YX=2.D0*A12+6.D0*A13*DY+
     *4.D0*A22*DX
      XF3X=6.D0*A30+24.D0*A40*DX+6.D0*A31*DY
      XF4X=24.D0*A40
      XF4Y=24.D0*A04
      XF2X2Y=4.D0*A22
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      A00=0.D0
      A10=0.D0
      A01=0.D0
      A20=0.D0
      A11=0.D0
      A02=0.D0
      A30=0.D0
      A21=0.D0
      A12=0.D0
      A03=0.D0
      A40=0.D0
      A31=0.D0
      A22=0.D0
      A13=0.D0
      A04=0.D0
C        I2=2 introduced to avoid compiler complainig when IZ=1...
      I2 = 2
      DO 502 M=1,NN
        IX=ICX+JX(M)
        IY=ICY+JY(M)
        TEMP = HC(I2,IX,IY,NT,IMAP)
        A00=A00+EE(1,M)*TEMP
        A10=A10+EE(2,M)*TEMP
        A01=A01+EE(3,M)*TEMP
        A20=A20+EE(4,M)*TEMP
        A11=A11+EE(5,M)*TEMP
        A02=A02+EE(6,M)*TEMP
        A30=A30+EE(7,M)*TEMP
        A21=A21+EE(8,M)*TEMP
        A12=A12+EE(9,M)*TEMP
        A03=A03+EE(10,M)*TEMP
        A40=A40+EE(11,M)*TEMP
        A31=A31+EE(12,M)*TEMP
        A22=A22+EE(13,M)*TEMP
        A13=A13+EE(14,M)*TEMP
502     A04=A04+EE(15,M)*TEMP

      YF0=A00+A10*DX+A01*DY+A20*DX2+A11*DX*DY+A02*DY2+
     *A30*DX3+A21*DY*DX2+A12*DX*DY2+A03*DY3+
     *A40*DX4+A31*DY*DX3+A22*DX*DX*DY2+A13*DX*DY3+A04*DY4
      YFY=A01+2.D0*A02*DY+A11*DX+3.D0*A03*DY2+A12*2.D0*DX*DY+A21*DX2
     *+4.D0*A04*DY3+3.D0*A13*DY*DY*DX+A22*2.D0*DY*DX2+A31*DX3
      YF2X=2.D0*A20+6.D0*A30*DX+2.D0*A21*DY+12.D0*A40*DX2+
     > 6.D0*A31*DX*DY+
     *2.D0*A22*DY2
      YF3X=6.D0*A30+24.D0*A40*DX+6.D0*A31*DY
      YF2Y=2.D0*A02+6.D0*A03*DY+2.D0*A12*DX+12.D0*A04*DY2+
     > 6.D0*A13*DX*DY+
     *2.D0*A22*DX2
      YF2XY=2.D0*A21+6.D0*A31*DX+
     *4.D0*A22*DY
      YF3Y=6.D0*A03+24.D0*A04*DY+6.D0*A13*DX
      YF4X=24.D0*A40
      YF4Y=24.D0*A04
      YF2X2Y=4.D0*A22
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCC              F                  CCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      B(1,1)=XF0+DZ*(ZFX+DZ*(-(XF2X+XF2Y)/2.D0
     *+DZ*((-ZF3X-ZF2YX)/6D0
     *+DZ*(XF4X+XF4Y+XF2X2Y*2.D0)/24.D0)))
      B(1,2)=YF0+DZ*(ZFY+DZ*(-(YF2X+YF2Y)/2.D0
     *+DZ*((-ZF3Y-ZF2XY)/6D0
     *+DZ*(YF4X+YF4Y+YF2X2Y*2.D0)/24.D0)))
      B(1,3)=ZF0+DZ*(-XFX-YFY+DZ*(-(ZF2X+ZF2Y)/2.D0
     *+DZ*((XF3X+YF3Y+XF2YX+YF2XY)/6D0
     *+DZ*(ZF4X+ZF4Y+ZF2X2Y*2.D0)/24.D0)))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCC             D2F / D  D        CCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DXX=XFX+DZ*(ZF2X+DZ*((-XF3X-YF2XY)/2.D0+DZ*(-ZF4X-ZF2X2Y)/6.D0))
      DYY=YFY+DZ*(ZF2Y+DZ*((-YF3Y-XF2YX)/2.D0+DZ*(-ZF4Y-ZF2X2Y)/6.D0))
      DZZ=-DXX-DYY
      DXY=XFY+DZ*(ZFYX+DZ*((-XF3Y-YF3X)/2.D0+DZ*(-ZF3XY-ZF3YX)/6.D0))
      DXZ=ZFX+DZ*(-XF2X-XF2Y+DZ*((-ZF3X-ZF2YX)/2.D0+
     *DZ*(XF4X+XF4Y+2.D0*XF2X2Y)/6.D0))
      DYZ=ZFY+DZ*(-YF2Y-YF2X+DZ*((-ZF3Y-ZF2XY)/2.D0+
     *DZ*(YF4Y+YF4X+2.D0*YF2X2Y)/6.D0))
      DB(1,1)=DXX
      DB(1,2)=DXY
      DB(1,3)=DXZ
      DB(2,1)=DXY
      DB(2,2)=DYY
      DB(2,3)=DYZ
      DB(3,1)=DXZ
      DB(3,2)=DYZ
      DB(3,3)=DZZ
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCC             D3F   / D  D  D        CCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DXXX=XF2X+DZ*(ZF3X+DZ*(-XF4X-XF2X2Y)/2.D0)
      DYYY=YF2Y+DZ*(ZF3Y+DZ*(-YF4Y-YF2X2Y)/2.D0)
      DXXY=YF2X+DZ*(ZF2XY+DZ*(-YF4X-YF2X2Y)/2.D0)
      DYYX=XF2Y+DZ*(ZF2YX+DZ*(-XF4Y-XF2X2Y)/2.D0)
      DXXZ=ZF2X+DZ*(-XF3X-YF2XY+DZ*(-ZF4X-ZF2X2Y)/2.D0)
      DYYZ=ZF2Y+DZ*(-YF3Y-XF2YX+DZ*(-ZF4Y-ZF2X2Y)/2.D0)
      DXYZ=ZFYX+DZ*(-YF3X-XF3Y+DZ*(-ZF3XY-ZF3YX)/2.D0)
      DZZZ=-DXXZ-DYYZ
      DZZX=-DXXX-DYYX
      DZZY=-DYYY-DXXY
      DDB(1,1,1)=DXXX
      DDB(1,1,2)=DXXY
      DDB(1,1,3)=DXXZ
      DDB(1,2,1)=DXXY
      DDB(1,2,2)=DYYX
      DDB(1,2,3)=DXYZ
      DDB(1,3,1)=DXXZ
      DDB(1,3,2)=DXYZ
      DDB(1,3,3)=DZZX
      DDB(2,1,1)=DXXY
      DDB(2,1,2)=DYYX
      DDB(2,1,3)=DXYZ
      DDB(2,2,1)=DYYX
      DDB(2,2,2)=DYYY
      DDB(2,2,3)=DYYZ
      DDB(2,3,1)=DXYZ
      DDB(2,3,2)=DYYZ
      DDB(2,3,3)=DZZY
      DDB(3,1,1)=DXXZ
      DDB(3,1,2)=DXYZ
      DDB(3,1,3)=DZZX
      DDB(3,2,1)=DXYZ
      DDB(3,2,2)=DYYZ
      DDB(3,2,3)=DZZY
      DDB(3,3,1)=DZZX
      DDB(3,3,2)=DZZY
      DDB(3,3,3)=DZZZ
      RETURN


 3001 CONTINUE

CC EXTRAPOLATION IS DONE FROM THE NT PLANE
C      NT=1

C          write(*,*) ' pavel, index, nt : ',index,nt

      IF(INDEX .EQ. 1) GOTO 302

      CALL SERV5
      CALL SERV32
      CALL SERV34
      INDEX=1
302   CONTINUE
      X0=A1
      Y0=R1
      Z0=Z1
      CALL JCENTR(X0,Y0,DX,DY,ICX,ICY,
     >                                KERK)
      DZ=Z0-Z(NT)

C          write(*,*) ' pavel, dz : ',Z0,Z(NT),NT,dz
C          write(*,*) ' pavel, x0,y0,z0  : ',x0,y0,z0

      DX2 = DX**2
      DY2 = DY**2

      A00=0.D0
      A10=0.D0
      A01=0.D0
      A20=0.D0
      A11=0.D0
      A02=0.D0
C        I3=3 introduced to avoid compiler complainig when IZ=1...
      I3 = 3
      DO 600 M=1,NN
        IX=ICX+JX(M)
        IY=ICY+JY(M)
        TEMP = HC(I3,IX,IY,NT,IMAP)
        A00=A00+EE(1,M)*TEMP
        A10=A10+EE(2,M)*TEMP
        A01=A01+EE(3,M)*TEMP
        A20=A20+EE(4,M)*TEMP
        A11=A11+EE(5,M)*TEMP
        A02=A02+EE(6,M)*TEMP
600   CONTINUE
      ZF0=A00+A10*DX+A01*DY+A20*DX2+A11*DX*DY+A02*DY2
      ZFY=A01+2.D0*A02*DY+A11*DX
      ZFX=A10+2.D0*A20*DX+A11*DY
      ZFYX=A11
      ZF2X=2.D0*A20
      ZF2Y=2.D0*A02
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      A00=0.D0
      A10=0.D0
      A01=0.D0
      A20=0.D0
      A11=0.D0
      A02=0.D0
      DO 601 M=1,NN
        IX=ICX+JX(M)
        IY=ICY+JY(M)
        TEMP = HC(1,IX,IY,NT,IMAP)
        A00=A00+EE(1,M)*TEMP
        A10=A10+EE(2,M)*TEMP
        A01=A01+EE(3,M)*TEMP
        A20=A20+EE(4,M)*TEMP
        A11=A11+EE(5,M)*TEMP
        A02=A02+EE(6,M)*TEMP
601   CONTINUE
      XF0=A00+A10*DX+A01*DY+A20*DX2+A11*DX*DY+A02*DY2
      XFY=A01+A11*DX+A02*DY*2.D0
      XFX=A10+2.D0*A20*DX+A11*DY
      XF2X=2.D0*A20
      XF2Y=2.D0*A02
C      PRINT 100,a1,r1,DX,DY, temp,XF0,A00
C100   FORMAT(1X,1p,7(G12.4,1X))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      A00=0.D0
      A10=0.D0
      A01=0.D0
      A20=0.D0
      A11=0.D0
      A02=0.D0
C        I2=2 introduced to avoid compiler complainig when IZ=1...
      I2 = 2
      DO 602 M=1,NN
        IX=ICX+JX(M)
        IY=ICY+JY(M)
        TEMP = HC(I2,IX,IY,NT,IMAP)
        A00=A00+EE(1,M)*TEMP
        A10=A10+EE(2,M)*TEMP
        A01=A01+EE(3,M)*TEMP
        A20=A20+EE(4,M)*TEMP
        A11=A11+EE(5,M)*TEMP
        A02=A02+EE(6,M)*TEMP
602     CONTINUE
      YF0=A00+A10*DX+A01*DY+A20*DX2+A11*DX*DY+A02*DY2
      YFY=A01+2.D0*A02*DY+A11*DX
      YF2X=2.D0*A20
      YF2Y=2.D0*A02
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCC              F                  CCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      B(1,1)=XF0+DZ*(ZFX+DZ*(-(XF2X+XF2Y)/2.D0))
      B(1,2)=YF0+DZ*(ZFY+DZ*(-(YF2X+YF2Y)/2.D0))
      B(1,3)=ZF0+DZ*(-XFX-YFY+DZ*(-(ZF2X+ZF2Y)/2.D0))
C      PRINT 101, '  PAVEL, B11  ',B(1,1)
C101   FORMAT(1X,A,(G12.4,1X))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCC             D2F / D  D        CCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DXX=XFX+DZ*ZF2X
      DYY=YFY+DZ*ZF2Y
      DZZ=-DXX-DYY
      DXY=XFY+DZ*ZFYX
      DXZ=ZFX+DZ*(-XF2X-XF2Y)
      DYZ=ZFY+DZ*(-YF2Y-YF2X)
      DB(1,1)=DXX
      DB(1,2)=DXY
      DB(1,3)=DXZ
      DB(2,1)=DXY
      DB(2,2)=DYY
      DB(2,3)=DYZ
      DB(3,1)=DXZ
      DB(3,2)=DYZ
      DB(3,3)=DZZ
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCC             D3F   / D  D  D        CCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DXXX=XF2X
      DYYY=YF2Y
      DXXY=YF2X
      DYYX=XF2Y
      DXXZ=ZF2X
      DYYZ=ZF2Y
      DXYZ=ZFYX
      DZZZ=-DXXZ-DYYZ
      DZZX=-DXXX-DYYX
      DZZY=-DYYY-DXXY
      DDB(1,1,1)=DXXX
      DDB(1,1,2)=DXXY
      DDB(1,1,3)=DXXZ
      DDB(1,2,1)=DXXY
      DDB(1,2,2)=DYYX
      DDB(1,2,3)=DXYZ
      DDB(1,3,1)=DXXZ
      DDB(1,3,2)=DXYZ
      DDB(1,3,3)=DZZX
      DDB(2,1,1)=DXXY
      DDB(2,1,2)=DYYX
      DDB(2,1,3)=DXYZ
      DDB(2,2,1)=DYYX
      DDB(2,2,2)=DYYY
      DDB(2,2,3)=DYYZ
      DDB(2,3,1)=DXYZ
      DDB(2,3,2)=DYYZ
      DDB(2,3,3)=DZZY
      DDB(3,1,1)=DXXZ
      DDB(3,1,2)=DXYZ
      DDB(3,1,3)=DZZX
      DDB(3,2,1)=DXYZ
      DDB(3,2,2)=DYYZ
      DDB(3,2,3)=DZZY
      DDB(3,3,1)=DZZX
      DDB(3,3,2)=DZZY
      DDB(3,3,3)=DZZZ
      RETURN

      ENTRY PAVELW(INDEXI,NTI)
      INDEX = INDEXI
C EXTRAPOLATION IS DONE FROM THE NT PLANE
      NT = NTI
      RETURN
      END
