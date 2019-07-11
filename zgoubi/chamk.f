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
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.CHAVE_2.H"     ! COMMON/CHAVE/ B(5,3),V(5,3),E(5,3)
      INCLUDE "C.DDBXYZ.H"     ! COMMON/DDBXYZ/ DB(3,3),DDB(3,3,3)
      INCLUDE "C.D3B_2.H"     ! COMMON/D3BXYZ/ D3BX(3,3,3), D3BY(3,3,3), D3BZ(3,3,3)
      INCLUDE "C.D4B.H"     ! COMMON/D4BXYZ/ D4BX(3,3,3,3) ,D4BY(3,3,3,3) ,D4BZ(3,3,3,3)

      INCLUDE "C.ORDRES.H"     ! COMMON/ORDRES/ KORD,IRD,IDS,IDB,IDE,IDZ

      DIMENSION BT(5,15)
      SAVE BT

      PARAMETER (I3=3)
      PARAMETER (MXC = 4)
      DIMENSION NMPT(MXC), NMPTI(MXC)

      PARAMETER (MXMAP=4)
      PARAMETER (MXAA2=MAX(24+MXC-1,20+10*MXMAP+9))
      DIMENSION UU(MXAA2), UUI(MXAA2)

      SAVE IOP, UU

      DATA NFIC / 1 /
      DATA NMPT / MXC * 1 /

      IF    (NFIC .EQ. 1) THEN
        CALL KSMAP(
     >           IMAP)
        CALL CHAMK1(A1,R1,Z1,IMAP,*999)

c         write(88,fmt='(1p,4e12.4,1x,2i4,1x,e12.4,a)')
c     >   A1,R1,A1N,R1N,i,imap,B(1,3),' A1,R1,A1N,R1N,imag, BT '

      ELSEIF(NFIC .GE. 2) THEN
C Case MOD=16. IMAP has been increazed by NFIC units in tosca, each field map stored in HC
        I = 1
        ID = 2  ! This is to be settled
        IRD=2

C Hyp UU : X_O2, Y_O2, theta, dY (magnet is first translated, then rotated, then Y-shifted)
        XO2 = UU(20+10*I+1)
        YO2 = UU(20+10*I+2)
        TTA = UU(20+10*I+3)
        DY =  UU(20+10*I+4)
C         if(xo2.ne.0.d0)   write(*,*) ' chamk ',i,xo2,yo2,tta,dy

        CT = COS(TTA)
        ST = SIN(TTA)
        A1N =  (A1-XO2)*CT + (R1-YO2)*ST
        R1N = -(A1-XO2)*ST + (R1-YO2)*CT - DY
        CALL CHAMK1(A1N,R1N,Z1,NMPT(I),*999)
        IOP = 1
        CALL ADPOL(ID,IOP,B,DB,DDB,D3BX,D3BY,D4BX,D4BY,BT)

c         write(88,fmt='(1p,4e12.4,1x,2i4,1x,e12.4,a)')
c     >   A1,R1,A1N,R1N,i,NMPT(I),B(1,3),' A1,R1,A1N,R1N,imag, BT '

        DO I = 2, NFIC-1
c               write(*,*) 'chamk NMPT(I) ',i,NMPT(I)
          XO2 = UU(20+10*I+1)
          YO2 = UU(20+10*I+2)
          TTA = UU(20+10*I+3)
          DY =  UU(20+10*I+4)
          CT = COS(TTA)
          ST = SIN(TTA)

C         if(xo2.ne.0.d0)   write(*,*) ' chamk ',i,xo2,yo2,tta,dy

          A1N =  (A1-XO2)*CT + (R1-YO2)*ST
          R1N = -(A1-XO2)*ST + (R1-YO2)*CT - DY
          CALL CHAMK1(A1N,R1N,Z1,NMPT(I),*999)
          CALL ADPOL(ID,IOP,B,DB,DDB,D3BX,D3BY,D4BX,D4BY,BT)

c         write(88,fmt='(1p,4e12.4,1x,2i4,1x,e12.4,a)')
c     >   A1,R1,A1N,R1N,i,NMPT(I),B(1,3),' A1,R1,A1N,R1N,imag, BT '

        ENDDO

        I = NFIC
c               write(*,*) 'chamk NMPT(I) ',i,NMPT(I)
        XO2 = UU(20+10*I+1)
        YO2 = UU(20+10*I+2)
        TTA = UU(20+10*I+3)
        DY =  UU(20+10*I+4)
C         if(xo2.ne.0.d0)   write(*,*) ' chamk ',i,xo2,yo2,tta,dy
        CT = COS(TTA)
        ST = SIN(TTA)
        A1N =  (A1-XO2)*CT + (R1-YO2)*ST
        R1N = -(A1-XO2)*ST + (R1-YO2)*CT - DY
        CALL CHAMK1(A1N,R1N,Z1,NMPT(I),*999)
        CALL ADPOL(ID,I3,B,DB,DDB,D3BX,D3BY,D4BX,D4BY,BT)

c         write(88,fmt='(1p,4e12.4,1x,2i4,1x,e12.4,a)')
c     >   A1,R1,A1N,R1N,i,NMPT(I),B(1,3),' A1,R1,A1N,R1N,imag, BT '

      ELSE
        CALL ENDJOB('Pgm chamk. No such possibility NFIC = ',NFIC)
      ENDIF

c               write(*,*) ' '

      RETURN

 999  RETURN 1

      ENTRY CHAMKW(NFICI,NMPTI,UUI)
      NFIC = NFICI
      NMPT = NMPTI
      UU = UUI
      RETURN

      END
