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
C  Brookhaven National Laboratory                    és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE INTPL3(A,R,Z,DA,DR,DZ,FTAB3,
     >                                       B,DB,DDB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ---------------------------------------------------------
C     Interpolation of field values in polar coord. 
C     --------------------------------------------------------
      DIMENSION FTAB3(3,3,3,3)
      DIMENSION B(5,3),DB(3,3),DDB(3,3,3)
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      INCLUDE "C.CHAMBR.H"     ! COMMON/CHAMBR/ LIMIT,IFORM,YLIM2,ZLIM2,SORT(MXT),FMAG,YCH,ZCH
 
      INCLUDE "C.INTEG.H"     ! COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      INCLUDE "C.MARK.H"     ! COMMON/MARK/ KART,KALC,KERK,KUASEX
C      LOGICAL ZSYM
      INCLUDE "C.TYPFLD.H"     ! COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
 
      DIMENSION A000(3),A100(3),A010(3),A001(3),A200(3),A020(3),
     >A002(3),A110(3),A101(3),A011(3)

      DO 415 L = 1,3
        A000(L)=0.D0
        A100(L)=0.D0
        A010(L)=0.D0
        A001(L)=0.D0
        A200(L)=0.D0
        A020(L)=0.D0
        A002(L)=0.D0
        A110(L)=0.D0
        A101(L)=0.D0
        A011(L)=0.D0
        DO 416 K=1,3
          KZ=K-2
          DO 416 J=1,3
            JR=2-J
            DO 416 I=1,3
              IA=I-2
              BIJK=FTAB3(L,I,J,K)     ! HC(L,IAC+IA,IRC+JR,IZC+KZ)
C              BMESH3(K,I,J) = BIJK
C                 write(88,*) FTAB3(L,I,J,K), L,IA,JR,KZ
              A000(L)=A000(L) + 
     >             DBLE(7-3*(IA*IA+JR*JR+KZ*KZ))/3.D0 *BIJK
              A100(L)=A100(L) +        IA        *BIJK
              A010(L)=A010(L) +        JR        *BIJK
              A001(L)=A001(L) +        KZ        *BIJK
              A200(L)=A200(L) +  DBLE(3*IA*IA-2)   *BIJK
              A020(L)=A020(L) +  DBLE(3*JR*JR-2)   *BIJK
              A002(L)=A002(L) +  DBLE(3*KZ*KZ-2)   *BIJK
              A110(L)=A110(L) +      IA*JR       *BIJK
              A101(L)=A101(L) +      IA*KZ       *BIJK
              A011(L)=A011(L) +      JR*KZ       *BIJK
 416    CONTINUE

C        CALL MAPLIM(*999, 27, BMESH3)

        A000(L)=A000(L)/( 9.D0      )*BRI
        A100(L)=A100(L)/(18.D0*DA   )*BRI
        A010(L)=A010(L)/(18.D0*DR   )*BRI
        A001(L)=A001(L)/(18.D0*DZ   )*BRI
        A200(L)=A200(L)/(18.D0*DA*DA)*BRI
        A020(L)=A020(L)/(18.D0*DR*DR)*BRI
        A002(L)=A002(L)/(18.D0*DZ*DZ)*BRI
        A110(L)=A110(L)/(12.D0*DA*DR)*BRI
        A101(L)=A101(L)/(12.D0*DA*DZ)*BRI
        A011(L)=A011(L)/(12.D0*DR*DZ)*BRI
 415  CONTINUE
 
      DO 417 L = 1,3
C           write(89,*) ' intpl3 b1L, L F : ',a000(l)/bri,l,a,r,z
        B(1,L)=A000(L)     
     >       + A100(L)*A   + A010(L)*R   + A001(L)*Z
     >       + A200(L)*A*A + A020(L)*R*R + A002(L)*Z*Z
     >       + A110(L)*A*R + A101(L)*A*Z + A011(L)*R*Z
C FM, Check using SIGMAPHI fieldmap
C                if(Z1.eq.0.D0) then    ! here I can see that A000 .ne.0 -> not correct
C                 write(89,fmt='(i4,6g14.6)') 
C     >               l,A000(L),b(1,1),b(1,2),b(1,3),db(1,1),db(1,2)
C                endif
        DB(1,L)  = A100(L) + 2.D0*A200(L)*A + A110(L)*R + A101(L)*Z
        DB(2,L)  = A010(L) + 2.D0*A020(L)*R + A110(L)*A + A011(L)*Z
        DB(3,L)  = A001(L) + 2.D0*A002(L)*Z + A101(L)*A + A011(L)*R
        DDB(1,1,L) = 2.D0*A200(L)
        DDB(1,2,L) = A110(L)
        DDB(2,1,L) = A110(L)
        DDB(2,2,L) = 2.D0*A020(L)
        DDB(1,3,L) = A101(L)
        DDB(3,1,L) = A101(L)
        DDB(2,3,L) = A011(L)
        DDB(3,2,L) = A011(L)
        DDB(3,3,L) = 2.D0*A002(L)
 417  CONTINUE

      RETURN
      END
