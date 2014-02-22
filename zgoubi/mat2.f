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
      SUBROUTINE MAT2(R,T,T3,T4)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(6,*) , T(6,6,*)
      DIMENSION  T3(5,*), T4(5,*)
C     ------------------------------------
C     Option  IORD = 2 :
C       Matrix ordre 1, coeff ordre 2, & 
C        coeffs > 2; Taylor series ordre 5
C     ------------------------------------
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
C      COMMON/DON/ A(09876,99),IQ(09876),IP(09876),NB,NOEL
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
 
      DIMENSION RPD(6,6), RMD(6,6), RPDO(6,6), RMDO(6,6) 
      SAVE RPD, RMD, DP

      CALL RAZ(R,6*6)
      R(5,5) = 1D0
      R(6,6) = 1D0
      CALL RAZ(T,5*6*6)
      S1 = 2.D0 * F(6,1)
 
C      ++++ Ordre  1  ,  &  ordres 2  to  4  uncoupled
C           +++ utilise  traj. 1-21
 
      DP = ( FO(1,18)-FO(1,19) )/( FO(1,18)+FO(1,19) )
      DP2 = DP*DP
      DO 10 J=2,5
        J1=J-1
        QP = F(J,18)
        QM = F(J,19)
        QPP= F(J,20)
        QMM= F(J,21)
        R(J1,6) =  ( 8.D0*(QP-QM) -(QPP-QMM) )/12.D0/DP
        T(J1,6,6)= (16.D0*(QP+QM) -(QPP+QMM) )/24.D0/DP2
        T3(J1,6)=-( 2.D0*(QP-QM) -(QPP-QMM) )/12.D0/DP2/DP
        T4(J1,6)=-( 4.D0*(QP+QM) -(QPP+QMM) )/24.D0/DP2/DP2
        DO 10 I=1,4
          II = 4*I - 2
          UO = FO(I+1,II)
          UO2 = UO*UO
          QP = F(J,II)
          QM = F(J,II+1)
          QPP= F(J,II+2)
          QMM= F(J,II+3)
          R(J1,I) =  ( 8.D0*(QP-QM) -(QPP-QMM) ) /12.D0 / UO
          T(J1,I,I)= (16.D0*(QP+QM) -(QPP+QMM) ) /24.D0 / UO2
          T3(J1,I)=-( 2.D0*(QP-QM) -(QPP-QMM) ) /12.D0 / UO2/UO
          T4(J1,I)=-( 4.D0*(QP+QM) -(QPP+QMM) ) /24.D0 /  UO2/UO2
          IF(J .EQ. 5) THEN
            R(5,I) = ( 8.D0*( F(6,II) - F(6,II+1) )
     >        - ( F(6,II+2) - F(6,II+3) )          ) /12.D0/UO
            T(5,I,I)=( 16.D0*( F(6,II) +F(6,II+1) - S1 )
     >        - ( F(6,II+2) + F(6,II+3) - S1 )     ) /24.D0/UO2
          ENDIF
 10   CONTINUE
      QP = F(6,18)
      QM = F(6,19)
      QPP= F(6,20)
      QMM= F(6,21)
      R(5,6)  = ( 8.D0*( QP-QM   ) - ( QPP-QMM   )  )/12.D0/DP
      T(5,6,6)= (16.D0*( QP+QM-S1) - ( QPP+QMM-S1)  )/24.D0/DP2
      T3(5,6)=-( 2.D0*( QP-QM   ) - ( QPP-QMM   )  )/12.D0/DP2/DP
      T4(5,6)=-( 4.D0*( QP+QM-S1) - ( QPP+QMM-S1)  )/24.D0/DP2/DP2
 
      CALL RAZ(RPD,4*4)
      DO 20 J=2,5
        J1=J-1
          UO = FO(2,34) - FO(2,37)
          RPD(J1,1) = (F(J,34) - F(J,37)) / UO
          UO = FO(3,46) - FO(3,49)
          RPD(J1,2) = (F(J,46) - F(J,49)) / UO
          UO = FO(4,54) - FO(4,57)
          RPD(J1,3) = (F(J,54) - F(J,57)) / UO
          UO = FO(5,58) - FO(5,61)
          RPD(J1,4) = (F(J,58) - F(J,61)) / UO
 20   CONTINUE
 
      CALL RAZ(RMD,4*4)
      DO 21 J=2,5
        J1=J-1
          UO = FO(2,35) - FO(2,36)
          RMD(J1,1) = (F(J,35) - F(J,36)) / UO
          UO = FO(3,47) - FO(3,48)
          RMD(J1,2) = (F(J,47) - F(J,48)) / UO
          UO = FO(4,55) - FO(4,56)
          RMD(J1,3) = (F(J,55) - F(J,56)) / UO
          UO = FO(5,59) - FO(5,60)
          RMD(J1,4) = (F(J,59) - F(J,60)) / UO
 21    CONTINUE
 
C      ++++ Ordre 2  coupled
C           +++ utilise  traj. >=  #22
 
      III=4
C     ++++ COUPLAGE ./YT
      II=22
      J1=1
      J2=2
      UO = FO(J1+1,II)
      VO = FO(J2+1,II)
C     ... J/J1.J2 : Y/YT  ... L/YT
      DO 12 J=1,5
        U=UO
        V=VO
        QPP = F(J+1,II)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPP=QPP-F(6,1)
        U=UO
        V=-VO
        QPM = F(J+1,II+1)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPM=QPM-F(6,1)
        U=-UO
        V=-VO
        QMM = F(J+1,II+2)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMM=QMM-F(6,1)
        U=-UO
        V=VO
        QMP = F(J+1,II+3)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMP=QMP-F(6,1)
        T(J,J1,J2)=(QPP+QMM-QPM-QMP)/(4.D0*UO*VO)   *0.5D0
 12   CONTINUE
 
C     ++++ COUPLAGE ./YZ
      II=II+III
      J1=1
      J2=3
      UO = FO(J1+1,II)
      VO = FO(J2+1,II)
C     ... J/J1.J2 : Y/YZ  ... L/YZ
      DO 13 J=1,5
        U=UO
        V=VO
        QPP = F(J+1,II)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPP=QPP-F(6,1)
        U=UO
        V=-VO
        QPM = F(J+1,II+1)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPM=QPM-F(6,1)
        U=-UO
        V=-VO
        QMM = F(J+1,II+2)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMM=QMM-F(6,1)
        U=-UO
        V=VO
        QMP = F(J+1,II+3)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMP=QMP-F(6,1)
        T(J,J1,J2)=(QPP+QMM-QPM-QMP)/(4.D0*UO*VO)    *0.5D0
 13   CONTINUE
 
C     ++++ COUPLAGE ./YP
      II=II+III
      J1=1
      J2=4
      UO = FO(J1+1,II)
      VO = FO(J2+1,II)
C     ... J/J1.J2 : Y/YP  ... L/YP
      DO 14 J=1,5
        U=UO
        V=VO
        QPP = F(J+1,II)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPP=QPP-F(6,1)
        U=UO
        V=-VO
        QPM = F(J+1,II+1)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPM=QPM-F(6,1)
        U=-UO
        V=-VO
        QMM = F(J+1,II+2)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMM=QMM-F(6,1)
        U=-UO
        V=VO
        QMP = F(J+1,II+3)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMP=QMP-F(6,1)
        T(J,J1,J2)=(QPP+QMM-QPM-QMP)/(4.D0*UO*VO)    *0.5D0
 14   CONTINUE
 
C     ++++ COUPLAGE ./YD
      II=II+III
      J1=1
      J2=6
      UO = FO(J1+1,II)
      VO = FO(   1,II)-FO(1,1)
C     ... J/J1.J2 : Y/YD  ... L/YD
      DO 16 J=1,5
        U=UO
        V=VO
        QPP = F(J+1,II)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPP=QPP-F(6,1)
        U=UO
        V=-VO
        QPM = F(J+1,II+1)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPM=QPM-F(6,1)
        U=-UO
        V=-VO
        QMM = F(J+1,II+2)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMM=QMM-F(6,1)
        U=-UO
        V=VO
        QMP = F(J+1,II+3)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMP=QMP-F(6,1)
        T(J,J1,J2)=(QPP+QMM-QPM-QMP)/(4.D0*UO*VO)    *0.5D0
 16   CONTINUE
 
C     ++++ COUPLAGE ./TZ
      II=II+III
      J1=2
      J2=3
      UO = FO(J1+1,II)
      VO = FO(J2+1,II)
C     ... J/J1.J2 : Y/TZ  ... L/TZ
      DO 23 J=1,5
        U=UO
        V=VO
        QPP = F(J+1,II)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPP=QPP-F(6,1)
        U=UO
        V=-VO
        QPM = F(J+1,II+1)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPM=QPM-F(6,1)
        U=-UO
        V=-VO
        QMM = F(J+1,II+2)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMM=QMM-F(6,1)
        U=-UO
        V=VO
        QMP = F(J+1,II+3)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMP=QMP-F(6,1)
        T(J,J1,J2)=(QPP+QMM-QPM-QMP)/(4.D0*UO*VO)    *0.5D0
 23   CONTINUE
 
C     ++++ COUPLAGE ./TP
      II=II+III
      J1=2
      J2=4
      UO = FO(J1+1,II)
      VO = FO(J2+1,II)
C     ... J/J1.J2 : Y/TP  ... L/TP
      DO 24 J=1,5
        U=UO
        V=VO
        QPP = F(J+1,II)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPP=QPP-F(6,1)
        U=UO
        V=-VO
        QPM = F(J+1,II+1)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPM=QPM-F(6,1)
        U=-UO
        V=-VO
        QMM = F(J+1,II+2)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMM=QMM-F(6,1)
        U=-UO
        V=VO
        QMP = F(J+1,II+3)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMP=QMP-F(6,1)
        T(J,J1,J2)=(QPP+QMM-QPM-QMP)/(4.D0*UO*VO)    *0.5D0
 24   CONTINUE
 
C     ++++ COUPLAGE ./TD
      II=II+III
      J1=2
      J2=6
      UO = FO(J1+1,II)
      VO = FO(   1,II)-FO(1,1)
C     ... J/J1.J2 : Y/TD  ... L/TD
      DO 26 J=1,5
        U=UO
        V=VO
        QPP = F(J+1,II)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPP=QPP-F(6,1)
        U=UO
        V=-VO
        QPM = F(J+1,II+1)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPM=QPM-F(6,1)
        U=-UO
        V=-VO
        QMM = F(J+1,II+2)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMM=QMM-F(6,1)
        U=-UO
        V=VO
        QMP = F(J+1,II+3)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMP=QMP-F(6,1)
        T(J,J1,J2)=(QPP+QMM-QPM-QMP)/(4.D0*UO*VO)    *0.5D0
 26   CONTINUE
 
C     ++++ COUPLAGE ./ZP
      II=II+III
      J1=3
      J2=4
      UO = FO(J1+1,II)
      VO = FO(J2+1,II)
C     ... J/J1.J2 : Y/ZP ...  L/ZP
      DO 34 J=1,5
        U=UO
        V=VO
        QPP = F(J+1,II)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPP=QPP-F(6,1)
        U=UO
        V=-VO
        QPM = F(J+1,II+1)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPM=QPM-F(6,1)
        U=-UO
        V=-VO
        QMM = F(J+1,II+2)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMM=QMM-F(6,1)
        U=-UO
        V=VO
        QMP = F(J+1,II+3)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMP=QMP-F(6,1)
        T(J,J1,J2)=(QPP+QMM-QPM-QMP)/(4.D0*UO*VO)    *0.5D0
 34   CONTINUE
 
C     ++++ COUPLAGE ./ZD
      II=II+III
      J1=3
      J2=6
      UO = FO(J1+1,II)
      VO = FO(   1,II)-F(1,1)
C     ... J/J1.J2 : Y/ZD  ... L/ZD
      DO 36 J=1,5
        U=UO
        V=VO
        QPP = F(J+1,II)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPP=QPP-F(6,1)
        U=UO
        V=-VO
        QPM = F(J+1,II+1)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPM=QPM-F(6,1)
        U=-UO
        V=-VO
        QMM = F(J+1,II+2)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMM=QMM-F(6,1)
        U=-UO
        V=VO
        QMP = F(J+1,II+3)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMP=QMP-F(6,1)
        T(J,J1,J2)=(QPP+QMM-QPM-QMP)/(4.D0*UO*VO)    *0.5D0
 36   CONTINUE
 
C     ++++ COUPLAGE ./PD
      II=II+III
      J1=4
      J2=6
      UO = FO(J1+1,II)
      VO = FO(   1,II)-F(1,1)
C     ... J/J1.J2 : Y/PD  ... L/PD
      DO 46 J=1,5
        U=UO
        V=VO
        QPP = F(J+1,II)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPP=QPP-F(6,1)
        U=UO
        V=-VO
        QPM = F(J+1,II+1)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPM=QPM-F(6,1)
        U=-UO
        V=-VO
        QMM = F(J+1,II+2)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMM=QMM-F(6,1)
        U=-UO
        V=VO
        QMP = F(J+1,II+3)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMP=QMP-F(6,1)
        T(J,J1,J2)=(QPP+QMM-QPM-QMP)/(4.D0*UO*VO)    *0.5D0
 46   CONTINUE
 
      RETURN

      ENTRY MAT2P(RPDO,DPO)
      DO 66 J=1,6
        DO 66 I=1,6
 66       RPDO(I,J) = RPD(I,J) 
      DPO = DP
      RETURN

      ENTRY MAT2M(RMDO,DPO)
      DO 67 J=1,6
        DO 67 I=1,6
 67        RMDO(I,J) = RMD(I,J) 
      DPO = DP
      RETURN

      END
