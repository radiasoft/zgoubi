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
      SUBROUTINE TUNES(R,F0,NMAIL,IERY,IERZ,OKPR,
     >                                           YNU,ZNU,CMUY,CMUZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(6,*), F0(6,*)
      LOGICAL OKPR

      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M

C      CHARACTER*15 TXTYNU, TXTZNU

      IERY=0
      IERZ=0
      CALL RAZ(F0, 6*6)

C      IF(OKPR) THEN
C        IF(NMAIL.LE.0) THEN
C          IF(NRES.GT.0) 
C     >    WRITE(NRES,*) '  NUMBER OF PERIODS = IFOC-10 = ERRONEOUS !'
C          RETURN
C        ELSE
C          IF(NRES.GT.0) WRITE(NRES,106) NMAIL
C 106      FORMAT(//,15X,' TWISS  parameters,  periodicity  of',
C     >           I4,'  is  assumed ',/
C     >             ,35X,' - UNCOUPLED -')
C        ENDIF
C      ENDIF
 
C  HORIZONTAL
      CMUY = .5D0 * (R(1,1)+R(2,2))
      IF(CMUY*CMUY .GE. 1.D0) THEN
        IERY=-1
        F0(1,1)=0.D0
        F0(1,2)=0.D0
        F0(2,1)=0.D0
        F0(2,2)=0.D0
        F0(1,6)=0.D0
        F0(2,6)=0.D0
        YNU = 0.D0
      ELSE

        COSMU=CMUY
CBETA style
        SINMU=SIGN( SQRT(1.0D0-COSMU*COSMU) , R(1,2) )
        YNU=ACOS(COSMU)*DBLE(NMAIL)/(2.D0 * PI) 
        IF (R(1,2) .LT. 0.0D0) YNU=DBLE(NMAIL)-YNU
        YNU=YNU-AINT(YNU)
C         write(*,*) ' beta style, sinmu, ynu = ',sinmu,ynu,nmail
CMAD style. Optimized precisison close to 1/2 ? 
        SINMU=SIGN(SQRT(-R(1,2)*R(2,1)-.25D0*(R(1,1)-R(2,2))**2),R(1,2))
        YNU = SIGN(ATAN2(SINMU,COSMU) /(2.D0 * PI) ,R(1,2))
        IF (R(1,2) .LT. 0.0D0) YNU=DBLE(NMAIL)+YNU
C         write(*,*) ' mad style, sinmu, ynu = ',sinmu,ynu

        BX0=R(1,2) / SINMU
        AX0=(R(1,1) - R(2,2)) / (2.D0*SINMU)
        F0(1,1)=BX0
        F0(1,2)=-AX0
        F0(2,1)=-AX0
        F0(2,2)=(1.D0+AX0*AX0)/BX0
        S2NU = 2.D0 - R(1,1) - R(2,2)
        F0(1,6) = ( (1.D0 - R(2,2))*R(1,6) + R(1,2)*R(2,6) ) / S2NU
        F0(2,6) = ( R(2,1)*R(1,6) + (1.D0 - R(1,1))*R(2,6) ) / S2NU

      ENDIF

C  VERTICAL
      CMUZ = .5D0 * (R(3,3)+R(4,4))
      IF(CMUZ*CMUZ .GE. 1.D0) THEN
        IERZ=-1
        F0(3,3)=0.D0
        F0(3,4)=0.D0
        F0(4,3)=0.D0
        F0(4,4)=0.D0
        F0(3,6)=0.D0
        F0(4,6)=0.D0
        ZNU = 0.D0
      ELSE
        COSMU = CMUZ
CBETA style
            ZNU=ACOS(COSMU)*DBLE(NMAIL)/(2.D0 * PI) 
            IF (R(3,4) .LT. 0.0D0) ZNU=DBLE(NMAIL)-ZNU
            ZNU=ZNU-AINT(ZNU)
            SINMU=SIGN( SQRT(1.0D0-COSMU*COSMU) , R(3,4) )
C         write(*,*) ' beta style, sinmu, znu = ',sinmu,znu,nmail
CMAD style. Optimized precisison close to 1/2 ? 
        SINMU=SIGN(SQRT(-R(3,4)*R(4,3)-.25D0*(R(3,3)-R(4,4))**2),R(3,4))
        ZNU = SIGN(ATAN2(SINMU,COSMU) /(2.D0 * PI) ,R(3,4))
        IF (R(3,4) .LT. 0.0D0) ZNU=DBLE(NMAIL)+ZNU
C         write(*,*) ' mad style, sinmu, znu = ',sinmu,znu
C         write(*,*) '                                tunes...' 

        BZ0=R(3,4)/SINMU
        AZ0=(R(3,3)-R(4,4))/2.D0/SINMU
        F0(3,3)=BZ0
        F0(3,4)=-AZ0
        F0(4,3)=-AZ0
        F0(4,4)=(1.D0+AZ0*AZ0)/BZ0
        S2NU = 2.D0 - R(3,3) - R(4,4)
        F0(3,6) = ( (1.D0 - R(4,4))*R(3,6) + R(3,4)*R(4,6) ) / S2NU
        F0(4,6) = ( R(4,3)*R(3,6) + (1.D0 - R(3,3))*R(4,6) ) / S2NU

      ENDIF

C      IF(OKPR) THEN
C        IF(NRES.GT.0) THEN
C          WRITE(NRES,103)
C 103      FORMAT(/,6X,
C     >    ' Beam  matrix  (beta/-alpha/-alpha/gamma)',
C     >    ' and  periodic  dispersion  (MKSA units)',/)
C          WRITE(NRES,104) (( F0(IA,IB) , IB=1,6) , IA=1,6)
C 104      FORMAT(6X,6F13.6)
C          WRITE(NRES,FMT='(/,35X,''Betatron  tunes'',/)') 
C          WRITE(TXTYNU,FMT='(A)') 'undefined'
C          WRITE(TXTZNU,FMT='(A)') 'undefined'

C         IF    (ABS(CMUY).LT.1.D0 .AND. ABS(CMUZ).LT.1.D0) THEN
C           WRITE(TXTYNU,FMT='(G15.8)') YNU
C           WRITE(TXTZNU,FMT='(G15.8)') ZNU
C         ELSEIF(ABS(CMUY).LT.1.D0 .OR. ABS(CMUZ).LT.1.D0) THEN
C           IF(CMUY*CMUY .LT. 1.D0) WRITE(TXTYNU,FMT='(G15.8)') YNU
C           IF(CMUZ*CMUZ .LT. 1.D0) WRITE(TXTZNU,FMT='(G15.8)') ZNU
C         ENDIF

C         WRITE(NRES,FMT='(15X,2(5X,A,A))') 
C     >              'NU_Y = ', TXTYNU, 'NU_Z = ', TXTZNU

C        ENDIF
C      ENDIF
      RETURN
      END
