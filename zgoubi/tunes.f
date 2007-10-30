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
      SUBROUTINE TUNES(R,F0,NMAIL,IERY,IERZ,OKPR,
     >                                           YNU,ZNU,CMUY,CMUZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(6,*), F0(6,*)
      LOGICAL OKPR

      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M

      CHARACTER*14 TXTYNU, TXTZNU

      IERY=0
      IERZ=0
      CALL RAZ(F0, 6*6)

      IF(OKPR) THEN
      IF(NMAIL.LE.0) THEN
        WRITE(NRES,*) '  NUMBER OF PERIODS = IFOC-10 = ERRONEOUS !'
        RETURN
      ELSE
        WRITE(NRES,106) NMAIL
 106    FORMAT(//,15X,' TWISS  parameters,  periodicity  of',
     >         I4,'  is  assumed :')
      ENDIF
      ENDIF
 
C  HORIZONTAL
      CMUY = .5D0 * (R(1,1)+R(2,2))
      IF(CMUY*CMUY .GE. 1.D0) THEN
        IERY=-1
C        WRITE(NRES,102) 'Y','NU_Y'
        F0(1,1)=0.D0
        F0(1,2)=0.D0
        F0(2,1)=0.D0
        F0(2,2)=0.D0
        F0(1,6)=0.D0
        F0(2,6)=0.D0
      ELSE

CBETA 
         COSMU=CMUY
            YNU=ACOS(COSMU)*NMAIL/(2.D0 * PI) 
            IF (R(1,2) .LT. 0.0D0) YNU=NMAIL-YNU
            YNU=YNU-AINT(YNU)
            SINMU=SIGN( SQRT(1.0D0-COSMU*COSMU) , R(1,2) )

CC        SINMU=SQRT(-R(1,2)*R(2,1)-.25D0*(R(1,1)-R(2,2))*(R(1,1)-R(2,2)))
CC        IF(R(1,2) .LT. 0.D0) SINMU = - SINMU
CC        YNU = ATAN2( SINMU , COSMU ) * NMAIL / (2.D0 * PI)  

        BX0=R(1,2) / SINMU
        AX0=(R(1,1) - R(2,2)) / (2.D0*SINMU)
        F0(1,1)=BX0
        F0(1,2)=-AX0
        F0(2,1)=-AX0
        F0(2,2)=(1+AX0*AX0)/BX0
        S2NU = 2.D0 - R(1,1) - R(2,2)
        F0(1,6) = ( (1.D0 - R(2,2))*R(1,6) + R(1,2)*R(2,6) ) / S2NU
        F0(2,6) = ( R(2,1)*R(1,6) + (1.D0 - R(1,1))*R(2,6) ) / S2NU

      ENDIF

C  VERTICAL
      CMUZ = .5D0 * (R(3,3)+R(4,4))
      IF(CMUZ*CMUZ .GE. 1.D0) THEN
        IERZ=-1
C        WRITE(NRES,102) 'Z', 'NU_Z'
        F0(3,3)=0.D0
        F0(3,4)=0.D0
        F0(4,3)=0.D0
        F0(4,4)=0.D0
        F0(3,6)=0.D0
        F0(4,6)=0.D0
      ELSE
CBETA 
        COSMU = CMUZ
            ZNU=ACOS(COSMU)*NMAIL/(2.D0 * PI) 
            IF (R(3,4) .LT. 0.0D0) ZNU=NMAIL-ZNU
            ZNU=ZNU-AINT(ZNU)
            SINMU=SIGN( SQRT(1.0D0-COSMU*COSMU) , R(3,4) )

C        SINMU=SQRT(-R(3,4)*R(4,3)-.25D0*(R(3,3)-R(4,4))*(R(3,3)-R(4,4)))
C        IF(R(3,4) .LT. 0.D0) SINMU = - SINMU
C        ZNU = ATAN2( SINMU , COSMU ) * NMAIL / (2.D0 * PI)  

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

      IF(OKPR) THEN

        WRITE(NRES,103)
 103    FORMAT(/,6X,
     >    ' Beam  matrix  (beta/-alpha/-alpha/gamma)',
     >    ' and  periodic  dispersion  (MKSA units)',/)
        WRITE(NRES,104) (( F0(IA,IB) , IB=1,6) , IA=1,6)
 104    FORMAT(6X,6F13.6)
        WRITE(NRES,FMT='(/,35X,''Betatron  tunes'',/)') 

       IF    (ABS(CMUY).LT.1.D0 .AND. ABS(CMUZ).LT.1.D0) THEN

C        WRITE(NRES,103)
C 103    FORMAT(/,6X,
C     >    ' Beam  matrix  (beta/-alpha/-alpha/gamma)',
C     >    ' and  periodic  dispersion  (MKSA units)',/)
C        WRITE(NRES,104) (( F0(IA,IB) , IB=1,6) , IA=1,6)
C 104    FORMAT(6X,6F13.6)
C        WRITE(NRES,FMT='(/,35X,''Betatron  tunes'',/)') 

         WRITE(TXTYNU,FMT='(G14.8)') YNU
         WRITE(TXTZNU,FMT='(G14.8)') ZNU

       ELSEIF(ABS(CMUY).LT.1.D0 .OR. ABS(CMUZ).LT.1.D0) THEN

C        WRITE(NRES,103)
C        WRITE(NRES,104) (( F0(IA,IB) , IB=1,6) , IA=1,6)
C        WRITE(NRES,FMT='(/,35X,''Betatron  tunes'',/)') 

         IF(CMUY*CMUY .LT. 1.D0) THEN
           WRITE(TXTYNU,FMT='(G14.8)') YNU
C          WRITE(NRES,100) 'NU_Y', YNU
C 100      FORMAT(35X,A,' = ',G14.8)
         ELSE
           WRITE(TXTYNU,FMT='(A14)') 'undefined'
         ENDIF
         IF(CMUZ*CMUZ .LT. 1.D0) THEN
           WRITE(TXTZNU,FMT='(G14.8)') ZNU
         ELSE
           WRITE(TXTZNU,FMT='(A14)') 'undefined'
         ENDIF
       ENDIF

       WRITE(NRES,FMT='(15X,2(5X,A,A14))') 
     >              'NU_Y = ', TXTYNU, 'NU_Z = ', TXTZNU

      ENDIF
      RETURN
      END
