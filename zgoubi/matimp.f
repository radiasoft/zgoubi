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
      SUBROUTINE MATIMP(R) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(6,*) , T(6,6,*)
      DIMENSION  T3(5,*)

      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
  
      DETY=R(1,1)*R(2,2)-R(1,2)*R(2,1)
      DETZ=R(3,3)*R(4,4)-R(3,4)*R(4,3)
      RIJ = R(2,2)
      IF(RIJ .EQ. 0.D0) RIJ = 1.D-10
      SFH = - R(1,2)/RIJ
      RIJ = R(4,4)
      IF(RIJ .EQ. 0.D0) RIJ = 1.D-10
      SFZ = - R(3,4)/RIJ
 
      I=1
      WRITE(NRES,103) I
 103  FORMAT(//,18X,'TRANSFER  MATRIX  ORDRE',I3,'  (MKSA units)',/)
      WRITE(NRES,104) (( R(IA,IB) , IB=1,6) , IA=1,6)
 104  FORMAT(6X,1P,6G16.6)
      WRITE(NRES,112) DETY-1.D0,DETZ-1.D0
112   FORMAT(/,10X,'DetY-1 = ',F18.10,',',4X,'DetZ-1 = ',F18.10)
      WRITE(NRES,FMT='(/,10X,''R12=0 at '',G12.4,'' m, '',7X, 
     >                       ''R34=0 at '',G12.4,'' m'')') SFH,SFZ

      CALL SYMPL(R)

      RETURN

      ENTRY MATIM2(R,T,T3)

C MODIFIED, FM, 04/97
C       ** CHANGE MATRICE TRIANGULAIRE EN CARREE SYMMETRIQUE/DIAG
        DO 11 IA=1,6
          DO 11 IB=1,6
            IC1=IB+1
            DO 11 IC=IC1,6
              T(IA,IC,IB)=T(IA,IB,IC)
 11     CONTINUE

      I=2
      WRITE(NRES,103) I
      DO 16 IA=1,6
        IF(IA.GT.1) WRITE(NRES,107)
 107    FORMAT(/)
        DO 16 IB=1,6
C          WRITE(NRES,108) ( IA,IC,IB, T(IA,IC,IB)  , IC=1,IB )
C MODIFIED, FM, 04/97
          WRITE(NRES,108) ( IA,IC,IB, T(IA,IC,IB)  , IC=1,6 )
 108      FORMAT( 6(I4,I2,I1,1P,G11.3) )
 16   CONTINUE
 
      CALL SYMPL2(R,T)
 
      WRITE(NRES,123) T3(1,1),T3(1,2),T3(1,3),T3(1,4)
 123  FORMAT(//,15X,'COEFFICIENTS  D''ORDRE  SUPERIEUR  ( MKSA ):'
     >,//,10X,' Y/Y3   ',5X,1P,G14.5
     > ,/,10X,' Y/T3   ',5X,   G14.5
     > ,/,10X,' Y/Z3   ',5X,   G14.5
     > ,/,10X,' Y/P3   ',5X,   G14.5,/)
      WRITE(NRES,124) T3(2,1),T3(2,2),T3(2,3),T3(2,4)
 124  FORMAT(
     >    10X,' T/Y3   ',5X,1P,G14.5
     > ,/,10X,' T/T3   ',5X,   G14.5
     > ,/,10X,' T/Z3   ',5X,   G14.5
     > ,/,10X,' T/P3   ',5X,   G14.5,/)
      WRITE(NRES,125) T3(3,1),T3(3,2),T3(3,3),T3(3,4)
 125  FORMAT(
     >    10X,' Z/Y3   ',5X,1P,G14.5
     > ,/,10X,' Z/T3   ',5X,   G14.5
     > ,/,10X,' Z/Z3   ',5X,   G14.5
     > ,/,10X,' Z/P3   ',5X,   G14.5,/)
      WRITE(NRES,126) T3(4,1),T3(4,2),T3(4,3),T3(4,4)
 126  FORMAT(
     >    10X,' P/Y3   ',5X,1P,G14.5
     > ,/,10X,' P/T3   ',5X,   G14.5
     > ,/,10X,' P/Z3   ',5X,   G14.5
     > ,/,10X,' P/P3   ',5X,   G14.5)
 
      CONTINUE
      WRITE(NRES,101) IEX(1),(F(J,1),J=1,7)
  101 FORMAT(' TRAJ 1 IEX,D,Y,T,Z,P,S,time :',I3,1P,5G12.4,2G17.5)
      RETURN
      END
