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
      SUBROUTINE REFER(IO,IORD,IFOC,IT1,IT2,IT3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     **********************************************
C     SETS THE REFERENCE FRAME
C     **********************************************
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      COMMON/CONST2/ ZERO, UN
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
 
      SAVE XI,YI,ALE,PATHL
      DATA PATHL / -9999.D0 / 
      GOTO (1,2) IO
 
 1    CONTINUE
C----- GOES TO NEW REFERENCE FRAME + computes new coordinates there
 
      IF    (IFOC .EQ. 0) THEN
        XI=0.D0
        YI=F(2,IT1)
        ALE=F(3,IT1)*.001D0
      ELSEIF(IFOC .EQ. 1) THEN
C------- RECHERCHE DES COORDONNEES DU POINT DE FOCALISATION
        IF    (IORD .EQ. 1) THEN
           CALL FOCAL1(IT1,IT2,IT3,XI,YI)
        ELSEIF(IORD .EQ. 2) THEN
           CALL FOCAL1(IT1,IT2,IT3,XI,YI)
        ENDIF
C------- LE SYSTEME DE REFERENCE POUR LE CALCUL DES COEFFICIENTS
C        DE TRANSFERT S'APPUIE SUR LA DIRECTION DE LA TRAJECTOIRE #1
        ALE=F(3,IT1)*.001D0
      ENDIF
      IF(IFOC.LE.1) THEN
        DO 8 I=1,IMAX
           CALL INITRA(I)
           CALL CHAREF(.FALSE.,XI,YI,ALE)
           CALL MAJTRA(I)
 8      CONTINUE
        IF(NRES .GT. 0) WRITE(NRES,100) XI,YI,ALE*DEG,ALE
 100    FORMAT(/,10X,' Frame for MATRIX calculation moved by :'
     >        ,/,10X,'  XC =',F9.3,' cm , YC =',F9.3,' cm ,   A ='
     >        ,F9.5,' deg  ( =',F9.6,' rad )',/)
      ENDIF
      PATHL = F(6,IT1)
      IF(NRES .GT. 0)WRITE(NRES,FMT='(/,1P,''Reference particle '',
     >''(#'',I4,''), path length :'',G16.8,'' cm'', 
     >''  relative momentum : '',G14.6)') IT1, F(6,IT1), F(1,IT1)
 
      RETURN
 
 2    CONTINUE
C----- COMES BACK TO OLD FRAME + OLD COORDINATES
 
      IF(IFOC.LE.1) THEN
         DO 81 I=1,IMAX
            CALL INITRA(I)
            CALL CHAREF(.FALSE.,ZERO,ZERO,-ALE)
            CALL CHAREF(.FALSE.,-XI,-YI,ZERO)
            CALL MAJTRA(I)
81       CONTINUE
      ENDIF
 
      RETURN

      ENTRY REFER1(
     >             PATHLO)
      PATHLO = PATHL
      RETURN
      
      END
