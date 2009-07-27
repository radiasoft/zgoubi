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
      SUBROUTINE IMPDEV
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      COMMON/CHAVE/ B(5,3),V(5,3),E(5,3)
      COMMON/CNTRLB/ XSTP,DIVB,ALAPL(4),RTN(3)
      COMMON/DDBXYZ/ DB(3,3),DDB(3,3,3)
      COMMON/D3BXYZ/ D3BX(3,3,3), D3BY(3,3,3), D3BZ(3,3,3)
      COMMON/D4BXYZ/ D4BX(3,3,3,3),D4BY(3,3,3,3),D4BZ(3,3,3,3)
      COMMON/DDEXYZ/ DE(3,3),DDE(3,3,3)
      COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      LOGICAL ZSYM
      COMMON/OPTION/ KFLD,MG,LC,ML,ZSYM
      COMMON/RIGID/ BORO,DPREF,DP,BR
      COMMON/TRAJ/ Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT
      COMMON/VITES/ U(6,3),DBR(6),DDT(6)
 
C------- CONDITIONS DE MAXWELL
        XSTP = XSTP + 1
        DIVB = DB(1,1) + DB(2,2) + DB(3,3)
        ALAPL(1) = DDB(1,1,1) + DDB(1,2,2) + DDB(1,3,3)
        ALAPL(2) = DDB(2,1,1) + DDB(2,2,2) + DDB(2,3,3)
        ALAPL(3) = DDB(3,1,1) + DDB(3,2,2) + DDB(3,3,3)
        RTN(1) = DB(2,3) - DB(3,2)
        RTN(2) = DB(3,1) - DB(1,3)
        RTN(3) = DB(1,2) - DB(2,1)

      IF(NRES .GT. 0) THEN
        IF(LST .EQ. 1) THEN
C--------- LST LE Champ SUR TRAJ. dans zgoubi.res
C          WRITE(*,105) IT,BR/BORO,Y,T,Z,P,X,SAR
          WRITE(NRES,105) IT,BR/BORO,Y,T,Z,P,X,SAR
 105      FORMAT(/,' SPGM   IMPDEV',5X,'TRAJ ',I3
     >    ,/,' D,Y,T,Z,P,X,S :',/,' ',1P,7E13.5)
 
          IF(KFLD .EQ. LC .OR. KFLD .EQ. ML) THEN
 
            WRITE(NRES,104) BR*E(1,1),BR*E(1,2),BR*E(1,3)
     >      ,BR*DE(1,1),BR*DE(2,1),BR*DE(3,1),BR*DE(2,2),BR*DE(3,2)
     >      ,BR*DDE(1,1,1),BR*DDE(2,1,1),BR*DDE(3,1,1)
     >      ,BR*DDE(2,2,1),BR*DDE(3,2,1),BR*DDE(3,3,1)
 104        FORMAT(' EX, EY, EZ (MV/cm):',1P,3E11.3
     >          ,/,' E'' :',5E11.3
     >          ,/,' E'''':',6E11.3)
          ENDIF
 
          IF(KFLD .EQ. MG .OR. KFLD .EQ. ML) THEN
 
C            WRITE(*,106) BR*B(1,1),BR*B(1,2),BR*B(1,3)
C     >      ,BR*DB(1,1),BR*DB(2,1),BR*DB(3,1),BR*DB(2,2),BR*DB(3,2)
C     >      ,BR*DDB(1,1,1),BR*DDB(2,1,1),BR*DDB(3,1,1)
C     >      ,BR*DDB(2,2,1),BR*DDB(3,2,1),BR*DDB(3,3,1)
            WRITE(NRES,106) BR*B(1,1),BR*B(1,2),BR*B(1,3)
     >      ,BR*DB(1,1),BR*DB(2,1),BR*DB(3,1),BR*DB(2,2),BR*DB(3,2)
     >      ,BR*DDB(1,1,1),BR*DDB(2,1,1),BR*DDB(3,1,1)
     >      ,BR*DDB(2,2,1),BR*DDB(3,2,1),BR*DDB(3,3,1)
 106        FORMAT(' BX, BY, BZ (kG)   :',1P,3E11.3
     >          ,/,' B'' :',5E11.3
     >          ,/,' B'''':',6E11.3)
          ENDIF
 
        ELSEIF(LST .EQ. 7) THEN
          WRITE(NRES,*)
          WRITE(NRES,105) IT,BR/BORO,Y,T,Z,P,X,SAR
          WRITE(NRES,109) ( J,  B(1,J),J=1,3)
 109      FORMAT('    B(J) ',3(1X,I1,1X,1P,E11.3))
          WRITE(NRES,100) (( I,J,  DB(I,J)  ,I=1,3),J=1,3)
 100      FORMAT('     DB(I,J) ',6(1X,2I1,1X,1P,E11.3))
          WRITE(NRES,101) ((( I,J,K,DDB(I,J,K) ,I=1,3),J=1,3),K=1,3)
 101      FORMAT('  DDB(I,J,K) ',6(1X,3I1,1X,1P,E11.3))
          WRITE(NRES,102)((( I,J,K,D3BX(I,J,K) ,I=1,3) ,J=1,3) ,K=1,3)
          WRITE(NRES,102)((( I,J,K,D3BY(I,J,K) ,I=1,3) ,J=1,3) ,K=1,3)
          WRITE(NRES,102)((( I,J,K,D3BZ(I,J,K) ,I=1,3) ,J=1,3) ,K=1,3)
 102      FORMAT('  D3B(I,J,K) ',6(1X,3I1,1X,1P,E11.3))
          WRITE(NRES,103) ((((D4BX(I,J,K,L),I=1,3),J=1,3),K=1,3),L=1,3)
          WRITE(NRES,103) ((((D4BY(I,J,K,L),I=1,3),J=1,3),K=1,3),L=1,3)
          WRITE(NRES,103) ((((D4BZ(I,J,K,L),I=1,3),J=1,3),K=1,3),L=1,3)
 103      FORMAT('  D4B(I,J,K,L) ',1P,9E11.3)

          WRITE(NRES,113) ( J,  E(1,J),J=1,3)
 113      FORMAT('    E(J) ',3(1X,I1,1X,1P,E11.3))
          WRITE(NRES,107) (( I,J,  DE(I,J)  ,I=1,3),J=1,3)
 107      FORMAT('     DE(I,J) ',6(1X,2I1,1X,1P,E11.3))
          WRITE(NRES,108) ((( I,J,K,DDE(I,J,K) ,I=1,3),J=1,3),K=1,3)
 108      FORMAT('  DDE(I,J,K) ',6(1X,3I1,1X,1P,E11.3))

        ELSEIF(LST .EQ. 8) THEN
          WRITE(NRES,110) ((U(I,J),I=1,6),J=1,3)
 110      FORMAT(' U(1-6,1-3) ',1P,9E11.3)
          WRITE(NRES,111) ((V(I,J),I=1,5),J=1,3)
 111      FORMAT(' V(1-5,1-3) ',1P,8E11.3)
          WRITE(NRES,112) ((B(I,J),I=1,5),J=1,3)
 112      FORMAT(' B(1-5,1-3) ',1P,8E11.3)

        ELSEIF(LST .EQ. 9) THEN

C--------- READ STEP #
          CALL INTEG4(NSTEP)
          K = LSTK()

C--------- Case  IL=2*10**n with n>1 (IL=20, 200, 2000...)
          IF(1+K*((NSTEP-1)/K) .NE. NSTEP) RETURN

          WRITE (88,fmt='(1p,20g12.4)') 
C                           BX     BY     BZ
     >    X,       Y,       B(1,1),B(1,2),B(1,3)
C         dBXx    dBXy    dBYy    dBXz    dBYz
     >    ,DB(1,1),DB(2,1),DB(2,2),DB(3,1),DB(3,2)
C         d2BXxx     d2BXxy     d2BXyy     d2BYyy  
     >    ,DDB(1,1,1),DDB(2,1,1),DDB(2,2,1),DDB(2,2,2)
C         d2BXxz     d2BXyz     d2BYyz  
c     >    ,DDB(3,1,1),DDB(3,2,1),DDB(3,2,2)
C         d3BXdy2z     d4BXdyz3      d3BZdz3  
c     >    ,D3BX(3,2,2),D4BX(3,3,3,2),D3BZ(3,3,3)
        ENDIF
      ENDIF
 
      RETURN
      END
