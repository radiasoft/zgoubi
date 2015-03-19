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
      SUBROUTINE IMPDEV
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      INCLUDE "C.CHAVE_2.H"     ! COMMON/CHAVE/ B(5,3),V(5,3),E(5,3)
      INCLUDE "C.CNTRLB.H"     ! COMMON/CNTRLB/ XSTP,DIVB,ALAPL(4),RTN(3)
      COMMON/DDBXYZ/ DB(3,3),DDB(3,3,3)
      INCLUDE "C.D3B_2.H"     ! COMMON/D3BXYZ/ D3BX(3,3,3), D3BY(3,3,3), D3BZ(3,3,3)
      INCLUDE "C.D4B.H"     ! COMMON/D4BXYZ/ D4BX(3,3,3,3) ,D4BY(3,3,3,3) ,D4BZ(3,3,3,3)
      COMMON/DDEXYZ/ DE(3,3),DDE(3,3,3)
      INCLUDE "C.INTEG.H"     ! COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
C      LOGICAL ZSYM
      INCLUDE "C.TYPFLD.H"     ! COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      INCLUDE "C.PTICUL.H"     ! COMMON/PTICUL/ AM,Q,G,TO
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      INCLUDE "C.SPIN.H"     ! COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)
      COMMON/TRAJ/ Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT
      COMMON/VITES/ U(6,3),DQBR(6),DDT(6)

      LOGICAL FIRST
      LOGICAL IDLUNI, OKIMP

      DATA FIRST, OKIMP / .TRUE., .FALSE. /
 
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
          WRITE(NRES,105) IT,QBR/(Q*BORO),Y,T,Z,P,X,SAR
 105      FORMAT(/,' SPGM   IMPDEV',5X,'TRAJ ',I3
     >    ,/,' D,Y,T,Z,P,X,S :',/,' ',1P,7E16.8)
 
          IF((KFLD .EQ. LC .OR. KFLD .EQ. ML).AND.BRI.NE.0D0) THEN
 
            WRITE(NRES,104) E(1,1)/BRI,E(1,2)/BRI,E(1,3)/BRI
     >      ,DE(1,1)/BRI,DE(2,1)/BRI,DE(3,1)/BRI,DE(2,2)/BRI,DE(3,2)/BRI
     >      ,DDE(1,1,1)/BRI,DDE(2,1,1)/BRI,DDE(3,1,1)/BRI
     >      ,DDE(2,2,1)/BRI,DDE(3,2,1)/BRI,DDE(3,3,1)/BRI
 104        FORMAT(' EX, EY, EZ (MV/cm):',1P,3E16.8
     >          ,/,' E'' :',5E16.8
     >          ,/,' E'''':',6E16.8)
          ENDIF
 
          IF((KFLD .EQ. MG .OR. KFLD .EQ. ML).AND.BRI.NE.0D0) THEN
 
            WRITE(NRES,106) B(1,1)/BRI,B(1,2)/BRI,B(1,3)/BRI
     >      ,DB(1,1)/BRI,DB(2,1)/BRI,DB(3,1)/BRI,DB(2,2)/BRI,DB(3,2)/BRI
     >      ,DDB(1,1,1)/BRI,DDB(2,1,1)/BRI,DDB(3,1,1)/BRI
     >      ,DDB(2,2,1)/BRI,DDB(3,2,1)/BRI,DDB(3,3,1)/BRI
 106        FORMAT(' BX, BY, BZ (kG)   :',1P,3E16.8
     >          ,/,' B'' :',5E16.8
     >          ,/,' B'''':',6E16.8)
          ENDIF
 
          WRITE(NRES,FMT='('' SF_X, _Y, _Z  |S| : '',1P,4(1X,E12.4))')
     >    (SF(J,IT), J=1,4)

        ELSEIF(LST .EQ. 7) THEN

          IF(FIRST) THEN

            FIRST = .FALSE. 
            IF(IDLUNI(
     >                LUN)) THEN
              OPEN(UNIT=LUN,FILE='zgoubi.impdev.out')
              OKIMP = .TRUE.
            ELSE
              OKIMP = .FALSE.
              WRITE(6,*)' SBR impdev. Could not opem zgoubi.impdev.out.'
            ENDIF
          
            IF(OKIMP) WRITE(LUN,FMT='(A,/,A)')
     >     '#1, 2, 3, 4, 5, 6, 7, 8-10,11-19, 20-46  , 47-127 , 128-370  
     >     ,371-373, 374-382, 383-409,  410,  411   ',
     >     '#D, Y, T, Z, P, X, S, B_J, DB_IJ, DDB_IJK, D3B_IJK, D4B_IJKL
     >     , E_J  ,  DE_IJ,   DDE_IJK,  IT,   BRho(kG.cm)  '
          ENDIF
 
          IF(OKIMP) THEN
            WRITE(LUN,109) QBR/(Q*BORO), Y, T, Z, P, X, SAR
     >      ,( B(1,J)/BRI,J=1,3)
     >      ,(( DB(I,J)/BRI  ,I=1,3),J=1,3)
     >      ,((( DDB(I,J,K)/BRI  ,I=1,3), J=1,3), K=1,3)
     >      ,((( D3BX(I,J,K)/BRI ,I=1,3) ,J=1,3) ,K=1,3)
     >      ,((( D3BY(I,J,K)/BRI ,I=1,3) ,J=1,3) ,K=1,3)
     >      ,((( D3BZ(I,J,K)/BRI ,I=1,3) ,J=1,3) ,K=1,3)
     >      ,((((D4BX(I,J,K,L)/BRI ,I=1,3),J=1,3),K=1,3),L=1,3)
     >      ,((((D4BY(I,J,K,L)/BRI ,I=1,3),J=1,3),K=1,3),L=1,3)
     >      ,((((D4BZ(I,J,K,L)/BRI ,I=1,3),J=1,3),K=1,3),L=1,3)
     >      ,( E(1,J)/BRI,J=1,3)
     >      ,(( DE(I,J)/BRI  ,I=1,3),J=1,3)
     >      ,((( DDE(I,J,K)/BRI ,I=1,3),J=1,3),K=1,3)
     >      ,IT,QBR/Q
 109        FORMAT( 1P,7E17.8
     >      , 3(1X,E13.4)
     >      , 9(1X,E13.4)
     >      ,27(1X,E13.4)
     >      ,27(1X,E13.4)
     >      ,27(1X,E13.4)
     >      ,27(1X,E13.4)
     >      ,81(1X,E13.4)
     >      ,81(1X,E13.4)
     >      ,81(1X,E13.4)
     >      , 3(1X,E13.4)
     >      , 9(1X,E13.4)
     >      ,27(1X,E13.4)
     >      ,   1X,I6
     >      ,  (1X,E13.4)
     >      )
          ENDIF
C          WRITE(NRES,*)
C          WRITE(NRES,105) IT,QBR/(Q*BORO),Y,T,Z,P,X,SAR
C          WRITE(NRES,109) ( J,  B(1,J),J=1,3)
C 109      FORMAT('    B(J) ',3(1X,I1,1X,1P,E11.3))
C          WRITE(NRES,100) (( I,J,  DB(I,J)  ,I=1,3),J=1,3)
C 100      FORMAT('     DB(I,J) ',6(1X,2I1,1X,1P,E11.3))
C          WRITE(NRES,101) ((( I,J,K,DDB(I,J,K) ,I=1,3),J=1,3),K=1,3)
C 101      FORMAT('  DDB(I,J,K) ',6(1X,3I1,1X,1P,E11.3))
C          WRITE(NRES,102)((( I,J,K,D3BX(I,J,K) ,I=1,3) ,J=1,3) ,K=1,3)
C          WRITE(NRES,102)((( I,J,K,D3BY(I,J,K) ,I=1,3) ,J=1,3) ,K=1,3)
C          WRITE(NRES,102)((( I,J,K,D3BZ(I,J,K) ,I=1,3) ,J=1,3) ,K=1,3)
C 102      FORMAT('  D3B(I,J,K) ',6(1X,3I1,1X,1P,E11.3))
C          WRITE(NRES,103) ((((D4BX(I,J,K,L),I=1,3),J=1,3),K=1,3),L=1,3)
C          WRITE(NRES,103) ((((D4BY(I,J,K,L),I=1,3),J=1,3),K=1,3),L=1,3)
C          WRITE(NRES,103) ((((D4BZ(I,J,K,L),I=1,3),J=1,3),K=1,3),L=1,3)
C 103      FORMAT('  D4B(I,J,K,L) ',1P,9E11.3)

C          WRITE(NRES,113) ( J,  E(1,J),J=1,3)
C 113      FORMAT('    E(J) ',3(1X,I1,1X,1P,E11.3))
C          WRITE(NRES,107) (( I,J,  DE(I,J)  ,I=1,3),J=1,3)
C 107      FORMAT('     DE(I,J) ',6(1X,2I1,1X,1P,E11.3))
C          WRITE(NRES,108) ((( I,J,K,DDE(I,J,K) ,I=1,3),J=1,3),K=1,3)
C 108      FORMAT('  DDE(I,J,K) ',6(1X,3I1,1X,1P,E11.3))

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
