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
      SUBROUTINE TWISS(READAT,*) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL READAT
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM

      CHARACTER*9   HMS
      
      DIMENSION T(6,6,6)
      SAVE      T
      DIMENSION TX3(5,6) , TX4(5,6)
      SAVE      TX3,       TX4

      DIMENSION REF(6)
      SAVE REF
      DIMENSION RREF(6,6), RPLUS(6,6), RMINUS(6,6) 
      SAVE RREF, RPLUS, RMINUS
      DIMENSION F0REF(6,6), F0P(6,6), F0M(6,6) 
      SAVE F0REF, F0P, F0M
      SAVE YNUREF, ZNUREF, YNUP, ZNUP, YNUM, ZNUM
      SAVE DNUYDP, DNUZDP, DNUYDY, DNUZDY, DNUZDZ
      SAVE UYREF, UYP

      SAVE IORD
       
      DATA KWRI6 / 1 /
      SAVE KWRI6

      KTWISS = A(NOEL,1)
      FACP = A(NOEL,2)
      FACA = A(NOEL,3)
      FAC1 = 1.D0/FACA

C--------- 2 more runs to get matrices around +/-dp/p chromatic closed orbits;
C          1 other run to get matrices with +/-dY amplitude 
C          1 other run to get matrices with +/-dZ amplitude 
          NRBLT = 4

      IF(KTWISS .EQ. 0) RETURN

      IF(KWRI6 .NE. 0) THEN
        CALL TIME2(HMS)
C        WRITE(6,FMT='(A,I7,3X,A)') ' IPASS=',IPASS,HMS
C        WRITE(6,FMT='(1H+,A1,A31,I7,A2,I7,A1,A11,A9,$)') CHAR(13),
C     >    ' # of pass achieved/Requested :',IPASS,' /', NRBLT+1,
C     >         ',','  at time :', HMS
C        IF(IPASS.EQ.1) WRITE(6,FMT='(/)')
C        CALL FLUSH(6)
      ENDIF
 
      IF( IPASS .LT. NRBLT ) THEN
        LUN=ABS(NRES) 
        IF(LUN.GT.0) 
     >     WRITE(LUN,100) IPASS, NRBLT+1
 100       FORMAT(/,30X,'  -----  TWISS procedure  -----',//,5X,
     >     'End of pass # ',I1,'/',I1,
     >     ' through the optical structure',/)
 
        IF(IPASS .EQ. 1) THEN
C--------- 2nd pass through structure will follow

C--------- Switch off print into zgoubi.res : 
          ANOEL2 = 0.1D0
          KWRIT = NINT(ANOEL2)
C--------- Switch on print to standard output :
          KWRI6=NINT(ANOEL2-KWRIT)*10

          IF(KWRIT .EQ. 0) THEN
C---------- inihibe WRITE si KWRIT=0 et plus de 2 passages
C            IF(NRES .GT. 0) NRES =-NRES
          ENDIF
          READAT = .FALSE.

          IF(KOBJ .EQ. 5) THEN
            IORD=1
          ELSEIF(KOBJ .EQ. 6) THEN
            IORD=2
          ENDIF 

          IF    (IORD .EQ. 1) THEN
            CALL REFER(1,1,0,1,4,5)
            CALL MAT1(RREF,T,1)
            CALL REFER(2,1,0,1,4,5)
          ELSEIF(IORD .EQ. 2) THEN
            CALL REFER(1,2,0,1,6,7)
            CALL MAT2(RREF,T,TX3,TX4)
            CALL REFER(2,2,0,1,6,7)
          ENDIF
          CALL MKSA(IORD,RREF,T,TX3,TX4)
          CALL MATIMP(RREF)
          CALL TUNES(RREF,F0REF,1,IERY,IERZ,.TRUE.,
     >                                          YNUREF,ZNUREF,CMUY,CMUZ)
 
          NUML = 1
          REF(1) = A(NUML,30)
          REF(2) = A(NUML,31)
          REF(3) = A(NUML,32)
          REF(4) = A(NUML,33)
          REF(5) = A(NUML,34)
          REF(6) = A(NUML,35)

C--------- Reset reference coordinates for OBJECT sampling : p -> p-dp
          NUML = 1
          FAP25 = FACP * A(NUML,25)
          A(NUML,35) =  REF(6) -  FAP25
          A(NUML,30) =  REF(1) -  F0REF(1,6) * FAP25
          A(NUML,31) =  REF(2) -  F0REF(2,6) * FAP25
          A(NUML,32) =  REF(3) -  F0REF(3,6) * FAP25
          A(NUML,33) =  REF(4) -  F0REF(4,6) * FAP25
          
          IPASS=IPASS+1
          NOEL=0 
          CALL SCUMS(0.D0)
          RETURN 1
 
        ELSEIF(IPASS .EQ. 2) THEN
C------- 3rd pass through structure will follow

        LUN=ABS(NRES)
        IF(LUN.GT.0) WRITE(LUN,100) IPASS
 
          IF    (IORD .EQ. 1) THEN
            CALL REFER(1,1,0,1,4,5)
            CALL MAT1(RMINUS,T,1)
            CALL REFER(2,1,0,1,4,5)
          ELSEIF(IORD .EQ. 2) THEN
            CALL REFER(1,2,0,1,6,7)
            CALL MAT2(RMINUS,T,TX3,TX4)
            CALL REFER(2,2,0,1,6,7)
          ENDIF
          CALL MKSA(IORD,RMINUS,T,TX3,TX4)
          CALL MATIMP(RMINUS)
          CALL TUNES(RMINUS,F0M,1,IERY,IERZ,.TRUE.,
     >                                             YNUM,ZNUM,CMUY,CMUZ)
 
C--------- Reset reference coordinates for OBJECT sampling : p -> p+dp
          NUML = 1
          FAP25 = FACP * A(NUML,25)
          A(NUML,35) =  REF(6) +  FAP25
          A(NUML,30) =  REF(1) +  F0REF(1,6) * FAP25
          A(NUML,31) =  REF(2) +  F0REF(2,6) * FAP25
          A(NUML,32) =  REF(3) +  F0REF(3,6) * FAP25
          A(NUML,33) =  REF(4) +  F0REF(4,6) * FAP25
          
          IPASS=IPASS+1
          NOEL=0 
          CALL SCUMS(0.D0)
          RETURN 1
 
        ELSEIF(IPASS .EQ. 3) THEN
C------- Chromatic calculations done

C--------- reactivate WRITE for printing results 
          IF(NRES.LT.0) NRES=-NRES

          IF    (IORD .EQ. 1) THEN
            CALL REFER(1,1,0,1,4,5)
            CALL MAT1(RPLUS,T,1)
            CALL REFER(2,1,0,1,4,5)
          ELSEIF(IORD .EQ. 2) THEN
            CALL REFER(1,2,0,1,6,7)
            CALL MAT2(RPLUS,T,TX3,TX4)
            CALL REFER(2,2,0,1,6,7)
          ENDIF
          CALL MKSA(IORD,RPLUS,T,TX3,TX4)
          CALL MATIMP(RPLUS)
          CALL TUNES(RPLUS,F0P,1,IERY,IERZ,.TRUE.,
     >                                            YNUP,ZNUP,CMUY,CMUZ)

C Momentum detuning
          NUML = 1
          DNUYDP = (YNUP-YNUM)/2.D0/A(NUML,25)
          DNUZDP = (ZNUP-ZNUM)/2.D0/A(NUML,25)
        
C--------- Amplitude detuning calculations follow
          NUML = 1
          A(NUML,35) =  REF(6)
          A(NUML,30) =  REF(1)
          A(NUML,31) =  REF(2)
          A(NUML,32) =  REF(3)
          A(NUML,33) =  REF(4)
          A(NUML,34) =  REF(5)
C--------- Reset reference coordinates for OBJECT sampling : y -> y+dy
          A(NUML,20) =  FACA * A(NUML,20)    
          A(NUML,21) =  FACA * A(NUML,21)    

          IPASS=IPASS+1
          NOEL=0 
          CALL SCUMS(0.D0)
          RETURN 1
 
        ENDIF

      ELSEIF(IPASS .EQ. NRBLT) THEN

          LUN=ABS(NRES)
          IF(LUN.GT.0) WRITE(LUN,100) IPASS
 
          IF    (IORD .EQ. 1) THEN
            CALL REFER(1,1,0,1,4,5)
            CALL MAT1(RPLUS,T,1)
            CALL REFER(2,1,0,1,4,5)
          ELSEIF(IORD .EQ. 2) THEN
            CALL REFER(1,2,0,1,6,7)
            CALL MAT2(RPLUS,T,TX3,TX4)
            CALL REFER(2,2,0,1,6,7)
          ENDIF
          CALL MKSA(IORD,RPLUS,T,TX3,TX4)
          CALL MATIMP(RPLUS)
          CALL TUNES(RPLUS,F0P,1,IERY,IERZ,.TRUE.,
     >                                            YNUP,ZNUP,CMUY,CMUZ)
 
C  Amplitude detuning, dY effects
          NUML = 1
          Y2 = FAC1*A(NUML,20)-A(NUML,30)
          YYP = Y2
          Y2 = Y2*Y2 
          YP2 = FAC1*A(NUML,21)-A(NUML,31)
          YYP = YYP*YP2
          YP2 = YP2*YP2 
          UYREF = F0REF(2,2)*Y2+2.D0*(-F0REF(2,1))*YYP+F0REF(1,1)*YP2 
          YY2 = A(NUML,20)-A(NUML,30)
          YYYP = YY2
          YY2 = YY2*YY2 
          YYP2 = A(NUML,21)-A(NUML,31)
          YYYP = YYYP*YYP2
          YYP2 = YYP2*YYP2 
          UYP = F0P(2,2)*YY2 + 2.D0*(-F0P(2,1))*YYYP + F0P(1,1)*YYP2 
          DNUYDY=(YNUP-YNUREF)/(UYP-UYREF)
          DNUZDY=(ZNUP-ZNUREF)/(UYP-UYREF)
        
C--------- Reset reference coordinates for OBJECT sampling : z -> z+dz
          NUML = 1
          A(NUML,20) =  FAC1 * A(NUML,20) 
          A(NUML,21) =  FAC1 * A(NUML,21) 
          A(NUML,22) =  FACA * A(NUML,22)
          A(NUML,23) =  FACA * A(NUML,23)
    
          IPASS=IPASS+1
          NOEL=0 
          CALL SCUMS(0.D0)
          RETURN 1

      ELSEIF(IPASS .GT. NRBLT) THEN
C------- Amplitude calculations done

C------- reactivate WRITE
        IF(NRES.LT.0) NRES=-NRES
        IF(NRES.GT.0) THEN
          WRITE(NRES,101) IPASS
 101      FORMAT(/,25X,' ****  End  of  ''TWISS''  procedure  ****',//
     >     ,5X,' There  has  been ',I10,
     >          '  pass  through  the  optical  structure ',/)
 
        ENDIF

        IF    (IORD .EQ. 1) THEN
          CALL REFER(1,1,0,1,4,5)
          CALL MAT1(RPLUS,T,1)
          CALL REFER(2,1,0,1,4,5)
        ELSEIF(IORD .EQ. 2) THEN
          CALL REFER(1,2,0,1,6,7)
          CALL MAT2(RPLUS,T,TX3,TX4)
          CALL REFER(2,2,0,1,6,7)
        ENDIF
        CALL MKSA(IORD,RPLUS,T,TX3,TX4)
        CALL MATIMP(RPLUS)
        CALL TUNES(RPLUS,F0P,1,IERY,IERZ,.TRUE.,
     >                                          YNUP,ZNUP,CMUY,CMUZ)

C Amplitude detuning, dZ effects
        NUML = 1
          Z2 = FAC1*A(NUML,22)-A(NUML,32)
          ZZP = Z2
          Z2 = Z2*Z2 
          ZP2 = FAC1*A(NUML,23)-A(NUML,33)
          ZZP = ZZP*ZP2
          ZP2 = ZP2*ZP2 
          UZREF = F0REF(4,4)*Z2+2.D0*(-F0REF(4,3))*ZZP+F0REF(3,3)*ZP2 
          ZZ2 = A(NUML,22)-A(NUML,32)
          ZZZP = ZZ2
          ZZ2 = ZZ2*ZZ2 
          ZZP2 = A(NUML,23)-A(NUML,33)
          ZZZP = ZZZP*ZZP2
          ZZP2 = ZZP2*ZZP2 
          UZP = F0P(4,4)*ZZ2 + 2.D0*(-F0P(4,3))*ZZZP + F0P(3,3)*ZZP2 

        DNUZDZ=(ZNUP-ZNUREF)/(UZP-UZREF)
        DNUYDZ=(YNUP-YNUREF)/(UZP-UZREF)
        
        WRITE(NRES,FMT='(/,34X,'' Chromaticities : '',//, 
     >  30X,''dNu_y / dp/p = '',G14.8,/, 
     >  30X,''dNu_z / dp/p = '',G14.8)') DNUYDP, DNUZDP

        WRITE(NRES,FMT='(/,38X,'' Amplitude  detunings : '',//, 
     >  42X,''/ dEps_y/pi       / dEps_z/pi'',/, 
     >  30X,''dNu_y'',7X,2(G14.8,3X),/, 
     >  30X,''dNu_z'',7X,2(G14.8,3X), //, 
     >  20X,''Nu_yRef = '',G14.8,'', Nu_zRef = '',G14.8, / 
     >  20X,''Nu_yP = '',G14.8,'',   Nu_zP = '',G14.8, / 
     >  20X,''Eps_yRef/pi = '',G14.8,'',   Eps_zRef/pi = '',G14.8, / 
     >  20X,''Eps_y+/pi = '',G14.8,'',   Eps_z+ = '',G14.8)')
     >      DNUYDY, DNUYDZ, DNUZDY, DNUZDZ, 
     >      YNUREF,ZNUREF,
     >      YNUP,  ZNUP,
     >      UYREF, UZREF,
     >      UYP,   UZP

      ENDIF
      RETURN
      END
