C  ZGOUBI, a program for computing the trajectories of charged particles
C  in electric and magnetic fields
C  Copyright (C) 1988-2007  Fran�ois M�ot
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
C  Fran�ois M�ot <fmeot@bnl.gov>
C  Brookhaven National Laboratory    
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  -------
      SUBROUTINE TWISS(LUN,
     >                 KOPTCS, READAT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL READAT
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      PARAMETER (LBLSIZ=10)
      CHARACTER(LBLSIZ) LABEL
      COMMON /LABEL/ LABEL(MXL,2)
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      COMMON/PTICUL/ AM,Q,G,TO
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      
      DIMENSION T(6,6,6)
      SAVE      T
      DIMENSION TX3(5,6) , TX4(5,6)
      SAVE      TX3,       TX4

      DIMENSION REF(6), PATHL(3)
      SAVE REF
      DIMENSION RREF(6,6), RPLUS(6,6), RMINUS(6,6) 
      SAVE RREF, RPLUS, RMINUS
      DIMENSION F0REF(6,6), F0P(6,6), F0M(6,6) 
      SAVE F0REF, F0P, F0M
      SAVE YNUREF, ZNUREF, YNUP, ZNUP, YNUM, ZNUM
      SAVE DNUYDP, DNUZDP, DNUYDY, DNUZDY, DNUZDZ
      SAVE UYREF, UYP

      SAVE IORD
      SAVE KWRI6

      SAVE ALPHA

C F2 contains seven 6-vectorss (2nd index), from ipass-6 (f2(1,*)) to ipass (f2(7,*))
      DIMENSION F2(7,6), XYS(6,12), KAUX(6)   !!, SM(6,6)
      SAVE F2
  
      PARAMETER (KSIZ=10)
      CHARACTER*(KSIZ)  KLOBJ
      LOGICAL IDLUNI
      INTEGER DEBSTR, FINSTR

      LOGICAL dolast
      save dolast
      
      dimension rturn(4,4)

      DATA KWRI6 / 1 /
      DATA dolast / .true. /

      KTW = INT(A(NOEL,1))
      KTW2 = INT( 10.D0*(A(NOEL,1) -KTW) )
      FACP = A(NOEL,2)
      FACA = A(NOEL,3)

      FAC1 = 1.D0/FACA

C KTW=1 : 1 pass, equivalent to matrix[1,11]
C KTW=2 : 2 more passes to get matrices around +/-dp/p chromatic closed orbits;
C KTW=3 : 1 more pass to get matrices with +/-dY amplitude and 
C        +1 more to get matrices with +/-dZ amplitude. 
C KTW .ge. 99 : compute linear functions from multi-turn tracking - installation not completed.

      IF(KTW.GE.99) GOTO 20
      IF(KTW.GE.1) GOTO 10
      RETURN

 10   CONTINUE
Compute optical functions, tunes, chromaticity, anharmonicities, from a few passes
C of 11 particles (based on MATRIX)
 
      if(ipass .eq. 4 .and. ktw .eq. 2) goto 222

      KLOBJ = 'OBJET'
      CALL GETNOL(KLOBJ,
     >                  NLOBJ)

      IF(IPASS .EQ. 1) THEN
C Compute periodic beta from first pass.
C 2nd pass through structure will follow iff KTW>1.

        dolast = .true.

        ISIGN = NRES/ABS(NRES)
        NRES = ISIGN*NRES
        WRITE(NRES,FMT='(/,25X,
     >    '' ****  End  of  pass #'',I1
     >  ,'' of TWISS  procedure  ****'',/)') IPASS

C------- Switch off print into zgoubi.res : 
        ANOEL2 = 0.1D0
        KWRIT = NINT(ANOEL2)
C------- Switch on print to standard output :
        KWRI6=NINT(ANOEL2-KWRIT)*10

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

        NRES = ISIGN*NRES

        IF(KTW.GE.2) THEN
  
          IF(KTW2 .EQ. 0) THEN
            IF(NRES .GT. 0) NRES =-NRES
          ENDIF

          CALL REFER1(
     >                PATHL(1)) 
C          NLOBJ = 1
          REF(1) = A(NLOBJ,30)
          REF(2) = A(NLOBJ,31)
          REF(3) = A(NLOBJ,32)
          REF(4) = A(NLOBJ,33)
          REF(5) = A(NLOBJ,34)
          REF(6) = A(NLOBJ,35)
C------- Reset reference coordinates for OBJECT sampling : p -> p-dp
C          NLOBJ = 1
          FAP25 = FACP * A(NLOBJ,25)
          A(NLOBJ,35) =  REF(6) -  FAP25
          A(NLOBJ,30) =  REF(1) -  1.D2*F0REF(1,6) * FAP25
          A(NLOBJ,31) =  REF(2) -  1.D3*F0REF(2,6) * FAP25
          A(NLOBJ,32) =  REF(3) -  1.D2*F0REF(3,6) * FAP25
          A(NLOBJ,33) =  REF(4) -  1.D3*F0REF(4,6) * FAP25
          
          IPASS=IPASS+1
          NOEL=0 
          CALL SCUMS(0.D0)

        ENDIF

        RETURN
 
      ELSEIF(IPASS .EQ. 2) THEN
C----- 3rd pass through structure will follow
        ISIGN = NRES/ABS(NRES)
        NRES = ISIGN*NRES
        IF(NRES.GT.0) WRITE(NRES,FMT='(/,25X,
     >    '' ****  End  of  pass #'',I1
     >  ,'' of TWISS  procedure  ****'',/)') IPASS

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
        CALL REFER1(
     >              PATHL(2)) 

        NRES = ISIGN*NRES
 
C------- Reset reference coordinates for OBJECT sampling : p -> p+dp
C        NLOBJ = 1
        FAP25 = FACP * A(NLOBJ,25)
        A(NLOBJ,35) =  REF(6) +  FAP25
        A(NLOBJ,30) =  REF(1) +  1.D2*F0REF(1,6) * FAP25
        A(NLOBJ,31) =  REF(2) +  1.D3*F0REF(2,6) * FAP25
        A(NLOBJ,32) =  REF(3) +  1.D2*F0REF(3,6) * FAP25
        A(NLOBJ,33) =  REF(4) +  1.D3*F0REF(4,6) * FAP25
          
        IPASS=IPASS+1
        NOEL=0 
        CALL SCUMS(0.D0)

        RETURN
 
      ELSEIF(IPASS .EQ. 3) THEN
C----- Chromatic tracking completed

        ISIGN = NRES/ABS(NRES)
        NRES = ISIGN*NRES
        IF(NRES.GT.0) WRITE(NRES,FMT='(/,25X,
     >    '' ****  End  of  pass #'',I1
     >  ,'' of TWISS  procedure  ****'',/)') IPASS

C------- reactivate WRITE for printing results 

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
        CALL REFER1(
     >              PATHL(3)) 

        NRES = ISIGN*NRES

C Momentum compaction
        ALPHA=( (PATHL(3)-PATHL(2))/PATHL(1) ) / (2.D0 * A(NLOBJ,25))

C Momentum detuning
C        NLOBJ = 1
        DNUYDP = (YNUP-YNUM)/2.D0/A(NLOBJ,25)
        DNUZDP = (ZNUP-ZNUM)/2.D0/A(NLOBJ,25)
        
        IF(KTW.GE.3) THEN
C------- Amplitude detuning tracking & calculations follow
C          NLOBJ = 1
          A(NLOBJ,35) =  REF(6)
          A(NLOBJ,30) =  REF(1)
          A(NLOBJ,31) =  REF(2)
          A(NLOBJ,32) =  REF(3)
          A(NLOBJ,33) =  REF(4)
          A(NLOBJ,34) =  REF(5)
C------- Reset reference coordinates for OBJECT sampling : y -> y+dy
          A(NLOBJ,20) =  FACA * A(NLOBJ,20)    
          A(NLOBJ,21) =  FACA * A(NLOBJ,21)    

          IPASS=IPASS+1
          NOEL=0 
          CALL SCUMS(0.D0)

          RETURN

        ELSE

          goto 222

        ENDIF

C        ENDIF

      ELSEIF(IPASS .EQ. 4) THEN

        ISIGN = NRES/ABS(NRES)
        NRES = ISIGN*NRES
        IF(NRES.GT.0) WRITE(NRES,FMT='(/,25X,
     >    '' ****  End  of  pass #'',I1
     >  ,'' of TWISS  procedure  ****'',/)') IPASS

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

          NRES = ISIGN*NRES
 
C  Amplitude detuning, dY effects
C          NLOBJ = 1
          Y2 = FAC1*A(NLOBJ,20)-A(NLOBJ,30)
          YYP = Y2
          Y2 = Y2*Y2 
          YP2 = FAC1*A(NLOBJ,21)-A(NLOBJ,31)
          YYP = YYP*YP2
          YP2 = YP2*YP2 
          UYREF = F0REF(2,2)/1.D2*Y2+2.D0*(-F0REF(2,1))*YYP+
     >             F0REF(1,1)*1.D2*YP2 
          YY2 = A(NLOBJ,20)-A(NLOBJ,30)
          YYYP = YY2
          YY2 = YY2*YY2 
          YYP2 = A(NLOBJ,21)-A(NLOBJ,31)
          YYYP = YYYP*YYP2
          YYP2 = YYP2*YYP2 
          UYP = F0P(2,2)/1.D2*YY2 + 2.D0*(-F0P(2,1))*YYYP + 
     >             F0P(1,1)*1.D2*YYP2 
          DNUYDY=(YNUP-YNUREF)/(UYP-UYREF)
          DNUZDY=(ZNUP-ZNUREF)/(UYP-UYREF)
        
C--------- Reset reference coordinates for OBJECT sampling : z -> z+dz
C          NLOBJ = 1
          A(NLOBJ,20) =  FAC1 * A(NLOBJ,20) 
          A(NLOBJ,21) =  FAC1 * A(NLOBJ,21) 
          A(NLOBJ,22) =  FACA * A(NLOBJ,22)
          A(NLOBJ,23) =  FACA * A(NLOBJ,23)
    
          IPASS=IPASS+1
          NOEL=0 
          CALL SCUMS(0.D0)
          RETURN 

C      ELSEIF(IPASS .GT. NRBLT) THEN
      ELSEIF(IPASS .EQ. 5) THEN
C------- Amplitude tracking completed

          ISIGN = NRES/ABS(NRES)
          NRES = ISIGN*NRES
          IF(NRES.GT.0) WRITE(NRES,FMT='(/,25X,
     >      '' ****  End  of  pass #'',I1
     >    ,'' of TWISS  procedure  ****'',/)') IPASS

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

          NRES = ISIGN*NRES

C Amplitude detuning, dZ effects
C        NLOBJ = 1
          Z2 = FAC1*A(NLOBJ,22)-A(NLOBJ,32)
          ZZP = Z2
          Z2 = Z2*Z2 
          ZP2 = FAC1*A(NLOBJ,23)-A(NLOBJ,33)
          ZZP = ZZP*ZP2
          ZP2 = ZP2*ZP2 
          UZREF = F0REF(4,4)/1.D2*Z2+2.D0*(-F0REF(4,3))*ZZP+
     >       F0REF(3,3)*1.D2*ZP2 
          ZZ2 = A(NLOBJ,22)-A(NLOBJ,32)
          ZZZP = ZZ2
          ZZ2 = ZZ2*ZZ2 
          ZZP2 = A(NLOBJ,23)-A(NLOBJ,33)
          ZZZP = ZZZP*ZZP2
          ZZP2 = ZZP2*ZZP2 
          UZP = F0P(4,4)/1.D2*ZZ2 + 2.D0*(-F0P(4,3))*ZZZP + 
     >          F0P(3,3)*1.D2*ZZP2 

             DNUZDZ=(ZNUP-ZNUREF)/(UZP-UZREF)
             DNUYDZ=(YNUP-YNUREF)/(UZP-UZREF)
        
      ENDIF

 222  continue

C-------------------------------------------------------------------
C-------------------------------------------------------------------
C Now make a last pass to get optical functions at all elements

      IF(dolast) THEN
        dolast = .false.

C----- So to print into zgoubi.OPTICS.out
        KOPTCS = 1

        ISIGN = NRES/ABS(NRES)
        NRES = ISIGN*NRES
        IF(NRES.GT.0) WRITE(NRES,FMT='(/,25X,
     >    '' ****  End  of  pass #'',I1
     >  ,'' of TWISS  procedure  ****'',/)') IPASS

        NRES = ISIGN*NRES
 
C------- Reset reference coordinates for OBJECT sampling : p -> p+dp
C        NLOBJ = 1
        FAP25 = FACP * A(NLOBJ,25)
        A(NLOBJ,35) =  REF(6) 
        A(NLOBJ,30) =  REF(1) 
        A(NLOBJ,31) =  REF(2) 
        A(NLOBJ,32) =  REF(3) 
        A(NLOBJ,33) =  REF(4) 
          
        IPASS=IPASS+1
        NOEL=0 
        CALL SCUMS(0.D0)

C        IF(LABEL(NOEL,1)(DEBSTR(LABEL(NOEL,1)):FINSTR(LABEL(NOEL,1))) 
C     >                                   .EQ.  'PRINT') THEN
C          IF(IDLUNI(
C     >              LUN)) THEN
C            OPEN(UNIT=LUN,FILE='zgoubi.TWISS.Out',ERR=96)
C          ELSE
C            GOTO 96
C          ENDIF

C--------- P0, AM  are  in  MEV/c, /c^2
          PREF = BORO*CL9*Q*DPREF
          if(am.le.1d-8) am = AMPROT
          Energy = sqrt(PREF*PREF + AM*AM)
          write(LUN,50)'@ NAME             %05s "TWISS"'
          write(LUN,50)'@ TYPE             %05s "TWISS"'
          write(LUN,50)'@ SEQUENCE         %04s "RING"'
          write(LUN,50)'@ PARTICLE         %00s ""'
          write(LUN,51)'@ MASS             %le', am/1.d3
          write(LUN,52)'@ CHARGE           %le', int(q)
          write(LUN,51)'@ ENERGY           %le', Energy/1.d3
          write(LUN,53)'@ PC               %le', PREF/1.D3,' 
     >                                             B.rho ',BORO/1.d3
          write(LUN,51)'@ GAMMA            %le', Energy/am
          write(LUN,50)'@ KBUNCH           %le                   1'
          write(LUN,50)'@ BCURRENT         %le                   0'
          write(LUN,50)'@ SIGE             %le                   0'
          write(LUN,50)'@ SIGT             %le                   1'
          write(LUN,50)'@ NPART            %le                   0'
          write(LUN,50)'@ EX               %le                   1'
          write(LUN,50)'@ EY               %le                   1'
          write(LUN,50)'@ ET               %le                   1'
          write(LUN,51)'@ LENGTH           %le', PathL(1)
          write(LUN,51)'@ ALFA             %le', Alpha
          write(LUN,50)'@ ORBIT5           %le                  -0'
          write(LUN,51)'@ GAMMATR          %le', sqrt(1.d0/Alpha)
          write(LUN,51)'@ Q1               %le', YNUREF + 8.d0
          write(LUN,51)'@ Q2               %le', ZNUREF + 8.d0
          write(LUN,51)'@ DQ1              %le', DNUYDP
          write(LUN,51)'@ DQ2              %le', DNUZDP
          write(LUN,51)'@ DXMAX            %le', 9999.
          write(LUN,51)'@ DYMAX            %le', 9999.
          write(LUN,50)'@ XCOMAX           %le                   0'
          write(LUN,50)'@ YCOMAX           %le                   0'
          write(LUN,51)'@ BETXMAX          %le', 9999.
          write(LUN,51)'@ BETYMAX          %le', 9999.
          write(LUN,50)'@ XCORMS           %le                   0'
          write(LUN,50)'@ YCORMS           %le                   0'
          write(LUN,51)'@ DXRMS            %le', 9999.
          write(LUN,51)'@ DYRMS            %le', 9999.
          write(LUN,50)'@ DELTAP           %le                   0'
          write(LUN,50)'@ SYNCH_1          %le                   0'
          write(LUN,50)'@ SYNCH_2          %le                   0'
          write(LUN,50)'@ SYNCH_3          %le                   0'
          write(LUN,50)'@ SYNCH_4          %le                   0'
          write(LUN,50)'@ SYNCH_5          %le                   0'
          write(LUN,50)'@ TITLE            %12s "Zgoubi model"'
          write(LUN,50)'@ ORIGIN           %12s "Zgoubi model"'
          write(LUN,50)'@ DATE             %08s "23/03/11"'
          write(LUN,50)'@ TIME             %08s "17.28.23"'
 50       format(a)
 51       format(a,G18.10)
 53       format(2(a,G18.10,1x))
 52       format(a,i6)

C          CLOSE(LUN)
C        ENDIF 

        RETURN
      ENDIF
C-------------------------------------------------------------------
C-------------------------------------------------------------------

      KOPTCS = 0

      IF(NRES.LT.0) NRES=-NRES
C------- reactivate READ in zgoubi.dat
        READAT = .TRUE.

        IF(NRES.GT.0) THEN
          WRITE(6,101) IPASS
          WRITE(NRES,101) IPASS
 101      FORMAT(/,T25,
     >   ' *********************************************************',/
     >   ,T25
     >   ,' **************  End  of  TWISS  procedure  **************',
     >   // ,5X,' There  has  been ',I10,
     >        '  pass  through  the  optical  structure ',/)

          WRITE(NRES,FMT='(/,34X,1P,'' Momentum compaction : '',//, 
     >    30X,''dL/L / dp/p = '',G15.8)') ALPHA
          WRITE(NRES,FMT='(5X,1P,''(dp = '',G13.6,5X  
     >    ,'' L(0)   = '',G14.4,'' cm, ''
     >    ,'' L(0)-L(-dp) = '',G14.4,'' cm, ''
     >    ,'' L(0)-L(+dp) = '',G14.4,'' cm) '' )') 
     >    A(1,25), pathl(1),(pathl(1)-pathl(2)),(pathl(1)-pathl(3))
          WRITE(NRES,FMT='(/,34X,1P,'' Transition gamma  = '',
     >    G15.8)') 1.d0/SQRT(ALPHA)

          WRITE(NRES,FMT='(/,34X,1P,'' Chromaticities : '',//, 
     >    30X,''dNu_y / dp/p = '',G15.8,/, 
     >    30X,''dNu_z / dp/p = '',G15.8)') DNUYDP, DNUZDP


          IF(KTW .GE.3) THEN
C             DNUZDZ=(ZNUP-ZNUREF)/(UZP-UZREF)
C             DNUYDZ=(YNUP-YNUREF)/(UZP-UZREF)
        
             WRITE(NRES,FMT='(/,38X,1P,'' Amplitude  detunings : '',//, 
     >      42X,''/ dEps_y/pi       / dEps_z/pi'',/, 
     >      30X,''dNu_y'',7X,2(G15.8,3X),/, 
     >      30X,''dNu_z'',7X,2(G15.8,3X), //, 
     >      20X,''Nu_y_Ref = '',G15.8,'', Nu_z_Ref = '',G15.8, / 
     >      20X,''Nu_y_+dp = '',G15.8,'',   Nu_z_+dp = '',G15.8, / 
     >      20X,''Eps_y_Ref/pi = '',G15.8,'',   Eps_z_Ref/pi = '',G15.8, / 
     >      20X,''Eps_y_+dA/pi = '',G15.8,'',   Eps_z_+dA = '',G15.8)')
     >      DNUYDY, DNUYDZ, DNUZDY, DNUZDZ, 
     >      YNUREF,ZNUREF,
     >      YNUP,  ZNUP,
     >      UYREF, UZREF,
     >      UYP,   UZP

          ENDIF
      ENDIF

      goto 97

 96   CONTINUE
      WRITE(ABS(NRES),FMT='(/,''SBR TWISS : '',
     >           ''Error open file zgoubi.TWISS.Out'')')
      WRITE(*        ,FMT='(/,''SBR TWISS : '',
     >           ''Error open file zgoubi.TWISS.Out'')')

 97   continue

      IPASS = 1
C      NRBLT = 0

      RETURN

 20   CONTINUE
Compute linear functions, from multiturn tracking. Number of the particle used for that is given by user. 
c  F2( KPM : 7-> 1 ) : from end of last pass to end 6 passes earlier
      KPM = MIN(IPASS,8)
      IF(KPM.LE.7) THEN
        DO IC = 1, 6
          iic1 = ic
          if(ic.eq.1) then
            iic = 6
          else
            iic = ic-1
          endif
          F2(KPM,IIC) = F(IIC1,NINT(A(NOEL,2)))
        ENDDO
      ELSE
        DO IC = 1, 6
          iic1 = ic
          if(ic.eq.1) then
            iic = 6
          else
            iic = ic-1
          endif
          DO KP = 1, KPM-1
            F2(KP,IIC) = F2(KP+1,IIC)
          ENDDO
          F2(7,IIC) = F(IIC1,NINT(A(NOEL,2)))
        ENDDO
      ENDIF

C      write(88,*) ' IPASS = ',IPASS
C      write(88,fmt='(1p,6e12.4,a)') ((f2(i,j),j=1,6),' twiss',i=7,1,-1)
C      write(88,fmt='(1p,10X,6e12.4)') (f(i,1),i=1,6)
C      write(88,*) ' -------------------- '
      
      IF(IPASS.LT.7) RETURN       

      do ic=1,6
        do i=1,6
          xys(i,ic) = f2(7-i,ic)
C          sm(i,ic) = f2(8-i,ic)
          XYS(I,IC+6) = f2(8-i,ic)
        enddo
      enddo

      write(88,*) ' IPASS = ',IPASS
      write(88,fmt='(1p,12e12.4)') ((xys(i,ic),ic=1,12),i=1,6)
C      write(88,fmt='(1p,10X,6e12.4)') ((sm(i,ic),ic=1,6),i=1,6)
      write(88,*) ' -------------------- '
      
      IER = 0
      call dlgau(6,6,6,XYS,KAUX,ier)

      write(88,*) ' IPASS = ',IPASS,IER
      write(88,fmt='(1p,6e12.4)') ((xys(i,ic),i=1,6),iC=1,6)
      write(88,*) ' +++++++++++++-------------------- '

      RETURN


      entry twiss1(
     >             rturn)
      DO J=1,4
        DO I=1,4
          rturn(I,J) = RREF(I,J)
        ENDDO
      ENDDO

      RETURN

      END
