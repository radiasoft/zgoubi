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
      SUBROUTINE TWISS(LUN,OKCPLD,
     >                            KOPTCS, READAT,KTW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL READAT, OKCPLD
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
      CHARACTER(KSIZ)  KLOBJ

      CHARACTER(KSIZ) KLEO

      LOGICAL IDLUNI
      INTEGER DEBSTR, FINSTR

      LOGICAL dolast
      save dolast
      
      dimension rturn(4,4)
      logical prdic, okorbt

      DIMENSION F0(6,6) 

      SAVE Q1, Q2, CC, OKORBT

      DATA KWRI6 / 1 /
      DATA DOLAST / .TRUE. /

      NMAIL = 1
      PRDIC = .TRUE.

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
C KTW2 = 1 is a request for orbit search prior to attacking twiss procedure

      IF(KTW.GE.99) GOTO 20  ! To be completed and tested
      IF(KTW.GE.1) GOTO 10
      RETURN

 10   CONTINUE
Compute optical functions, tunes, chromaticity, anharmonicities, from a few passes
C of 11 particles (based on MATRIX)

      okorbt = .not. ktw2 .eq. 1 
      if(okorbt) Then
        if(nres.gt.0) then
          write(nres,*)
          write(nres,*) ' Closed orbit search was not requested, '
     >    //' particle 1 is assumed on closed orbit.'
          write(nres,*)
        endif
      else
        if(nres.gt.0) then
          write(nres,*)
          write(nres,*) ' Closed orbit search was requested, closed or'
     >    //'bit coordinates first searched and assigned to particle 1.'
          write(nres,*)
        endif
      endif

      if(.not. okorbt) call tworbt

 
      if(ipass .eq. 4 .and. ktw .eq. 2) goto 222

      KLOBJ = 'OBJET'
      CALL GETNOL(KLOBJ,
     >                  NLOBJ)      

      IF(IPASS .EQ. 1) THEN
C Compute periodic beta from first pass.
C 2nd pass through structure will follow iff KTW>1.

        DOLAST = .TRUE.

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
          CALL MAT1(1,
     >                RREF,T)
          CALL REFER(2,1,0,1,4,5)
        ELSEIF(IORD .EQ. 2) THEN
          CALL REFER(1,2,0,1,6,7)
          CALL MAT2(RREF,T,TX3,TX4)
          CALL REFER(2,2,0,1,6,7)
        ENDIF

        CALL MKSA(IORD,RREF,T,TX3,TX4)
C        CALL MATIMP(RREF)
        CALL TUNES(RREF,F0REF,NMAIL,IERY,IERZ,.TRUE.,
     >                                          YNUREF,ZNUREF,CMUY,CMUZ)
        CALL MATIMP(RREF,F0REF,YNUREF,ZNUREF,CMUY,CMUZ,NMAIL,PRDIC,1)

        CALL BEAMAT(RREF,PRDIC,OKCPLD,
     >                                 F0,Q1,Q2,CC)        
            
c             write(*,*) q1, q2, cc, f0
c                read(*,*)
        NRES = ISIGN*NRES

c          call ZGKLE(iq(noel), 
c     >                             kleo)
c           write(*,*) ' 1 twiss noel iq(noel) ',noel, kleo 
c          call ZGKLE(iq(noel-1), 
c     >                             kleo)
c           write(*,*) ' noel -1 ',kleo 
c          call ZGKLE(iq(noel+1), 
c     >                             kleo)
c          call ZGKLE(iq(noel+1), 
c     >                             kleo)
c           write(*,*) ' noel +1 ',kleo 
c          call ZGKLE(iq(noel+2), 
c     >                             kleo)
c           write(*,*) ' noel +2 ',kleo 

        IF(KTW.GE.2) THEN
  
c          IF(KTW2 .EQ. 0) THEN
c            IF(NRES .GT. 0) NRES =-NRES
c          ENDIF

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
          CALL MAT1(1,
     >                RMINUS,T)
          CALL REFER(2,1,0,1,4,5)
        ELSEIF(IORD .EQ. 2) THEN
          CALL REFER(1,2,0,1,6,7)
          CALL MAT2(RMINUS,T,TX3,TX4)
          CALL REFER(2,2,0,1,6,7)
        ENDIF
        CALL MKSA(IORD,RMINUS,T,TX3,TX4)
C        CALL MATIMP(RMINUS)
        CALL TUNES(RMINUS,F0M,NMAIL,IERY,IERZ,.TRUE.,
     >                                           YNUM,ZNUM,CMUY,CMUZ)
        CALL MATIMP(RMINUS,F0M,YNUM,ZNUM,CMUY,CMUZ,NMAIL,PRDIC,1)
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

c          call ZGKLE(iq(noel), 
c     >                             kleo)
c           write(*,*) ' 2 twiss noel iq(noel) ',noel, kleo 
c          call ZGKLE(iq(noel-1), 
c     >                             kleo)
c           write(*,*) ' noel -1 ',kleo 
c          call ZGKLE(iq(noel+1), 
c     >                             kleo)
c           write(*,*) ' noel +1 ',kleo 
c          call ZGKLE(iq(noel+2), 
c     >                             kleo)
c           write(*,*) ' noel +2 ',kleo 

   
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
          CALL MAT1(1,
     >                RPLUS,T)
          CALL REFER(2,1,0,1,4,5)
        ELSEIF(IORD .EQ. 2) THEN
          CALL REFER(1,2,0,1,6,7)
          CALL MAT2(RPLUS,T,TX3,TX4)
          CALL REFER(2,2,0,1,6,7)
        ENDIF
        CALL MKSA(IORD,RPLUS,T,TX3,TX4)
C        CALL MATIMP(RPLUS)
        CALL TUNES(RPLUS,F0P,NMAIL,IERY,IERZ,.TRUE.,
     >                                          YNUP,ZNUP,CMUY,CMUZ)
        CALL MATIMP(RPLUS,F0P,YNUP,ZNUP,CMUY,CMUZ,NMAIL,PRDIC,1)
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

c          call ZGKLE(iq(noel), 
c     >                             kleo)
c           write(*,*) ' 3 twiss noel iq(noel) ',noel, kleo 
c          call ZGKLE(iq(noel-1), 
c     >                             kleo)
c           write(*,*) ' noel -1 ',kleo 
c          call ZGKLE(iq(noel+1), 
c     >                             kleo)
c           write(*,*) ' noel +1 ',kleo 
c          call ZGKLE(iq(noel+2), 
c     >                             kleo)
c           write(*,*) ' noel +2 ',kleo 

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
            CALL MAT1(1,
     >                  RPLUS,T)
            CALL REFER(2,1,0,1,4,5)
          ELSEIF(IORD .EQ. 2) THEN
            CALL REFER(1,2,0,1,6,7)
            CALL MAT2(RPLUS,T,TX3,TX4)
            CALL REFER(2,2,0,1,6,7)
          ENDIF
          CALL MKSA(IORD,RPLUS,T,TX3,TX4)
C          CALL MATIMP(RPLUS)
          CALL TUNES(RPLUS,F0P,NMAIL,IERY,IERZ,.TRUE.,
     >                                            YNUP,ZNUP,CMUY,CMUZ)
          CALL MATIMP(RPLUS,F0P,YNUP,ZNUP,CMUY,CMUZ,NMAIL,PRDIC,1)

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
            CALL MAT1(1,
     >                  RPLUS,T)
            CALL REFER(2,1,0,1,4,5)
          ELSEIF(IORD .EQ. 2) THEN
            CALL REFER(1,2,0,1,6,7)
            CALL MAT2(RPLUS,T,TX3,TX4)
            CALL REFER(2,2,0,1,6,7)
          ENDIF
          CALL MKSA(IORD,RPLUS,T,TX3,TX4)
C          CALL MATIMP(RPLUS)
          CALL TUNES(RPLUS,F0P,NMAIL,IERY,IERZ,.TRUE.,
     >                                          YNUP,ZNUP,CMUY,CMUZ)
          CALL MATIMP(RPLUS,F0P,YNUP,ZNUP,CMUY,CMUZ,NMAIL,PRDIC,1)

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

      IF(DOLAST) THEN
        DOLAST = .FALSE.

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
          write(LUN,51)'@ |C|              %le', CC
          write(LUN,51)'@ Q1               %le', Q1
          write(LUN,51)'@ Q2               %le', Q2
          write(LUN,50)'@ SYNCH_4          %le                   0'
          write(LUN,50)'@ SYNCH_5          %le                   0'
          write(LUN,50)'@ TITLE            %12s "Zgoubi model"'
          write(LUN,50)'@ ORIGIN           %12s "twiss.f"'
          write(LUN,50)'@ DATE             %08s "  "'
          write(LUN,50)'@ TIME             %08s "  "'
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
          WRITE(NRES,FMT='(5X,1P,''(dp = '',E13.6,5X  
     >    ,'' L(0)   = '',E14.5,'' cm, ''
     >    ,'' L(0)-L(-dp) = '',E14.5,'' cm, ''
     >    ,'' L(0)-L(+dp) = '',E14.5,'' cm) '' )') 
     >    A(1,25), pathl(1),(pathl(1)-pathl(2)),(pathl(1)-pathl(3))
          WRITE(NRES,FMT='(/,34X,1P,'' Transition gamma  = '',
     >    E15.8)') 1.d0/SQRT(ALPHA)

          WRITE(NRES,FMT='(/,34X,1P,'' Chromaticities : '',//, 
     >    30X,''dNu_y / dp/p = '',E15.8,/, 
     >    30X,''dNu_z / dp/p = '',E15.8)') DNUYDP, DNUZDP


          IF(KTW .GE.3) THEN
C             DNUZDZ=(ZNUP-ZNUREF)/(UZP-UZREF)
C             DNUYDZ=(YNUP-YNUREF)/(UZP-UZREF)
        
             WRITE(NRES,FMT='(/,38X,1P,'' Amplitude  detunings : '',//, 
     >      42X,''/ dEps_y/pi       / dEps_z/pi'',/, 
     >      30X,''dNu_y'',7X,2(E15.8,3X),/, 
     >      30X,''dNu_z'',7X,2(E15.8,3X), //, 
     >      20X,''Nu_y_Ref = '',E15.8,'', Nu_z_Ref = '',E15.8, / 
     >      20X,''Nu_y_+dp = '',E15.8,'',   Nu_z_+dp = '',E15.8, / 
     >      20X,''Eps_y_Ref/pi = '',E15.8,'',   Eps_z_Ref/pi = '',E15.8, / 
     >      20X,''Eps_y_+dA/pi = '',E15.8,'',   Eps_z_+dA = '',E15.8)')
     >      DNUYDY, DNUYDZ, DNUZDY, DNUZDZ, 
     >      YNUREF,ZNUREF,
     >      YNUP,  ZNUP,
     >      UYREF, UZREF,
     >      UYP,   UZP

          ENDIF
      ENDIF

      IPASS = 1

      RETURN

C-------------------------------------------------------------------
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

      IF(IPASS.LT.7) RETURN       

      DO IC=1,6
        DO I=1,6
          XYS(I,IC) = F2(7-I,IC)
C          SM(I,IC) = F2(8-I,IC)
          XYS(I,IC+6) = F2(8-I,IC)
        ENDDO
      ENDDO

C      WRITE(88,*) ' IPASS = ',IPASS
C      WRITE(88,FMT='(1P,12E12.4)') ((XYS(I,IC),IC=1,12),I=1,6)
CC      WRITE(88,FMT='(1P,10X,6E12.4)') ((SM(I,IC),IC=1,6),I=1,6)
C      WRITE(88,*) ' -------------------- '
      
      IER = 0
      CALL DLGAU(6,6,6,XYS,KAUX,IER)

C      WRITE(88,*) ' IPASS = ',IPASS,IER
C      WRITE(88,FMT='(1P,6E12.4)') ((XYS(I,IC),I=1,6),IC=1,6)
C      WRITE(88,*) ' +++++++++++++-------------------- '

      RETURN


      ENTRY TWISS1(
     >             RTURN)
      DO J=1,4
        DO I=1,4
          RTURN(I,J) = RREF(I,J)
        ENDDO
      ENDDO

      RETURN

      END
