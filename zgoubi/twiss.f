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
      SUBROUTINE TWISS(
     >                 READAT) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL READAT
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      
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
  
      DATA KWRI6 / 1 /

c         WRITE(*,*) ' TWISS ',ipass, nrblt

      KTW = NINT(A(NOEL,1))
      FACP = A(NOEL,2)
      FACA = A(NOEL,3)
      FAC1 = 1.D0/FACA

C KTW=1 : 1 pass, equivalent to matrix[1,11]
C KTW=2 : 2 more passes to get matrices around +/-dp/p chromatic closed orbits;
C KTW=3 : 1 more pass to get matrices with +/-dY amplitude and an additional 
C           one to get matrices with +/-dZ amplitude. 

C      IF    (KTW .EQ. 1) THEN
C        NRBLT = 0
C      ELSEIF(KTW .EQ. 2) THEN
C        NRBLT = 2
C      ELSEIF(KTW .EQ. 3) THEN
C        NRBLT = 4
C      ENDIF

      IF(KTW.GE.99) GOTO 20
      IF(KTW.GE.1) GOTO 10
      RETURN

 10   CONTINUE
Compute optical functions, tunes, chromaticity, anharmonicities, from a few passes
C of 11 particles (based on MATRIX)
C      IF(IPASS .LT. NRBLT ) THEN
C        LUN=ABS(NRES) 
C        IF(LUN.GT.0) 
C     >     WRITE(LUN,100) IPASS  !, NRBLT+1
C          WRITE(*,*) IPASS  ,' NRBLT+1'
C 100       FORMAT(/,30X,'  -----  TWISS procedure  -----',//,5X,'End '
C     >     ,'of pass # ',I1,/)
C     >     ,'of pass # ',I1,'/',I1,' through the optical structure',/)
 
        IF(IPASS .EQ. 1) THEN
C Compute periodic beta from first pass.
C 2nd pass through structure will follow iff KTW>1.

C--------- Switch off print into zgoubi.res : 
          ANOEL2 = 0.1D0
          KWRIT = NINT(ANOEL2)
C--------- Switch on print to standard output :
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

          IF(KTW.GE.2) THEN
  
            IF(NRES .GT. 0) NRES =-NRES

            CALL REFER1(
     >                  PATHL(1)) 
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
            A(NUML,30) =  REF(1) -  1.D2*F0REF(1,6) * FAP25
            A(NUML,31) =  REF(2) -  1.D3*F0REF(2,6) * FAP25
            A(NUML,32) =  REF(3) -  1.D2*F0REF(3,6) * FAP25
            A(NUML,33) =  REF(4) -  1.D3*F0REF(4,6) * FAP25
          
            IPASS=IPASS+1
            NOEL=0 
            CALL SCUMS(0.D0)

          ENDIF

          RETURN
 
        ELSEIF(IPASS .EQ. 2) THEN
C------- 3rd pass through structure will follow

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
     >                PATHL(2)) 
 
C--------- Reset reference coordinates for OBJECT sampling : p -> p+dp
          NUML = 1
          FAP25 = FACP * A(NUML,25)
          A(NUML,35) =  REF(6) +  FAP25
          A(NUML,30) =  REF(1) +  1.D2*F0REF(1,6) * FAP25
          A(NUML,31) =  REF(2) +  1.D3*F0REF(2,6) * FAP25
          A(NUML,32) =  REF(3) +  1.D2*F0REF(3,6) * FAP25
          A(NUML,33) =  REF(4) +  1.D3*F0REF(4,6) * FAP25
          
          IPASS=IPASS+1
          NOEL=0 
          CALL SCUMS(0.D0)

          RETURN
 
        ELSEIF(IPASS .EQ. 3) THEN
C------- Chromatic tracking completed

C--------- reactivate WRITE for printing results 
C          IF(NRES.LT.0) NRES=-NRES

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
          CALL REFER1(
     >                PATHL(3)) 

C Momentum compaction
          ALPHA = ( (PATHL(3)-PATHL(2))/PATHL(1) ) / (2.D0 * A(NUML,25))

C Momentum detuning
          NUML = 1
          DNUYDP = (YNUP-YNUM)/2.D0/A(NUML,25)
          DNUZDP = (ZNUP-ZNUM)/2.D0/A(NUML,25)
        
          IF(KTW.GE.3) THEN
C--------- Amplitude detuning tracking & calculations follow
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

            RETURN
 
          ENDIF

C        ENDIF

C      ELSEIF(IPASS .EQ. NRBLT) THEN
        ELSEIF(IPASS .EQ. 4) THEN

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
          UYREF = F0REF(2,2)/1.D2*Y2+2.D0*(-F0REF(2,1))*YYP+
     >             F0REF(1,1)*1.D2*YP2 
          YY2 = A(NUML,20)-A(NUML,30)
          YYYP = YY2
          YY2 = YY2*YY2 
          YYP2 = A(NUML,21)-A(NUML,31)
          YYYP = YYYP*YYP2
          YYP2 = YYP2*YYP2 
          UYP = F0P(2,2)/1.D2*YY2 + 2.D0*(-F0P(2,1))*YYYP + 
     >             F0P(1,1)*1.D2*YYP2 
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
          RETURN 

C      ELSEIF(IPASS .GT. NRBLT) THEN
        ELSEIF(IPASS .EQ. 5) THEN
C------- Amplitude tracking completed

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
        UZREF = F0REF(4,4)/1.D2*Z2+2.D0*(-F0REF(4,3))*ZZP+
     >       F0REF(3,3)*1.D2*ZP2 
        ZZ2 = A(NUML,22)-A(NUML,32)
        ZZZP = ZZ2
        ZZ2 = ZZ2*ZZ2 
        ZZP2 = A(NUML,23)-A(NUML,33)
        ZZZP = ZZZP*ZZP2
        ZZP2 = ZZP2*ZZP2 
        UZP = F0P(4,4)/1.D2*ZZ2 + 2.D0*(-F0P(4,3))*ZZZP + 
     >          F0P(3,3)*1.D2*ZZP2 

      ENDIF

      IF(NRES.LT.0) NRES=-NRES
C----- reactivate READ in zgoubi.dat
      READAT = .TRUE.

      IF(NRES.LT.0) NRES=-NRES
      IF(NRES.GT.0) THEN
        WRITE(*,101) IPASS
        WRITE(NRES,101) IPASS
 101    FORMAT(/,25X,' ****  End  of  ''TWISS''  procedure  ****',//
     >   ,5X,' There  has  been ',I10,
     >        '  pass  through  the  optical  structure ',/)

        WRITE(NRES,FMT='(/,34X,1P,'' Momentum compaction : '',//, 
     >  30X,''dL/L / dp/p = '',G15.8)') ALPHA
        WRITE(NRES,FMT='(5X,1P,''(dp = '',G13.6,5X  
     >  ,'' L(0)   = '',G14.4,'' cm, ''
     >  ,'' L(0)-L(-dp) = '',G14.4,'' cm, ''
     >  ,'' L(0)-L(+dp) = '',G14.4,'' cm) '' )') 
     >  A(1,25), pathl(1),(pathl(1)-pathl(2)),(pathl(1)-pathl(3))
        WRITE(NRES,FMT='(/,34X,1P,'' Transition gamma  = '',
     >  G15.8)') 1.d0/SQRT(ALPHA)

        WRITE(NRES,FMT='(/,34X,1P,'' Chromaticities : '',//, 
     >  30X,''dNu_y / dp/p = '',G15.8,/, 
     >  30X,''dNu_z / dp/p = '',G15.8)') DNUYDP, DNUZDP


        IF(KTW .GE.3) THEN
          DNUZDZ=(ZNUP-ZNUREF)/(UZP-UZREF)
          DNUYDZ=(YNUP-YNUREF)/(UZP-UZREF)
        
          WRITE(NRES,FMT='(/,38X,1P,'' Amplitude  detunings : '',//, 
     >    42X,''/ dEps_y/pi       / dEps_z/pi'',/, 
     >    30X,''dNu_y'',7X,2(G15.8,3X),/, 
     >    30X,''dNu_z'',7X,2(G15.8,3X), //, 
     >    20X,''Nu_y_Ref = '',G15.8,'', Nu_z_Ref = '',G15.8, / 
     >    20X,''Nu_y_+dp = '',G15.8,'',   Nu_z_+dp = '',G15.8, / 
     >    20X,''Eps_y_Ref/pi = '',G15.8,'',   Eps_z_Ref/pi = '',G15.8, / 
     >    20X,''Eps_y_+dA/pi = '',G15.8,'',   Eps_z_+dA = '',G15.8)')
     >    DNUYDY, DNUYDZ, DNUZDY, DNUZDZ, 
     >    YNUREF,ZNUREF,
     >    YNUP,  ZNUP,
     >    UYREF, UZREF,
     >    UYP,   UZP

        ENDIF
      ENDIF

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
      END
