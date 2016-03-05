      SUBROUTINE MAT3RD(LUN)
C     ------------------------------------------------------
C     Computation of transport coefficients up to 3rd order. 
C     From Nick Tsoupas, RAYTRACE, Oct 2015
C     ------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      INCLUDE "C.OBJET.H"     ! COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),IREP(MXT),AMQLU,PABSLU

      DIMENSION XO(MXT), YO(MXT), VXO(MXT), VYO(MXT)
C      DIMENSION XO(MXT), YO(MXT), ZO(MXT), VXO(MXT), VYO(MXT)
C     1, VZO(MXT), RTL(MXT),RLL(MXT)
      DIMENSION XI(MXT), YI(MXT), VXI(MXT), VYI(MXT)
C      DIMENSION XI(MXT), YI(MXT), ZI(MXT), VXI(MXT), VYI(MXT)
C     1, VZI(MXT),DELP(MXT)
      COMMON/ABERR/  R(6,6),T2(5,6,6),W3(5,6,6,6),Z4(5,6,6,6,6)
      DIMENSION X(MXT)
      DIMENSION PL_D(MXT)
C
      save IFLAGCO, ITERCO
      data IFLAGCO / 3 / ! Indicates that closed orbit search in RAYTRACE converged
      data ITERCO / 0 /

      do j = 1, mxj 
        do it = 1, imax

          XO(IT) = f(2,it)*1d-2
          YO(IT) = f(4,it)*1d-2
c          ZO(IT) = f(6,it)
          VXO(IT) = f(3,it)*1d-3
          VYO(IT) = f(5,it)*1d-3
c          VZO(IT) = f(2,it)
c          RTL(IT = f(7,it)*1d-6
c          RLL(IT) = f(6,it)*1d-2
          PL_D(IT) = f(6,it)*1d-2

          XI(IT) = fO(2,it)*1d-2
          YI(IT) = fO(4,it)*1d-2
c          ZI(IT) = fO(6,it)
          VXI(IT) = fO(3,it)*1d-3
          VYI(IT) = fO(5,it)*1d-3
c          VZI(IT) = fO(2,it)
        enddo        
      enddo

	N6=6        
	IF(IFLAGCO.EQ.1) N6=2
      DO 21  I1= 1,N6
      DO 21  I2= 1,N6
      R(I1,I2) = 0.          
	IF(IFLAGCO.EQ.1) GO TO 2111
      DO 21 I3= 1,5
      T2(I3,I1,I2) = 0.0
      DO 21 I4=1,6
      W3(I3,I4,I2,I1)=0.0
      DO 21 I5=1,6
      Z4(I3,I5,I4,I2,I1)=0.0
21    CONTINUE
2111  CONTINUE 
C****
      		XM= FO(2,5)*1d-2  !XI(5)
      		TM= FO(3,9)*1d-3  !VXI(9)
      		YM= FO(4,13)*1d-2  !YI(13)
      		PM= FO(5,17)*1d-3  !VYI(17)
      		DM= FO(1,21) -1.d0  !DELP(21)

      IF(IMAX.EQ.22) GO TO 1000

      X2M=XM*XM
      XTM=XM*TM
      XYM=XM*YM
      XPM=XM*PM
      XDM=XM*DM
C*
      T2M=TM*TM
      TYM=TM*YM
      TPM=TM*PM
      TDM=TM*DM
C*
      Y2M=YM*YM
      YPM=YM*PM
      YDM=YM*DM
C*
      P2M=PM*PM
      PDM=PM*DM
C*
      D2M=DM*DM
C*
      IF(IMAX.EQ.82) GO TO 1000
C*
      X3M=X2M*XM
C*
      X2TM=X2M*TM
      X2YM=X2M*YM
      X2PM=X2M*PM
      X2DM=X2M*DM
C*
      XT2M=XM*T2M
      XTYM=XTM*YM
      XTPM=XTM*PM
      XTDM=XTM*DM
      XY2M=XM*Y2M
      XYPM=XM*YM*PM
      XYDM=XYM*DM
C*
      XP2M=XM*P2M
      XPDM=XPM*DM
C*
      XD2M=XM*D2M
C*
      T3M =T2M*TM
C*
      T2YM=T2M*YM
      T2PM=T2M*PM
      T2DM=T2M*DM
C*
      TY2M=TM*Y2M
      TYPM=TM*YPM
      TYDM=TYM*DM
C*
      TP2M=TM*P2M
      TPDM=TPM*DM
C*
      TP2M=TM*P2M
      TPDM=TM*PDM
C*
      TD2M=TM*D2M
C*
      Y3M=Y2M*YM
C*
      Y2PM=Y2M*PM
      Y2DM=Y2M*DM
C*
      YP2M=YM*P2M
      YPDM=YPM*DM
C*
      YD2M=YM*D2M
C*
      P3M =P2M*PM
C*
      P2DM=P2M*DM
      PD2M=PM*D2M
C*
      D3M=D2M*DM
C*
      X4M=X3M*XM
      T4M=T3M*TM
      Y4M=Y3M*YM
      P4M=P3M*PM
      D4M=D3M*DM
C*
      X2T2M=X2M*T2M
      X2Y2M=X2M*Y2M
      X2P2M=X2M*P2M
      X2D2M=X2M*D2M
C*
      T2Y2M=T2M*Y2M
      T2P2M=T2M*P2M
      T2D2M=T2M*D2M
C*
      Y2P2M=Y2M*P2M
      Y2D2M=Y2M*D2M
C*
      P2D2M=P2M*D2M
C*
1000  CONTINUE
C****
C**** CALCULATE COEFFICIENTS
C****
C***
C      call MAT3ca(XO,YO,ZO,VXO,VYO,VZO,RTL,RLL,PL_D,XM,TM,YM,PM,DM,j)
      call MAT3ca(XO,YO,VXO,VYO,PL_D,XM,TM,YM,PM,DM,j)

        IF(IMAX.EQ.22) RETURN
C*
C*     SECOND ORDER ABERRATION COEFFICIENTS.
C*****
C*
      DO 100 L=1,5
C**
      NRC=82
       CALL MAT3TR(L,IMAX,X,XO,YO,VXO,VYO,PL_D)
C       CALL MAT3TR(L,IMAX,X,XO,YO,ZO,VXO,VYO,VZO,PL_D,RLL)
C*
      T2(L,1,1)= ( 81.0*(X(3)+X(4))+(X(5)+X(6)) )/(16.0*X2M)
C*
      T2(L,1,2)=( 81.0*(X(23)+X(24))-(X(25)+X(28))
     1           -81.0*(X(3)+X(4))-81.0*(X(7)+X(8))
     1           +(X(5)+X(6))+(X(9)+X(10)) )/(16.0*XTM)
      T2(L,2,2)=( 81.0*(X(7)+X(8))-(X(9)+X(10)) )/(16.0*T2M)
C*
      T2(L,1,3)=( 81.0*(X(29)+X(30))-(X(31)+X(34))
     1           -81.0*(X(3)+X(4))-81.0*(X(11)+X(12))
     1           +(X(5)+X(6))+(X(13)+X(14)) )/(16.0*XYM)
      T2(L,2,3)=( 81.0*(X(47)+X(48))-(X(49)+X(52))
     1           -81.0*(X(7)+X(8))-81.0*(X(11)+X(12))
     1           +(X(9)+X(10))+(X(13)+X(14)) )/(16.0*TYM)
      T2(L,3,3)=( 81.0*(X(11)+X(12))-(X(13)+X(14)) )/(16.0*Y2M)
C*
      T2(L,1,4)=( 81.0*(X(35)+X(36))-(X(37)+X(40))
     1           -81.0*(X(3)+X(4))-81.0*(X(15)+X(16))
     1           +(X(5)+X(6))+(X(17)+X(18)) )/(16.0*XPM)
      T2(L,2,4)=( 81.0*(X(53)+X(54))-(X(55)+X(58))
     1           -81.0*(X(7)+X(8))-81.0*(X(15)+X(16))
     1           +(X(9)+X(10))+(X(17)+X(18)) )/(16.0*TPM)
      T2(L,3,4)=( 81.0*(X(65)+X(66))-(X(67)+X(70))
     1           -81.0*(X(11)+X(12))-81.0*(X(15)+X(16))
     1           +(X(13)+X(14))+(X(17)+X(18)) )/(16.0*YPM)
      T2(L,4,4)=( 81.0*(X(15)+X(16))-(X(17)+X(18)) )/(16.0*P2M)
C*
      T2(L,1,6)=( 81.0*(X(41)+X(42))-(X(43)+X(46))
     1           -81.0*(X(3)+X(4))-81.0*(X(19)+X(20))
     1           +(X(5)+X(6))+(X(21)+X(22)) )/(16.0*XDM)
      T2(L,2,6)=( 81.0*(X(59)+X(60))-(X(61)+X(64))
     1           -81.0*(X(7)+X(8))-81.0*(X(19)+X(20))
     1           +(X(9)+X(10))+(X(21)+X(22)) )/(16.0*TDM)
      T2(L,3,6)=( 81.0*(X(71)+X(72))-(X(73)+X(76))
     1           -81.0*(X(11)+X(12))-81.0*(X(19)+X(20))
     1           +(X(13)+X(14))+(X(21)+X(22)) )/(16.0*YDM)
      T2(L,4,6)=( 81.0*(X(77)+X(78))-(X(79)+X(82))
     1           -81.0*(X(15)+X(16))-81.0*(X(19)+X(20))
     1           +(X(17)+X(18))+(X(21)+X(22)) )/(16.0*PDM)
      T2(L,6,6)=( 81.0*(X(19)+X(20))-(X(21)+X(22)) )/(16.0*D2M)
C*
100   CONTINUE
C*
      WRITE(LUN,120)
  120 FORMAT(/,40X, '  *2nd order TRANSPORT coeffs.* (units m, rad)')
      DO 24 I1= 1,5
      DO 25 I2= 1,6
      WRITE(LUN,121) ( I1,I3,I2, T2(I1,I3,I2), I3=1,I2 )
C      write(*,121) ( I1,I3,I2, T2(I1,I3,I2), I3=1,I2 )
  121 FORMAT(9X,6(I4,I2,I1, 1PE11.3)  )
   25 CONTINUE
      WRITE(LUN,122)
  122 FORMAT( 1X  )
   24 CONTINUE
C*
       IF(IMAX.EQ.82) RETURN
C*
C*    CALCULATE THIRD ORDER TERMS
C*
      DO 200 L=1,5
C*
      NRC=102
      CALL MAT3TR(L,IMAX,X,XO,YO,VXO,VYO,PL_D)
C      CALL MAT3TR(L,IMAX,X,XO,YO,ZO,VXO,VYO,VZO,PL_D,RLL)
C*
      W3(L,1,1,1) = ( -27.0*(X(3)-X(4))+9.0*(X(5)-X(6)) )/(16.0*X3M)
C*
      W3(L,1,1,2) = ( (X(25)+X(26))-(X(27)+X(28))-2.0*(X(9)-X(10)) )
     1               /(4.0*X2TM)
C*
      W3(L,1,2,2) = ( (X(25)-X(26))+(X(27)-X(28))-2.0*(X(5)-X(6)) )
     1               /(4.0*XT2M)
      W3(L,2,2,2) = ( -27.0*(X(7)-X(8))+9.0*(X(9)-X(10)) )/(16.0*T3M)
C*
      W3(L,1,1,3) = ( (X(31)+X(32))-(X(33)+X(34))-2.0*(X(13)-X(14)) )
     1               /(4.0*X2YM)
C*
      W3(L,1,2,3)=((X(83)-X(84))+(X(5)-X(6))+(X(9)-X(10))+(X(13)-X(14))
     1          -(X(25)-X(28))-(X(31)-X(34))-(X(49)-X(52)) )/(2.0*XTYM)
      W3(L,2,2,3) = ( (X(49)+X(50))-(X(51)+X(52))-2.0*(X(13)-X(14)) )
     1               /(4.0*T2YM)
C*
      W3(L,1,3,3) = ( (X(31)-X(32))+(X(33)-X(34))-2.0*(X(5)-X(6)) )
     1               /(4.0*XY2M)
      W3(L,2,3,3) = ( (X(49)-X(50))+(X(51)-X(52))-2.0*(X(9)-X(10)) )
     1               /(4.0*TY2M)
      W3(L,3,3,3) =(-27.0*(X(11)-X(12))+9.0*(X(13)-X(14)) )/(16.0*Y3M)
C*
      W3(L,1,1,4) = ( (X(37)+X(38))-(X(39)+X(40))-2.0*(X(17)-X(18)) )
     1               /(4.0*X2PM)
C*
      W3(L,1,2,4)=((X(85)-X(86))+(X(5)-X(6))+(X(9)-X(10))+(X(17)-X(18))
     1           -(X(25)-X(28))-(X(37)-X(40))-(X(55)-X(58)) )/(2.0*XTPM)
      W3(L,2,2,4) = ( (X(55)+X(56))-(X(57)+X(58))-2.0*(X(17)-X(18)) )
     1               /(4.0*T2PM)
C*
      W3(L,1,3,4)=((X(89)-X(90))+(X(5)-X(6))+(X(13)-X(14))+(X(17)-X(18))
     1           -(X(31)-X(34))-(X(37)-X(40))-(X(67)-X(70)))/(2.0*XYPM)
      W3(L,2,3,4)=((X(95)-X(96))+(X(9)-X(10))+(X(13)-X(14))+X(17)-X(18)
     1           -(X(49)-X(52))-(X(55)-X(58))-(X(67)-X(70)))/(2.0*TYPM)
      W3(L,3,3,4)=( (X(67)+X(68))-(X(69)+X(70))-2.0*(X(17)-X(18)))
     1             /(4.0*Y2PM)
C*
      W3(L,1,4,4)=( (X(37)-X(38))+(X(39)-X(40))-2.0*(X(5)-X(6)))
     1             /(4.0*XP2M)
      W3(L,2,4,4)=( (X(55)-X(56))+(X(57)-X(58))-2.0*(X(9)-X(10)) )
     1             /(4.0*TP2M)
      W3(L,3,4,4)=( (X(67)-X(68))+(X(69)-X(70))-2.0*(X(13)-X(14)) )
     1             /(2.0*YP2M)
      W3(L,4,4,4)=( -27.0*(X(15)-X(16))+9.0*(X(17)-X(18)) )
     1             /(16.0*P3M)
C*
      W3(L,1,1,6)=( (X(43)+X(44))-(X(45)+X(46))-2.0*(X(21)-X(22)) )
     1             /(4.0*X2DM)
C*
      W3(L,1,2,6)=((X(87)-X(88))+(X(5)-X(6))+(X(9)-X(10))+(X(21)-X(22))
     1           -(X(25)-X(28))-(X(43)-X(46))-(X(61)-X(64)))/(2.0*XTDM)
      W3(L,2,2,6)=( (X(61)+X(62))-(X(63)+X(64))-2.0*(X(21)-X(22)))
     1             /(4.0*T2DM)
C*
      W3(L,1,3,6)=((X(91)-X(92))+(X(5)-X(6))+(X(13)-X(14))+(X(21)-X(22))
     1           -(X(31)-X(34))-(X(43)-X(46))-(X(73)-X(76)))/(2.0*XYDM)
      W3(L,2,3,6)=((X(97)-X(98))+X(9)-X(10)+(X(13)-X(14))+(X(21)-X(22))
     1           -(X(49)-X(52))-(X(61)-X(64))-(X(73)-X(76)))/(2.0*TYDM)
      W3(L,3,3,6)=( (X(73)+X(74))-(X(75)+X(76))-2.0*(X(21)-X(22)) )
     1             /(4.0*Y2DM)
C*
      W3(L,1,4,6)=((X(93)-X(94))+(X(5)-X(6))+(X(17)-X(18))+(X(21)-X(22))
     1           -(X(37)-X(40))-(X(43)-X(46))-(X(79)-X(82)))/(2.0*XPDM)
      W3(L,2,4,6)=((X(99)-X(100))+(X(9)-X(10))+(X(17)-X(18))+X(21)-X(22)
     1           -(X(55)-X(58))-(X(61)-X(64))-(X(79)-X(82)))/(2.0*TPDM)
      W3(L,3,4,6)=((X(101)-X(102))+(X(13)-X(14))+X(17)-X(18)+X(21)-X(22)
     1           -(X(67)-X(70))-(X(73)-X(76))-(X(79)-X(82)))/(2.0*YPDM)
      W3(L,4,4,6)=( (X(79)+X(80))-(X(81)+X(82))-2.0*(X(21)-X(22)) )
     1             /(4.0*P2DM)
C*
      W3(L,1,6,6)=( (X(43)-X(44))+(X(45)-X(46))-2.0*(X(5)-X(6)) )
     1             /(4.0*XD2M)
      W3(L,2,6,6)=( (X(61)-X(62))+(X(63)-X(64))-2.0*(X(9)-X(10)) )
     1             /(4.0*TD2M)
      W3(L,3,6,6)=( (X(73)-X(74))+(X(75)-X(76))-2.0*(X(13)-X(14)) )
     1             /(4.0*YD2M)
      W3(L,4,6,6)=( (X(79)-X(80))+(X(81)-X(82))-2.0*(X(17)-X(18)) )
     1             /(4.0*PD2M)
      W3(L,6,6,6)=( -27.0*(X(19)-X(20))+9.0*(X(21)-X(22)) )
     1             /(16.0*D3M)
C*
200   CONTINUE
C*
C*    PRINT THIRD ORDER TRANSFORM
C*
      WRITE(LUN,500)
C	write(*,500)
500   FORMAT(/,40X,' *3rd order TRANSPORT coeffs.* (units m, rad)')
C*
      DO 300 I1=1,5
502   FORMAT(1X)
      DO 299 I4=1,6
      DO 298 I3=1,I4
C*
       WRITE(LUN,501) (I1,I2,I3,I4, W3(I1,I2,I3,I4), I2=1,I3)
C       write(*,501) (I1,I2,I3,I4, W3(I1,I2,I3,I4), I2=1,I3)
501   FORMAT(8x,6(I4,I2,I1,I1, 1PE11.3) )
298   CONTINUE
299   CONTINUE
       WRITE(LUN,502)
300   CONTINUE
C*
C*    COMPUTE FEW FORTH ORDER COEFFICIENTS
C*
      DO 700 L=1,5
C*
      NRC=102
      CALL MAT3TR(L,IMAX,X,XO,YO,VXO,VYO,PL_D)
C      CALL MAT3TR(L,IMAX,X,XO,YO,ZO,VXO,VYO,VZO,PL_D,RLL)
C*
      Z4(L,1,1,1,1)=(-81.0*(X(3)+X(4))+9.0*(X(5)+X(6)) )
     1                /(16.0*X4M)
      Z4(L,2,2,2,2)=(-81.0*(X(7)+X(8))+9.0*(X(9)+X(10)) )
     1               /(16.0*T4M)
      Z4(L,3,3,3,3)=(-81.0*(X(11)+X(12))+9.0*(X(13)+X(14)))
     1               /(16.0*Y4M)
      Z4(L,4,4,4,4)=(-81.0*(X(15)+X(16))+9.0*(X(17)+X(18)))
     1               /(16.0*P4M)
      Z4(L,6,6,6,6)=(-81.0*(X(19)+X(20))+9.0*(X(21)+X(22)))
     1               /(16.0*D4M)
C*
      Z4(L,1,1,2,2)=( (X(25)+X(26))+(X(27)+X(28))
     1            -2.0*(X(5)+X(6))-2.0*(X(9)+X(10)) )/(4.0*X2T2M)
      Z4(L,1,1,3,3)=( (X(31)+X(32))+(X(33)+X(34))
     1            -2.0*(X(5)+X(6))-2.0*(X(13)+X(14)) )/(4.0*X2Y2M)
      Z4(L,1,1,4,4)=( (X(37)+X(38))+(X(39)+X(40))
     1            -2.0*(X(5)+X(6))-2.0*(X(17)+X(18)) )/(4.0*X2P2M)
      Z4(L,1,1,6,6)=( (X(43)+X(44))+(X(45)+X(46))
     1            -2.0*(X(5)+X(6))-2.0*(X(21)+X(22)) )/(4.0*X2D2M)
C*
      Z4(L,2,2,3,3)=( (X(49)+X(50))+(X(51)+X(52))
     1            -2.0*(X(9)+X(10))-2.0*(X(13)+X(14)))/(4.0*T2Y2M)
      Z4(L,2,2,4,4)=( (X(55)+X(56))+(X(57)+X(58))
     1            -2.0*(X(9)+X(10))-2.0*(X(17)+X(18)))/(4.0*T2P2M)
      Z4(L,2,2,6,6)=( (X(61)+X(62))+(X(63)+X(64))
     1            -2.0*(X(9)+X(10))-2.0*(X(21)+X(22)))/(4.0*T2D2M)
C*
      Z4(L,3,3,4,4)=( (X(67)+X(68))+(X(69)+X(70))
     1           -2.0*(X(13)+X(14))-2.0*(X(17)+X(18)))/(4.0*Y2P2M)
      Z4(L,3,3,6,6)=( (X(73)+X(74))+(X(75)+X(76))
     1           -2.0*(X(13)+X(14))-2.0*(X(21)+X(22)))/(4.0*Y2D2M)
C*
      Z4(L,4,4,6,6)=( (X(79)+X(80))+(X(81)+X(82))
     1           -2.0*(X(17)+X(18))-2.0*(X(21)+X(22)))/(4.0*P2D2M)
C*
700   CONTINUE
C*
C***  PRINT FEW OF THE FOURTH ORDER COEFFICIENTS
C*
      WRITE(LUN,800)
800   FORMAT( /,46X,'  * Few 4th order TRANSPORT coefficients * ' )
C*
      DO 900 I1=1,5
      WRITE(LUN,502)
      DO 900 I3=1,6
C*
      WRITE(LUN,801) (I1,I2,I2,I3,I3, Z4(I1,I2,I2,I3,I3), I2=1,I3 )
801   FORMAT( 6(I4,I2,I1,I1,I1,1PE11.3) )
900   CONTINUE
C*
      RETURN
      END
C*
