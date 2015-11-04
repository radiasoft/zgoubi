        SUBROUTINE RAY102(p)
C****
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       dimension p(*)
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      INCLUDE "C.OBJET.H"     ! COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT

      DIMENSION XI(MXT), YI(MXT), VXI(MXT), VYI(MXT),DELP(MXT)
       save tmin, pmin
      save xmax, tmax, ymax, pmax, dmax

          data tmin, pmin / 1.d0, 1.d0 /

        NR = 102
        xmax = p(2)
        tmax = p(3)
        ymax = p(4)
        pmax = p(5)
        dmax = p(1)

        DO 1 I=1,NR
        XI(I)=0.d0
        YI(I)=0.d0
        VXI(I)=0.d0
        VYI(I)=0.d0
        DELP(I)=0.d0
1       CONTINUE        
        IF (TMIN.EQ.0.) TMIN=1.0d0
        IF (PMIN.EQ.0.) PMIN=1.0d0
        XMAX3 = XMAX/3.0
        YMAX3 = YMAX/3.0
        PMAX3 = PMAX/3.0
        DMAX3 = DMAX/3.0
        TMAX2 = TMAX/2.0
        TMAX3 = TMAX/3.0
        PMAX2 = PMAX/2.0
        PMAX3 = PMAX/3.0
        IF (NR.EQ.2) GO TO 2
        IF (NR.EQ.6) GO TO 2
        IF (NR.EQ.14) GO TO 2
        IF (NR.GT.14) GO TO 3
C        WRITE(NRES,100) NR
C        CALL EXIT
2       VXI(2)=TMIN
        VYI(2)=PMIN
        IF (NR.EQ.2) GO TO 5
        VXI(3)=TMAX2
        VXI(4)=-TMAX2
        VXI(5)=TMAX
        VXI(6)=-TMAX
        IF (NR.EQ.6) GO TO 5
        VYI(7)=PMAX2
        VXI(8)=TMAX2
        VYI(8)=PMAX2
        VXI(9)=-TMAX2
        VYI(9)=PMAX2
        VXI(10)=TMAX
        VYI(10)=PMAX2
        VXI(11)=-TMAX
        VYI(11)=PMAX2
        VYI(12)=PMAX
        VXI(13)=TMAX2
        VYI(13)=PMAX
        VXI(14)=-TMAX2
        VYI(14)=PMAX
C****
C****
C****
    5   DO 4 I=1,NR
        XI(I) = XMAX
        YI(I) = YMAX
    4   DELP(I) = DMAX
        RETURN
C****
C****
C****
3       VXI(2)=TMIN
        VYI(2)=PMIN
c*
        XI(3)=XMAX3
        XI(4)=-XMAX3
        XI(5)=XMAX
        XI(6)=-XMAX
C*
        VXI(7)=TMAX3
        VXI(8)=-TMAX3
        VXI(9)=TMAX
        VXI(10)=-TMAX
C*
        YI(11)=YMAX3
        YI(12)=-YMAX3
        YI(13)=YMAX
        YI(14)=-YMAX
C*
        VYI(15)=PMAX3
        VYI(16)=-PMAX3
        VYI(17)=PMAX
        VYI(18)=-PMAX
C*
        DELP(19)=DMAX3
        DELP(20)=-DMAX3
        DELP(21)=DMAX
        DELP(22)=-DMAX
C*
C*
        XI(23)=XMAX3
        VXI(23)=TMAX3
        XI(24)=-XMAX3
        VXI(24)=-TMAX3
        XI(25)=XMAX
        VXI(25)=TMAX
        XI(26)=-XMAX
        VXI(26)=TMAX
        XI(27)=XMAX
        VXI(27)=-TMAX
        XI(28)=-XMAX
        VXI(28)=-TMAX
C*
        XI(29)=XMAX3
        YI(29)=YMAX3
        XI(30)=-XMAX3
        YI(30)=-YMAX3
        XI(31)=XMAX
        YI(31)=YMAX
        XI(32)=-XMAX
        YI(32)= YMAX
        XI(33)= XMAX
        YI(33)=-YMAX
        XI(34)=-XMAX
        YI(34)=-YMAX
C*
        XI(35)= XMAX3
        VYI(35)=PMAX3
        XI(36) = -XMAX3
        VYI(36)= -PMAX3
        XI(37) = XMAX
        VYI(37)= PMAX
        XI(38) =-XMAX
        VYI(38)= PMAX
        XI(39) = XMAX
        VYI(39)=-PMAX
        XI(40) =-XMAX
        VYI(40)=-PMAX
C*
        XI(41) = XMAX3
        DELP(41)= DMAX3
        XI(42) =-XMAX3
        DELP(42)=-DMAX3
        XI(43) = XMAX
        DELP(43)= DMAX
        XI(44) =-XMAX
        DELP(44)= DMAX
        XI(45) = XMAX
        DELP(45)=-DMAX
        XI(46) =-XMAX
        DELP(46)=-DMAX
C*
        VXI(47)= TMAX3
        YI(47)= YMAX3
        VXI(48)=-TMAX3
        YI(48)=-YMAX3
        VXI(49)= TMAX
        YI(49)= YMAX
        VXI(50)=-TMAX
        YI(50)= YMAX
        VXI(51)= TMAX
        YI(51)=-YMAX
        VXI(52)=-TMAX
        YI(52)=-YMAX
C*
        VXI(53)= TMAX3
        VYI(53)= PMAX3
        VXI(54)=-TMAX3
        VYI(54)=-PMAX3
        VXI(55)= TMAX
        VYI(55)= PMAX
        VXI(56)=-TMAX
        VYI(56)= PMAX
        VXI(57)= TMAX
        VYI(57)=-PMAX
        VXI(58)=-TMAX
        VYI(58)=-PMAX
C*
        VXI(59) = TMAX3
        DELP(59)= DMAX3
        VXI(60) =-TMAX3
        DELP(60)=-DMAX3
        VXI(61) = TMAX
        DELP(61)= DMAX
        VXI(62) =-TMAX
        DELP(62)= DMAX
        VXI(63) = TMAX
        DELP(63)=-DMAX
        VXI(64) =-TMAX
        DELP(64)=-DMAX
C*
        YI(65) = YMAX3
        VYI(65)= PMAX3
        YI(66) =-YMAX3
        VYI(66)=-PMAX3
        YI(67) = YMAX
        VYI(67)= PMAX
        YI(68) =-YMAX
        VYI(68)= PMAX
        YI(69) = YMAX
        VYI(69)=-PMAX
        YI(70) =-YMAX
        VYI(70)=-PMAX
C*
        YI(71) = YMAX3
        DELP(71)=DMAX3
        YI(72)  =-YMAX3
        DELP(72)=-DMAX3
        YI(73)  = YMAX
        DELP(73)= DMAX
        YI(74)  =-YMAX
        DELP(74)= DMAX
        YI(75)  = YMAX
        DELP(75)=-DMAX
        YI(76)  =-YMAX
        DELP(76)=-DMAX
C*
        VYI(77) = PMAX3
        DELP(77)= DMAX3
        VYI(78) =-PMAX3
        DELP(78)=-DMAX3
        VYI(79) = PMAX
        DELP(79)= DMAX
        VYI(80) =-PMAX
        DELP(80)= DMAX
        VYI(81) = PMAX
        DELP(81)=-DMAX
        VYI(82) =-PMAX
        DELP(82)=-DMAX
C*
        XI(83) = XMAX
        VXI(83)= TMAX
        YI(83) = YMAX
        XI(84) =-XMAX
        VXI(84)=-TMAX
        YI(84) =-YMAX
C*
        XI(85) = XMAX
        VXI(85)= TMAX
        VYI(85)= PMAX
        XI(86) =-XMAX
        VXI(86)=-TMAX
        VYI(86)=-PMAX
C*
        XI(87) = XMAX
        VXI(87)= TMAX
        DELP(87)= DMAX
        XI(88)  =-XMAX
        VXI(88) =-TMAX
        DELP(88)=-DMAX
C*
        XI(89)  = XMAX
        YI(89)  = YMAX
        VYI(89) = PMAX
        XI(90)  =-XMAX
        YI(90)  =-YMAX
        VYI(90) =-PMAX
C*
        XI(91)  = XMAX
        YI(91)  = YMAX
        DELP(91)= DMAX
        XI(92)  =-XMAX
        YI(92)  =-YMAX
        DELP(92)=-DMAX
C*
        XI(93) = XMAX
        VYI(93) = PMAX
        DELP(93)= DMAX
        XI(94)  =-XMAX
        VYI(94) =-PMAX
        DELP(94)=-DMAX
C*
        VXI(95) = TMAX
        YI(95)  = YMAX
        VYI(95) = PMAX
        VXI(96) =-TMAX
        YI(96)  =-YMAX
        VYI(96) =-PMAX
C*
        VXI(97) = TMAX
        YI(97)  = YMAX
        DELP(97)= DMAX
        VXI(98) =-TMAX
        YI(98)  =-YMAX
        DELP(98)=-DMAX
C*
        VXI(99) = TMAX
        VYI(99) = PMAX
        DELP(99)= DMAX
        VXI(100) =-TMAX
        VYI(100) =-PMAX
        DELP(100)=-DMAX
C*
        YI(101)  = YMAX
        VYI(101) = PMAX
        DELP(101)= DMAX
        YI(102)  =-YMAX
        VYI(102) =-PMAX
        DELP(102)=-DMAX
C*

        do i = 1, nr
            fo(2,i) = xi(i) 
            fo(3,i) = vxi(i) 
            fo(4,i) = yi(i) 
            fo(5,i) = vyi(i) 
            fo(6,i) = 0.d0
            fo(1,i) = delp(i)
        enddo

      RETURN
      END
