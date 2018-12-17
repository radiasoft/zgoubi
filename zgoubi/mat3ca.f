
      SUBROUTINE MAT3CA(
     >XO,YO,VXO,VYO,PL_D,XM,TM,YM,PM,DM,J)
C     >XO,YO,ZO,VXO,VYO,VZO,RTL,RLL,PL_D,XM,TM,YM,PM,DM,J)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      DIMENSION XO(*),YO(*),ZO(*),VXO(*),VYO(*)
C     >, VZO(*),RTL(*),RLL(*),PL_D(*)
      DIMENSION XO(*),YO(*),VXO(*),VYO(*),PL_D(*)
C****
      INCLUDE "C.CDF.H"    ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.ABERR.H"  ! COMMON/ABERR/ R(6,6),T2(5,6,6),W3(5,6,6,6),Z4(5,6,6,6,6)
      INCLUDE "C.AV_STD.H" ! COMMON/AV_STD/ RFO(51,6,6)

        DIMENSION RM(6,6)
c****
c****
      R(1,1) =  ( 27.0*( XO(3)-XO(4) )-( XO(5)-XO(6) ) )/(16.0*XM)
      R(1,2) =  ( 27.0*( XO(7)-XO(8) )-( XO(9)-XO(10) ))/(16.0*TM)
C	IFLAGCO=1 CLOSED ORBIT IS NOT FOUND YET THUS CALCULATE
C	ONLY R11,R12,R21,R22
      IF(IFLAGCO.EQ.1) GO TO 101
      R(1,3) =  ( 27.0*(XO(11)-XO(12))-(XO(13)-XO(14) ))/(16.0*YM)
      R(1,4) =  ( 27.0*(XO(15)-XO(16))-(XO(17)-XO(18) ))/(16.0*PM)
      R(1,6) =  ( 27.0*(XO(19)-XO(20))-(XO(21)-XO(22) ))/(16.0*DM)
C*
101   CONTINUE
      R(2,1) =  (27.0*(VXO(3)-VXO(4) )-(VXO(5)-VXO(6) ) )/(16.0*XM)
      R(2,2) =  (27.0*(VXO(7)-VXO(8) )-(VXO(9)-VXO(10) ))/(16.0*TM)
      IF(IFLAGCO.EQ.1) RETURN
      R(2,3) = (27.0*(VXO(11)-VXO(12))-(VXO(13)-VXO(14) ))/(16.0*YM)
      R(2,4) = (27.0*(VXO(15)-VXO(16))-(VXO(17)-VXO(18) ))/(16.0*PM)
      R(2,6) = ( 27.0*(VXO(19)-VXO(20))-(VXO(21)-VXO(22) ))/(16.0*DM)
C*
      R(3,1) =  ( 27.0*( YO(3)-YO(4) )-( YO(5)-YO(6) ) )/(16.0*XM)
      R(3,2) =  ( 27.0*( YO(7)-YO(8) )-( YO(9)-YO(10) ))/(16.0*TM)
      R(3,3) =  ( 27.0*(YO(11)-YO(12))-(YO(13)-YO(14) ))/(16.0*YM)
      R(3,4) =  ( 27.0*(YO(15)-YO(16))-(YO(17)-YO(18) ))/(16.0*PM)
      R(3,6) =  ( 27.0*(YO(19)-YO(20))-(YO(21)-YO(22) ))/(16.0*DM)
C*
      R(4,1) =  (27.0*(VYO(3)-VYO(4) )-(VYO(5)-VYO(6) ) )/(16.0*XM)
      R(4,2) =  (27.0*(VYO(7)-VYO(8) )-(VYO(9)-VYO(10) ))/(16.0*TM)
      R(4,3) =  (27.0*(VYO(11)-VYO(12))-(VYO(13)-VYO(14) ))/(16.0*YM)
      R(4,4) =  (27.0*(VYO(15)-VYO(16))-(VYO(17)-VYO(18) ))/(16.0*PM)
      R(4,6) =  (27.0*(VYO(19)-VYO(20))-(VYO(21)-VYO(22) ))/(16.0*DM)
C*
      R(5,1) =  (27.0*(PL_D(3)-PL_D(4) )-(PL_D(5)-PL_D(6) ) )/(16.0*XM)
      R(5,2) =  (27.0*(PL_D(7)-PL_D(8) )-(PL_D(9)-PL_D(10) ))/(16.0*TM)
      R(5,3) = (27.0*(PL_D(11)-PL_D(12))-(PL_D(13)-PL_D(14) ))/(16.0*YM)
      R(5,4) = (27.0*(PL_D(15)-PL_D(16))-(PL_D(17)-PL_D(18) ))/(16.0*PM)
      R(5,5) =  1.0
      R(5,6) = (27.0*(PL_D(19)-PL_D(20))-(PL_D(21)-PL_D(22) ))/(16.0*DM)
C*
      R( 6,1 )  =  0.
      R( 6,2 )  =  0.
      R( 6,3 )  =  0.
      R( 6,4 )  =  0.
      R( 6,5 )  =  0.
      R( 6,6 )  =  1.
C****
c     To Calculate Average of matrix elements and stand dev.
      Do i1=1,6
         Do i2=1,6
         rfo(j,i1,i2)=R(i1,i2)
         rm(i1,i2)=R(i1,i2)
         enddo
      enddo
C****
      WRITE(NRES,22)
   22 FORMAT(/,40X, ' *1st order transform* (Units m, rad)'  )
      DO IR=1,6
        WRITE(NRES,221)  ( R(IR, IJ), IJ=1,6)
      ENDDO
  221 FORMAT( 12X, 6F10.5  )
C****
      RETURN
      END
