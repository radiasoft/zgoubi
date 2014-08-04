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
      SUBROUTINE READCO(NL,
     >                     KART,LET,YZXB,NDX,*,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ----------------------------------------------------
C     Look for and read coordinates, etc. of particle # NT
C     ----------------------------------------------------
      CHARACTER(1) LET
      INCLUDE 'MXVAR.H'
      DIMENSION YZXB(MXVAR),NDX(5)

      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      COMMON/LABCO/ ORIG(MXL,6) 
      COMMON/LUN/ NDAT,NRES,NPLT,NFAI,NMAP,NSPN
      INCLUDE 'MAXNTR.H'          
      COMMON/TRACKM/COOR(NTRMAX,9),NPTS,NPTR
      INCLUDE 'MAXCOO.H'
      COMMON/UNITS/ UNIT(MXJ-1) 
      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB

      PARAMETER (MXS=4)
      DIMENSION FO(MXJ),F(MXJ),SI(MXS),SF(MXS)
      INCLUDE 'MAXTRA.H'
      DIMENSION SX(MXT),SY(MXT),SZ(MXT)
      
      PARAMETER (LBLSIZ=10)
      CHARACTER(LBLSIZ) LBL1, LBL2

      PARAMETER (KSIZ=10)
      CHARACTER(KSIZ) KLEY
      CHARACTER(1) TX1

      LOGICAL BINARY,BINAR,OKKP,OKKT,OKKL
      DIMENSION XVEC(3), YVEC(3), ZVEC(3)

      CHARACTER(1) KLET, KLETO, KLETI

      SAVE KP1, KP2, KP3, BINARY
      SAVE KL1, KL2
      SAVE KT1, KT2
      SAVE KKEX, KLET

      SAVE MOD, RFR, RFR2 
      SAVE NOEL1, NOC
      SAVE XVEC, YVEC, ZVEC

      DATA MOD / 0 /
      DATA RFR, RFR2 / 0.D0, 0.D0 /

      DATA KP1, KP2, KP3 / 1, 99999999, 1 /
      DATA KT1, KT2 / 1, MXT /
      DATA KL1, KL2 / 1, 999999 /
      DATA KKEX, KLET / 1, '*' / 

      DATA NOEL1, NOC / -1, 0 /
      data xvec / 1.d0, 0.d0, 0.d0 /
      data yvec / 0.d0, 1.d0, 0.d0 /
      data zvec / 0.d0, 0.d0, 1.d0 /

      KT3 = 1

      IF(NL .EQ. NSPN) THEN
C--------- read in zgoubi.spn type storage file

          IMAX = 0
          IF(BINARY) THEN
 111         CONTINUE
            READ(NL,ERR=99,END=10) 
     >      KEX,(SI(J),J=1,4),(SF(J),J=1,4),F6I
     >      ,ENERG,IT,IMAX,IPASS,NOEL,KLEY,LBL1,LBL2,LET

            IF(.NOT. OKKT(KT1,KT2,KT3,IT,KEX,LET,
     >                                       IEND)) GOTO 111
            IF(.NOT. OKKP(KP1,KP2,KP3,IPASS,
     >                                       IEND)) THEN
              IF(IEND.EQ.0) THEN
                GOTO 111
              ELSEIF(IEND.EQ.1) THEN
                GOTO 91
              ENDIF
            ENDIF
            IF(.NOT. OKKL(KL1,KL2,NOEL,
     >                                       IEND)) GOTO 111

          ELSE
 1          READ(NL,101,ERR=99,END=10) 
     >      KEX,(SI(J),J=1,4),(SF(J),J=1,4),F6I
     >      ,ENERG,IT,IMAX,IPASS,NOEL
     >      ,TX1,KLEY,TX1,TX1,LBL1,TX1,TX1,LBL2,TX1,TX1,LET,TX1
            INCLUDE 'FRMSPN.H'

            IF(.NOT. OKKT(KT1,KT2,KT3,IT,KEX,LET,
     >                                       IEND)) GOTO 1        
            IF(.NOT. OKKP(KP1,KP2,KP3,IPASS,
     >                                       IEND)) THEN
              IF(IEND.EQ.0) THEN
                GOTO 1
              ELSEIF(IEND.EQ.1) THEN
                GOTO 91
              ENDIF
            ENDIF
            IF(.NOT. OKKL(KL1,KL2,NOEL,
     >                                       IEND)) GOTO 1

          ENDIF
        
        YZXB(6) = F6I* UNIT(5)     

        YZXB(21) = SF(1)
        YZXB(22) = SF(2)
        YZXB(23) = SF(3)
        YZXB(24) = SF(4)

        IF(IPASS.EQ. 1) THEN
          SX(IT) = 0.D0
          SY(IT) = 0.D0
          SZ(IT) = 0.D0
        ENDIF
        SX(IT) = SX(IT)+SF(1)
        SY(IT) = SY(IT)+SF(2)
        SZ(IT) = SZ(IT)+SF(3)
        YZXB(25) = SX(IT)/IPASS
        YZXB(26) = SY(IT)/IPASS
        YZXB(27) = SZ(IT)/IPASS

        YZXB(20) = ENERG
        YZXB(39) = IPASS 
        YZXB(57) = NOEL

      ELSE

        IF(NL .EQ. NFAI) THEN
C--------- read in zgoubi.fai type storage file

          IMAX = 0
          IF(BINARY) THEN
 222        CONTINUE
            READ(NL,ERR=99,END=10) 
     >      KEX,(FO(J),J=1,7),
     >      (F(J),J=1,7), 
     >      (SI(J),J=1,4),(SF(J),J=1,4),
     >      ENEKI, ENERG, 
     >      IT, IREP, SORT, AMQ1,AMQ2,AMQ3,AMQ4,AMQ5, RET, DPR, PS,
     >      BORO, IPASS, NOEL ,KLEY,LBL1,LBL2,LET

            IF(.NOT. OKKT(KT1,KT2,KT3,IT,KEX,LET,
     >                                IEND)) GOTO 222
            IF(.NOT. OKKP(KP1,KP2,KP3,IPASS,
     >                                IEND))  THEN
C              IF(IEND.EQ.0) THEN
                                             GOTO 222
C              ELSEIF(IEND.EQ.1) THEN
C                GOTO 91
C              ENDIF
            ENDIF
            IF(.NOT. OKKL(KL1,KL2,NOEL,
     >                                IEND)) GOTO 222
            IF(IEND.EQ.1) GOTO 91

          ELSE
 21         READ(NL,110,ERR=99,END=10)
     >      KEX,(FO(J),J=1,7),
     >      (F(J),J=1,7), 
     >      (SI(J),J=1,4),(SF(J),J=1,4),
     >      ENEKI, ENERG, 
     >      IT, IREP, SORT, AMQ1,AMQ2,AMQ3,AMQ4,AMQ5, RET, DPR,  PS,
     >      BORO, IPASS,NOEL, 
     >      TX1,KLEY,TX1,TX1,LBL1,TX1,TX1,LBL2,TX1,TX1,LET,TX1

            INCLUDE "FRMFAI.H"

            IF(.NOT. OKKT(KT1,KT2,KT3,IT,KEX,LET,
     >                                IEND)) GOTO 21
            IF(.NOT. OKKP(KP1,KP2,KP3,IPASS,
     >                                IEND))  THEN
C              IF(IEND.EQ.0) THEN
                                             GOTO 21
C              ELSEIF(IEND.EQ.1) THEN
C                GOTO 91
C              ENDIF
            ENDIF
            IF(.NOT. OKKL(KL1,KL2,NOEL,
     >                                IEND)) GOTO 21
            IF(IEND.EQ.1) GOTO 91

          ENDIF

        ELSEIF(NL .EQ. NPLT) THEN
C--------- read in zgoubi.plt type storage file

          IMAX = 0
          IF(BINARY) THEN
 232         CONTINUE
            READ(NL,ERR=99,END=10) 
     >      KEX,(FO(J),J=1,7),
     >      (F(J),J=1,7), BTI, DS, 
     >      KART, IT, IREP, SORT, XX, BX, BY, BZ, RET, DPR, PS,
     >      (SI(J),J=1,4),(SF(J),J=1,4),
     >      EX, EY, EZ, BORO, IPASS, NOEL, KLEY,LBL1,LBL2,LET
C            IF(LM .NE. -1) THEN
C              IF(LM .NE. NOEL) GOTO 232
C            ENDIF

            IF(.NOT. OKKT(KT1,KT2,KT3,IT,KEX,LET,
     >                                IEND)) GOTO 232
            IF(.NOT. OKKP(KP1,KP2,KP3,IPASS,
     >                                IEND))  THEN
C              IF(IEND.EQ.0) THEN
                                             GOTO 232
C              ELSEIF(IEND.EQ.1) THEN
C                GOTO 91
C              ENDIF
            ENDIF
            IF(.NOT. OKKL(KL1,KL2,NOEL,
     >                                IEND)) GOTO 232
            IF(IEND.EQ.1) GOTO 91

          ELSE
 31         READ(NL,100,ERR=99,END=10)
     >      KEX,(FO(J),J=1,MXJ),
     >      (F(J),J=1,MXJ), BTI, DS,
     >      KART, IT, IREP, SORT, XX, BX, BY, BZ, RET, DPR, PS,
     >      (SI(J),J=1,4),(SF(J),J=1,4),
     >      EX, EY, EZ, BORO, IPASS, NOEL,
     >      TX1,KLEY,TX1,TX1,LBL1,TX1,TX1,LBL2,TX1,TX1,LET,TX1
            INCLUDE "FRMPLT.H"
CCCCCCCCCCC           if(it.eq.1) yref = f(2)

C            IF(LM .NE. -1) THEN
C              IF(LM .NE. NOEL) GOTO 31
C            ENDIF

            IF(.NOT. OKKT(KT1,KT2,KT3,IT,KEX,LET,
     >                                IEND)) GOTO 31
            IF(.NOT. OKKP(KP1,KP2,KP3,IPASS,
     >                                IEND))  THEN
C              IF(IEND.EQ.0) THEN
                                             GOTO 31
C              ELSEIF(IEND.EQ.1) THEN
C                GOTO 91
C              ENDIF
            ENDIF
            IF(.NOT. OKKL(KL1,KL2,NOEL,
     >                                IEND)) GOTO 31
            IF(IEND.EQ.1) GOTO 91

          ENDIF
        ENDIF        !NL = NFAI, NPLT

C------- dp/p
        J = 1
        JU = 6
        YZXB(J)   =   F(J)   * UNIT(JU)     
C        YZXB(J)   =   1.D0 + F(J)  
        YZXB(J+10) =  FO(J)   * UNIT(JU)        ! dp/p_initial

        DO J=2,MXJ
          JU = J-1
C------- J=2,7 : Y, T, Z, P, S, Time
          YZXB(J)   =  F(J)   * UNIT(JU)     
CCCCCCC          if(j.eq.2) YZXB(J) = (f(j)-yref) * UNIT(JU)
C------- J=2,7 : Y_0, ..., Time_0
          YZXB(J+10) = FO(J)  * UNIT(JU) 
         ENDDO

C------- KART=1 : Cartesian coordinates, X is current x-coordinate (normally 
C        ranging in [XI,XF] as defined in quasex.f)
C        KART=2 : Cylindrical coordinates, X is current angle (normally 
C        ranging in [XI,XF] as defined in aimant.f)
C        write(*,*) '  * * * * * * * * ',xx
        YZXB(8) = XX

        IF(KART .EQ. 1) THEN
          YZXB(8) = YZXB(8) * UNIT(5)
        ELSE
C Rustine  RACCAM pour plot avec ffag-spi
          if(noel.ne.noel1) then
            if(kley .eq. 'FFAG-SPI') noc = noc+1
            noel1 = noel
          endif
          nbCell = 10
          pnCell = 4.d0 * atan(1.d0) / DBLE(nbCell)
c            write(*,*) '  noc, yzxb(8) ', noc, yzxb(8)
          YZXB(8) = YZXB(8) + pnCell * (2.d0*noc -1.d0) 
C          YZXB(2) = YZXB(2) +   DY * UNIT(5) 
        ENDIF

C         step size :
        YZXB(9) = DS       * UNIT(5)
C         r = sqrt(y^2+z^2) :
        YZXB(10) = SQRT(YZXB(2)*YZXB(2) + YZXB(4)*YZXB(4))
        YZXB(18) = RET
C------- (p_ps)/ps
        YZXB(19) = DPR            
C-------- momentum
C        YZXB(19) = BORO * (1.D0+F(1))*0.299792458D0   
        YZXB(20) = ENEKI
C Spin stories
        YZXB(21) = SF(1)
        YZXB(22) = SF(2)    
        YZXB(23) = SF(3)    
        YZXB(24) = SQRT(SF(1)**2 + SF(2)**2 + SF(3)**2)
C        YZXB(24) = SF(4)
C         convert B from kG to T
        YZXB(25) = acos(vscal(SI,SF,3)) * deg
        YZXB(26) = acos(vscal(YVEC,SF,3)) * deg

        YZXB(30) = BX      * .1D0
        YZXB(31) = BY      * .1D0
        YZXB(32) = BZ      * .1D0
        YZXB(33) = SQRT(BX*BX + BY*BY +  BZ*BZ) * .1D0

        YZXB(34) = EX * 1.d8    ! YZXB(34, 35, 36) in V/m
        YZXB(35) = EY * 1.d8 
        YZXB(36) = EZ * 1.d8
        YZXB(37) = SQRT(EX*EX + EY*EY +  EZ*EZ) * 1.d8

C AMAG is magnyfying factor for plotting of element synoptic and trajectories
C        CALL INSY1(
C     >             AMAG)        
C        PI2 = 0.5D0 * PI     !2.D0*ATAN(1.D0)
        YINL = F(2)* UNIT(1)  
        ZINL = F(4)* UNIT(3)  
C FM, Dec. 05       XINL = XX* UNIT(5) - ORIG(NOEL,5)
        IF    (KART .EQ. 1) THEN
          XINL = XX* UNIT(5) + ORIG(NOEL,5)
        ELSEIF(KART .EQ. 2) THEN
C--------- Plot from tracking in polar frame
          XINL = XX
          X = XINL
          Y = YINL
          IF( (KX .EQ. 48 .AND. KY .EQ. 42) ) THEN
C----------- Plot in Lab X-Y
            IF    (MOD.NE.0) THEN
C             mod=22 for 150MeV FFAG field map
C             mod=20 for RACCAM spiral field map
C             write(*,*) mod,rfr,'  ploter'
              IF(MOD.NE.20) THEN
                Y = Y + RFR
                TEMP = X
                X = Y * SIN(TEMP) 
                Y = Y * COS(TEMP) - RFR2
              ENDIF
            ELSEIF(MOD.EQ.0) THEN
C Example : SPES3 using DIPOLE-M
C           FFAG-SPI
              Y = Y + RFR
              TEMP = X 
C-------------------------
C Rustine  RACCAM pour plot avec ffag-spi
              if(noel.ne.noel1) then
                if(kley .eq. 'FFAG-SPI') noc = noc+1
                noel1 = noel
              endif
              nbCell = 10
              pnCell = 4.d0 * atan(1.d0) / DBLE(nbCell)
C               write(*,*) '  noc, yzxb(8) ', noc, yzxb(8)
              temp = temp + pnCell * (2.d0*noc -1.d0) 
C-------------------------
              X = Y * SIN(TEMP) 
              Y = Y * COS(TEMP)  - RFR2
            ENDIF
            XINL = X - ORIG(NOEL,5)
            YINL = Y
          ENDIF
        ENDIF
        YZXB(44) = ZINL 
        PHI = ORIG(NOEL,6)  
        CT = COS(ORIG(NOEL,4)+PHI) 
        ST = SIN(ORIG(NOEL,4)+PHI)
        YZXB(48) = ( XINL*CT - YINL*ST) + ORIG(NOEL,1)
        YZXB(42) = ( XINL*ST + YINL*CT) + ORIG(NOEL,2)

        CONTINUE

      ENDIF ! NL = NSPN, NFAI, NPLT

C      Location about where particle was lost
      YZXB(38) = SORT * 1.D-2
      YZXB(39) = IPASS 
      YZXB(57) = NOEL
      YZXB(58) = IT
      YZXB(59) = ENERG / AMQ1 *AMQ3   ! G.gamma
      YZXB(60) = SQRT(YZXB(2)**2 + YZXB(4)**2)
      alf = 0.d0 
      bta = 26.608d0
      YZXB(61) = atan2(alf*F(2)/1.d2+bta*F(3)/1.d3,F(2)/1.d2)
      alf = 0.d0 
      bta = 11.302651d0
      YZXB(62) = atan2(alf*F(4)/1.d2+bta*F(5)/1.d3,F(4)/1.d2)

C        write(77,*) KP1,KP2,KP3,IPASS,' readco'


C- For RACCAM design --------------------------
c       nCell = 8
c        pi = 4.d0 *atan(1.d0)
c      if(noel.ne.noel1) noc = noc+1
c      YZXB(62) = (Y+RFR) * SIN(XX + 2.d0 * pi / nCell * DBLE(noc-1))
c      YZXB(68) = (Y+RFR) * COS(XX + 2.d0 * pi / nCell * DBLE(noc-1))
C-----------------------------------------------

      NDX(1)=KEX
      NDX(2)=IT
      NDX(3)=IREP
      NDX(4)=IMAX
      NDX(5)=NOEL

      RETURN

 91   CONTINUE
      write(*,*) ipass,' AT 91,   readco'
      write(*,*) ipass,' AT 91,   readco'
      write(*,*) ipass,' AT 91,   readco'
      write(*,*) ipass,' AT 91,   readco'
      write(*,*) ipass,' AT 91,   readco'
      write(*,*) ipass,' AT 91,   readco'
      RETURN 1

C------------------ Pass # KP1 to KP2, ipass-modulo KP3
      ENTRY READC1(
     >              KP1O,KP2O,KP3O)
C------- Read pass #  KP1 to KP2 step KP3
      KP1O=KP1
      KP2O=KP2
      KP3O=KP3
      RETURN
C--
      ENTRY READC2(LN)
C------- Write pass #,  KP1 to KP2, ipass-modulo KP3
 12   WRITE(6,FMT='(''  Option status is now : KP1='',I9
     >  ,'',   KP2='',I9,'', ipass-modulo ='',I9)') KP1, KP2, KP3
        WRITE(6,FMT='(''    Expected data : '',
     >  /,10X,'' KP1>0,  KP2>=KP1 : will plot in range [KP1,KP2]''
     >  ,'',  ipass-modulo  KP3>0''
C     >  /,10X,''(ii) KP1=-1, KP2 > 0 : will plot all ipass-modulo KP2 ''
     >  )')
        WRITE(6,FMT='(/
     >  ,'' Enter desired values KP1>0, KP2>KP1, modulo-KP3>0''
     >  )')
        READ(LN,*,ERR=12) KP1NEW, KP2NEW, KP3NEW
        IF(KP1NEW.LE.0 .OR. KP2NEW.LE.0 .OR. KP3NEW.LE.0) GOTO 12
      GOTO 11
C--
      ENTRY READC2B(KP1W,KP2W,KP3W)
C------- Pass KP1 to pass KP2, KP-modulo
        KP1NEW=KP1W
        KP2NEW=KP2W
        KP3NEW=KP3W
 11     CONTINUE
        IF(KP1NEW.NE.0) KP1=KP1NEW
        IF(KP2NEW.NE.0) KP2=KP2NEW
        IF(KP3NEW.NE.0) KP3=KP3NEW
C          IF(KP1.NE.-1) THEN
            CALL PLOT31        ! OKS set to .FALSE.
            WRITE(6,*)
            WRITE(6,FMT='(
     >      '' Warning : although possibly requested (Menu-7/3/11),'')')
            WRITE(6,FMT='(''      reset of S coordinate upon '')')
            WRITE(6,FMT='(''      multiturn plot is NOW inhibited'')')
C            WRITE(6,FMT='(''      (only compatible with KP1=-1) '')')
C          ENDIF
      RETURN
C----------------------------

C------------------ Element #, KL1 to KL2
      ENTRY READC3(
     >             KL1O,KL2O)
C------- Read lmnt #,  KL1 to KL2
      KL1O=KL1
      KL2O=KL2
      RETURN
C--
      ENTRY READC4(LN)
        KLA = KL1
        KLB = KL2
C------- Write lmnt #,  KL1 to KL2
        WRITE(6,FMT='('' Observation now at elements  KL1='',I6,
     >   ''to   KL2='',I6)') KL1, KL2
        WRITE(6,FMT='(''     Available options for elment are : '',
     >   /,10X,'' 0 < KL1 < KL2  : will plot for elements '',
     >   ''in range [KL1,KL2] (assuming savings there))'',
     >   /,10X,'' KL1=-1, KL2 > 0 : will plot every KL2 other element ''
     >   ,''(assuming savings consistent'')')
        WRITE(6,FMT='(/,
     >        '' Enter desired values KL1, KL2  : '')')
        READ(LN,fmt=*,ERR=33) KL1NEW, KL2NEW
      GOTO 32
 33     KL1NEW = KLA
        KL2NEW = KLB
      GOTO 32
C--
      ENTRY READC4B(KL1W,KL2W)
C------- Lmnt  KL1 to lmnt KL2
        KL1NEW=KL1W
        KL2NEW=KL2W
 32     CONTINUE
        IF(KL1NEW.NE.0) KL1=KL1NEW
        IF(KL2NEW.NE.0) KL2=KL2NEW
      RETURN
C----------------------------

C------------------ Traj #, KT1 to traj KT2
      ENTRY READC5(
     >             KT1O,KT2O)
C------- Read traj.,  KT1 to KT2
        KT1O=KT1
        KT2O=KT2
      RETURN
C--
      ENTRY READC6(LN)
C------- Specify traj. #,  KT1 to KT2
 50     WRITE(6,FMT='(''  Option status is now : KT1='',I6,
     >   '',   KT2='',I6)') KT1,KT2
        WRITE(6,FMT='(''     Available options are : '',
     >  /,9X,'' KT1, KT2 > 0  : will treat range [KT1,KT2]'', 
     >  /,9X,'' KT1=-1, KT2 > 0 '',
     >                    '' : will treat every KT2 other particle '')')
        WRITE(6,FMT='(/,
     >        '' Enter desired values KT1, KT2  ( 0 0 to exit ) : '')')
        READ(LN,*,ERR=50) KT1NEW, KT2NEW
        IF(KT2NEW.GT.MXT) THEN
          WRITE(6,FMT='(9X,'' N2  cannot exceed '',I6)') MXT
          GOTO 50
        ENDIF
        IF(KT1NEW.GT.0) THEN
          IF(KT2NEW.LT.KT1NEW) GOTO 50
        ELSEIF(KT1NEW.NE.-1) THEN
          GOTO 50
        ENDIF
      GOTO 51
C--
      ENTRY READC6B(KT1W,KT2W)
C------- Traj KT1 to traj KT2
        KT1NEW=KT1W
        KT2NEW=KT2W
 51     CONTINUE
        IF(KT1NEW.NE.0) KT1=KT1NEW
        IF(KT2NEW.NE.0) KT2=KT2NEW
      RETURN
C----------------------------

      ENTRY READC7(BINAR)
      BINAR=BINARY
      RETURN

      ENTRY READC8(BINAR)
      BINARY=BINAR
      RETURN

      ENTRY READC9(
     >             KKEXO,KLETO)
        KKEXO=KKEX
        KLETO=KLET
      RETURN

      ENTRY READCA(KKEXI,KLETI)
        KKEX=KKEXI
        KLET=KLETI
      RETURN      

      ENTRY READCC(MODI,RFRI,RFR2I)
        MOD = MODI
        RFR = RFRI
        RFR2 = RFR2I
      RETURN

 10   continue
        noc = 0
        noel1 = 0
        RETURN 1
 99   continue
        noc = 0
        noel1 = 0
        RETURN 2      

      END
