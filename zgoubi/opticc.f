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
      SUBROUTINE OPTICC(NOEL,PRDIC,OKCPLD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL OKCPLD, PRDIC

      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      DIMENSION R(6,6), F0(6,6), AKL(3)
      PARAMETER (KSIZ=10)
      CHARACTER(KSIZ) KLE
      LOGICAL OKLNO, OKLNOI
      SAVE OKLNO
      SAVE LNOPT
      SAVE DXMA, DYMA, XMA, YMA, BTXMA, BTYMA,
     >              XM, YM, DXM, DYM,
     >              XM2, YM2, DXM2, DYM2, DLTP, NC
      DIMENSION RSAV(6,6), RLOC(6,6), NW(6), BW(6)
      SAVE RSAV
      PARAMETER (PI = 4.D0 *(ATAN(1.D0)))

      DATA OKLNO / .FALSE. /

      DATA DXMA, DYMA / -1.D10, -1.D10 /
      DATA DXMI, DYMI / 1.D10, 1.D10 /
      DATA XMA, YMA / -1.D10, -1.D10 /
      DATA XMI, YMI / 1.D10, 1.D10 /
      DATA BTXMA, BTYMA / -1.D0, -1.D0 /
      DATA BTXMI, BTYMI / 99999.D0, 99999.D0 /
      DATA XM, YM / 0.D0, 0.D0 /
      DATA XM2, YM2 / 0.D0, 0.D0 /
      DATA DXM, DYM / 0.D0, 0.D0 /
      DATA DXM2, DYM2 / 0.D0, 0.D0 /
      DATA DLTP / 0.D0 /

      DATA NC / 0 /
C      DATA ((RSAV(I,J),I=1,6),J=1,6) / 36*0.D0 /
      DATA (RSAV(I,I),I=1,6) / 6*1.D0 /

      IF(NRES.GT.0)
     >WRITE(NRES,FMT='(/,''---------------------------------------------
     > Local OPTICS monitoring : '')')

      IER = 0
      IORD = 1
      IFOC = 0
      KWR = 0
      CALL MATRIC(IORD,IFOC,KWR,OKCPLD,
     >                                 IER)
      IF(IER.NE.0) GOTO 99
      CALL MATRI1(
     >            R)
c        WRITE(NRES,*)
c        WRITE(NRES,*)  '  R_SAV : '
c        WRITE(NRES,104) (( Rsav(IA,IB) , IB=1,6) , IA=1,6)

      CALL ALAIN(6,6,RSAV,NW,BW,IER)
      RLOC = MATMUL(R,RSAV)

      IF(NRES.GT.0) THEN
        WRITE(NRES,113)
 113    FORMAT(//,18X,'TRANSFER  MATRIX  OF  LAST  ELEMENT',
     >  '  (MKSA units)',/)
        WRITE(NRES,104) (( RLOC(IA,IB) , IB=1,6) , IA=1,6)
 104    FORMAT(6X,1P,6G16.6)
        WRITE(NRES,112) (RLOC(1,1)*RLOC(2,2)-RLOC(1,2)*RLOC(2,1))-1.D0,
     >  (RLOC(3,3)*RLOC(4,4)-RLOC(3,4)*RLOC(4,3))-1.D0
112     FORMAT(/,10X,'DetY-1 = ',F18.10,',',4X,'DetZ-1 = ',F18.10)
      ENDIF
      DO J = 1, 6
        DO I = 1, 6
          RSAV(I,J) = R(I,J)
        ENDDO
      ENDDO

      CALL REFER3(
     >            XI,YI,ALE,PATHL,ZE,PE)
      CALL BEAMAT(R,PRDIC,OKCPLD,
     >                           F0,PHY,PHZ,CSTRN,RPRM)
      NC = NC + 1
      DX = F0(1,6)
      DY = F0(3,6)
      IF(F0(1,6) .GT. DXMA) DXMA = DX
      IF(F0(3,6) .GT. DYMA) DYMA = DY
      IF(YI .GT. XMA) XMA = YI
      IF(ZE .GT. YMA) YMA = ZE
      IF(F0(1,1) .GT. BTXMA) BTXMA = F0(1,1)
      IF(F0(3,3) .GT. BTYMA) BTYMA = F0(3,3)
      IF(F0(1,6) .LT. DXMI) DXMI = DX
      IF(F0(3,6) .LT. DYMI) DYMI = DY
      IF(YI .LT. XMI) XMI = YI
      IF(ZE .LT. YMI) YMI = ZE
      IF(F0(1,1) .LT. BTXMI) BTXMI = F0(1,1)
      IF(F0(3,3) .LT. BTYMI) BTYMI = F0(3,3)
      XM = XM + YI
      YM = YM + ZE
      XM2 = XM2 + YI*YI
      YM2 = YM2 + ZE*ZE
      DXM = DXM + DX
      DYM = DYM + DY
      DXM2 = DXM2 + DX*DX
      DYM2 = DYM2 + DY*DY

      CALL BEAIMP(F0,PHY,PHZ)    ! print to zgoubi.res

      CALL ZGKLEY(KLE)
      IF(KLE.EQ.'AGSMM') THEN
         CALL AGSMKL(
     >        AL, AK1, AK2, AK3)
         AKL(1) = AK1 * AL *1.D-2
         AKL(2) = AK2 * AL *1.D-2
         AKL(3) = AK3 * AL *1.D-2
         DEV = 2.D0 * PI /240.
      ELSEIF(KLE.EQ.'AGSQUAD') THEN
         CALL AGSQKL(
     >        AL, AK1)
         AKL(1) = 0.D0
         AKL(2) = AK1 * AL *1.D-1
         DEV = 0.D0
      ELSEIF(KLE.EQ.'MULTIPOL') THEN
         CALL MULTKL(
     >        AL, AK1, AK2, AK3, DEVI)
         AKL(1) = AK1 * AL *1.D-2
         AKL(2) = AK2 * AL *1.D-2
         AKL(3) = AK3 * AL *1.D-2
         DEV = DEVI
      ELSEIF(KLE.EQ.'BEND') THEN
         CALL BENKL(
     >              DEVI)
         AKL = 0.D0
         DEV = DEVI   ! rad
      ELSE
         AKL = 0.D0
         DEV = 0.D0
      ENDIF

      IF(OKLNO)
     > CALL OPTIMP(LNOPT,NOEL,F0,PHY,PHZ,AKL,CSTRN,RPRM,R,
     > AL,DEV,                     ! print to zgoubi.OPTICS.out (OPTICS keyword)
     >     PP0)                    ! or to zgoubi.TWISS.out (TWISS keyword)

      RETURN

 99   CONTINUE
      IF(NRES.GT.0) WRITE(NRES,*) 'Pgm opticc.  Exit upon IER=-1.'
     >//' Matrix cannot be computed. Hint : check OBJET/KOBJ=5 or 6.'
      RETURN

      ENTRY OPTIC1(
     >              DXMAO, DYMAO, DXMIO, DYMIO,
     >              XMAO, YMAO, XMIO, YMIO,
     >              BTXMAO, BTYMAO, BTXMIO, BTYMIO,
     >              XMO, YMO, DXMO, DYMO,
     >              XM2O, YM2O, DXM2O, DYM2O, DLTPO, NCO)
      DXMAO = DXMA
      DYMAO = DYMA
      XMAO = XMA
      YMAO = YMA
      BTXMAO = BTXMA
      BTYMAO = BTYMA
      DXMIO = DXMI
      DYMIO = DYMI
      XMIO = XMI
      YMIO = YMI
      BTXMIO = BTXMI
      BTYMIO = BTYMI
      XMO = XM
      YMO = YM
      DXMO = DXM
      DYMO = DYM
      XM2O = XM2
      YM2O = YM2
      DXM2O = DXM2
      DYM2O = DYM2
      DLTPO = PP0
      NCO = NC
      RETURN

      ENTRY OPTIC2(OKLNOI,LNOPTI)
      OKLNO = OKLNOI
      LNOPT = LNOPTI
      RETURN

      END
