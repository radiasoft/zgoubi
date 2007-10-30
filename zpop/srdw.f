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
      SUBROUTINE SRDW(NLOG,KSC,KPL,NL,LM,NT,IT,OX,Q,AM,FE,HZ,YNRM,
     >  MNU,OKECH)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION OX(*)
      LOGICAL OKECH

      COMMON/CDF/ IES,IORDRE,LCHA,LIST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN
      COMMON/CONST/ CL,PI,DPI,RAD,DEG,QE,AH
      INCLUDE 'MXSTEP.H'
      PARAMETER (MSAM=1000)
      COMMON/SRTAB/ EFD(MXSTEP,3),D3W(MSAM,9),NPTS
      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB

      DIMENSION DW(MSAM,6), SDW(6), SYDX(4)

      CHARACTER XVA(3)*10, YVA(3)*15
      DATA XVA / 'omga',    'WaveL.' , 'omga'/
      DATA YVA / 'dW/dNu.dO', 'dW/dWaveL.dO', 'dN/dt.(do/o).dO' /
      CHARACTER XDI(3)*10, YDI(3)*15, SY(3)*12
      DATA XDI / '(keV)',    '(Mu-m)' , '(keV)'/
      DATA YDI / '(J/Hz.srd)', '(J/m.srd)' , '(Phot/s.BW.srd)' /
      DATA SY / ' (J/srd)' ,' (J/srd)', ' (Phot/srd)'/

      WRITE(6,*) ' SBR SRDW     Busy, calculate electric field'
      WRITE(6,*) '               May take a while...'

      CALL SREF(3,NL,LM,NT,IT,OX, 1, 2,Q,AM,FE, 0, 0,
     >                                        GAM,R,NOC,NRD,*79)

 79   CONTINUE
      WRITE(6,*) ' SBR SRDW     Busy, calculate spectrum'

C------- Calculate DW_1 = omga (KeV)
C                  DW_2,3 = dW_y,z/dNu.dOmga (J/Hz.srd)
      CALL SRDWC(ANU1,MNU,HZ,R,
     >                         DW)

      DO 1 IN=1,MNU
        D3W(IN,2) = 0.D0
        D3W(IN,3) = 0.D0
 1    CONTINUE

C--------- Store spectrum
C          Calculate energy density : SDW_2,3 = dW_x,y/dOmga  (J/srd)
      DOM = 1.D0
      DPH = 1.D0
      DPS = 1.D0
      NY = 1
      NZ = 1
      CALL SRDWST(MNU,NY,NZ,DW,HZ,DOM,DPH,DPS,
     >                                        D3W,SDW)

      IF(.NOT. OKECH) THEN

        WRITE(6,*) '               Busy, calculate scales'
        CALL SRDWSC(KSC,KPL,1, 2, 3,MNU,DW,YNRM,
     >                                          OKECH,IER)
        OKECH = IER .GE. 0

      ENDIF

      IF(OKECH) THEN

        WRITE(6,*) '               Busy, plotting'
        CALL SRAX(OKECH,KSC,XMI,XMA,YMI,YMA,XDI(KSC),YDI(KSC))
        CALL SRDWPL(NLOG,KSC,KPL,LM,NT,1, 2, 3,MNU,DW,YNRM,
     >                                                     SYDX)
        CALL SRPLV(NLOG,KSC,LM,NT,MNU,SYDX,HZ,
     >    XVA(KSC),XDI(KSC),YVA(KSC),YDI(KSC),SY(KSC),YNRM)

      ELSE
        WRITE(6,*) ' Scale error detected in SBR SRDWSC'
        WRITE(6,*) ' Could not plot E-field'
      ENDIF

      RETURN
      END
