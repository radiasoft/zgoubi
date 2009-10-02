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
      SUBROUTINE SRDWM(NLOG,NL,LM,KPL,OX,WF,Q,AM,FE,HZ,YNRM,
     >  NOMFIC,OKECH,OKOPN,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION OX(*), WF(*)
      CHARACTER*(*) NOMFIC
      LOGICAL OKECH, OKOPN

      COMMON/CONST/ CL,PI,DPI,RAD,DEG,QE,AH
      PARAMETER (MSAM=1000)

      LOGICAL INPECH, OKBIN
      INCLUDE 'FILPLT.H'
C      DATA KPL / 23 /

      ANU1 = WF(1)
      ANU2 = WF(2)
      MNU  = WF(3)

      GOTO 21

 20   CONTINUE      
      CALL FBGTXT
      WRITE(6,FMT='(/,''  Press RETURN for more'')') 
      READ(5,200,ERR=20) REP 
 200  FORMAT(A1)

 21   CONTINUE
      CALL FBGTXT
C      CALL HOMCLR

      WRITE(6,100) ANU1, ANU2, MNU, (OX(I), I=1, 3)
 100  FORMAT(//,3X,60('-'),//,20X,' MENU  - Spectrum :' ,//
     1 ,5X,'  1   Frequency range and sampling:',/
     > ,5X,'        Frq1, Frq2:',1P,2G12.4,' (keV),  Nfreq =',I6,/
     2 ,5X,'  2    Observer position w.r.t magnet frame:',/
     > ,5X,'         X, Y, Z = ',3G12.4,' m',/
     3 ,5X,'  3    Plot  options',/
     4 ,5X,'  4    Hand scales ',/
     5 ,5X,'  5    Automatic scales ',/
     7 ,5X,'       **    Available PLOTS : **',/
     7 ,5X,'  71     dW/dNu.dOmega  v.s.  omega',/
     7 ,5X,'  72     dW/dWaveL.dOmga  v.s.  WaveL',/
     7 ,5X,'  73     dN/dt.(domga/omga).dOmga  v.s.  omega',/
     8 ,5X,'  8    Print screen',/
     9 ,5X,'  9    EXIT  THIS  MENU     ',/
     2 ,5X,' 12    ERASE  DISPLAY    ',/
     2 ,3X,60('-'),//)

      WRITE(6,FMT='(''$  Option  number : '')')
      READ(5,FMT='(I2)',ERR=21) IOPT
      IF(IOPT .GE. 71 .AND. IOPT .LE. 73) THEN
        IOP  = IOPT - 70
        IOPT = 7
      ELSE
        IOP = 0
      ENDIF        
      GOTO ( 1, 2, 3, 4, 5,21, 7, 8,99,21,21,12) IOPT  
      GOTO 21

 1    WRITE(6,*)
      WRITE(6,*) ' 1micro-m = 1.24E-3 keV '
      WRITE(6,*) ' Visible range: 0.8     ->  0.4micro-m'
      WRITE(6,*) '              = 1.55E-3 ->  3.1e-3 keV'
      WRITE(6,*) ' Give Nu-min, Nu-max (keV) '
      READ(5,*,ERR=1) ANU1, ANU2
      IF(ANU1 .GE. ANU2) GOTO 1
 11   WRITE(6,*) ' Give number of samples (<',MSAM,') '
      READ(5,*,ERR=11) MNU
      IF(MNU .GT. MSAM) GOTO 11
      WF(1) = ANU1
      WF(2) = ANU2
      WF(3) = MNU
      GOTO 21

 2    CONTINUE
 22   WRITE(6,*) 
      WRITE(6,*) ' Observer position  X,Y,Z (m):'
      READ(5,*,ERR=22) (OX(I),I=1,3)
C      OKECH = .FALSE.
      GOTO 21

 3    CONTINUE
        CALL PLTOPT(
     >              LM,OKECH,OKBIN,KOPT)
        IF(KOPT.EQ.8) THEN
          CALL READC5(KT1,KT2)
          IF(KT2.NE.KT1) THEN   
            WRITE(6,*)
     >      ' *** WARNING : SR utility will plot trajectory #KT1= ',
     >                                                KT1,' only'
            CALL READC6B(KT1,KT1)
          ENDIF
          CALL XYZBR(NL)
          GOTO 20
        ELSEIF(KOPT.EQ.88) THEN
 31       WRITE(6,*) ' Plot Ey, Ez, both (23,2,3)'
          READ(5,*,ERR=31) KPL 
          IF(KPL .NE. 2 .AND. KPL .NE. 3 .AND. KPL .NE. 23) KPL=23
        ENDIF
      GOTO 21 

 4    CONTINUE
        OKECH=INPECH()
      GOTO 21
        
 5    CONTINUE
        IF(.NOT. OKOPN) CALL OPNDEF(NPLT,FILPLT,NL,
     >                                             NOMFIC,OKOPN) 
        IF(.NOT. OKOPN) RETURN 1
      GOTO 20
         
 7    CONTINUE

        IF(.NOT. OKOPN) CALL OPNDEF(NPLT,FILPLT,NL,
     >                                             NOMFIC,OKOPN) 
        IF(.NOT. OKOPN) RETURN 1

        KSC = IOP
        CALL SRDNU(KSC,WF)

        GOTO (71,72,73) IOP

  71      CONTINUE
C----------- dW/dNu.dOmga (J/Hz.srd) v.s. omga(keV)
            CALL SRDW(NLOG,KSC,KPL,OX,Q,AM,FE,HZ,YNRM,
     >        MNU,OKECH)
          GOTO 20

  72      CONTINUE
C----------- dW/dWaveL.dOmga (J/m.srd) v.s. WaveL(micro-m)
            IF(ANU1*ANU2 .EQ. 0.D0) THEN
              WRITE(6,*) ' Infinite Wavelength in scales definition,'
              WRITE(6,*) ' unacceptable - Change it'
            ELSE
              CALL SRDW(NLOG,KSC,KPL,OX,Q,AM,FE,HZ,YNRM,
     >          MNU,OKECH)
            ENDIF
          GOTO 20

  73      CONTINUE
C----------- dN/dt.(domga/omga).dOmga (phot/s.BW.srd) 
C                        v.s. omga(keV)
            CALL SRDW(NLOG,KSC,KPL,OX,Q,AM,FE,HZ,YNRM,
     >        MNU,OKECH)

      GOTO 20

 8    CONTINUE
C        CALL MENVCF
        CALL SAVPLT
      GOTO 21

 12   CONTINUE
        CALL CLSCR
        OKECH=.FALSE.
      GOTO 21

 99   RETURN
      END
