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
      SUBROUTINE SRINM(NLOG,NL,LM,KPL,OX,WF,WO,Q,AM,FE,HZ,YNRM,
     >  NOMFIC,OKOPN,GNUFIL,STORE,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION OX(*), WF(*), WO(*)
      CHARACTER*(*) NOMFIC, GNUFIL
      LOGICAL OKOPN, STORE
C     -------------------------------------------------------
C     Menu for integral calculations: dP/do, dP/dPhidPsi, etc
C     -------------------------------------------------------
      COMMON/CONST/ CL,PI,DPI,RAD,DEG,QE,AH
      LOGICAL OKECH, OKVAR, OKBIN
      COMMON/ECHL/OKECH, OKVAR, OKBIN
      INCLUDE 'MXVAR.H'
      CHARACTER KVAR(MXVAR)*7, KPOL(2)*9, KDIM(MXVAR)*7 
      COMMON/INPVR/ KVAR, KPOL, KDIM
      COMMON/LUN/NDAT,NRES,NPLT,NFAI,NMAP,NSPN

      PARAMETER (MSAM=1000)

      LOGICAL INPECH
      CHARACTER REP
      LOGICAL CHANGE

      CHARACTER * (9)   DMY, HMSI, HMSF

      SAVE KSC, LUO
      INCLUDE 'FILPLT.H'

C      CALL LINTYP(1)

      GOTO 21

 20   CONTINUE      
      CALL FBGTXT
      WRITE(6,FMT='(/,''  Press RETURN for more'')') 
      READ(5,200,ERR=20) REP 
 200  FORMAT(A1)

 21   CONTINUE
      CALL FBGTXT
C      CALL HOMCLR

      II = WO(5)
      II = 2 * ( II/2) + 1
      WO(5) = II
      II = WO(6)
      II = 2 * ( II/2) + 1
      WO(6) = II

      WRITE(6,100) WF(1), WF(2), INT(WF(3)), (OX(I),I=1,3),
     >WO(2), WO(3), INT(WO(5)), INT(WO(6))
 100  FORMAT(//,3X,60('-'),//,20X,' Integrals. MENU:' ,//
     1,/,5X,' 1   Frequency range and sampling:'
     >,/,5X,'       Frq1, Frq2:',1P,2G12.4,' (keV),  Nfreq =',I6
     2,/,5X,' 21  Observation window:'
     >,/,5X,'       Position X,Y,Z   : ',3G12.4,' (M)'
     >,/,5X,'       Widths WY,WZ(+/-): ',2G12.4,' (M)'
     >,/,5X,' 22    Sampling NY,NZ   : ',2I6
     3,/,5X,' 3     Plot options'
     4,/,5X,' 4     Hand scales'
     6,/,5X,'     **    Available PLOTS : **'
     7,/,5X,' 71    dW/dNu  v.s.  omga'
     7,/,5X,' 72    dW/dWaveL  v.s.  Wavelength'
     7,/,5X,' 73    dN/dt.(domga/omga)  v.s.  omga'
     7,/,5X,' 74    dW/dPhi  v.s.  Phi'
     7,/,5X,' 75    dW/dPsi  v.s.  Psi'
     8,/,5X,' 8   Print screen'
     9,/,5X,' 9   EXIT  THIS  MENU'
     2,/,5X,'12   ERASE  DISPLAY'
     >,/,3X,60('-'),//)

      WRITE(6,FMT='(''$  Option  number : '')')
      READ(5,FMT='(I2)',ERR=21) IOPT
      IF(IOPT .GE. 21 .AND. IOPT .LE. 22) THEN
        IOP  = IOPT - 20
        IOPT = 2
      ELSEIF(IOPT .GE. 71 .AND. IOPT .LE. 75) THEN
        IOP  = IOPT - 70
        IOPT = 7
      ELSE
        IOP = 0   
      ENDIF        
      GOTO ( 1, 2, 3, 4,21,21, 7, 8, 9,21,21,12) IOPT  
      GOTO 21    

 1    CONTINUE
 11     WRITE(6,*) 
        WRITE(6,*) ' Min and max frequencies (keV):'
        READ(5,*,ERR=11) WF(1), WF(2)
        IF(WF(1) .GE. WF(2)) GOTO 11
 19     WRITE(6,*) ' Sampling (<',MSAM,'):'
        READ(5,*,ERR=19) WF(3)
        IF(WF(3) .GT. MSAM) GOTO 19
        OKECH = .FALSE.
      GOTO 21

 2    CONTINUE
        GOTO(22,23) IOP
 
 22       WRITE(6,*) 
          WRITE(6,*) ' Window position X, Y, Z (m):'
          READ(5,*,ERR=22) OX(1), OX(2), OX(3)
          WRITE(6,*) ' Window width  WY, WZ (+/-, m):'
          READ(5,*,ERR=22) WO(2), WO(3)
          IF(IOP .EQ. 1) GOTO 21

 23       CONTINUE
          WRITE(6,*) ' Sampling NY, NZ (reasonable odd integers):'
          READ(5,*,ERR=23) WO(5), WO(6)
          OKECH = .FALSE.
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
        ELSEIF(KOPT.EQ.88) THEN
          IP = 0
          WRITE(6,*) ' Plot intermediate graphs E(t) ? (N/Y)'
          READ(5,FMT='(A1)',ERR=99) REP
          IF(REP .EQ. 'Y' .OR. REP .EQ. 'y') IP = 1
          WRITE(6,*) ' Store d2W/dNudOmga(Phi,Psi) in GNUFILe ? (N/Y)'
          READ(5,FMT='(A1)',ERR=99) REP
C          IF(REP .EQ. 'Y' .OR. REP .EQ. 'y') STORE = .TRUE.
          STORE = REP .EQ. 'Y' .OR. REP .EQ. 'y'

          IF(STORE) THEN
            CALL TIME2(HMSI)
            CALL OPNGNU(GNUFIL,
     >                         STORE,LUO)
          ENDIF

 32       WRITE(6,*) ' Plot Ey, Ez, both (23,2,3)'
          READ(5,*,ERR=32) KPL 
          IF(KPL .NE. 2 .AND. KPL .NE. 3 .AND. KPL .NE. 23) KPL=23
        ENDIF
      GOTO 20
        
 4    CONTINUE
        OKECH=INPECH()
      GOTO 21
                
 7    CONTINUE

 79     WRITE(6,FMT='(/,''  New calculation (N/Y) ?'')') 
        READ(5,200,ERR=79) REP 
C        IF(REP .EQ. 'Y' .OR. REP .EQ. 'y') CHANGE = .TRUE.
        CHANGE = REP .EQ. 'Y' .OR. REP .EQ. 'y'

C------- KSC = 1 : dW/dNu (J/Hz) v.s. omga (keV)
C------- KSC = 2 : dW/dWaveL v.s. Wavelength
C------- KSC = 3 : dN/dt.(domga/omga) (phot/s.BW)  v.s.  omga(keV)
C------- KSC = 4 : dW/dPhi (J/rad) v.s. Gamma.Psi (vertical angle, rad)
C------- KSC = 5 : dW/dPsi (J/rad) v.s. Gamma.Psi (vertical angle, rad)
        KSC = IOP
        IF(IOP .EQ. 0) KSC = 1
        
        IF(KSC .EQ. 2) THEN
          IF(WF(1)*WF(2) .EQ. 0.D0) THEN
            WRITE(6,*) ' Infinite Wavelength in scales definition,'
            WRITE(6,*) '      unacceptable - change it'
          ENDIF
        ENDIF

        CALL SRDNU(KSC,WF)

        IF(CHANGE) THEN

          IF(.NOT. OKOPN) CALL OPNDEF(NPLT,FILPLT,NL,
     >                                               NOMFIC,OKOPN) 
          IF(.NOT. OKOPN) RETURN 1

          CALL PRTIME(1,
     >                  DMY,HMSI,HMSF)
          IF(STORE) THEN
              CALL OPNGNU(GNUFIL, 
     >                           STORE,LUO)
              WRITE(6,FMT=
     >        '('' SR integrals stored in '',A,'',  unit #'',I3,/)') 
     >        GNUFIL, LUO

              WRITE(LUO,101) WF(1), WF(2), INT(WF(3)), (OX(I),I=1,3),
     >          WO(2), WO(3), INT(WO(5)), INT(WO(6))
 101          FORMAT(
     1             '#  Frequency range and sampling:'
     >          ,/,'#    Frq1, Frq2:',1P,2G12.4,' (keV),  Nfreq =',I6
     2          ,/,'#  Observation window:'
     >          ,/,'#    Position X,Y,Z   : ',3G12.4,' (M)'
     >          ,/,'#    Widths WY,WZ(+/-): ',2G12.4,' (M)'
     >          ,/,'#    Sampling NY,NZ   : ',2I6,/,'#') 
              WRITE(LUO,FMT=
     >  '('' Phi, Psi, dWy,z/dOmga (*YNRM)'',/,''# YNRM = '',1P,G14.6)')
     >        YNRM
          ELSE
            WRITE(6,FMT='('' SR integrals will not be stored'')')
          ENDIF      

          CALL SRINC
     >      (NLOG,KSC,KPL,IP,OX,WF,WO,Q,AM,FE,HZ,STORE,LUO,
     >                                                   YNRM,NOC,NRD)
          CALL PRTIME(2,
     >                  DMY,HMSI,HMSF)
          CHANGE = .FALSE.
        ENDIF

        CALL SRINPL(NLOG,KSC,KPL,IP,OX,WF,WO,Q,AM,FE,HZ,
     >                                                 YNRM,NOC,NRD)
        IF( STORE ) THEN
          WRITE(6,FMT='(/,'' SR integrals stored in '',A)') GNUFIL
        ELSE
          WRITE(6,FMT='(/,'' SR integrals not stored '',A)')
        ENDIF
      GOTO 20

 8    CONTINUE
C        CALL MENVCF
        CALL SAVPLT
      GOTO 21

 9    CONTINUE
C        OKECH = .FALSE.
      GOTO 99

 12   CONTINUE
        CALL CLSCR
        OKECH=.FALSE.
      GOTO 21

 99   RETURN

      ENTRY SRINM1(STORE)
      RETURN
      END
