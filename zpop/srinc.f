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
      SUBROUTINE SRINC
     >  (NLOG,KSC,KPL,IP,LM,NT,OX,WF,WO,Q,AM,FE,HZ,STORE,LUO,
     >                                                 YNRM,NOC,NRD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION OX(*), WF(*), WO(*)
      LOGICAL STORE

      COMMON/CDF/ IES,IORDRE,LCHA,LIST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN
      COMMON/CONST/ CL,PI,DPI,RAD,DEG,QE,AH
      LOGICAL OKECH, OKVAR, OKBIN
      COMMON/ECHL/OKECH, OKVAR, OKBIN
C      PARAMETER (MXSTEP=90000,MSAM=1000)
      INCLUDE 'MXSTEP.H'
      PARAMETER (MSAM=1000)
      COMMON/SRTAB/ EFD(MXSTEP,3),D3W(MSAM,9),NPTS
      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB

      DIMENSION XO(3)
      SAVE XO
      DIMENSION DW(MSAM,6), SDW(6), SYDX(4)

      CHARACTER XVA(5)*10, YVA(5)*15
      CHARACTER XDI(5)*10, YDI(5)*15, SY(5)*12

      DIMENSION ANU(MSAM), DNU(MSAM)
      COMMON/FREQ/ ANU, DNU

      SAVE MNU,NZ,DZ,DW

      CHARACTER * (9)   DMY, HMS

      DATA XVA / 'omga',    'WaveL.' , 'omga', 'Phi', 'Psi'/
      DATA YVA / 
     >  'dW/dNu', 'dW/dWaveL.', 'dN/dt.(do/o)', 'dW/dPhi', 'dW/dPsi'/
      DATA XDI / '(keV)',    '(Mu-m)' , '(keV)', 2*'(rad)'/
      DATA YDI / '(J/Hz)', '(J/m)', '(Phot/s.BW.srd)', 2*'(J/rad)'/
      DATA SY / 2*'(J in wndow)',  '(phot./wndw)', 2*'(J in wndow)' /
      DATA ONE / 1.D0 /

      WRITE(6,*) ' SBR SRINC    Busy, calculating power over window'
      
      MNU = WF(3)
      DO 1 IN=1,MNU
        D3W(IN,1) = ANU(IN)
        D3W(IN,2) = 0.D0
        D3W(IN,3) = 0.D0
 1    CONTINUE

      XO(1) = OX(1)
      NY = WO(5)
      NY2 = (NY+1)/2
      NZ = WO(6)
      NZ2 = (NZ+1)/2
      DY = 2.D0 * WO(2) / WO(5)
      DZ = 2.D0 * WO(3) / WO(6)

      DO 2 JY=1, NY
        D3W(JY,5) = 0.D0
        D3W(JY,6) = 0.D0
 2    CONTINUE
      DO 3 KZ=1, NZ
        D3W(KZ,8) = 0.D0
        D3W(KZ,9) = 0.D0
 3    CONTINUE

      IF(KSC .EQ. 4) THEN
        MX = NZ
      ELSE
        MX = MNU
      ENDIF

      IW = 0
      XO(3) = OX(3) - DZ*(NZ/2) - DZ
      DO 20 KZ=1, NZ
        XO(3) = XO(3) + DZ
        XO(2) = OX(2) - DY*(NY/2) - DY
        RPS = 1.D0
        IF(KZ.EQ.1 .OR. KZ.EQ.NZ) RPS=0.5D0
        DO 21 JY=1, NY
          IW = IW + 1
          XO(2) = XO(2) + DY

          PH = ATAN2( XO(2), XO(1) )
          DPH = ATAN2( (XO(2) + 0.5D0*DY), XO(1) )
     >        - ATAN2( (XO(2) - 0.5D0*DY), XO(1) )
          RXY = SQRT(XO(1)*XO(1) + XO(2)*XO(2))
          PS = ATAN2( XO(3), RXY )
          DPS = ATAN2( (XO(3) + 0.5D0*DZ), RXY )
     >        - ATAN2( (XO(3) - 0.5D0*DZ), RXY ) 

C--------- Solid angle dOmga (rad*rad)
          DPS = DPS * RPS
          IF(JY.EQ.1 .OR. JY.EQ.NY) DPH = DPH/2.D0
C          DOM = COS( PS ) * DPH * DPS
          DOM =  DPH * DPS

          IF(KZ .EQ. NZ2) D3W(JY,4) = PH
          IF(JY .EQ. NY2) D3W(KZ,7) = PS

          CALL DATE2(DMY)                     
          CALL TIME2(HMS)

          IF(IP .EQ. 0) THEN
            CALL SREF(3,XO, 1, 2,Q,AM,FE, 0, 0,
     >                                              GAM,R,NOC,NRD,*79)

          ELSEIF(IP .EQ. 1) THEN
C----------- Intermediate plot of E-field - at each fragment of the window

            CALL SREF(1,XO, 1, 2,Q,AM,FE, 0, 0,
     >                                              GAM,R,NOC,NRD,*179)
 179        CONTINUE
            IF(XMI.LT.XMA .AND. YMI.LT.YMA) THEN
C              CALL TXTFBG
              CALL TRAXES(XMI,XMA,YMI,YMA,-1) 
            ELSE
              WRITE(6,*) 
     >         ' Scale min-max problem ; cannot plot E-field !'
            ENDIF
            CALL SREF(2,XO, 1, 2,Q,AM,FE, 0, 0,
     >                                              GAM,R,NOC,NRD,*79)
          ENDIF

C--------- Calculate DW_1 = omga (KeV)
C                    DW_2,3 = dW_y,z/dNu.dOmga (J/Hz.srd)
 79       CALL SRDWC(ANU1,MNU,HZ,R,
     >                             DW)

C--------- Store sum over fragments  -> D3W_2,3(Nu) = dW_x,y/dNu  (J/Hz),
C          calculate energy in present frgmnt -> SDW_2,3 = W_x,y  (J)
C          and store energy at prsnt Phi -> D3W_5,6(Phi) = dW_x,y/dPhi (J/rad)
C                    energy at prsnt Psi -> D3W_8,9(Psi) = dW_x,y/dPsi (J/rad)
          CALL SRDWST(MNU,JY,KZ,DW,HZ,DOM,DPH,DPS,
     >                                            D3W,SDW)

         IF(JY .EQ. NY2) THEN
            WRITE(6,*) '---------------------------------------'
            WRITE(6,FMT='(A1,A9,2X,A9,$)') CHAR(13), DMY, HMS
            WRITE(6,FMT=
     >'('' Energy in fragmnt # '',I5,''  at y,z='',1P,2G12.4,'' m :'')')
     >      IW, XO(2),XO(3)
C----------- On HP stations only :
            WRITE( *,*) JY,KZ,' Py =',SDW(2),';  Pz =',SDW(3),' (J)'
            CALL FLUSH(6)

            WRITE(NLOG,*) '---------------------------------------'
            WRITE(NLOG,FMT='(A9,2X,A9)') DMY, HMS
            WRITE(NLOG,FMT=
     >'('' Energy in fragmnt # '',I5,''  at y,z='',1P,2G12.4,'' m :'')')
     >      IW, XO(2),XO(3)
C----------- On HP stations only :
            WRITE(NLOG,*) JY,KZ,' Py =',SDW(2),';  Pz =',SDW(3),' (J)'
            CALL FLUSH(NLOG)
          ENDIF

          IF(STORE) THEN
C----------- Store W_y,z in present fragment 
            YN = YNRM/DOM
            WRITE(LUO,FMT='(1P,5G14.6)')
     >        PH,PS,SDW(2)*YN, SDW(3)*YN, (SDW(2)+SDW(3))*YN
            IF(JY .EQ. NY) WRITE(LUO,*)
          ENDIF

          IF(IP .EQ. 1) THEN
C----------- Intermediate plot of dW_y,z/dNu.dOmga - at each fragment

            CALL SRDWSC(  0,KPL,1, 2, 3,MNU,DW,YNRM,
     >                                             OKECH,IER)
            
            IF(IER .GE. 0) THEN
              CALL CLSCR

              CALL SRAX(.TRUE.,  0,XMI,XMA,YMI,YMA,XDI(KSC),YDI(KSC))
              CALL SRDWPL(NLOG,  0,KPL,1, 2, 3,MNU,DW,YNRM,
     >                                                          SYDX)
              CALL SRPLV(NLOG,  0,LM,NT,MNU,SYDX,HZ,
     >          XVA(KSC),XDI(KSC),YVA(KSC),YDI(KSC),SY(KSC),YNRM)
            ELSE
              WRITE(6,*) ' Scale error detected in SBR SRDWSC'
              WRITE(6,*) ' Could not plot E-field'
            ENDIF
          ENDIF
          
 21     CONTINUE
 20   CONTINUE

      CLOSE(LUO)
      RETURN

      ENTRY SRINPL(NLOG,KSC,KPL,IP,LM,NT,OX,WF,WO,Q,AM,FE,HZ,
     >                                                 YNRM,NOC,NRD)

      IF(KSC .EQ. 4) THEN
C------- dW/dPhi(Phi)
        J1 = 4
        J2 = 5
        J3 = 6
        MX = NY
C        DX = DY
        FAC = 1.D0
      ELSEIF(KSC .EQ. 5) THEN
C------- dW/dPsi(Psi)
        J1 = 7
        J2 = 8
        J3 = 9
        MX = NZ
C        DX = DZ
        FAC = 1.D0
      ELSE
C------- dW/dNu(omga)
        J1 = 1
        J2 = 2
        J3 = 3
        MX = MNU
C        DX = DNU
        FAC = HZ
      ENDIF

      IF(.NOT. OKECH) THEN
        WRITE(6,*) '               Busy, calculate scales'
        CALL SRDWSC(KSC,KPL,J1,J2,J3,MX,D3W,YNRM,
     >                                        OKECH,IER)

        IF(IER .GE. 0) THEN
          WRITE(6,*) '               Busy, plotting'
          CALL SRAX(.TRUE.,KSC,XMI,XMA,YMI,YMA,XDI(KSC),YDI(KSC))
        ELSE
          WRITE(6,*) ' Scale error detected in SBR SRDWSC'
          WRITE(6,*) ' Could not plot '
        ENDIF
      ENDIF

      CALL SRDWPL(NLOG,KSC,KPL,J1,J2,J3,MX,D3W,YNRM,
     >                                                    SYDX)
      CALL SRPLV(NLOG,KSC,LM,NT,MX,SYDX,FAC,
     >            XVA(KSC),XDI(KSC),YVA(KSC),YDI(KSC),SY(KSC),YNRM)

      RETURN
      END
