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
      SUBROUTINE CAVITE(*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      COMMON/PTICUL/ AM,Q,G,TO
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      INCLUDE 'MXFS.H'
      COMMON/SCAL/SCL(MXF,MXS),TIM(MXF,MXS),NTIM(MXF),KSCL
      COMMON/SYNCH/ PH(MXT), DPR(MXT), PS
      COMMON/UNITS/ UNIT(MXJ)
 
      DIMENSION WF1(MXT), PHAS(MXT)
      SAVE WF1, PHAS

      CHARACTER*9 SKAV(8)
      DIMENSION DTI0(MXT)
      SAVE DTI0

      LOGICAL OKOPEN, OKIMP, IDLUNI
      SAVE LUN, OKOPEN, OKIMP

      SAVE TIOLD, PHIOLD

      SAVE PP0

      DATA WF1, PHAS / MXT*0.D0, MXT*0.D0 /
      DATA SKAV/'** OFF **','OPTION 1 ','OPTION 2 ','OPTION 3 ', 
     >   'OPTION 4 ', 'OPTION 5 ', '   FFAG  ', ' Isochro.' /

      data dti0 / mxt*0.D0 /

      DATA LUN / 99 /
      DATA OKOPEN, OKIMP /.FALSE., .FALSE. /

      DATA TIOLD, PHIOLD /  2 * 0.D0 /
      DATA PP0 / 1.D0 / 

      KCAV = NINT(A(NOEL,1))
      OKIMP = (NINT(10.D0*A(NOEL,1)) - 10*KCAV) .EQ. 1
      IF(OKIMP) THEN 
        IF(.NOT.OKOPEN) THEN
          IF(IDLUNI(
     >              LUN)) THEN
            OPEN(UNIT=LUN,FILE='zgoubi.CAVITE.Out',
     >                     FORM='FORMATTED',ERR=99, IOSTAT=IOS)
          ELSE
            GOTO 99
          ENDIF
          IF(IOS.NE.0) GOTO 99
          OKOPEN = .TRUE.
        ENDIF
      ENDIF

      IF(NRES.GT.0) WRITE(NRES,100) SKAV(KCAV+1)
 100  FORMAT(15X,' Accelerating  cavity',/,16X,A,/)
      IF(KCAV .EQ. 0) RETURN
 
      AN10 = A(NOEL,10)
      AN11 = A(NOEL,11)
      VLT= A(NOEL,20)
      AN21= A(NOEL,21)
      PHS= AN21
 
      IF(NRES.GT.0) THEN
        IF(Q*AM .EQ. 0.D0) THEN
          WRITE(NRES,106)
 106      FORMAT(//,15X,' SVP  INDIQUER  M  ET  Q  DES  PROJECTILES !'
     >         ,/,15X,' - UTILISER  LE  MOT-CLE  ''PARTICUL''',/)
          RETURN 1
        ENDIF
      ENDIF

C----- P0, AM  are  in  MEV/c, /c^2
      P0 = BORO*CL9*Q
      AM2= AM*AM
      QV = VLT*Q *1.D-6
 

      GOTO(10,20,30,40,50,60,70,80) KCAV
      IF(NRES .GT. 0) 
     >         WRITE(NRES,*) ' SBR  CAVITE:  ERRONEOUS  INPUT  DATA'
      RETURN 1
  

 10   CONTINUE
C Ph_s is computed from rigidity law as specified using SCALING/CAVITE
      ORBL = AN10
      HARM = AN11
C----- PARTICULE SYNCHRONE, ENTREE DE LA CAVITE
      PS = P0*SCALER(IPASS,NOEL,
     >                          DTA1)
      BTS = PS/SQRT(PS*PS+AM2)
      DTS = ORBL / ( CL * BTS)
      OMRF = 2.D0*PI*HARM / DTS
      WS = PS/BTS
      FREV = HARM/DTS
C----- PARTICULE SYNCHRONE, SORTIE DE LA CAVITE
      PS = P0*SCALER(IPASS+1,NOEL,
     >                            DTA1)
      PP0 = PS / P0
      DWS = SQRT(PS*PS+AM2) - WS
      WS = WS + DWS
      PHS=ASIN(DWS/QV)
      GOTO 1
 
 20   CONTINUE
      ORBL = AN10
      HARM = AN11
C----- PARTICULE SYNCHRONE, ENTREE DE LA CAVITE
      IF(IPASS .EQ. 1) PS = P0
      BTS = PS/SQRT(PS*PS+AM2)
      DTS = ORBL / ( CL * BTS)
      OMRF = 2.D0*PI*HARM / DTS
      WKS = PS/BTS - AM
      FREV = HARM/DTS
C----- PARTICULE SYNCHRONE, SORTIE DE LA CAVITE
      DWS = QV*SIN(PHS)
      WKS = WKS + DWS
      PS = SQRT(WKS*(WKS+2.D0*AM))
      PP0 = PS / P0
      WS = WKS + AM
!      WS = SQRT(PS*PS+AM2)
      GOTO 1
 
 30   CONTINUE
      HARM = AN11
      CALL SCUMR(
     >           DUM,SCUM,TCUM) 
      OMRF = 2.D0*PI*HARM 
      DWS = QV*SIN(PHS)
 
      IF(NRES.GT.0) THEN
        WRITE(NRES,130) 
     >  OMRF/(2.D0*PI),PHS,HARM,DWS,SCUM*UNIT(5),TCUM+AN10, AM, Q*QE
 130    FORMAT(
     >  /,20X,'Cavity  frequency                  =',1P,E15.4,' Hz',
     >  /,20X,'Synchronous  phase                 =',  E15.4,' rad',
     >  /,20X,'Harmonic                           =',  E15.4,' ',
     >  /,20X,'Synchronous energy  gain           =',  E15.4,' MeV',
     >  /,20X,'Cumulated  distance  from  origin  =',  E15.4,' m',
     >  /,20X,'Synchronous  time                  =',  E15.4,' s',
     >  /,20X,'Particle mass           =',  E15.5,' MeV/c2',
     >  /,20X,'         charge         =',  E15.6,' C')
      ENDIF
 
      DO 31 I=1,IMAX
        IF(IEX(I) .LT. -1) GOTO 31
        TTA = F(3,I)*.001D0
        PHI = F(5,I)*.001D0
        P = P0*F(1,I)
        PCP = P*COS(PHI)
        PX = PCP * COS(TTA)
        PY = PCP * SIN(TTA)
        PZ = P * SIN(PHI)
 
        AM2 = AMQ(1,I)*AMQ(1,I)
        ENRG = SQRT(P*P+AM2)
        WF1(I) = ENRG-AMQ(1,I)
        BTA = P/ENRG
        TI = F(7,I)*UNIT(7)
        DPHI = (TI-(TCUM+AN10))*OMRF
        DPHI = DPHI - INT(DPHI/(2.D0*PI))*2.D0*PI-PI
c        IF(DPHI .GT.  PI) DPHI =DPHI -2.D0*PI
c        IF(DPHI .LT. -PI) DPHI =DPHI +2.D0*PI
        PH(I)=DPHI-PHS
        WF = WF1(I) + QV*SIN(PH(I))
        WF1(I) = WF
        P = SQRT(WF*(WF + 2.D0*AMQ(1,I)))
        PX=SQRT( P*P -PY*PY-PZ*PZ)
 
        DPR(I)=WF-QV*SIN(PHS)
        F(1,I) = P/P0
        F(3,I) = ATAN2(PY,PX)*1000.D0
        F(5,I) = ATAN2(PZ,SQRT(PX*PX+PY*PY))*1000.D0
 31   CONTINUE
      GOTO 88
 
 40   CONTINUE
      ORBL = AN10
      HARM = AN11
C------ For FNAL p-driver, Nov. 2000. Also ok with Saturne (~/saturne/sat_cav4 cases)
C     ... PARTICULE SYNCHRONE, ENTREE DE LA CAVITE
      PS = SCALPS()
      WS = SQRT(PS*PS+AM2)
      BTS = PS/WS
      DTS = ORBL / ( CL * BTS)
      OMRF = 2.D0*PI*HARM / DTS
      FREV = HARM /DTS
C     ... PARTICULE SYNCHRONE, SORTIE DE LA CAVITE
C-------- Watch out ! This DTS calcul. assumes CAVITE IS THE LAST optical lmnt in zgoubi.dat
      PS = SCALDP(DTS,
     >                TIME)    
      DWS = SQRT(PS*PS+AM2) - WS
      WS = WS + DWS
      PHS=ASIN(DWS/QV)
      GOTO 1
 
 50   CONTINUE
      ORBL = AN10
      HARM = AN11
C------ For FNAL p-driver, Nov. 2000. Vrf is read from a file
C     ... PARTICULE SYNCHRONE, ENTREE DE LA CAVITE
      PS = SCALPS()
      WS = SQRT(PS*PS+AM2)
      BTS = PS/WS
      DTS = ORBL / ( CL * BTS)
      OMRF = 2.D0*PI*HARM / DTS
      FREV = HARM/DTS
C     ... PARTICULE SYNCHRONE, SORTIE DE LA CAVITE
C-------- Watch out ! This DTS calcul. assumes CAVITE is at the end of THE LAST optical lmnt
      PS = SCALDP(DTS,
     >                TIME)    
      DWS = SQRT(PS*PS+AM2) - WS
      WS = WS + DWS
      QV = VALQV(TIME)*Q *1.D-6
      PHS=ASIN(DWS/QV)
      GOTO 1
 
 60   CONTINUE
C----- FFAG acceleration. 
C      Single cavity, assumed located at end of zgoubi.dat     
C WS0 is the synchronous energy at start
      WS0 = AN11
C      HN = AN12
      HN = 1.D0
      FREV0 = AN10/HN
      PS0 = SQRT((WS0+AM)**2 - AM2)
      BTS0 = PS0/SQRT(PS0*PS0+AM2)
C     ... Conditions at cavity entrance 
C        write(*,*) frev0,hn,SCALER(IPASS,NOEL,DTA1,DTA2,DTA3),ipass,noel
      OMRF = 2.D0*PI*FREV0*HN
C     ... Synchronous conditions at cavity exit
      DWS = QV * SIN(PHS)
      WS = WS0 + DBLE(IPASS) * DWS
      PS = SQRT((WS+AM)**2 - AM2)
      BTS = PS/SQRT(PS*PS+AM2)

c        write(*,*) ' sbr cavite  goto 1 '
c         goto 1

      IF(NRES.GT.0) THEN
        GTRNUS = SQRT(ABS(QV*AN11*COS(PHS) / (2.D0*PI*WS)))
        ACCDP=  SQRT(QV/(AN11*WS)) * 
     >  SQRT(ABS(-(2.D0*COS(PHS)/PI+(2.D0*PHS/PI-1.D0)*SIN(PHS))))
        WRITE(NRES,126) HN, QV/(Q *1.D-6), PHS, OMRF/(2.D0*PI),  
     >       QV*SIN(PHS),COS(PHS),   GTRNUS, ACCDP
 126    FORMAT(1P, 
     >       /,20X,'RF  harmonic          =',G15.4,
     >       /,20X,'Peak  voltage         =',G15.4,' V',
     >       /,20X,'Synchronous  phase    =',G15.4,' rd',
     >       /,20X,'RF  frequency         =',G15.4,' Hz',
     >       /,20X,'qV.SIN(Phi_s)         =',G15.4,' MeV',
     >       /,20X,'cos(Phi_s)            =',G15.4,' ',
     >       /,20X,'Nu_s/sqrt(alpha)      =',G15.4,'  ',
     >       /,20X,'dp-acc*sqrt(alpha)    =',G15.4,'  ')
      ENDIF

      DO 63 I=1,IMAX 

        IF(IEX(I) .LT. -1) GOTO 63

        TTA = F(3,I)*.001D0
        PHI = F(5,I)*.001D0
        P = P0*F(1,I)
        PCP = P*COS(PHI)
        PX = PCP * COS(TTA)
        PY = PCP * SIN(TTA)
        PZ = P * SIN(PHI)
 
        AM2 = AMQ(1,I)*AMQ(1,I)
        ENRG = SQRT(P*P+AM2)
        WF1(I) = ENRG-AMQ(1,I)
        BTA = P/ENRG
C F(7,I) is time in mu_s. of course, TI is in s
        TI = F(7,I) * UNIT(7) 
C        PHI = OMRF * TI + PHS
        DUM = SCALE4(F(7,I),WF1(I))
        PHI = SCALER(IPASS,NOEL,
     >                          coTime) + PHS

c        ekbef = sqrt(p*p + am*am) -am
c        zero = 0.
c          write(88,*) ipass, zero, phi, ti, ekbef

C Phase, in [-pi,pi] interval
        PHI = PHI - INT(PHI/(2.D0*PI)) * 2.D0*PI 
        IF    (PHI .GT.  PI) THEN
          PHI =PHI - 2.D0*PI
        ELSEIF(PHI .LT. -PI) THEN
          PHI =PHI + 2.D0*PI
        ENDIF

        DWF =  QV * SIN(PHI)
        PH(I) = PHI
        WF1(I) = WF1(I) + DWF
        WF = WF1(I)

        P = SQRT(WF*(WF + 2.D0*AMQ(1,I)))
        DPR(I)=(P-PS)/PS
        PX=SQRT( P*P -PY*PY-PZ*PZ)

c        ekaft = sqrt(p*p + am*am) -am
c        qvhat = (ekaft - ekbef)/SIN(PHI)
c        write(89,fmt='(1p, 5e24.16)') 
c     >        ekbef, TI-TIOLD, coTime, PHI-PHIOLD, 
c     >                 2.d0*PI*(TI-TIOLD)/(PHI-PHIOLD)
        TIOLD = TI
        PHIOLD = PHI


        F(1,I) = P/P0
        F(3,I) = ATAN(PY/PX)*1000.D0
        F(5,I) = ATAN(PZ/SQRT(PX*PX+PY*PY))*1000.D0

        IF(OKIMP) THEN
          scala = SCALER(IPASS,NOEL,
     >                              DTA1)
          WRITE(LUN,FMT='(1P,I6,1x,5G14.6,1x,I6,A)') 
     >            ipass, omrf/(2.*pi),
C     >            scala, 
C     >            PHI-PHS,
     >            PHI,
     >            TI,
     >             wf,
C     >              QV*SIN(PH(I))/(Q*1.D-6), 
     >               P-PS, I, 
     >            ' ipass freq phi-phs ti qV*sin p-ps ITraj'
        ENDIF
 63   CONTINUE
      GOTO 88
 
 
 70   CONTINUE
C Bucketless acceleration. Ok for non-scaling FFAG
C Used for muon, EMMA
      CALL SCUMR(
     >           DUM,SCUM,TCUM) 
C Orbit length between 2 cavities, RF freq., phase of 1st cavity (ph0=0 is 
C at V(t)=0)
      ORBL = AN10
      FCAV = AN11
      PH0 = AN21
      PS = P0
      BTS = PS/SQRT(PS*PS+AM2)
      DTS = ORBL / ( CL * BTS)
      HARM = DTS * FCAV
      OMRF = 2.D0 * PI * FCAV
C      WS = PS / BTS
C      TS = TS + DTS
 
      IF(NRES.GT.0) THEN
        WRITE(NRES,170) FCAV,HARM,QV,BORO,DTS,
     >                    SCUM*UNIT(5),TCUM*UNIT(7),AM,Q*QE
 170    FORMAT(
     >  /,20X,'Cavity  frequency                 =',1P,G15.6,' Hz',
     >  /,20X,'Harmonic                          =',   G15.6,' ',
     >  /,20X,'Max energy  gain                  =',   G15.6,' MeV',
     >  /,20X,'TOF for BRho_ref (',G10.2,') is',       G15.6,' s',
     >  /,20X,'Cumulated distance since origin   =',   G15.6,' m',
     >  /,20X,'Cumulated   TOF      "     "      =',   G15.6,' s',
     >  /,20X,'Particle mass                     =',   G15.6,' MeV/c2',
     >  /,20X,'         charge                   =',   G15.6,' C')
      ENDIF
 
      DO 71 I=1,IMAX
        IF(IEX(I) .LT. -1) GOTO 71
        TTA = F(3,I)*.001D0
        PHI = F(5,I)*.001D0
        P = P0*F(1,I)
        PCP = P*COS(PHI)
        PX = PCP * COS(TTA)
        PY = PCP * SIN(TTA)
        PZ = P * SIN(PHI)
 
        AM2 = AMQ(1,I) * AMQ(1,I)
        ENRG = SQRT(P*P + AM2)
        WF1(I) = ENRG - AMQ(1,I)
        BTA = P / ENRG
C F(7,I) is time in mu_s. of course, TI is in s
        TI = F(7,I) * UNIT(7) 
        PHI = OMRF * TI + PH0
          
C Phase, in [-pi,pi] interval
        PHI = PHI - INT(PHI/(2.D0*PI)) * 2.D0*PI 
        IF    (PHI .GT.  PI) THEN
          PH(I) =PHI - 2.D0*PI
        ELSEIF(PHI .LT. -PI) THEN
          PH(I) =PHI + 2.D0*PI
        ELSE
          PH(I) =PHI 
        ENDIF

        DWF =  QV * SIN(PH(I))
C------- Rustine etude ffag muon
        IF(OMRF.LE.0.D0) DWF = QV
C------------------------------
        WF1(I) = WF1(I) + DWF
        WF = WF1(I)

C Kin. energy, MeV
        DPR(I)=WF

        P = SQRT(WF*(WF + 2.D0*AMQ(1,I)))
        PX=SQRT( P*P -PY*PY-PZ*PZ)
        F(1,I) = P / P0
        F(3,I) = ATAN2(PY,PX) / UNIT(2)
        F(5,I) = ATAN2(PZ,SQRT(PX*PX+PY*PY)) / UNIT(4)

        IF(OKIMP) 
     >          WRITE(LUN,FMT='(1P,4G14.6,2I6)') PH(I),DPR(I),
     >            TI,QV*SIN(PH(I))/(Q*1.D-6), I , IPASS

 71   CONTINUE

      GOTO 88
 

 80   CONTINUE
C For HNJ acceleration developements 
      CALL SCUMR(
     >           DUM,SCUM,TCUM) 
C Orbit length between 2 cavities, RF freq., phase of 1st cavity (ph0=0 is 
C at V(t)=0)
      ORBL = AN10
      FCAV = AN11
      PH0 = AN21
      PS = P0
      BTS = PS/SQRT(PS*PS+AM2)
      DTS = ORBL / ( CL * BTS)
      HARM = DTS * FCAV
      OMRF = 2.D0 * PI * FCAV
C      WS = PS / BTS
C      TS = TS + DTS
 
      IF(NRES.GT.0) THEN
        WRITE(NRES,180) FCAV,HARM,QV,BORO,DTS,
     >                    SCUM*UNIT(5),TCUM*UNIT(7),AM,Q*QE
 180    FORMAT(
     >  /,20X,'Cavity  frequency                 =',1P,G15.6,' Hz',
     >  /,20X,'Harmonic                          =',   G15.6,' ',
     >  /,20X,'Max energy  gain                  =',   G15.6,' MeV',
     >  /,20X,'TOF for BRho_ref (',G10.2,') is',       G15.6,' s',
     >  /,20X,'Cumulated distance since origin   =',   G15.6,' m',
     >  /,20X,'Cumulated   TOF      "     "      =',   G15.6,' s',
     >  /,20X,'Particle mass                     =',   G15.6,' MeV/c2',
     >  /,20X,'         charge                   =',   G15.6,' C')
      ENDIF
 
      DO 81 I=1,IMAX
        IF(IEX(I) .LT. -1) GOTO 81
        TTA = F(3,I)*.001D0
        PHI = F(5,I)*.001D0
        P = P0*F(1,I)
        PCP = P*COS(PHI)
        PX = PCP * COS(TTA)
        PY = PCP * SIN(TTA)
        PZ = P * SIN(PHI)
 
        AM2 = AMQ(1,I) * AMQ(1,I)
        ENRG = SQRT(P*P + AM2)
        WF1(I) = ENRG - AMQ(1,I)
        BTA = P / ENRG
C F(7,I) is time in mu_s. of course, TI is in s
        TI = F(7,I) * UNIT(7) 
C        PHI = OMRF * (TI - TS)
        alpha = an10
        ddphi = omrf * ti * alpha * (p-p0)/p0
        PHI = OMRF * TI + PH0  - ddphi
C         write(*,*) ' sbr cav phi, dphi ',phi,  ddphi

C Phase, in [-pi,pi] interval
        PHI = PHI - INT(PHI/(2.D0*PI)) * 2.D0*PI 
        IF    (PHI .GT.  PI) THEN
          PH(I) =PHI - 2.D0*PI
        ELSEIF(PHI .LT. -PI) THEN
          PH(I) =PHI + 2.D0*PI
        ELSE
          PH(I) =PHI 
        ENDIF

        DWF =  QV * SIN(PH(I))
C------- Rustine etude ffag muon
        IF(OMRF.LE.0.D0) DWF = QV
C------------------------------
        WF1(I) = WF1(I) + DWF
        WF = WF1(I)

C Kin. energy, MeV
        DPR(I)=WF

        P = SQRT(WF*(WF + 2.D0*AMQ(1,I)))
        F(1,I) = P / P0
        PX=SQRT( P*P -PY*PY-PZ*PZ)
        F(3,I) = ATAN2(PY,PX) / UNIT(2)
        F(5,I) = ATAN2(PZ,SQRT(PX*PX+PY*PY)) / UNIT(4)

        IF(OKIMP) 
     >          WRITE(LUN,FMT='(1P,4G14.6,2I6)') PH(I),DPR(I),
     >            TI,QV*SIN(PH(I))/(Q*1.D-6), I , IPASS

 81   CONTINUE

      GOTO 88
 
 1    CONTINUE

      IF(NRES.GT.0) THEN
        GTRNUS = SQRT(ABS(QV*AN11*COS(PHS) / (2.D0*PI*WS)))
        ACCDP=  SQRT(QV/(AN11*WS)) * 
     >  SQRT(ABS(-(2.D0*COS(PHS)/PI+(2.D0*PHS/PI-1.D0)*SIN(PHS))))
        WRITE(NRES,120) ORBL, AN11, QV/(Q *1.D-6), FREV, PHS, DTS, 
     >       QV*SIN(PHS),COS(PHS),   GTRNUS, ACCDP
 120    FORMAT(1P, 
     >       /,20X,'Orbit  length         =',G15.4,' m',
     >       /,20X,'RF  harmonic          =',G15.4,
     >       /,20X,'Peak  voltage         =',G15.4,' V',
     >       /,20X,'RF  frequency         =',G15.4,' Hz',
     >       /,20X,'Synchronous  phase    =',G15.4,' rd',
     >       /,20X,'Isochronous  time     =',  G15.4,' s',
     >       /,20X,'qV.SIN(Phi_s)         =',G15.4,' MeV',
     >       /,20X,'cos(Phi_s)            =',G15.4,' ',
     >       /,20X,'Nu_s/sqrt(alpha)      =',G15.4,'  ',
     >       /,20X,'dp-acc*sqrt(alpha)    =',G15.4,'  ')

C        IF(KCAV .EQ. 1) WRITE(NRES,199) SCALER(IPASS+1,NOEL,DTA1,DTA2,DTA3)
C 199    FORMAT(/,20X,'Post acceleration SCALING factor is ',1P,G16.8)
C            write(33,*) time, QV/(Q/QE *1.D-6)/1.6D6, QV*SIN(PHS)/1.3, 
C     >        ps/16, AN11/DTS/(53e6), GTRNUS/40/0.11
      ENDIF

C----- Initial conditions of  the IMAX particles
      IF(IPASS .EQ. 1) THEN

        DO 2 I=1,IMAX
          F(6,I)=F(6,I)-FO(6,I)
          P = P0*F(1,I)
          AM2 = AMQ(1,I)*AMQ(1,I)
          WF1(I) = SQRT(P*P+AM2)-AMQ(1,I)
          PHAS(I) = PHS
C Introduction dti0 juin 05 pour tenir compte de la phase initiale
C A verifier...
C          bta = p / sqrt(p*p + am2*am2)
C          dti0(i) = -fo(6,I)*.01D0/(bta * cl)
 2      CONTINUE

      ENDIF
 
      DO 3 I=1,IMAX 

        IF(IEX(I) .LT. -1) GOTO 3

        TTA = F(3,I)*.001D0
        PHI = F(5,I)*.001D0
        P = P0*F(1,I)
        PCP = P*COS(PHI)
        PX = PCP * COS(TTA)
        PY = PCP * SIN(TTA)
        PZ = P * SIN(PHI)
 
        AM2 = AMQ(1,I)*AMQ(1,I)
        ENRG = SQRT(P*P+AM2)
        WF1(I) = ENRG-AMQ(1,I)
        BTA = P/ENRG
C At all pass#,  F6i is the previous-turn path length, DTI is the time it 
C   took since the last passage in CAVITE 
        DTI = F(6,I)*.01D0 / (BTA*CL)
C rajouté en juin 05, pour tenir compte de s à l'origine:
           dti = dti + dti0(i)

        PHAS(I) = PHAS(I) + (DTI-DTS)*OMRF
CCC        PHAS(I) = PHAS(I) + PHS
        IF(PHAS(I) .GT.  PI) PHAS(I) =PHAS(I) -2.D0*PI
        IF(PHAS(I) .LT. -PI) PHAS(I) =PHAS(I) +2.D0*PI
        WF = WF1(I) + QV*SIN(PHAS(I))

C          write(*,fmt='(/,1p,3e12.4,2i5)') 
C     >             dti,F(7,I),F(6,I),i,ipass
C          write(*,fmt='(/,1p,4e12.4,2i5)') 
C     >             dti*ipass,, F(7,I),F(6,I),i,ipass

CCC        WF = WF1(I) + QV*(SIN(PHAS(I)) - SIN(PHS))
        WF1(I) = WF
        P = SQRT(WF*(WF + 2.D0*AMQ(1,I)))
        PX=SQRT( P*P -PY*PY-PZ*PZ)
 
C        DPR(I)=P-PS
        DPR(I)=(P-PS)/PS
        PH(I)=PHAS(I)

        F(1,I) = P/P0
        F(3,I) = ATAN(PY/PX)*1000.D0
        F(5,I) = ATAN(PZ/SQRT(PX*PX+PY*PY))*1000.D0
 
        F(6,I)=0.D0 

        IF(OKIMP) 
     >    WRITE(LUN,FMT='(1P,6(E14.6,1X),2(I6,1X))') PH(I),PHS,P-PS,
     >    OMRF,TI,QV*SIN(PH(I))/(Q*1.D-6), I , IPASS

 3    CONTINUE
      GOTO 88
 
 88   CONTINUE
      RETURN

 99   IF(NRES .GT. 0) WRITE(NRES,101) 'CAVITE',
     >                          ' OPEN zgoubi.CAVITE.Out'
 101  FORMAT(/,' ******  SBR ',A,' : ERROR',A,/)
      RETURN

      ENTRY CAVIT1(
     >             SCALO)
      SCALO = PP0
      END
