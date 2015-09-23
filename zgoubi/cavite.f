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
C  USA
C  -------
      SUBROUTINE CAVITE(*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      PARAMETER (ISZTA=80)
      CHARACTER(ISZTA) TA
      PARAMETER (MXTA=45)
      INCLUDE "C.DONT.H"     ! COMMON/DONT/ TA(MXL,MXTA)
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
      INCLUDE "C.OBJET.H"     ! COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      INCLUDE "C.PTICUL.H"     ! COMMON/PTICUL/ AM,Q,G,TO
      INCLUDE "C.REBELO.H"   ! COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      INCLUDE 'MXFS.H'
      INCLUDE 'MXSCL.H'
      INCLUDE "C.SCAL.H"     ! COMMON/SCAL/ SCL(MXF,MXS,MXSCL),TIM(MXF,MXS),NTIM(MXF),KSCL
      INCLUDE "C.SCALP.H"     ! COMMON/SCALP/ VPA(MXF,MXP),JPA(MXF,MXP)
      INCLUDE "C.SYNCH.H"     ! COMMON/SYNCH/ PH(MXT), DPR(MXT), PS
      INCLUDE "C.UNITS.H"     ! COMMON/UNITS/ UNIT(MXJ)
 
      DIMENSION WF1(MXT), PHAS(MXT)
      SAVE WF1, PHAS

      SAVE SYNCH_TIME
      DIMENSION KFM(MXSCL)
      
      CHARACTER(10) SKAV(12)
      DIMENSION DTI0(MXT)
      SAVE DTI0

      LOGICAL OKOPEN, OKIMP, IDLUNI
      SAVE LUN, OKOPEN, OKIMP

      SAVE TIOLD, PHIOLD

      SAVE DWS, PHS_PREV
!     SR loss
      LOGICAL SRLOSS
      PARAMETER (SQRT2 = SQRT(2.D0),SQRT8 = SQRT(8.D0), I0=0.D0)

      CHARACTER(60) TYPCH(5)
      INTEGER DEBSTR, FINSTR

      DATA WF1, PHAS / MXT*0.D0, MXT*0.D0 /
      DATA SKAV /'** OFF **','OPTION 1 ','OPTION 2 ','OPTION 3 ', 
     >   'OPTION 4 ', 'OPTION 5 ', '   FFAG  ', 'Isochron.',
     >  ' ' , ' ' , 'eRHIClinac', 'SR + accel.' /

      DATA DTI0 / MXT*0.D0 /

      DATA LUN / 99 /
      DATA OKOPEN, OKIMP /.FALSE., .FALSE. /

      DATA TIOLD, PHIOLD /  2 * 0.D0 /
      DATA DWS / 0.D0 / 

      DATA TYPCH / 'no motion damping (det(M) forced to 1)',
     >'DE/E<<1 aproximation, no motion damping (det(M) forced to 1)',
     >'cavity is inhibited, transport is identity matrix',
     >'DE/E<<1 aproximation, regular transport including damping',
     >'regular transport, including damped motion' /

      DUM = SCALER(1, NOEL, 
     >                     DUM)

      SRLOSS = .FALSE.
      dwsr = 0.D0

      CALL SCALE9(
     >            KFM )
      DO IFM = 1, MXSCL
        IF(KFM(IFM) .LE. 0) THEN
          GOTO 121
        ELSE
          IF(KFM(IFM).GT.MXD .OR. KFM(IFM).GT.MXF) 
     >    CALL ENDJOB('Pgm cavite. Exceed array size, KFM = ',KFM(IFM))
        ENDIF
        DO I= 1 , JPA(KFM(IFM),MXP)
          A(NOEL,JPA(KFM(ifm),I)) = VPA(KFM(IFM) ,I)
        ENDDO
      ENDDO
 121  CONTINUE

      KCAV = NINT(A(NOEL,1))
      IF(IPASS .EQ. 1) THEN 
        OKIMP = TA(NOEL,1) .EQ. 'PRINT'
        IF(OKIMP) THEN 
          IF(.NOT.OKOPEN) THEN
            IF(IDLUNI(
     >              LUN)) THEN
              OPEN(UNIT=LUN,FILE='zgoubi.CAVITE.Out',
     >                     FORM='FORMATTED',ERR=99, IOSTAT=IOS)
            ELSE
              OKIMP = .FALSE.
              GOTO 99
            ENDIF
            IF(IOS.NE.0) GOTO 99
            OKOPEN = .TRUE.
          ENDIF
        ENDIF
      ENDIF

      IF(NRES.GT.0) THEN
        WRITE(NRES,100) SKAV(KCAV+1)
 100    FORMAT(15X,' Accelerating cavity. Type is :',3X,A,/)
        IF(OKIMP) WRITE(NRES,FMT='(15X,
     >  '' Cavite parameters saved in zgoubi.CAVITE.Out'',/)') 
      ENDIF

      IF(KCAV .EQ. 0) RETURN
 
      AN10 = A(NOEL,10)
      AN11 = A(NOEL,11)
      VLT= A(NOEL,20)
      AN21= A(NOEL,21)
      AN22= A(NOEL,22)
      PHS= AN21
 
      IF(Q*AM .EQ. 0.D0) CALL ENDJOB('Pgm cavite. Give mass & charge'
     >//'of particles. Use ''PARTICUL'' keyword',-99)

C----- P0, AM  are  in  MEV/c, /c^2
      P0 = BORO*CL9*Q
      AM2= AM*AM
      QV = VLT*Q *1.D-6
 

      GOTO(10,20,30,40,50,60,70,80,77,110,21) KCAV
      CALL ENDJOB(' Sbr cavite : No such option KCAV =',KCAV)
  
C------------------------------------------- 
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
      DWS = SQRT(PS*PS+AM2) - WS
      WS = WS + DWS
      PHS=ASIN(DWS/QV)
      GOTO 1
 
C------------------------------------------- 
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
C--- Case SR loss in storage ring (no acceleration)
      CALL SRLOS3(
     >            SRLOSS)
      IF(SRLOSS) WKS = WKS - DWS
      PS = SQRT(WKS*(WKS+2.D0*AM))
      WS = WKS + AM
      GOTO 1

C------------------------------------------- 
C Works like 20 if no SR, or if SR in storage ring. Diff is allows SR with acceleration. 
C It uses ph_s to accelerate (not for compensation of SR unlike 20) (so, ph_s=0 
C in storage mode); ph_s is corrected for SR compensation : 
C requires A(noel,22) = theoretical SR loss at first pass, then, SR loss assumes ~gamma^4 
C dependence at subsequent passes. 
 21   CONTINUE
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
      WS = WKS + AM
C--- Case SR loss in storage ring (no acceleration)
      CALL SRLOS3(
     >            SRLOSS)
      IF(SRLOSS) THEN
        gg4 = ( ws / (ws-dws)   )**4
        dwsr = a(noel,22) * gg4         
      ENDIF
      GOTO 1

C------------------------------------------- 
 30   CONTINUE
      HARM = AN11
      CALL SCUMR(
     >           DUM,SCUM,TCUM) 
      OMRF = 2.D0*PI*HARM 
      DWS = QV*SIN(PHS)
C----- PARTICULE SYNCHRONE, ENTREE DE LA CAVITE
      IF(IPASS .EQ. 1) PS = P0
      BTS = PS/SQRT(PS*PS+AM2)
      WKS = PS/BTS - AM
C----- PARTICULE SYNCHRONE, SORTIE DE LA CAVITE
      WKS = WKS + DWS
C--- Case SR loss
      CALL SRLOS3(
     >            SRLOSS)
      IF(SRLOSS) WKS = WKS - DWS
      PS = SQRT(WKS*(WKS+2.D0*AM))
 
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
C        PH(I)=DPHI-PHS
        PH(I) = PHS
        WF = WF1(I) + QV*SIN(PH(I))
        WF1(I) = WF
        P = SQRT(WF*(WF + 2.D0*AMQ(1,I)))
        PX=SQRT( P*P -PY*PY-PZ*PZ)


C      write(88,*) ' cavite ', QV*SIN(PH(I)), dphi,PH(I)*deg
 
        DPR(I)=WF-QV*SIN(PHS)
        F(1,I) = P/P0
        F(3,I) = ATAN2(PY,PX)*1000.D0
        F(5,I) = ATAN2(PZ,SQRT(PX*PX+PY*PY))*1000.D0
 31   CONTINUE
      GOTO 88
 
C------------------------------------------- 
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
 
C------------------------------------------- 
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
 
C------------------------------------------- 
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
      OMRF = 2.D0*PI*FREV0*HN
C     ... Synchronous conditions at cavity exit
      DWS = QV * SIN(PHS)
      WS = WS0 + DBLE(IPASS) * DWS
      PS = SQRT((WS+AM)**2 - AM2)
      BTS = PS/SQRT(PS*PS+AM2)

      IF(NRES.GT.0) THEN
        GTRNUS = SQRT(ABS(QV*AN11*COS(PHS) / (2.D0*PI*WS)))
        ACCDP=  SQRT(QV/(AN11*WS)) * 
     >  SQRT(ABS(-(2.D0*COS(PHS)/PI+(2.D0*PHS/PI-1.D0)*SIN(PHS))))
        WRITE(NRES,126) HN, QV/(Q *1.D-6), PHS, OMRF/(2.D0*PI),  
     >       QV*SIN(PHS),COS(PHS),   GTRNUS, ACCDP
 126    FORMAT(1P, 
     >       /,20X,'RF  harmonic          =',E15.4,
     >       /,20X,'Peak  voltage         =',E15.4,' V',
     >       /,20X,'Synchronous  phase    =',E15.4,' rd',
     >       /,20X,'RF  frequency         =',E15.4,' Hz',
     >       /,20X,'qV.SIN(Phi_s)         =',E15.4,' MeV',
     >       /,20X,'cos(Phi_s)            =',E15.4,' ',
     >       /,20X,'Nu_s/sqrt(alpha)      =',E15.4,'  ',
     >       /,20X,'dp-acc*sqrt(alpha)    =',E15.4,'  ')
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
C F(7,I) is time in mu_s. Of course, TI is in s
        TI = F(7,I) * UNIT(7) 
C        PHI = OMRF * TI + PHS
        DUM = SCALE4(F(7,I),WF1(I))
        PHI = SCALER(IPASS,NOEL,
     >                          coTime) + PHS

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

        TIOLD = TI
        PHIOLD = PHI


        F(1,I) = P/P0
        F(3,I) = ATAN(PY/PX)*1000.D0
        F(5,I) = ATAN(PZ/SQRT(PX*PX+PY*PY))*1000.D0

        IF(OKIMP) THEN
          SCALA = SCALER(IPASS,NOEL,
     >                              DTA1)
          WRITE(LUN,FMT='(1P,I6,1x,5G14.6,1x,I6,A)') 
     >            IPASS, OMRF/(2.*PI),
C     >            scala, 
C     >            PHI-PHS,
     >            PHI,
     >            TI,
     >             WF,
C     >              QV*SIN(PH(I))/(Q*1.D-6), 
     >               P-PS, I, 
     >            '       ! ipass freq phi ti WF p-ps ITraj'
C     >            ' ipass freq phi-phs ti qV*sin p-ps ITraj'
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
     >  /,20X,'Cavity  frequency                 =',1P,E15.6,' Hz',
     >  /,20X,'Harmonic                          =',   E15.6,' ',
     >  /,20X,'Max energy  gain                  =',   E15.6,' MeV',
     >  /,20X,'TOF for BRho_ref (',G10.2,') is',       E15.6,' s',
     >  /,20X,'Cumulated distance since origin   =',   E15.6,' m',
     >  /,20X,'Cumulated   TOF      "     "      =',   E15.6,' s',
     >  /,20X,'Particle mass                     =',   E15.6,' MeV/c2',
     >  /,20X,'         charge                   =',   E15.6,' C')
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
     >     WRITE(LUN,FMT='(1P,5e14.6,2I6)') PH(I),DPR(I),
     >     TI, ti-dble(ipass)*dts, QV*SIN(PH(I))/(Q*1.D-6), I , IPASS

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
     >  /,20X,'Cavity  frequency                 =',1P,E15.6,' Hz',
     >  /,20X,'Harmonic                          =',   E15.6,' ',
     >  /,20X,'Max energy  gain                  =',   E15.6,' MeV',
     >  /,20X,'TOF for BRho_ref (',G10.2,') is',       E15.6,' s',
     >  /,20X,'Cumulated distance since origin   =',   E15.6,' m',
     >  /,20X,'Cumulated   TOF      "     "      =',   E15.6,' s',
     >  /,20X,'Particle mass                     =',   E15.6,' MeV/c2',
     >  /,20X,'         charge                   =',   E15.6,' C')
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

      synch_time = synch_time + 1.d0/FREV*AN11


      IF(NRES.GT.0) THEN
        GTRNUS = SQRT(ABS(QV*AN11*COS(PHS) / (2.D0*PI*WS)))
        ACCDP=  SQRT(QV/(AN11*WS)) * 
     >  SQRT(ABS(-(2.D0*COS(PHS)/PI+(2.D0*PHS/PI-1.D0)*SIN(PHS))))
        DGDT = QV*SIN(PHS)/(ORBL/CL)/AM 
        WRITE(NRES,220) ORBL, AN11, QV/(Q *1.D-6), FREV, PHS, DTS, 
     >       QV*SIN(PHS),COS(PHS),GTRNUS, ACCDP,DGDT,
     >         QV/(Q *1.D-6)*SIN(PHS)/ORBL
 220    FORMAT(1P, 
     >       /,20X,'Orbit  length           =',E15.4,' m',
     >       /,20X,'RF  harmonic            =',E15.4,
     >       /,20X,'Peak  voltage           =',E15.4,' V',
     >       /,20X,'RF  frequency           =',E15.4,' Hz',
     >       /,20X,'Synchronous  phase      =',E15.4,' rd',
     >       /,20X,'Isochronous  time       =',  E15.4,' s',
     >       /,20X,'qV.SIN(Phi_s)           =',E15.4,' MeV',
     >       /,20X,'cos(Phi_s)              =',E15.4,' ',
     >       /,20X,'Nu_s/sqrt(alpha)        =',E15.4,'  ',
     >       /,20X,'dp-acc*sqrt(alpha)      =',E15.4,'  '
     >       /,20X,'dgamma/dt               =',E15.4,' /s '
     >       /,20X,'rho*dB/dt               =',E15.4,' T.m/s '
     >       )
        IF(SRLOSS) WRITE(NRES,FMT='(1P,
     >       /,20X,''SR loss compensation    ='',E19.8,'' MeV'')')
     >       DWSR

C        IF(KCAV .EQ. 1) WRITE(NRES,199) SCALER(IPASS+1,NOEL,DTA1,DTA2,DTA3)
C 199    FORMAT(/,20X,'Post acceleration SCALING factor is ',1P,G16.8)
C            write(33,*) time, QV/(Q/QE *1.D-6)/1.6D6, QV*SIN(PHS)/1.3, 
C     >        ps/16, AN11/DTS/(53e6), GTRNUS/40/0.11
      ENDIF

cC----- Initial conditions of  the IMAX particles
      IF ( PHS .NE. PHS_prev .AND. IPASS .GT. 1 ) then
        DO I=1,IMAX
           PHAS(I) = PHAS(I) + PHS - PHS_prev
        ENDDO
      ENDIF
      PHS_prev = PHS
 
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
C At all pass#,  F6i is the previous-turn path length (see below : F(6,I) set to 0), 
C DTI is the time it took since the last passage in CAVITE 
        DTI = F(6,I)*.01D0 / (BTA*CL)
        IF(IPASS .EQ. 1) THEN
           PHAS(I) = PHS + (qv/abs(qv))*(DTI-DTS)*OMRF
        ELSE
           PHAS(I) = PHAS(I) + (qv/abs(qv))*(DTI-DTS)*OMRF
        ENDIF
        IF(PHAS(I) .GT.  PI) PHAS(I) =PHAS(I) -2.D0*PI
        IF(PHAS(I) .LT. -PI) PHAS(I) =PHAS(I) +2.D0*PI
        WF = WF1(I) + QV*SIN(PHAS(I))
        if(srloss) wf = wf + dwsr
        WF1(I) = WF
        P = SQRT(WF*(WF + 2.D0*AMQ(1,I)))
        PX=SQRT( P*P -PY*PY-PZ*PZ)
 

C        DPR(I)=P-PS
        DPR(I)=(P-PS)/PS
        PH(I)=PHAS(I) 
C        BLAG=(PHAS(I)-PHS)/OMRF
C        BLNG=BLAG*(BTA*CL)
C        PH(I)=BLAG
     

        F(1,I) = P/P0
        F(3,I) = ATAN(PY/PX)*1000.D0
        F(5,I) = ATAN(PZ/SQRT(PX*PX+PY*PY))*1000.D0
 
        F(6,I)=0.D0 

        IF(OKIMP) 
     >  WRITE(LUN,FMT='(1P,7(E14.6,1X),2(I6,1X),10(1X,E14.6),A)') 
     >  PH(I),PHS,P-PS,OMRF,DTI,DTS,QV*SIN(PH(I))/(Q*1.D-6), I,IPASS,
     >  ORBL, HARM, BTA,BTS, OMRF,FREV,DWS,WKS,PS,WS,
     >  ' PH(I),PHS,P-PS,OMRF,DTI,DTS,QV*SIN(PH(I))/(Q*1.D-6),I,IPASS'
     >  //' ORBL, HARM, BTA,BTS, OMRF,FREV,DWS,WKS,PS,WS'

 3    CONTINUE
      GOTO 88

 110  CONTINUE
Cavite is modeled by Chambers style matrix, so accounting for transverse focusing.
C After J.Rosenzweig, L.Serafini, Phys Rev E Vo. 49, Num 2, 1994.
C Source code moved from BETA on Sept. 2015. Origin of phase is on >0 crest.
C Orbit length between 2 cavities, RF freq., phase of 1st cavity (ph0=0 is at V(t)=0)
      CAVL = AN10*1.d2     ! cavL in cm
      CALL SCUMW(0.5D0*CAVL/UNIT(5))
      CALL SCUMR(
     >           DUM,SCUM,TCUM) 
      FCAV = AN11          ! RF freq. in Hz
      PH0 = AN21  ! RF phase reference. Normally zero for e-linac
      PS = P0
      BTS = PS/SQRT(PS*PS+AM2)
      HARM = 1.d0
      DTS = SCUM*UNIT(5) / ( CL * BTS) 
      OMRF = 2.D0 * PI * FCAV
      IDMP = NINT(AN22)
 
      IF(NRES.GT.0) THEN
        WRITE(NRES,200) IDMP,
     >  TYPCH(IDMP+3)(DEBSTR(TYPCH(IDMP+3)):FINSTR(TYPCH(IDMP+3))),
     >  FCAV,cavl,QV,BORO,DTS,
     >  SCUM*UNIT(5),TCUM*UNIT(7),AM,Q*QE
 200    FORMAT(1P,
     >  /,15X,'CHAMBERS  CAVITY  STYLE',
     >  /,15X,'Transport option : ',I2,'  (',A,')',
     >  /,20X,'Cavity  frequency                     =',E15.6,' Hz',
     >  /,20X,'        length                        =',E15.6,' m',
     >  /,20X,'Max energy  gain                      =',E15.6,' MeV',
     >  /,20X,'TOF for BRho_ref (',G10.2,')          =',E15.6,' s',
     >  /,20X,'Cumulated distance at cavity center   =',E15.6,' m',
     >  /,20X,'Cumulated   TOF      "         "      =',E15.6,' s',
     >  /,20X,'Particle mass                         =',E15.6,' MeV/c2',
     >  /,20X,'         charge                       =',E15.6,' C')
      ENDIF

      DO I=1,IMAX
        IF(IEX(I) .GE. -1) THEN 
          P = P0*F(1,I) 
          AM2 = AMQ(1,I) * AMQ(1,I)
          ENRG = SQRT(P*P + AM2)
          WF1(I) = ENRG - AMQ(1,I)
          BTA = P / ENRG

          TI = F(7,I) * UNIT(7) 
          DSAR2=0.5D0*CAVL /(COS(F(3,I)*1.D-3)*COS(F(5,I)*1.D-3))
          TI = TI + dsar2 / (bta*cl)

C F(7,I) is time in mu_s. of course, TI is in s
          PHI = OMRF * TI + PH0   
C Phase, in [-pi,pi] interval

          WRITE(*,*) ' cav ',PHI, INT(PHI/(2.D0*PI))
     >               ,phi-pi,phi+pi
          PHI = PHI - INT(PHI/(2.D0*PI)) * 2.D0*PI 
          IF    (PHI .GT.  PI) THEN
            PH(I) =PHI - 2.D0*PI
          ELSEIF(PHI .LT. -PI) THEN
            PH(I) =PHI + 2.D0*PI
          ELSE
            PH(I) =PHI 
          ENDIF

          write(*,*) 'cavite ph(i)', phi, ph(i)
         

          COSRF=COS(PH(I))
          DWF=QV*COSRF
          WI = WF1(I) 
          WF1(I) = WF1(I) + DWF
          WF = WF1(I)
C Kin. energy, MeV
          DPR(I)=WF
          P = SQRT(WF*(WF + 2.D0*AMQ(1,I)))
        
          IF     (DWF.EQ.0.D0 .OR. IDMP.EQ.0.) THEN
C        CAVITY + DRIFT
            V11= 1.D0
            V12= 0.D0
            V21= 0.D0                               
            V22= 1.D0
          ELSE IF(IDMP.EQ.1) THEN
C        CHAMBERS CAVITY WITH DE/E<<1 APROXIMATION Det(M)#1
            FAC=SQRT(WI/WF)
            V11= FAC
            V12= CAVL*FAC
            V21= 0.D0                               
            V22= FAC
          ELSE IF(IDMP.EQ.-1) THEN
C        CHAMBERS CAVITY WITH DE/E<<1 APROXIMATION Det(M)=1
            FAC=SQRT(WI/WF)
            V11= FAC
            V12= CAVL*FAC
            V21= 0.D0     
            V22= FAC
            DWFT=DSQRT(V11*V22-V21*V12)
            V11=V11/DWFT
            V12=V12/DWFT
            V22=V22/DWFT
            V21=V21/DWFT
          ELSE IF(IDMP.EQ.2) THEN
C        CHAMBERS CAVITY Det(M)#1
            EFEI=WF/WI
            FAC=DLOG(EFEI)/SQRT8/COSRF
            COSFAC=COS(FAC)
            SINFAC=SIN(FAC)
            RAP=SQRT8*CAVL/DWF
            V11= COSFAC-SQRT2*SINFAC*COSRF
            V12=    SINFAC*RAP*WI*COSRF
            V21=-SINFAC/RAP/WF*(2*COSRF+1.D0/COSRF)
            V22=(COSFAC+SQRT2*SINFAC*COSRF)/EFEI
          ELSE IF(IDMP.EQ.-2) THEN
C        CHAMBERS CAVITY Det(M)=1
            EFEI=WF/WI 
            FAC=DLOG(EFEI)/SQRT8/COSRF
            COSFAC=COS(FAC)
            SINFAC=SIN(FAC)
            RAP=SQRT8*CAVL/DWF
            V11= COSFAC-SQRT2*SINFAC*COSRF
            V12=    SINFAC*RAP*WI*COSRF
            V21=-SINFAC/RAP/WF*(2*COSRF+1./COSRF)
            V22=(COSFAC+SQRT2*SINFAC*COSRF)/EFEI
            DWFT=DSQRT(V11*V22-V21*V12)
            V11=V11/DWFT
            V12=V12/DWFT
            V22=V22/DWFT
            V21=V21/DWFT
          ENDIF 

          write(*,*) 'cavite , WI, WF ', WI, WF 
                 read(*,*)

          F(1,I) = P / P0
          F(2,I) = (v11 * F(2,I)*.01D0 + v12 * F(3,I)*.001D0)*1.D2
          F(3,I) = (v21 * F(2,I)*.01D0 + v22 * F(3,I)*.001D0)*1.D3
          F(4,I) = (v11 * F(4,I)*.01D0 + v12 * F(5,I)*.001D0)*1.D2
          F(5,I) = (v21 * F(4,I)*.01D0 + v22 * F(5,I)*.001D0)*1.D3
          F(6,I) = F(6,I) + 2.D0 * dsar2/unit(5)
          F(7,I) = F(7,I) + 2.D0* dsar2 / (bta*cl) / unit(7)
          CALL SCUMW(0.5D0*CAVL/UNIT(5))

          IF(OKIMP) 
     >    WRITE(LUN,FMT='(1P,5e14.6,2I6,e14.6)') PH(I),DPR(I),
     >    TI, wi, wf, I , IPASS, DWF

        ENDIF   
      ENDDO

      CALL SCUMW(0.5D0*CAVL/UNIT(5))

      GOTO 88

 
 88   CONTINUE
      DPREF = PS / P0
      RETURN

 99   CONTINUE
      IF(NRES .GT. 0) WRITE(NRES,101) 'CAVITE',
     >                          ' OPEN zgoubi.CAVITE.Out'
 101  FORMAT(/,' ******  SBR ',A,' : ERROR',A,/)
      RETURN

 77   RETURN

      ENTRY CAVIT1(
     >             DWSO)
      DWSO = DWS
      RETURN

      ENTRY CAVIT3(
     >     A_synch_time)
      A_synch_time = synch_time
      RETURN 

      END
