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
C  Upton, NY, 11973,  USA
C  -------
      SUBROUTINE CAVITE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
C      PARAMETER (LNTA=132) ; CHARACTER(LNTA) TA
C      PARAMETER (MXTA=45)
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

      SAVE SYNCT
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
      PARAMETER (SQRT2 = SQRT(2.D0),SQRT8 = SQRT(8.D0))

      CHARACTER(60) TYPCH(5)
      INTEGER DEBSTR, FINSTR

      PARAMETER (CG=8.846D-14)  ! m/MeV^3

      LOGICAL STRCON
      PARAMETER (MXH=5)
      DIMENSION HRM(MXH), VHRM(MXH)

      LOGICAL CEBAF
      LOGICAL CRNLSY
      LOGICAL ERCRCS
 
      SAVE OMGA, BT0

      DATA WF1, PHAS / MXT*0.D0, MXT*0.D0 /
      DATA SKAV /'** OFF **','OPTION 1 ','OPTION 2 ','OPTION 3 ', 
     >   'OPTION 4 ', 'OPTION 5 ', '   FFAG  ', 'Isochron.',
     >  ' ' , ' ' , 'RLA', 'SR +accel.' /

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

      DATA CEBAF, CRNLSY, ERCRCS / .FALSE., .FALSE., .FALSE. /
      DATA SYNCT / 0.D0 /

      DUM = SCALER(1, NOEL, 
     >                     DUM)

      SRLOSS = .FALSE.
      U0 = 0.D0
      DWSR = 0.D0
      NBH = 1

      DUM = SCALE9(
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
        OKIMP = STRCON(TA(NOEL,1),'PRINT',
     >                                    IS)
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
        WRITE(NRES,102) SKAV(KCAV+1)
 102    FORMAT(15X,' Accelerating cavity. Type is :',3X,A,/)
        IF(OKIMP) WRITE(NRES,FMT='(15X,
     >  '' Cavite parameters saved in zgoubi.CAVITE.Out'',/)') 
      ENDIF

      IF(KCAV .EQ. 0) RETURN
 
      AN10 = A(NOEL,10)
      AN11 = A(NOEL,11)
      AN20= A(NOEL,20)
C      VLT= A(NOEL,20)
      AN21= A(NOEL,21)
      AN22= A(NOEL,22)
      AN23= A(NOEL,23)
      PHS= AN21
 
      IF(Q*AM .EQ. 0.D0) CALL ENDJOB('Pgm cavite. Give mass & charge'
     >//'of particles. Use ''PARTICUL'' keyword',-99)

C----- P0, AM  are  in  MEV/c, /c^2
      P0 = BORO*CL9*Q
      AM2= AM*AM
C      QV = VLT*Q *1.D-6
      QV = AN20 *Q *1.D-6
 

      GOTO(10,20,30,40,50,60,70,80,999,100,110) KCAV
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
      FRF = HARM/DTS
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
      NBH = NINT(A(NOEL,19))
      IF(NBH .LE. 0 .OR. NBH .GT. 2) 
     >CALL ENDJOB('Pgm cavite.  No such possibility NBH=',NBH)
      PHS= A(NOEL,20+NBH)
      IF(NBH .GT. 1) THEN
        DO I = 1, NBH
          HRM(I) = A(NOEL,10+I)
          VHRM(I) = A(NOEL,19+I)
        ENDDO
        AK = VHRM(2) / VHRM(1)
        RN = HRM(2) / HRM(1)
        PHN = ATAN(TAN(PHS)/RN) / RN
      ENDIF
C----- PARTICULE SYNCHRONE, ENTREE DE LA CAVITE
      IF(IPASS .EQ. 1) PS = P0
      BTS = PS/SQRT(PS*PS+AM2)
      DTS = ORBL / ( CL * BTS)
      OMRF = 2.D0*PI*HARM / DTS
      WKS = PS/BTS - AM
      FRF = HARM/DTS
C----- PARTICULE SYNCHRONE, SORTIE DE LA CAVITE
      DWS = QV*SIN(PHS)
      IF(NBH .EQ. 2) DWS = DWS + QV*AK*SIN(RN*PHN)
      WKS = WKS + DWS
C--- Case SR loss in storage ring (no acceleration). In that case PS is constant, DWS expected to be equal to SR loss.
      CALL SRLOS3(
     >            SRLOSS)
      IF(SRLOSS) WKS = WKS - DWS
      PS = SQRT(WKS*(WKS+2.D0*AM))
      WS = WKS + AM
      GOTO 1

C------------------------------------------- 
C Works like CAVITE/2 (tag 20) if no SR, or if SR in storage ring. Diff is allows SR with acceleration. 
C It uses ph_s to accelerate (not for compensation of SR unlike 20) (so, ph_s=0 
C in storage mode); ph_s is corrected for SR compensation : 
C requires A(noel,22) = theoretical SR loss at first pass, then, SR loss assumes ~gamma^4 
C dependence at subsequent passes. 
 110  CONTINUE
      CRNLSY = STRCON(TA(NOEL,1), 'CornellSynch',
     >                                         IS)
      ERCRCS = STRCON(TA(NOEL,1), 'eRHIC_RCS',
     >                                        IS)
      ORBL = AN10
      HARM = AN11
      RADIUS = AN12  ! ***************
C----- PARTICULE SYNCHRONE, ENTREE DE LA CAVITE
      IF(IPASS .EQ. 1) THEN
        STS = 0.D0
        IF(CRNLSY) OMGA = 2.D0*PI*60.D0
        PS = P0
        BT0 = P0/SQRT(P0*P0+AM2)
      ENDIF
      BTS = PS/SQRT(PS*PS+AM2)
      DTS = ORBL / ( CL * BTS)
      STS = STS + DTS
      FRF = HARM/DTS
      OMRF = 2.D0*PI*FRF
      WKS = PS/BTS - AM
C----- PARTICULE SYNCHRONE, SORTIE DE LA CAVITE
      CALL SRLOS3(
     >            SRLOSS)

      IF(SRLOSS) THEN 
        IF(ERCRCS) THEN
C        SI2 = 0.0568860943517   !  =I2=2pi/rho
C        U0 = CG * (PS/BTS)**4 / (2.D0*PI) * SI2
C        U0 = 0.51552876 /(12008.775542289)**4 *  (PS/BTS)**4 
C normalisation de U0=11.071834e-3 keV, a E=9.8932E+02MeV
c         U0 = 11.071834e-3 * ( (PS/BTS) / 9.8932E+02)**4           ! C_gamma=88.46276*1e-6 m/GeV^3


Case Vahid's eRHIC RCS
          radius = 283.860202518    ! bend radius (m), assumed isofield lattice
          u0 = 88.46276*((ps*1d-3)/bts)**4/radius *1d-3  ! u0 (MeV)  (elctrn with bta~1 : 88.463*E[GeV]^4/rho[m]*(Ang/2pi))
        ELSEIF(CRNLSY) THEN
          U0 = (A(NOEL,22)*1d-6) *( (ps/bts) / (p0/bt0) )**4  ! AN22(eV)=Energy loss at first pass. U0 in MeV 
********************
        ELSE
          U0 = 0.D0
        ENDIF
      ENDIF
      WKS = WKS -U0

      IF    (ERCRCS) THEN
C        QV = AN20 *Q *1.D-6 *(1.d0 + 1.8d-3* dble(ipass)) ! peak voltage. Varies from 15 to 120MV in 4000 turns 
c        QV = Q * (AN20 *1.D-6  + dble(ipass) * 76.14d-3) ! 1970 turns:   (final ^V - initial ^V) / number of turns = (180-30)/1970 
        QV = Q * (AN20 *1.D-6  + dble(ipass) * 375.14d-3) ! 400 turns:   = (180-30)/400 
C                                                     ! 120MV justified by SR loss at 20GeV being 50MV
C                                                     ! and willing sin(ph_s) to stay ~1/2
C      qvfrac = AN20 *Q *1.D-6 /2.d0       ! energy gain per turn
        qvfrac = AN20 *Q *1.D-6 * sin(2.8d0)       ! energy gain per turn
        phs = asin((qvfrac + u0)/qv)   ! |U0| is the energy loss per turn
        phs = pi - phs  
        DWS = QV*SIN(PHS)

      ELSEIF(CRNLSY) THEN
        QV = Q *(4.4D0* SIN(OMGA*STS) + 8.8D0*(SIN(OMGA*STS/2.D0) )**8)  ! MeV
C        PHS = PI/2.D0
        PHS = AN21
        DWS = QV*SIN(PHS) - U0

      ENDIF

c            write(*,*) ' cavite srloss ',srloss,u0,dws,qv

      WKS = WKS + DWS
      PS = SQRT(WKS*(WKS+2.D0*AM))
      WS = WKS + AM
C--- Case SR loss in storage ring or booster
c      CALL SRLOS3(
c     >            SRLOSS)
c      IF(SRLOSS) THEN
c        gg4 = ( ws / (ws-dws)   )**4
c        dwsr = a(noel,22) * gg4         
c      ENDIF
C      GOTO 1

C      SYNCT = SYNCT + 1.D0/FRF*AN11
C      SYNCT = STS

      IF(NRES.GT.0) THEN
        GTRNUS = SQRT(ABS(QV*AN11*COS(PHS) / (2.D0*PI*WS)))
        ACCDP=  SQRT(QV/(AN11*WS)) * 
     >  SQRT(ABS(-(2.D0*COS(PHS)/PI+(2.D0*PHS/PI-1.D0)*SIN(PHS))))
        DGDT = QV*SIN(PHS)/(ORBL/CL)/AM 
        WRITE(NRES,220) ORBL, AN11, QV/(Q *1.D-6), FRF, PHS, DTS, 
     >       QV*SIN(PHS),COS(PHS),GTRNUS, ACCDP,DGDT,
     >         QV/(Q *1.D-6)*SIN(PHS)/ORBL,U0
      ENDIF

      DO I=1,IMAX 

       IF(IEX(I) .GT. 0) THEN

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
c        if(srloss) wf = wf + dwsr
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
     >  WRITE(LUN,FMT='(1P,7(E14.6,1X),2(I6,1X),12(1X,E14.6),A)') 
     >  PH(I),PHS,P-PS,OMRF,DTI,DTS,QV*SIN(PH(I))/(Q*1.D-6), I,IPASS,
     >  ORBL, HARM, BTA,BTS, OMRF,FRF,DWS,WKS,PS,WS,U0,QV,
     >  ' PH(I),PHS,P-PS,OMRF,DTI,DTS,QV*SIN(PH(I))/(Q*1.D-6),I,IPASS'
     >  //' ORBL, HARM, BTA,BTS, OMRF,FRF,DWS,WKS,PS,WS,U0,QV'

       ENDIF

      ENDDO
      GOTO 88


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
      FRF = HARM /DTS
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
      FRF = HARM/DTS
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
      FRF0 = AN10/HN
      PS0 = SQRT((WS0+AM)**2 - AM2)
      BTS0 = PS0/SQRT(PS0*PS0+AM2)
C     ... Conditions at cavity entrance 
      OMRF = 2.D0*PI*FRF0*HN
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
C Bucketless acceleration. Ok for non-scaling FFAG and cyclotron
C Used for muon, EMMA, cyclotron
      CALL SCUMR(
     >           DUM,SCUM,TCUM) 
C Orbit length between 2 cavities, RF freq., phase of 1st cavity (ph0=0 is 
C at V(t)=0)
      HARM = AN10
      IF(HARM.LE.0.D0) HARM=1.D0   ! for compatibility w/ older version where AN10 was unused
      FCAV = AN11
      TREF = HARM / FCAV  ! sec. Synchronous time.
      PH0 = AN21
      PS = P0
      BTS = PS/SQRT(PS*PS+AM2)
      OMRF = 2.D0 * PI * FCAV
C      WS = PS / BTS
C      TS = TS + DTS
 
      IF(NRES.GT.0) THEN
        WRITE(NRES,170) FCAV,HARM,QV,
     >                    SCUM*UNIT(5),TCUM*UNIT(7),AM,Q*QE
 170    FORMAT(
     >  /,20X,'Cavity  frequency                 =',1P,E15.6,' Hz',
     >  /,20X,'Harmonic                          =',   E15.6,' ',
     >  /,20X,'Max energy  gain                  =',   E15.6,' MeV',
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

C Particle phase wrt. RF :
        DTI = TI - DBLE(IPASS) *TREF
        PHI = OMRF * DTI + PH0

C        PHI = OMRF * TI + PH0
C Phase, in [-pi,pi] interval
C        PHI = PHI - INT(PHI/(2.D0*PI)) * 2.D0*PI 
c        IF    (PHI .GT.  PI) THEN
c          PH(I) =PHI - 2.D0*PI
c        ELSEIF(PHI .LT. -PI) THEN
c          PH(I) =PHI + 2.D0*PI
c        ELSE
          PH(I) = modulo(PHI, 2.d0*pi)
c        ENDIF

C        DWF =  QV * SIN(PH(I))
        DWF =  QV * SIN(PHI)
C------- Rustine etude ffag muon
C        IF(OMRF.LE.0.D0) DWF = QV
C------------------------------

C Kin. energy, MeV
        WF1(I) = WF1(I) + DWF
        WF = WF1(I)
        DPR(I)=WF

           write(*,*) ' cavite p, px av : ',p,px

        P = SQRT(WF*(WF + 2.D0*AMQ(1,I)))
        PX=SQRT( P*P -PY*PY-PZ*PZ)

           write(*,*) ' cavite p, px ap : ',p,px
           write(*,*) ' '
        F(1,I) = P / P0
        F(3,I) = ATAN2(PY,PX) / UNIT(2)
        F(5,I) = ATAN2(PZ,SQRT(PX*PX+PY*PY)) / UNIT(4)

        IF(OKIMP) 
     >  WRITE(LUN,FMT='(1P,4(e14.6,1x),2(I6,1x),7(e14.6,1x),a)') 
     >  PHI,DWF,TI, SIN(PHI), I , IPASS
     >  ,phi/(2.d0*pi),omrf,omrf*ti,wf,ph(i),phi+ph0,DTI
     >  ,' phi, dwf, t, sin(ph+ph0),i,ipass,ph(i)/2pi,omrf,omrf*t,wf,'
     >  //'ph(i),phi+ph0'
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

      SYNCT = SYNCT + 1.D0/FRF*AN11

      IF(NRES.GT.0) THEN
        GTRNUS = SQRT(ABS(QV*AN11*COS(PHS) / (2.D0*PI*WS)))
        ACCDP=  SQRT(QV/(AN11*WS)) * 
     >  SQRT(ABS(-(2.D0*COS(PHS)/PI+(2.D0*PHS/PI-1.D0)*SIN(PHS))))
        DGDT = QV*SIN(PHS)/(ORBL/CL)/AM 
        WRITE(NRES,220) ORBL, AN11, QV/(Q *1.D-6), FRF, PHS, DTS, 
     >       QV*SIN(PHS),COS(PHS),GTRNUS, ACCDP,DGDT,
     >         QV/(Q *1.D-6)*SIN(PHS)/ORBL,U0
 220    FORMAT(1P, 
     >       /,20X,'Orbit  length           =',E19.8,' m',
     >       /,20X,'RF  harmonic            =',E19.8,
     >       /,20X,'Peak  voltage           =',E19.8,' V',
     >       /,20X,'RF  frequency           =',E19.8,' Hz',
     >       /,20X,'Synchronous  phase      =',E19.8,' rd',
     >       /,20X,'Isochronous  time       =',E19.8,' s',
     >       /,20X,'qV.sin(phi_s)           =',E19.8,' MeV',
     >       /,20X,'cos(phi_s)              =',E19.8,' ',
     >       /,20X,'Nu_s/sqrt(alpha)        =',E19.8,'  ',
     >       /,20X,'dp-acc*sqrt(alpha)      =',E19.8,'  '
     >       /,20X,'dgamma/dt               =',E19.8,' /s ',
     >       /,20X,'rho*dB/dt               =',E19.8,' T.m/s ',
     >       /,20X,'SR loss, this pass      =',E19.8,' MeV ',
     >       /)

        IF(NBH .EQ. 2) THEN
          WRITE(NRES,221) NINT(HRM(1)),NINT(HRM(2)),VHRM(1),VHRM(2),
     >    FRF,RN*FRF,NINT(RN),AK,ATAN(TAN(PHS)/RN)/RN,
     >    QV*(SIN(PHS)+AK*SIN(RN*PHN))
 221      FORMAT(1P, 
     >       /,20X,'h1, h2                   :    ',I0,',  ',I0,
     >       /,20X,'Peak  voltage  1, 2      :',E19.8,',  ',E19.8,' V',
     >       /,20X,'RF  frequency  1, 2      :',E19.8,',  ',E19.8,' Hz',
     >       /,20X,'N = h2 / h1              =    ',I0,
     >       /,20X,'k = V2 / V1              =',E19.8,
     >       /,20X,'phi_n = atan(tan(phi_s)/n)/n    =',E19.8,' rd',
     >       /,20X,'qV1.(sin(phi_s)+k.sin(n.phi_n)) =',E19.8,' MeV',
     >       /)
        ENDIF

        IF(SRLOSS) WRITE(NRES,FMT='(1P,
     >  20X,''SR loss compensation    ='',E19.8,'' MeV'',/)') DWSR


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
        IF(NBH .EQ. 2) WF = WF + QV*AK*SIN(RN*(PHAS(I) - PHS + PHN))

        IF(SRLOSS) WF = WF + DWSR

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
     >  WRITE(LUN,FMT='(1P,7(E14.6,1X),2(I6,1X),11(1X,E14.6),A)') 
     >  PH(I),PHS,P-PS,OMRF,DTI,DTS,QV*SIN(PH(I))/(Q*1.D-6), I,IPASS,
     >  ORBL, HARM, BTA,BTS, OMRF,FRF,DWS,WKS,PS,WS,U0,
     >  ' PH(I),PHS,P-PS,OMRF,DTI,DTS,QV*SIN(PH(I))/(Q*1.D-6),I,IPASS'
     >  //' ORBL, HARM, BTA,BTS, OMRF,FRF,DWS,WKS,PS,WS,U0'

 3    CONTINUE
      GOTO 88

 100  CONTINUE
Cavite is modeled by Chambers style matrix, so accounting for transverse focusing.
C After J.Rosenzweig, L.Serafini, Phys Rev E Vo. 49, Num 2, 1994.
C Source code moved from BETA on Sept. 2015. Origin of phase is on >0 crest.
C Orbit length between 2 cavities, RF freq., phase of 1st cavity (ph0=0 is at V(t)=0)
      CAVM = AN10          ! cavLength /m
      FCAV = AN11          ! RF freq. in Hz
      CAVL = CAVM*1.D2     ! cavLength /cm
      PH0 = AN21   -PI/2.D0        ! RF phase.  '-PI/2.D0' is just for convention
      IDMP = NINT(AN22) ! Chambers matrix options
      PHREF = AN23  ! Used in updating of DPREF

      CALL SCUMW(0.5D0*CAVL)
      CALL SCUMR(
     >           DUM,SCUM,TCUM) 
      PS = P0 * DPREF
      BTS = PS/SQRT(PS*PS+AM2)
      HARM = 1.D0
      OMRF = 2.D0 * PI * FCAV
      WS0 = SQRT(PS**2 +AM**2) - AM

      DWS = QV * COS(PH0)
      WSF = WS0 +  DWS
      PSF = SQRT((WSF+AM)**2 - AM2)

      IF(NRES.GT.0) THEN
        WRITE(NRES,200) IDMP,
     >  TYPCH(IDMP+3)(DEBSTR(TYPCH(IDMP+3)):FINSTR(TYPCH(IDMP+3))),
     >  FCAV,CAVM,PH0,DWS,WSF/WS0,BORO*DPREF,DPREF,BORO*PSF/P0,PSF/P0,
     >  SCUM*UNIT(5),TCUM,AM,Q*QE
 200    FORMAT(1P,
     >  / ,15X,'CHAMBERS  CAVITY  STYLE',
     >  //,15X,'Transport option : ',I2,'  (',A,')',
     >  / ,20X,'Cavity  frequency                     =',E15.6,' Hz',
     >  / ,20X,'        length                        =',E15.6,' m',
     >  / ,20X,'        RF phase phi_0                =',E15.6,' rad',
     >  / ,20X,'Synch energy gain qV.cos(phi_0)       =',E15.6,' MeV',
     >  / ,20X,'WF/WI                                 =',E15.6,' ',
     >  //,20X,'BRho_ref in (dp_ref in)               =',E14.6,' kG.cm',
     >                                                 3x,'(',E14.6,')',
     >  / ,20X,'BRho_ref out (dp_ref out)             =',E14.6,' kG.cm',
     >                                                 3x,'(',E14.6,')',
     >  / ,20X,'Cumulated distance at cavity center   =',E15.6,' m',
     >  / ,20X,'Cumulated   TOF      "         "      =',E15.6,' s',
     >  //,20X,'Particle mass                         =',E15.6,' MeV/c2'
     >  ,/ ,20X,'         charge                       =',E15.6,' C')
      ENDIF

      TIAV = 0.D0
      II = 0
      DO I=1,IMAX
        IF(IEX(I) .GE. -1) THEN 
          II = II + 1
          P = P0*F(1,I) 
          AM2 = AMQ(1,I) * AMQ(1,I)
          ENRG = SQRT(P*P + AM2)
          BTA = P / ENRG
Compute particle time at center of cavity 
C          TI = F(7,I) * UNIT(7) 
          DSAR2=0.5D0*CAVL /(COS(F(3,I)*1.D-3)*COS(F(5,I)*1.D-3))
          F7I = F(7,I) + (dsar2*unit(5)) / (bta*cl) / unit(7) 
C F(7,I) is time in mu_s. Of course, TI is in s
          TI = F7I * UNIT(7) 
          TIAV = TIAV + TI
        ENDIF   
      ENDDO
      TIAV = TIAV / DBLE(II)

      DO I=1,IMAX
        IF(IEX(I) .GE. -1) THEN 
          P = P0*F(1,I) 
          AM2 = AMQ(1,I) * AMQ(1,I)
          ENRG = SQRT(P*P + AM2)
          WF1(I) = ENRG - AMQ(1,I)
          BTA = P / ENRG

Compute particle time at center of cavity 
          DSAR2=0.5D0*CAVL /(COS(F(3,I)*1.D-3)*COS(F(5,I)*1.D-3))
          F(6,I) = F(6,I) + DSAR2
          F(7,I) = F(7,I) + (dsar2*unit(5)) / (bta*cl) / unit(7) 
C FM Dec 2015 - wrong units          TI = TI + dsar2 / (bta*cl)
C F(7,I) is time in mu_s. Of course, TI is in s
          TI = F(7,I) * UNIT(7) 

C Relative time to bunch centroid
          TI = TI-TIAV
          PHI = OMRF * TI + PH0   
          PHIAV = PHIAV + PHI
          PH(I) =PHI 
          DWF=QV*COS(PHI)
          CPH = COS(PHI)

CCCCCCCCCCCCCCCCCCCCCCCCCCC
c      CEBAF = STRCON(TA(NOEL,1), 'CEBAF',
c     >                                   IS)
c     >                                        IS)
C tests cebaf
c          IF(CEBAF) THEN 
c            IF    (PHI .GT.  PI) THEN
c              PH(I) =PHI - 2.D0*PI
c            ELSEIF(PHI .LT. -PI) THEN
c              PH(I) =PHI + 2.D0*PI
c            ELSE
c              PH(I) =PHI 
c            ENDIF
c          ENDIF
ccccccccccccccccccccccccccc

          WI = WF1(I) 
          WF1(I) = WF1(I) + DWF
          WF = WF1(I)
C Kin. energy, MeV
          DPR(I)=WF
          P = SQRT(WF*(WF + 2.D0*AMQ(1,I)))
          
          IF     (DWF.EQ.0.D0 .OR. IDMP.EQ.0.) THEN
C        CAVITY + DRIFT
            V11= 1.D0
            V12= CAVM
            V21= 0.D0                               
            V22= 1.D0
          ELSE IF(IDMP.EQ.1) THEN
C        CHAMBERS CAVITY WITH DE/E<<1 APROXIMATION Det(M)#1
            FAC=SQRT(WI/WF)
            V11= FAC
            V12= CAVM*FAC
            V21= 0.D0                               
            V22= FAC
          ELSE IF(IDMP.EQ.-1) THEN
C        CHAMBERS CAVITY WITH DE/E<<1 APROXIMATION Det(M)=1
            FAC=SQRT(WI/WF)
            V11= FAC
            V12= CAVM*FAC
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
            FAC=DLOG(EFEI)/SQRT8/CPH
            COSFAC=COS(FAC)
            SINFAC=SIN(FAC)
            RAP=SQRT8*CAVM/DWF
            V11= COSFAC-SQRT2*SINFAC*CPH
            V12=    SINFAC*RAP*WI*CPH
            V21=-SINFAC/RAP/WF*(2*CPH+1.D0/CPH)
            V22=(COSFAC+SQRT2*SINFAC*CPH)/EFEI
          ELSE IF(IDMP.EQ.-2) THEN
C        CHAMBERS CAVITY Det(M)=1
            EFEI=WF/WI 
            FAC=DLOG(EFEI)/SQRT8/CPH
            COSFAC=COS(FAC)
            SINFAC=SIN(FAC)
            RAP=SQRT8*CAVM/DWF
            V11= COSFAC-SQRT2*SINFAC*CPH
            V12=    SINFAC*RAP*WI*CPH
            V21=-SINFAC/RAP/WF*(2*CPH+1./CPH)
            V22=(COSFAC+SQRT2*SINFAC*CPH)/EFEI
            DWFT=DSQRT(V11*V22-V21*V12)
            V11=V11/DWFT
            V12=V12/DWFT
            V22=V22/DWFT
            V21=V21/DWFT
          ENDIF 

          F(1,I) = P / P0
          F(2,I) = (v11 * F(2,I)*.01D0 + v12 * F(3,I)*.001D0)*1.D2
          F(3,I) = (v21 * F(2,I)*.01D0 + v22 * F(3,I)*.001D0)*1.D3
          F(4,I) = (v11 * F(4,I)*.01D0 + v12 * F(5,I)*.001D0)*1.D2
          F(5,I) = (v21 * F(4,I)*.01D0 + v22 * F(5,I)*.001D0)*1.D3
          DSAR2=0.5D0*CAVL /(COS(F(3,I)*1.D-3)*COS(F(5,I)*1.D-3))           
          F(6,I) = F(6,I) + DSAR2
          BTA = P / SQRT(P*P + AM2)
          F(7,I) = F(7,I) + (dsar2*unit(5)) / (bta*cl) / unit(7) 

          IF(OKIMP) 
     >    WRITE(LUN,FMT='(1P,5e14.6,2I6,e14.6)') PH(I),DPR(I),
     >    TI, WI, WF, I , IPASS, DWF
          
c          ii = ii + 1

        ENDIF   
      ENDDO
      PHIAV = PHIAV / DBLE(II)

      PS = PSF
      CALL SCUMW(0.5D0*CAVL)

      IF(NRES.GT.0) THEN
C        TIAV = TIAV / DBLE(II)
C        PHIAV = PHIAV / DBLE(II)
        WRITE(NRES,fmt='(1P,
     >   /,20X,''Averaged over the '',I0,'' particles : '',
     >  /,25X,''- <arrival time> at cavity      = '',2E15.6,'' mu_s'',
     >  /,25X,''- and resulting <phase>         = '',E15.6,
     >  /,25X,''- resulting qV.cos(<phase>)     = '',E15.6,'' MeV''
     >  )') ii, Tiav, tiav2, phiav, QV*COS(PHiav)
      ENDIF

      GOTO 88
 
 88   CONTINUE
      DPREF = PS / P0

      RETURN

 99   CONTINUE
      IF(NRES .GT. 0) WRITE(NRES,101) 'CAVITE',
     >                          ' OPEN zgoubi.CAVITE.Out'
 101  FORMAT(/,' ******  SBR ',A,' : ERROR',A,/)
      RETURN

 999  RETURN   

      ENTRY CAVIT1(
     >             DWSO)
      DWSO = DWS
      RETURN

      ENTRY CAVIT3(
     >             SYNCTO)
      SYNCTO = SYNCT
      RETURN 

      END 
