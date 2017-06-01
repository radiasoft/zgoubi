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
C  Upton, LI, NY, 11973, USA
C  -------
      FUNCTION SCALER(IPASS,NOEL, 
     >                           D1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      PARAMETER (LBLSIZ=20)
      INCLUDE 'MXLD.H'
      CHARACTER(LBLSIZ) LABEL
      INCLUDE "C.LABEL.H"     ! COMMON/LABEL/ LABEL(MXL,2)
      INCLUDE "C.PTICUL.H"     ! COMMON/PTICUL/ AM,Q,G,TO
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      INCLUDE 'MXFS.H'
      INCLUDE 'MXSCL.H'
      INCLUDE "C.SCAL.H"     ! COMMON/SCAL/ SCL(MXF,MXS,MXSCL),TIM(MXF,MXS),NTIM(MXF),KSCL
      INCLUDE "C.SCALP.H"     ! COMMON/SCALP/ VPA(MXF,MXP),JPA(MXF,MXP)
      PARAMETER (KSIZ=10)
      CHARACTER(KSIZ) FAM
      CHARACTER(LBLSIZ) LBF
      INCLUDE "C.SCALT.H"     ! COMMON/SCALT/ FAM(MXF),LBF(MXF,MLF)
      INCLUDE "MAXTRA.H"
      INCLUDE "C.SPIN.H"     ! COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)

      CHARACTER(KSIZ) KLEY

      LOGICAL EMPTY, STRCON, SCALEX

      SAVE TIME, ICTIM
      SAVE ARGDP, FACDP, PSYN
      SAVE TEMP

      PARAMETER (ND=70000)
      DIMENSION XM(ND), YM(ND), XMI(ND), YMI(ND)
      DIMENSION DAT1(ND), DAT2(ND), DAT3(ND), 
     >          DAT1I(ND),DAT2I(ND),DAT3I(ND)
      SAVE XM, YM, NFRQ, DAT1, DAT2, DAT3
      SAVE OCLOCK, EKIN

      DIMENSION MODSCL(MXF)

      DIMENSION SCL2I(MXF,MXD),TIM2I(MXF,MXD),NTIM2I(MXF)
      DIMENSION SCL2(MXF,MXD), TIM2(MXF,MXD), NTIM2(MXF)
      SAVE SCL2, TIM2, NTIM2

      DIMENSION KFM(MXSCL), KFMO(MXSCL)
      SAVE KFM 
      INTEGER DEBSTR, FINSTR
      LOGICAL OK3

      LOGICAL IDLUNI, OK, OKOPN, OKPRT, OKPRTI
      LOGICAL STRWLD

      PARAMETER (NLB=1)
      CHARACTER(LBLSIZ) LBLST(NLB)

      SAVE OKPRT, OKOPN, LPRT

      DATA TIME, ICTIM / 0.D0 , 0/
      DATA TEMP / 1.D0 /

      DATA XM, YM / ND*0.D0, ND*0.D0/
      DATA OK3 / .TRUE. /
      DATA OKPRT / .FALSE. /

      OK3 = .TRUE.
      SCALER = 1.D0
      CALL SCALI5(
     >            MODSCL,NFAM)
C      IF(NFAM .GT. MXF ) CALL ENDJOB('Pgm scaler. Need NFAM < ',MXF)
      CALL ZGKLEY(
     >            KLEY)
C----- Looks whether current kley is registered for scaling (in FAM(KF), when declared in 'SCALING'). 
C        Looks for possible limitation due to LABEL[s] associated with FAM(KF). 
      DO I = 1, MXSCL
        KFM(I) = -99
      ENDDO

      KF1 = 1
      IFM = 0

 3    CONTINUE

      DO KF = KF1, NFAM

        IF(KLEY .EQ. FAM(KF)) THEN
C--------- CURRENT KLEY RECORDED FOR SCALING 
        
          IF( .NOT. EMPTY(LBF(KF,1)) ) THEN
C------------ LBF(,1-MLF) has to match current KLEY's label 1 in part or in full

            KL = 1
            DOWHILE(.NOT. EMPTY(LBF(KF,KL)) .AND. KL.LE.MLF)  

C              IF(STRCON(LBF(KF,KL),'*',
C     >                                 IS)) THEN
C                LBFA = DEBSTR(LBF(KF,KL))
C                LBFB = FINSTR(LBF(KF,KL))
C                LLBF = LBFB-LBFA+1
C                LABA = DEBSTR(LABEL(NOEL,1))
C                LABB = FINSTR(LABEL(NOEL,1))
C                LLAB = LABB-LABA+1
C
CC               ... either LBF ends with '*' ...
C                IF(  LLBF-1 .GE. 1 ) THEN
C                  IF(  LABEL(NOEL,1)(1:LLBF-1)
C     >                      .EQ. LBF(KF,KL)(1:LLBF-1)) GOTO 2
C                ELSE
CC Dec 2015. To be checked before release
CCC Means that LBF='*'. Scaling applies to all keywords
CC                  GOTO 2
C                ENDIF


                  lblst(1) = LBF(KF,KL)

c              write(*,*) ' scaler ',lblst(1) ,LABEL(NOEL,1),
c     >           STRWLD(NLB,LBLST,LABEL(NOEL,1),IS),is,
c     >             LABEL(NOEL,2),
c     >           STRWLD(NLB,LBLST,LABEL(NOEL,2),IS),is


              IF  (STRWLD(NLB,LBLST,LABEL(NOEL,1),
     >                                          IS) 
     >        .OR. STRWLD(NLB,LBLST,LABEL(NOEL,2),
     >                                          IS)) THEN
              
                 jj = 1
                 IF(STRWLD(NLB,LBLST,LABEL(NOEL,2),
     >                                          IS))  jj = 2

c                 write(*,*) ' scaler ',jj,LABEL(NOEL,1),
c     >             LABEL(NOEL,2)
C                         read(*,*)
                  ok3 = .false.

                LBFA = DEBSTR(LBF(KF,KL))
                LBFB = FINSTR(LBF(KF,KL))
                LLBF = LBFB-LBFA+1
                LABA = DEBSTR(LABEL(NOEL,jj))
                LABB = FINSTR(LABEL(NOEL,jj))
                LLAB = LABB-LABA+1
                IF( (LLAB-LLBF+2).GT.0) THEN                 !yann to protect -1 in LABEL table
                   IF(LABEL(NOEL,jj)(LLAB-LLBF+2:LLAB)
     >                  .EQ.LBF(KF,KL)(2:LBFB)) GOTO 2
                ENDIF

                GOTO 2

              ELSE
C               ... or it as the right label...
                IF(LABEL(NOEL,1).EQ. LBF(KF,KL)) THEN
C Yann, Oct 2014 : commented to allow for multiple scaling "lines/families" to point to the same element.
C                  ok3 = .false.

                  GOTO 2
                ELSE
           
                ENDIF

              ENDIF

              KL = KL + 1
            ENDDO

          ELSE

            GOTO 2

          ENDIF
        ENDIF
       
      ENDDO

      GOTO 88

 2    CONTINUE


C      IF(IFM .EQ. MXSCL) CALL ENDJOB('Pgm scaler. Need IFM < ',MXSCL)
      IFM = IFM + 1
      KFM(IFM) = KF

      KTI = NTIM(KF)

      IF(KTI .NE. 0) THEN

        IF(KTI .GT. 0) THEN

          DO 1 I=1,KTI
 
            IF    (MODSCL(KF) .LT. 10) THEN

              IT1=NINT(TIM(KF,I))
              IF(I .LT. KTI ) THEN
                I2 = I+1
                IT2 = NINT(TIM(KF,I2))
              ELSE
                I2 = I
                IT2=IT1
              ENDIF  

C       IF(I2 .GT. MXS) CALL ENDJOB('Pgm scaler. Need I2 < ',MXS)

              IF( IPASS .GE. IT1 .AND. IPASS .LE. IT2 ) THEN
                SCALER = SCL(KF,I,1)
                IF(IT2 .NE. IT1) SCALER= SCALER+ (SCL(KF,I2,1)-SCALER )*
     >             DBLE( IPASS - IT1 ) / DBLE(IT2 - IT1)
C FM 08/99
C     >              DBLE( IPASS - IT1 ) / (1.D0+ IT2 - IT1 )

C FM 17/15
C                GOTO 88

              ENDIF

            ELSEIF((MODSCL(KF) .EQ. 10) .OR.
     >             (MODSCL(KF) .EQ. 11)) THEN

              XT1=TIM(KF,I)
              IF(I .LT. KTI ) THEN
                I2 = I+1
                XT2 = TIM(KF,I2)
              ELSE
                I2 = I
                XT2=XT1
              ENDIF  
              BRO = BORO * DPREF

C              write(66,*) kf, bro, dpref, xt1,xt2
              IF( BRO .GE. XT1 .AND. BRO .LE. XT2 ) THEN

C        IF(I .GT. MXS) CALL ENDJOB('Pgm scaler. Need I2 < ',MXS)
                SCALER = SCL(KF,I,1)
                IF(XT2 .NE. XT1) SCALER= SCALER+ (SCL(KF,I2,1)- SCALER)*
     >                 ( BRO - XT1 ) / (XT2 - XT1)
               
                IF(MODSCL(KF) .EQ. 11) THEN

                  SCAL2 = 1.D0

                  IIT = NTIM2(KF)
                  DO IT = 1, IIT
                    YT1=TIM2(KF,IT)
                    IF(IT .LT. IIT ) THEN
                      I2 = IT+1
                      YT2 = TIM2(KF,I2)
                    ELSE
                      I2 = IT
                      YT2=YT1
                    ENDIF  
                    BRO = BORO * DPREF

c                      write(66,*) ' scaler SCL2(KF,IT) ',SCL2(KF,IT)
c                      write(66,*) ' scaler ',BRO,YT1 ,YT2 

C         IF(IT .GT. MXD) CALL ENDJOB('Pgm scaler. Need IT < ',MXD)

                    IF( BRO .GE. YT1 .AND. BRO .LE. YT2 ) THEN
                      SCAL2 = SCL2(KF,IT)
                      IF(YT2 .NE. YT1) SCAL2= SCAL2+ 
     >                  (SCL2(KF,I2)- SCAL2) * (BRO - YT1) / (YT2 - YT1)

                      GOTO 11
                    ENDIF
                  ENDDO

 11               CONTINUE

                  SCALER = SCALER * SCAL2

                ENDIF

                GOTO 88
              ENDIF
            ELSEIF(MODSCL(KF) .EQ. 12  ) THEN
               SCALER = 1.D0

            ELSEIF(MODSCL(KF) .EQ. 13  ) THEN
               XT1=TIM(KF,I)
               IF(I .LT. KTI ) THEN
                  I2 = I+1
                  XT2 = TIM(KF,I2)
               ELSE
                  I2 = I
                  XT2=XT1
               ENDIF  
               CALL CAVIT3(
     >                     STIME)
            
             
               IF( STIME .GE. XT1 .AND. STIME .LE. XT2 ) THEN
                  SCALER = SCL(KF,I,1)

                  IF(XT2 .NE. XT1) SCALER= SCALER+ (SCL(KF,I2,1)
     >                 - SCALER)*( STIME - XT1 ) / (XT2 - XT1)
               GOTO 88
               ENDIF


            ENDIF

 1        CONTINUE

        ELSEIF(KTI .EQ. -1) THEN

c          call cavit1(
c     >                PP0,GAMMA,dWs)
c          scaler = SCL(KF,1) * pp0
          SCALER = SCL(KF,1,1) * DPREF
        
c              write(*,*) 'scaler kti=-1. kf, scal ',KF,SCALER
C                   read(*,*)
    
        ELSEIF(KTI .EQ. -2) THEN
C--------- Field law for scaling FFAG, LPSC, Sept. 2007
c          xv = ipass
c          scaler = CUBSPL(xm,ym,xv,nd,nfrq)/ym(1)
          XV = OCLOCK
C                         time   phase
          SCALER = CUBSPL(XM,    YM,    XV,ND,NFRQ)
C                         
          IF(NFRQ.GT.ND) CALL ENDJOB(' SBR SCALER, too many data '
     >    //'sent to cubspl. Max allowed is ND = ',ND)

          XV = EKIN
C                          ekin freq
          COTIME1 = CUBSPL(DAT3,DAT2,XV,ND,NFRQ)
          COTIME = 1.D0/CUBSPL(DAT3,DAT2,XV,ND,NFRQ)
          D1 = COTIME

C        ELSEIF(KTI .EQ. -60) THEN
CC--------- AGS dipoles, K1 and K2 laws
C          CALL MULTP2(.TRUE.)
C          CALL CAVIT1(
C     >                PP0,GAMMA,DWS)
CC          print *, PP0,GAMMA,DWS
CC            print *, PP0,BORO,CL9,Q,PP0*BORO*CL9/1.D3 
C          CALL AGSKS(PP0*BORO*CL9/1.D3)
C          scaler = SCL(KF,1) * pp0

        ELSEIF(KTI .EQ. -77) THEN
C-------- Field law protn driver, FNAL, Nov.2000
          IF( IPASS .GE. TIM(KF,1) .AND. IPASS .LE. TIM(KF,2)  ) THEN
            BRMI = SCL(KF,1,1)
            BRMA = SCL(KF,2,1)
            BREF = SCL(KF,3,1)
            FREP = SCL(KF,4,1)
            OMGAT = 2.D0*PI*FREP*TIME
            
C            TEMP = 0.5D0*( (BRMA+BRMI) - (BRMA-BRMI)* ( COS(OMGAT) 
C     >       - 0.25D0*SIN(2.D0*OMGAT) ) ) / BREF 
C            SCALER=TEMP 
            SCALER=TEMP * BRMI
C------- CONV converts from Brho (kG.cm) to momentum (MeV/c)
            CONV=1.D3 * CL9*Q
            FACDP = PI*FREP*(BRMA-BRMI) * CONV
            ARGDP = 2.D0*PI*FREP
C------------ PSYN(MeV/c) = mom. of synchronous particle
            IF(IPASS .EQ. TIM(KF,1)) PSYN = BRMI * CONV
C            WRITE(NLOG,*) TIME, SCALER, p/p0, IPASS, NOEL,
C     >            ' SBR SCALER :   TIME, SCALER, p/p0, IPASS, NOEL'
          ENDIF

        ELSEIF(KTI .EQ. -89) THEN
C-------- Field law AC dipole for Mei, Delta-Airlines, 2nd Oct. 2009
          PHAS = SCL(KF,2,1)
          Q1   = SCL(KF,3,1)
          Q2   = SCL(KF,4,1)
          PP   = SCL(KF,5,1)
          BSC   = SCL(KF,1,1)  
          RAMPN = TIM(KF,1)
          FLATN = TIM(KF,2)
          DOWNN = TIM(KF,3)
          DBLIP = DBLE(IPASS)
          IF    (IPASS .LE. RAMPN+FLATN+DOWNN) THEN
            IF    (IPASS .LE. RAMPN) THEN
              SCALER = DBLIP/RAMPN
              QN = Q1
            ELSEIF(IPASS .GT. RAMPN .AND. IPASS .LE. RAMPN+FLATN) THEN
              SCALER = 1.D0
              QN = Q1+(Q2-Q1) * (DBLIP-RAMPN)/FLATN
            ELSEIF(IPASS .GT. RAMPN+FLATN .AND. 
     >                              IPASS .LE. RAMPN+FLATN+DOWNN) THEN
              SCALER = (RAMPN+FLATN+DOWNN-DBLIP)/DOWNN
              QN = Q2
            ENDIF
            SCALER = BSC * SCALER  * COS(2.D0*PP*DBLIP*QN + PHAS)
          ELSEIF(IPASS .GT. RAMPN+FLATN+DOWNN) THEN
            SCALER = 0.D0
          ENDIF
C          write(88,*) ' scaler ',IPASS,scaler,NINT(RAMPN+FLATN+DOWNN)

        ELSEIF(KTI .EQ. -88) THEN
C-------- Field law AC dipole for Mei, Delta-Airlines, 2nd Oct. 2009
C-------- removed 2pi from input, added constant scaler SCLMX, jorg
          PHAS = SCL(KF,1,1)
          Q1   = SCL(KF,2,1)
          Q2   = SCL(KF,3,1)
          SCLMX = SCL(KF,4,1)
          RININ = TIM(KF,1)
          RAMPN = TIM(KF,2)
          FLATN = TIM(KF,3)
          DOWNN = TIM(KF,4)
          acturns = DBLE(IPASS)-RININ
          rampstart = RININ
          sweepstart = RININ+RAMPN
          sweepend = sweepstart+FLATN
          downend = sweepend+DOWNN
c          write(88,*) phas, q1, q2, SCLMX, RININ, RAMPN,
c     &      FLATN,DOWNN,acturns,rampstart,sweepstart,sweepend,downend

          IF (IPASS .LT. RININ .or. IPASS .gt. downend) THEN
             SCAL = 0.d0
             QN=0
             iwhere=1
             totalphasep=0
          ELSE
             IF  (IPASS .LE. sweepstart) THEN
                   SCAL = acturns/RAMPN
                   QN = Q1
                   totalphasep=QN*acturns
c             write(77,*) acturns,RAMPN,QN,acturns
             iwhere=2
             ELSEIF  (IPASS .LE. sweepend) THEN
                   SCAL = 1.D0
                   totalturn=IPASS
                   sweepturns = acturns-RAMPN
                   sweeprate = (Q2-Q1)/FLATN
                   QN = Q1+ sweeprate * sweepturns
                   totalphasep=Q1*acturns+
     +               sweeprate*sweepturns*sweepturns/2.d0
             iwhere=3
             ELSE
                   SCAL = (RAMPN+FLATN+DOWNN-acturns)/DOWNN
                   QN = Q2
                   totalphasep=QN*(acturns-RAMPN-FLATN)+Q1*(RAMPN+FLATN)
             iwhere=4
             ENDIF
 
c                write(77,*) SCAL,2.D0*PI,totalphasep, PHAS, SCLMX
             SCAL = SCAL * COS(2.D0*PI*totalphasep + PHAS) * SCLMX
c                write(77,*) SCAL,ipass,iwhere

          ENDIF
c         scal is used for debugging -jkl
          SCALER = scal

c           IF (IPASS .ge. RININ ) then
c          write(*,*) IPASS,QN,totalphasep,scaler
c             read(*,*)
c           endif

c         write(77,*) ipass,noel,scaler,qn,totalphasep,PHAS,Q1,Q2,
c     &         SCLMX,RININ,RAMPN,FLATN,
c     &      DOWNN,acturns,rampstart,sweepstart,sweepend,downend, iwhere
c         write(77,*) ipass,noel,scal,qn,totalphasep,iwhere,phas



C        ELSEIF(KTI .EQ. -88) THEN
CC-------- Field law AC dipole for Mei, Delta-Airlines, 2nd Oct. 2009
C          PHAS = SCL(KF,1,1)
C          Q1   = SCL(KF,2,1)
C          Q2   = SCL(KF,3,1)
C          PP   = SCL(KF,4,1)
C          RAMPN = TIM(KF,1)
C          FLATN = TIM(KF,2)
C          DOWNN = TIM(KF,3)
C          DBLIP = DBLE(IPASS)
C          IF    (IPASS .LE. RAMPN+FLATN+DOWNN) THEN
C            IF    (IPASS .LE. RAMPN) THEN
C              SCALER = DBLIP/RAMPN
C              QN = Q1
C            ELSEIF(IPASS .GT. RAMPN .AND. IPASS .LE. RAMPN+FLATN) THEN
C              SCALER = 1.D0
C              QN = Q1+(Q2-Q1) * (DBLIP-RAMPN)/FLATN
C            ELSEIF(IPASS .GT. RAMPN+FLATN .AND. 
C     >                              IPASS .LE. RAMPN+FLATN+DOWNN) THEN
C              SCALER = (RAMPN+FLATN+DOWNN-DBLIP)/DOWNN
C              QN = Q2
C            ENDIF
C            SCALER = SCALER * COS(2.D0*PP*DBLIP*QN + PHAS)
C          ELSEIF(IPASS .GT. RAMPN+FLATN+DOWNN) THEN
C            SCALER = 1.D0
C          ENDIF
C          write(88,*) ' scaler ',IPASS,scaler,NINT(RAMPN+FLATN+DOWNN)

        ELSEIF(KTI .EQ. -87) THEN
C AGS Q-Jump quads
c            call cavit1(
c     >                  PP0,GAMMA,dWs)
            CALL CAVIT1(
     >                  DWS)
            IF(AM.LE.0.D0) CALL ENDJOB(' SBR scaler : need mass >',0)
            P2 = DPREF * BORO*CL9*Q
            GAMMA = SQRT(P2 + AM*AM)/AM
            GG = G*GAMMA
            IF(GG .GT. TIM(KF,1)) THEN
C             Q-JUMP STARTED 
              SWITCH=1.D0
              DN = TIM(KF,2)
              DTRN = TIM(KF,3)
              DGG = G * DWS/AM
              DINT = G*GAMMA - INT(G*GAMMA)
              TRMP1 = DN - (0.5D0*DTRN)*DGG
              IF    (DINT .LE. TRMP1) THEN
                FAC=0.D0
              ELSE
                TRMP2 = DN + (0.5D0*DTRN)*DGG
                IF(DINT.GT.TRMP1 .AND. DINT .LT. TRMP2) THEN
                  FAC = (DINT-TRMP1) / (TRMP2-TRMP1)
                ELSE
                  TRMP3 = (1.D0-DN) - (0.5D0*DTRN)*DGG
                  IF(DINT.LT.TRMP3) THEN
                    FAC = 1.D0
                  ELSE
                    TRMP4 = (1.D0-DN) + (0.5D0*DTRN)*DGG
                    IF(DINT.LT.TRMP4) THEN
                      FAC = (TRMP4-DINT) / (TRMP4-TRMP3)
                    ELSE
                      FAC = 0.D0  
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ELSE
              SWITCH=0.D0
            ENDIF           
            
            IF(SWITCH.NE.0.D0) THEN
C              SCALER = FAC * SCL(KF,1) * PP0
              SCALER = FAC * SCL(KF,1,1) * DPREF
            ELSE
              SCALER = 0.D0
            ENDIF
             SZ = SF(3,1)

        ELSE
          STOP 'FCTN SCALER :  invalid input data NTIM(1)'
        ENDIF
      ELSE
        IF(NRES .GT. 0) WRITE(NRES,*) ' FCT SCALER : timing is empty'
      ENDIF

      KF1 = KFM(IFM) + 1
      IF(OK3) GOTO 3

 88   CONTINUE

      IF(OKPRT) WRITE(LPRT,*) IPASS,SCALER,NOEL,KLEY,LABEL(NOEL,1),
     >LABEL(NOEL,2)

c                 write(*,*) ' scaler ',jj,LABEL(NOEL,1),
c     >             LABEL(NOEL,2), scaler

      RETURN

      ENTRY SCALDP(DTIM,
     >                  TIMOUT)
      TIME = TIME + DTIM
      TIMOUT = TIME
      ICTIM=ICTIM+1
      PSYN = PSYN + FACDP * SIN(ARGDP*TIME) * DTIM
      TEMP = (PSYN)/( BRMI * CONV)
      SCALDP = PSYN
      RETURN

      ENTRY SCALPS()
      SCALPS = PSYN
      RETURN

      ENTRY SCALE2(SCL2I,TIM2I,NTIM2I,JF)
      NTIM2(JF) = NTIM2I(JF)
      DO JT = 1, NTIM2(JF)
         SCL2(JF,JT) = SCL2I(JF,JT) 
         TIM2(JF,JT) = TIM2I(JF,JT) 
      ENDDO
      SCALE2 = 0.D0
      RETURN

      ENTRY SCALE4(OCLOCI,EKINI)
      OCLOCK = OCLOCI
      EKIN = EKINI 
      SCALE4 = 0.D0
      RETURN

      ENTRY SCALE6(XMI,YMI,DAT1I,DAT2I,DAT3I,NFRQI)
C               OCLOCK PHI TURN#  FREQ  EKIN 
      NFRQ = NFRQI
      IF(NFRQ.GT.ND) CALL ENDJOB(' SBR SCALER, too many data '
     >//'sent to cubspl. Max allowed is ND = ',ND)
      DO 6 I=1,NFRQ
        XM(I) = XMI(I)
        YM(I) = YMI(I)
        DAT1(I) = DAT1I(I)
        DAT2(I) = DAT2I(I)
        DAT3(I) = DAT3I(I)
 6    CONTINUE
      SCALE6 = 0.D0
      RETURN

C      ENTRY SCALE8(FREV0I,E0I,AKI)
C      FREV0 = FREV0I
C      E0 = E0I
C      AK = AKI
C      SCALE8 = 0.D0 
C      RETURN

      ENTRY SCALE9(
     >             KFMO)
      DO J = 1, MXSCL
        KFMO(J) = KFM(J)
        IF((KFMO(J).GE.0)) THEN
         IF ((NTIM(KFMO(J)).NE.-1 )) THEN
           BRO = BORO*DPREF
           DO I=1, JPA(KFM(J),MXP) 
              VPA(KFM(J),I) = SPLINT(TIM(KFM(J),1:NTIM(KFM(J))),
     >        SCL(KFM(J),1:NTIM(KFM(J)),I),NTIM(KFM(J)),BRO)
           ENDDO
         ENDIF   
        ENDIF
      ENDDO
      SCALE9 = 0.D0 
      RETURN

      ENTRY SCALEX(OKPRTI)
      OKPRT = OKPRTI
      IF(OKPRT) THEN 
        IF(.NOT.OKOPN) THEN
          OK =(IDLUNI(
     >                LPRT)) 
          OPEN(UNIT=LPRT,FILE='zgoubi.SCALING.Out',
     >                   ERR=66,IOSTAT=IOS)
          WRITE(LPRT,*) '#  IPASS, SCALER, lmnt #, KLEY, '
     >    //'LABEL(NOEL,1), LABEL(NOEL,2)'
          OKOPN = .TRUE.
        ENDIF
      ELSE
        IF(OKOPN) THEN
          CLOSE(LPRT)
          OKOPN = .FALSE.
        ENDIF
      ENDIF
      GOTO 67
 66   CONTINUE
      IF(NRES.GT.0) WRITE(NRES,FMT='(A)') 
     >'                *** WARNING, pgm scaler : ' 
     >//' could not open zgoubi.SCALING.out. Will skip printing.'
      OKPRT = .FALSE.
 67   CONTINUE
      IF(NRES.GT.0) THEN
        IF(OKPRT) THEN 
          WRITE(NRES,FMT='(/,15X,''PRINT option is ON, '')') 
          WRITE(NRES,FMT='(  15X,
     >    ''SCALING parameters will be logged in zgoubi.SCALING.Out.'',
     >    /)') 
        ELSE
          WRITE(NRES,FMT='(/,15X,''PRINT option is OFF.'')') 
        ENDIF
      ENDIF
      SCALEX = OKPRT
      RETURN

      END
