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
C  Upton, LI, NY, 11973
C  -------
      FUNCTION SCALER(IPASS,NOEL, 
     >                           D1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      PARAMETER (LBLSIZ=10)
      INCLUDE 'MXLD.H'
      CHARACTER(LBLSIZ) LABEL
      COMMON /LABEL/ LABEL(MXL,2)
      COMMON/PTICUL/ AM,Q,G,TO
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      INCLUDE 'MXFS.H'
      INCLUDE 'MXSCL.H'
      COMMON/SCAL/ SCL(MXF,MXS,MXSCL),TIM(MXF,MXS),NTIM(MXF),KSCL
      COMMON/SCALP/ VPA(MXF,MXP),JPA(MXF,MXP)
      PARAMETER (KSIZ=10)
      CHARACTER(KSIZ) FAM
      CHARACTER(LBLSIZ) LBF
      COMMON/SCALT/ FAM(MXF),LBF(MXF,MLF)
      INCLUDE "MAXTRA.H"
      COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)

      CHARACTER(KSIZ) KLEY

      LOGICAL EMPTY, STRCON

      SAVE TIME, ICTIM
      SAVE ARGDP, FACDP, PSYN
      SAVE TEMP

      PARAMETER (ND=200000)
      DIMENSION XM(ND), YM(ND), XMI(ND), YMI(ND)
      dimension dat1(ND), dat2(ND), dat3(ND), 
     >          dat1i(ND),dat2i(ND),dat3i(ND)
      SAVE XM, YM, NFRQ, dat1, dat2, dat3
      SAVE OCLOCK, ekin

      DIMENSION MODSCL(MXF)

      DIMENSION SCL2I(MXF,MXD),TIM2I(MXF,MXD),NTIM2I(MXF)
      DIMENSION SCL2(MXF,MXD), TIM2(MXF,MXD), NTIM2(MXF)
      SAVE SCL2, TIM2, NTIM2

      DIMENSION KFM(MXSCL), KFMO(MXSCL)
      SAVE KFM 
      INTEGER DEBSTR, FINSTR
      logical ok3

      DATA TIME, ICTIM / 0.D0 , 0/
      DATA TEMP / 1.D0 /

      DATA XM, YM / ND*0.D0, ND*0.D0/
      data ok3 / .true. /

      DOUBLE PRECISION SPLINT

      ok3 = .true.
      SCALER = 1.D0
      CALL SCALI5(
     >            MODSCL,NFAM)
      CALL ZGKLEY(
     >            KLEY)
C----- Looks whether current kley is registered for scaling (in FAM(KF), when declared in 'SCALING'). 
C        Looks for possible limitation due to LABEL[s] associated with FAM(KF). 
      do i = 1, MXSCL
        KFM(i) = -99
      enddo

      KF1 = 1
      IFM = 0

 3    continue

      DO KF = KF1, NFAM

        IF(KLEY .EQ. FAM(KF)) THEN
C--------- Current KLEY recorded for scaling 
          
          IF( .NOT. EMPTY(LBF(KF,1)) ) THEN
C------------ Current KLEY will undergo scaling if...

             
            DO  KL=1,MLF
              IF(STRCON(LBF(KF,KL),'*',
     >                                 IS)) THEN
                LBFA = DEBSTR(LBF(KF,KL))
                LBFB = FINSTR(LBF(KF,KL))
                LLBF = LBFB-LBFA+1
                LABA = DEBSTR(LABEL(NOEL,1))
                LABB = FINSTR(LABEL(NOEL,1))
                LLAB = LABB-LABA+1

c               if(noel.eq.684 .or. noel .eq. 759 ) then
c                write(*,*) 'scaler  NOEL, IS, ',NOEL,IS,
c     >        LABEL(NOEL,1)(1:LLBF-1),' ' , LBF(KF,KL)(1:LLBF-1)
c                write(*,*) 'scaler  ',IS,
c     >        LABEL(NOEL,1)(LLAB-LLBF+2:LLAB),' ', LBF(KF,KL)(2:LBFB)
c              write(*,*) 
c     >        (LABEL(NOEL,1)(1:LLBF-1) .EQ. LBF(KF,KL)(1:LLBF-1))
c              write(*,*) 
c     >        (LABEL(NOEL,1)(LLAB-LLBF+2:LLAB).EQ.LBF(KF,KL)(2:LBFB))
c                    read(*,*)
c                 endif

C               ... either LBF ends with '*' ...
                IF(  LABEL(NOEL,1)(1:LLBF-1) .EQ. LBF(KF,KL)(1:LLBF-1))
     >               goto 2
                IF( (LLAB-LLBF+2).GT.0) then !yann to protect -1 in LABEL table
                   IF(LABEL(NOEL,1)(LLAB-LLBF+2:LLAB)
     >                  .EQ.LBF(KF,KL)(2:LBFB)) 
     >                  GOTO 2
                ENDIF

              ELSE
C               ... or it as the right label...
                IF(LABEL(NOEL,1).EQ. LBF(KF,KL)) THEN
                  ok3 = .false.
                  GOTO 2
c      write(88,*) ' scaler lab?        '
c     >                         ,LABEL(NOEL,1), LBF(KF,KL),noel,kf,kl
                ENDIF
              ENDIF
            ENDDO
          ELSE
C------------ ...or if it has no label at all
c      write(88,*) ' scaler Nolab ',LABEL(NOEL,1), LBF(KF,KL),noel,kf,kl
            GOTO 2
          ENDIF
        ENDIF

      ENDDO

      GOTO 99

 2    CONTINUE

c      write(88,*) ' scaler >2 ',LABEL(NOEL,1), LBF(KF,KL),noel,kf,kl

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

              IF( IPASS .GE. IT1 .AND. IPASS .LE. IT2 ) THEN
                SCALER = SCL(KF,I,1)
                IF(IT2 .NE. IT1) SCALER= SCALER+ (SCL(KF,I2,1)-SCALER )*
     >             DBLE( IPASS - IT1 ) / DBLE(IT2 - IT1)
C FM 08/99
C     >              DBLE( IPASS - IT1 ) / (1.D0+ IT2 - IT1 )
                GOTO 99

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

c                    if(kf.eq.1)
cC                     if(noel.eq.17)
c     >               write(66,*) ' scaler ',scaler, 
c     >               scal2,scaler*scal2, BRO,yt1,yt2,kf,it
c                     stop ' SBR scaler -  TESTS'

                GOTO 99
              ENDIF
            ELSEIF(MODSCL(KF) .EQ. 12  ) THEN
               SCALER = 1.d0

            ELSEIF(MODSCL(KF) .EQ. 13  ) THEN
               XT1=TIM(KF,I)
               IF(I .LT. KTI ) THEN
                  I2 = I+1
                  XT2 = TIM(KF,I2)
               ELSE
                  I2 = I
                  XT2=XT1
               ENDIF  
               CALL CAVIT3(STIME)
            
             
               IF( STIME .GE. XT1 .AND. STIME .LE. XT2 ) THEN
                  SCALER = SCL(KF,I,1)



                  IF(XT2 .NE. XT1) SCALER= SCALER+ (SCL(KF,I2,1)
     >                 - SCALER)*( STIME - XT1 ) / (XT2 - XT1)
               GOTO 99
               ENDIF


            ENDIF

 1        CONTINUE

        ELSEIF(KTI .EQ. -1) THEN

c          call cavit1(
c     >                PP0,GAMMA,dWs)
c          scaler = SCL(KF,1) * pp0
          scaler = SCL(KF,1,1) * dpref

        ELSEIF(KTI .EQ. -2) THEN
C--------- Field law for scaling FFAG, LPSC, Sept. 2007
c          xv = ipass
c          scaler = CUBSPL(xm,ym,xv,nd,nfrq)/ym(1)
          xv = OCLOCK
C                         time   phase
          scaler = CUBSPL(xm,    ym,    xv,nd,nfrq)
C                         
          xv = ekin
C                          ekin freq
          coTime1 = CUBSPL(dat3,dat2,xv,nd,nfrq)
          coTime = 1.d0/CUBSPL(dat3,dat2,xv,nd,nfrq)
          D1 = coTime
C        write(*,*) ' scaler ',xv, cotime1, coTime

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

        ELSEIF(KTI .EQ. -88) THEN
C-------- Field law AC dipole for Mei, Delta-Airlines, 2nd Oct. 2009
          PHAS = SCL(KF,1,1)
          Q1   = SCL(KF,2,1)
          Q2   = SCL(KF,3,1)
          PP   = SCL(KF,4,1)
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
            SCALER = SCALER * COS(2.D0*PP*DBLIP*QN + PHAS)
          ELSEIF(IPASS .GT. RAMPN+FLATN+DOWNN) THEN
            SCALER = 1.D0
          ENDIF
C          write(88,*) ' scaler ',IPASS,scaler,NINT(RAMPN+FLATN+DOWNN)

        ELSEIF(KTI .EQ. -87) THEN
C AGS Q-Jump quads
c            call cavit1(
c     >                  PP0,GAMMA,dWs)
            call cavit1(
     >                  dWs)
            if(am.le.0.d0) call endjob(' SBR scaler : need mass >',0)
            P2 = dpref * BORO*CL9*Q
            GAMMA = SQRT(P2 + AM*AM)/AM
            gg = G*GAMMA
            if(GG .GT. TIM(KF,1)) THEN
C             Q-jump started 
              switch=1.d0
              dN = TIM(KF,2)
              dTrn = TIM(KF,3)
              dgg = G * dWs/AM
              dint = G*GAMMA - INT(G*GAMMA)
              trmp1 = dN - (0.5d0*dTrn)*dgg
              if    (dint .le. trmp1) then
                fac=0.d0
              else
                trmp2 = dN + (0.5d0*dTrn)*dgg
                if(dint.gt.trmp1 .and. dint .lt. trmp2) then
                  fac = (dint-trmp1) / (trmp2-trmp1)
                else
                  trmp3 = (1.d0-dN) - (0.5d0*dTrn)*dgg
                  if(dint.lt.trmp3) then
                    fac = 1.d0
                  else
                    trmp4 = (1.d0-dN) + (0.5d0*dTrn)*dgg
                    if(dint.lt.trmp4) then
                      fac = (trmp4-dint) / (trmp4-trmp3)
                    else
                      fac = 0.d0  
                    endif
                  endif
                endif
              endif
            else
              switch=0.d0
            endif           
            
            if(switch.ne.0.d0) then
C              scaler = fac * SCL(KF,1) * pp0
              scaler = fac * SCL(KF,1,1) * dpref
            else
              scaler = 0.d0
            endif

C             write(77,fmt='(i6,8(1x,F10.4),a)')
C     >       ipass,fac,gg,dint,trmp1,trmp2,trmp3,trmp4,dws,' scaler'
             SZ = SF(3,1)
c             write(77,fmt='(i6,8(1x,F10.4),a)')
c     >       ipass,fac,gg,SZ,dint,trmp1,trmp2,trmp3,trmp4,' scaler'

        ELSE
          STOP 'FCTN SCALER :  invalid input data NTIM(1)'
        ENDIF
      ELSE
        IF(NRES .GT. 0) WRITE(NRES,*) ' FCT SCALER : timing is empty'
      ENDIF

      KF1 = KFM(ifm) + 1
      if(ok3) goto 3

 99   CONTINUE


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
      do JT = 1, NTIM2(JF)
         SCL2(JF,JT) = SCL2I(JF,JT) 
         TIM2(JF,JT) = TIM2I(JF,JT) 
      enddo

c               if(jf.eq.1 .or. jf.eq.2)   then
cc               if(jf.eq.1 )   then
c                 write(*,*) ' SCALE2 / if,iit : ',Jf,iit
c                 write(*,*) '  (SCL2(jF,IT),it=1,iit) : '
c                 write(*,*) (SCL2(jF,IT),it=1,iit)
c                 write(*,*) ' (TIM2(jF,IT),it=1,iit) : '
c                 write(*,*) (TIM2(jF,IT),it=1,iit)
c                 write(*,*) ' '
c              endif

      RETURN

      ENTRY SCALE4(OCLOCI,ekinI)
      OCLOCK = OCLOCI
      ekin = ekinI 
      SCALE4 = 0.D0
      RETURN

      ENTRY SCALE6(xmi,ymi,dat1i,dat2i,dat3i,nfrqi)
C               oclock phi turn#  freq  Ekin 
      nfrq = nfrqi
      do 6 i=1,nfrq
        xm(i) = xmi(i)
        ym(i) = ymi(i)
        dat1(i) = dat1i(i)
        dat2(i) = dat2i(i)
        dat3(i) = dat3i(i)
C          write(*,*) i, dat2(i), dat3(i)
 6    continue
      SCALE6 = 0.D0
      RETURN

c      ENTRY SCALE8(FREV0I,E0I,AKI)
c      FREV0 = FREV0I
c      E0 = E0I
c      AK = AKI
c      SCALE8 = 0.D0 
c      RETURN

      ENTRY SCALE9(
     >             KFMO)
      do j = 1, MXSCL
        KFMO(J) = KFM(J)
        IF((KFMO(j).GE.0)) then
        IF ((NTIM(KFMO(J)).NE.-1 )) then
           BRO = BORO*DPREF
           DO I=1, JPA(KFM(j),MXP) 
              VPA(KFM(j),I) = splint(TIM(KFM(j),1:NTIM(KFM(j))),
     >             SCL(KFM(j),1:NTIM(KFM(j)),I),NTIM(KFM(j)),BRO)
           ENDDO
        ENDIF   
        ENDIF
       enddo
      END
