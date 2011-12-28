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
      COMMON/PTICUL/ AM,Q,G,TO
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      INCLUDE 'MXFS.H'
      COMMON/SCAL/SCL(MXF,MXS),TIM(MXF,MXS),NTIM(MXF),KSCL
      INCLUDE 'MXLD.H'
      PARAMETER (LBLSIZ=8)
      PARAMETER (KSIZ=10)
      CHARACTER FAM*(KSIZ),LBF*(LBLSIZ),KLEY*(KSIZ),LABEL*(LBLSIZ)
      COMMON/SCALT/ FAM(MXF),LBF(MXF,2),KLEY,LABEL(MXL,2)
      INCLUDE "MAXTRA.H"
      COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)

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

      DATA TIME, ICTIM / 0.D0 , 0/
      DATA TEMP / 1.D0 /

      DATA XM, YM / ND*0.D0, ND*0.D0/

      SCALER = 1.D0

C----- Looks whether current kley is registered for scaling (in FAM(KF)). 
C        Looks for possible limitation due to LABEL[s] associated with FAM(KF). 
      DO 10 KF = 1, MXF
        IF(KLEY .EQ. FAM(KF)) THEN
C--------- Current KLEY recorded for scaling 

          IF( .NOT. EMPTY(LBF(KF,1)) ) THEN
C------------ Current KLEY will undergo scaling...
            DO  KL=1,MXLF
              IF(STRCON(LBF(KF,KL),'*',
     >                                 IS)) THEN
C               ... if either LBF ends with '*' ...
                IF(LABEL(NOEL,1)(1:IS-1) .EQ. LBF(KF,KL)(1:IS-1)) GOTO 2
              ELSE
C               ... or it as the right label...
                IF(LABEL(NOEL,1).EQ. LBF(KF,KL)) GOTO 2
              ENDIF
            ENDDO
          ELSE
C------------ ...or it has no label at all
            GOTO 2
          ENDIF
        ENDIF
 10   CONTINUE

      RETURN

 2    CONTINUE

      KTI = NTIM(KF)
C         print*, kti, ' scaler '
      IF(KTI .NE. 0) THEN

        IF(KTI .GT. 0) THEN

          DO 1 I=1,KTI
            IT1=TIM(KF,I)
            IF(I .LT. KTI ) THEN
              I2 = I+1
              IT2 = TIM(KF,I2)
            ELSE
              I2 = I
              IT2=IT1
            ENDIF  
 
            IF( IPASS .GE. IT1 .AND. IPASS .LE. IT2 ) THEN
              SCALER = SCL(KF,I)
              IF(IT2 .NE. IT1) SCALER = SCALER + (SCL(KF,I2) - SCALER )*
     >           DBLE( IPASS - IT1 ) / DBLE(IT2 - IT1)
C FM 08/99
C     >            DBLE( IPASS - IT1 ) / (1.D0+ IT2 - IT1 )
              RETURN
            ENDIF
 1        CONTINUE

        ELSEIF(KTI .EQ. -1) THEN

          call cavit1(
     >                PP0,GAMMA,dWs)
          scaler = SCL(KF,1) * pp0

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
            BRMI = SCL(KF,1)
            BRMA = SCL(KF,2)
            BREF = SCL(KF,3)
            FREP = SCL(KF,4)
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
          PHAS = SCL(KF,1)
          Q1   = SCL(KF,2)
          Q2   = SCL(KF,3)
          PP   = SCL(KF,4)
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
          write(88,*) ' scaler ',IPASS,scaler,NINT(RAMPN+FLATN+DOWNN)

        ELSEIF(KTI .EQ. -87) THEN
C AGS Q-Jump quads
            call cavit1(
     >                  PP0,GAMMA,dWs)
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
              scaler = fac * SCL(KF,1) * pp0
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

      ENTRY SCALE8(FREV0I,E0I,AKI)
      FREV0 = FREV0I
      E0 = E0I
      AK = AKI
      SCALE8 = 0.D0 
      RETURN

      END
