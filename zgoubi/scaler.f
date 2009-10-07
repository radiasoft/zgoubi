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
      FUNCTION SCALER(IPASS,NOEL, 
     >                           D1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      COMMON/PTICUL/ AM,Q,G,TO
      INCLUDE 'MXFS.H'
      COMMON/SCAL/SCL(MXF,MXS),TIM(MXF,MXS),NTIM(MXF),KSCL
      INCLUDE 'MXLD.H'
      CHARACTER FAM*8,LBF*8,KLEY*10,LABEL*8
      COMMON/SCALT/ FAM(MXF),LBF(MXF,2),KLEY,LABEL(MXL,2)

      LOGICAL EMPTY

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
C------------ Current KLEY will undergo scaling if either it as the right label...
            DO 11 KL=1,MXLF
 11           IF(LABEL(NOEL,1).EQ. LBF(KF,KL)) GOTO 2
          ELSE
C------------ ...or it has no label at all
            GOTO 2
          ENDIF
        ENDIF
 10   CONTINUE

      RETURN

 2    CONTINUE

      KTI = NTIM(KF)

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
     >                  PP0)
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
          IF( IPASS .GE. TIM(KF,1) .AND. IPASS .LE. TIM(KF,2)  ) THEN
          ENDIF
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
      DP = FACDP * SIN(ARGDP*TIME) * DTIM
      PSYN = PSYN + DP
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
