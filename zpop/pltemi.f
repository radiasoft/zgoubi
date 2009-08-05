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
      SUBROUTINE PLTEMI(NL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ---------------------------
C     Compute and plot emittances
C     ---------------------------
      COMMON/CONST/ CL,PI,DPI,RAD,DEG,QE,AH
      LOGICAL OKECH, OKVAR, OKBIN
      COMMON/ECHL/OKECH, OKVAR, OKBIN
      COMMON/LUN/NDAT,NRES,NPLT,NFAI,NMAP,NSPN
      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB
      CHARACTER LET
      INCLUDE 'MXVAR.H'
      DIMENSION YZXB(MXVAR),NDX(5)
      INCLUDE 'MXLD.H'
      PARAMETER(MSC=6,MSC2=9)
      DIMENSION NPART(MXL)
      DIMENSION SC(MSC,MXL),SC2(MSC2,MXL),SCUM(MXL),EsPi(3,MXL)
      DIMENSION GAMMA(3,MXL),ALPHA(3,MXL),BETA(3,MXL)
      DIMENSION DEM(3,MXL), DEMI(3), DEMA(3), JEL(MXL)
      LOGICAL OKOPN
      PARAMETER(ZERO=0.D0)
      CHARACTER REP
      LOGICAL CHANGE,NEWP,NEWL, BINARY,EMPTY
      INCLUDE 'FILFAI.H'
      CHARACTER*20 FNAME, FNAM
      SAVE FNAME, BINARY
      INTEGER DEBSTR,FINSTR
      DATA FNAME / FILFAI /
      DATA MXPAS0, MXLM0 / 99999 , 99999/

 79     WRITE(6,FMT='(/,''  New calculation (N/Y) ?'')')
        READ(5,FMT='(A)',ERR=79) REP
        CHANGE = REP .EQ. 'Y' .OR. REP .EQ. 'y'
        IF(.NOT.CHANGE) GOTO 70

 80     WRITE(6,FMT='(/,''Give name of file to open - default is '',A)')  
     >                                      FNAME
        READ(5,FMT='(A)',ERR=80) FNAM
        IF(EMPTY(FNAM)) FNAM=FNAME
        FNAME=FNAM(DEBSTR(FNAM):FINSTR(FNAM))
        BINARY=FNAME(1:2).EQ.'B_' .OR. FNAME(1:2).EQ. 'b_'
        CLOSE(NL)
        CALL OPNDEF(NFAI,FNAME,NL,
     >                            FNAM,OKOPN)

 69     WRITE(6,FMT='(/,
     >    ''  Max pass, max lmnt (default is'',2I6,'') ?'')')
     >               MXPAS0, MXLM0
        READ(5,*,ERR=69,END=68) MXPAS, MXLM
 68     CONTINUE
        IF(MXLM.LE.0) MXLM=MXLM0
        IF(MXPAS.LE.0) MXPAS=MXPAS0
          WRITE(6,*) '   MXPAS, MXLM = ',MXPAS, MXLM    
      NRD=0
      NPASS = 0
      NELMAX = 0
      CALL RAZI(NPART,MXL)
      CALL RAZ(SC,MSC*MXL)
      WRITE(6,*)
      WRITE(6,*) ' Reading into ',FNAME
      WRITE(6,*) '        Calculating first order momenta'
      WRITE(6,*)
 44   CONTINUE
      NRD=NRD+1
C----- Read next coordinate
      CALL READCO(NL,
     >                  KART,LET,YZXB,NDX,*50,*59)
C     >                     KART,LET,YZXB,NDX,*50,*44)
C      IF(NDX(1) .LT. -1) GOTO 44

C      NEWP=NPASS.NE.NINT(YZXB(20))
C      NPASS=NINT(YZXB(20))
      NEWP=NPASS.NE.NINT(YZXB(39))
      NPASS=NINT(YZXB(39))
      IF(NEWP) THEN
C------- Detects start of a pass ; pass counter
        WRITE(6,*) ' Reading coordinates at PASS number ',NPASS
C------- Reset lmnt counter at start of this pass
        IF(NPASS.GT.MXPAS) GOTO 50 
        NEL=1
        IF(NEL.GT.NELMAX) NELMAX=NEL
      ELSE
        NEWL=NOEL.NE.NDX(5)
C--------- Gets into a new element
        IF(NEWL) THEN
          NEL=NEL+1
          IF(NEL.GT.NELMAX) NELMAX=NEL
        ENDIF
      ENDIF
      NOEL=NDX(5)

      X =YZXB(2)
      XP=YZXB(3)
      Y =YZXB(4)
      YP=YZXB(5)
C----- YZXB(7)=time or length, YZXB(1)=dp/p
      Z=YZXB(7)
      ZP=YZXB(1)

      NPART(NEL)=NPART(NEL)+1
C----- SC(1-6) -> x,x',y,y', dt,dp, current
      SC(1,NEL)=SC(1,NEL)+X
      SC(2,NEL)=SC(2,NEL)+XP
      SC(3,NEL)=SC(3,NEL)+Y
      SC(4,NEL)=SC(4,NEL)+YP
      SC(5,NEL)=0.D0          ! SC(5,NEL)+Z
      SC(6,NEL)=SC(6,NEL)+ZP

      GOTO 44          
C     ----------------------------------------------
 59   CONTINUE
      CALL FBGTXT
      WRITE(6,*) ' SBR PLTEMI: ERROR DURING READ OF EVENT #',NRD

 50   CONTINUE
      NEL=NELMAX
      WRITE(6,*) '         Data reading ended, at pass # ',NPASS
      WRITE(6,*) '                        now computing means' 
      DO 51 IEL=1, NEL
          XPART=FLOAT(NPART(IEL))
          DO 51 JC=1,MSC
            SC(JC,IEL)=SC(JC,IEL)/XPART
 51   CONTINUE

      CALL REWIN2(NL,*99)

      NRD=0
      NPASS = 0
      CALL RAZI(NPART,MXL)
      CALL RAZ(SC2,MSC2*MXL)
      CALL RAZ(SCUM,MXL)

      WRITE(6,*)
      WRITE(6,*) ' Reading into ',FNAME
      WRITE(6,*) '        Calculating second order momenta'
      WRITE(6,*)
 45   CONTINUE
      NRD=NRD+1
C----- Read next coordinate
      CALL READCO(NL,
     >                  KART,LET,YZXB,NDX,*10,*19)
C      IF(NDX(1) .LT. -1) GOTO 44

C      NEWP=NPASS.NE.NINT(YZXB(20))
C      NPASS=NINT(YZXB(20))
      NEWP=NPASS.NE.NINT(YZXB(39))
      NPASS=NINT(YZXB(39))
      IF(NEWP) THEN
C------- Detects start of a pass ; pass counter
        WRITE(6,*) ' Reading coordinates at PASS number ',NPASS
C------- Reset lmnt counter at start of this pass
        IF(NPASS.GT.MXPAS) GOTO 10 
        NEL=1
      ELSE
        NEWL=NOEL.NE.NDX(5)
C--------- Gets into a new element
        IF(NEWL) NEL=NEL+1
      ENDIF
      NOEL=NDX(5)

      X =YZXB(2)
      XP=YZXB(3)
      Y =YZXB(4)
      YP=YZXB(5)
C----- YZXB(7)=time, YZXB(1)=dp/p
      Z =YZXB(7)
      ZP=YZXB(1)

      NPART(NEL)=NPART(NEL)+1
      SC2(1,NEL)=SC2(1,NEL)+(X -SC(1,NEL))**2
      SC2(2,NEL)=SC2(2,NEL)+(XP-SC(2,NEL))**2
      SC2(3,NEL)=SC2(3,NEL)+(Y -SC(3,NEL))**2
      SC2(4,NEL)=SC2(4,NEL)+(YP-SC(4,NEL))**2
      SC2(5,NEL)=0.D0   ! SC2(5,NEL)+(Z -SC(5,NEL))**2
      SC2(6,NEL)=SC2(6,NEL)+(ZP-SC(6,NEL))**2
      SC2(7,NEL)=SC2(7,NEL)+(X -SC(1,NEL))*(XP-SC(2,NEL))
      SC2(8,NEL)=SC2(8,NEL)+(Y -SC(3,NEL))*(YP-SC(4,NEL))
      SC2(9,NEL)=0.D0   ! SC2(9,NEL)+(Z -SC(5,NEL))*(ZP-SC(6,NEL))

      SCUM(NEL)=SCUM(NEL)+YZXB(6)
      JEL(NEL)=NDX(5)

      GOTO 45
C     ------------------------------------------------

 19   CONTINUE
      CALL FBGTXT
      WRITE(6,*) ' SBR PLTEMI: ERROR DURING READ OF EVENT #',NRD

 10   CONTINUE
      NEL=NELMAX
      WRITE(6,*) '         Data reading ended, at pass # ',NPASS
      WRITE(6,*) '                       now ending calculations.' 
      DO 11 IEL=1,NEL
          XPART=NPART(IEL)
          SCUM(IEL)=SCUM(IEL)/XPART
          DO 11 JC=1,MSC2
            SC2(JC,IEL)=SC2(JC,IEL)/XPART
 11       CONTINUE

      CALL RAZ(DEMA,3)
      CALL RAZ(DEMI,3)
      WRITE(6,FMT='(/,A,/)')  ' Emittance growth : '
      DO 28 IEL=1,NEL
        DO 29 JC=1,MSC,2
          JJ=(JC+1)/2
C---------- Emittance/4
          SQ2=4.D0*(SC2(JC,IEL)*SC2(JC+1,IEL)-SC2(JJ+6,IEL)**2 )
          IF(SQ2.LT.0.D0) THEN
            WRITE(6,*) ' Warning : imaginary emittance at iel,jc=',
     $            IEL,JJ
            WRITE(6,*) '        Calculation skipped...'
          ELSE
            SQ=SQRT(SQ2)
            EsPi(JJ,IEL)= SQ/2.D0
            IF(SQ.EQ.0.D0) THEN
              WRITE(6,FMT='(I2,A)') JJ,'  ZERO EMITTANCE'
              WRITE(6,*) '      Cannot compute beta functions...'
            ELSE
              BETA(JJ,IEL) =2.D0* SC2(JC,IEL)  /SQ
              GAMMA(JJ,IEL)=2.D0* SC2(JC+1,IEL)/SQ
              ALPHA(JJ,IEL)=2.D0* SC2(JJ+6,IEL)/SQ
              WRITE(6,FMT='(A,I4,A,I3,A,1P,2G11.3,2G12.4)') ' Lmnt # ',
     >         JEL(IEL),'  ',JJ,'  Beta, Alpha, Eps/pi = ',BETA(JJ,IEL),
     >              ALPHA(JJ,IEL),EsPi(JJ,IEL),EsPi(JJ,IEL)-EsPi(JJ,1)
C     >                ALPHA(JJ,IEL),EsPi(JJ,IEL),EsPi(JJ,IEL)-EsPi0(JJ)
            ENDIF
C---------- Emtitance growth
C            DEM(JJ,IEL)=EsPi(JJ,IEL)-EsPi0(JJ)
            DEM(JJ,IEL)=EsPi(JJ,IEL)-EsPi(JJ,1)
            IF(DEM(JJ,IEL).LT.0.D0) THEN
              WRITE(6,FMT='(A,I4,2X,I3,A)') ' Lmnt # ',JEL(IEL),
     >         JJ,' negative growth -> Will be forced to zero...'
              DEM(JJ,IEL)=0.D0
            ENDIF
            IF(DEM(JJ,IEL).GT.DEMA(JJ)) DEMA(JJ)=DEM(JJ,IEL)
            IF(DEM(JJ,IEL).LT.DEMI(JJ)) DEMI(JJ)=DEM(JJ,IEL)
          ENDIF
 29     CONTINUE
         
        WRITE(78,FMT='(1P,G10.2,I4,6G11.3,2G12.4,G11.3)')
C Col #     1          2         3         4       5  
C-------    s                   <X>       <Y>     <dp>
     >    SCUM(IEL),JEL(IEL),SC(1,IEL),SC(3,IEL),SC(6,IEL), 
C Col #                   6  7                              8
C--------             env X, Z                         sigma_{dE/E}
     >   (SQRT(EsPi(JJJ,IEL)*BETA(JJJ,IEL)),JJJ=1,2),SQRT(SC2(6,IEL)),
C Col #           9 10             11
C-------         Ex Ez            dEx
     >   (EsPi(JJJ,IEL),JJJ=1,2),DEM(1,IEL)
 28   CONTINUE

 70   CONTINUE
      TEMP =DEMA(1)
      IF(DEMA(2).GT.TEMP) TEMP=DEMA(2)
      IF(DEMA(3).GT.TEMP) TEMP=DEMA(3)
      IF(.NOT.OKECH) THEN
        CALL CLSCR
        CALL TRAXES(ZERO,SCUM(NEL),ZERO,TEMP,2)
        OKECH=.TRUE.
      ENDIF
      IL=-1
      DO 34 JJ=1,2
      IL=IL+2
      CALL LINTYP(IL)
        CALL VECTPL(ZERO,ZERO,4)
        DO 34 IEL=1,NEL
          CALL VECTPL(SCUM(IEL),DEM(JJ,IEL),2)
 34     CONTINUE
      CALL FBGTXT
      WRITE(6,*)
      WRITE(6,FMT='(2(I5,A))') npass,' pass',nel,' elements'
      WRITE(6,FMT='(1p,2g12.4)') SCUM(NEL),TEMP

 99   CONTINUE
      CLOSE(NL)
      CLOSE(88)
      RETURN
      END
