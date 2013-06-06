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
C  -------
      SUBROUTINE OBJ3(KOBJ2,BORO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     **************************************
C     CONSTITUTION DE L'OBJET INITIAL KOBJ=3
C     **************************************
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      COMMON/CHAMBR/ LIMIT,IFORM,YLIM2,ZLIM2,SORT(MXT),FMAG,BMAX
     > ,YCH,ZCH
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IIP(MXL),NB,NOEL
      CHARACTER*80 TA
      PARAMETER (MXTA=45)
      COMMON/DONT/ TA(MXL,MXTA)
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      CHARACTER LET
      COMMON/FAISCT/ LET(MXT)
      CHARACTER  KAR(41)
      COMMON/KAR/ KAR
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      COMMON/PTICUL/ AAM,Q,G,TO
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      COMMON/SYNCH/ RET(MXT), DPR(MXT),PS

      CHARACTER*80 NOMFIC
      INTEGER DEBSTR,FINSTR
      LOGICAL IDLUNI
      LOGICAL BINARI, BINARY
      CHARACTER*11 FRMT
      LOGICAL OKKT, OKKP
      CHARACTER TDUMX*10
      CHARACTER TDUM8*8

      CHARACTER LETI, LETAG
      CHARACTER*130 TXT

C----- Reset particle counter
      IF(IPASS.EQ.1) CALL CNTRST

      KT1 = NINT(A(NOEL,20)      )
      KT2 = NINT(A(NOEL,21)      )
      KTSTP = A(NOEL,22)
      KP1 = NINT(A(NOEL,30))
      KP2 = NINT(A(NOEL,31)) 
      KP3 = NINT(A(NOEL,32))
      YFAC = A(NOEL,40)
      TFAC = A(NOEL,41)
      ZFAC = A(NOEL,42)
      PFAC = A(NOEL,43)
      SFAC = A(NOEL,44)
      DPFAC= A(NOEL,45)
      TIFAC= A(NOEL,46)
      LETAG=TA(NOEL,1)(1:1)
      YREF = A(NOEL,50)
      TREF = A(NOEL,51)
      ZREF = A(NOEL,52)
      PREF = A(NOEL,53)
      SREF = A(NOEL,54)
      DREF= A(NOEL,55)
      TIREF= A(NOEL,56)
      IF(KOBJ2.EQ.0) DREF=1.D0
      INITC = NINT(A(NOEL,60))
      NOMFIC=TA(NOEL,2)
      

      NOMFIC=NOMFIC(DEBSTR(NOMFIC):FINSTR(NOMFIC))

      BINARY=BINARI(NOMFIC,
     >                     IDUM)
      FRMT='FORMATTED'
      IF(BINARY) FRMT='UNFORMATTED'
      IF(IPASS .EQ. 1) THEN
        IF(IDLUNI(
     >            NL)) THEN
          IF(NRES.GT.0) WRITE(NRES,FMT='(/,''   Opening input file  '',
     >    A,''  Unit # : '',I3)') 
     >    NOMFIC(DEBSTR(NOMFIC):FINSTR(NOMFIC)),NL
          OPEN(UNIT=NL,FILE=NOMFIC,STATUS='OLD',FORM=FRMT,ERR=96, 
     >                                          IOSTAT=IOS)
            IF(IOS.NE.0) GOTO 96
            IF(NRES.GT.0) WRITE(NRES,FMT='(/,''   Reading  initial'',2X
     >      ,''conditions  in  '',A)') NOMFIC
            CALL HEADER(NL,NRES,4,BINARY,
     >                                   *999)
        ELSE
          GOTO 96
        ENDIF
      ENDIF

C----- Flag for lmnt number - not used (for the moment...)
      LM = -1
C----- Traj. counter
      I = 0
      IT1 = 0
      IKAR = 0
      IPASSR = 0
      IEND = 0

C      AMQLU(3) = (BINARY.OR.KOBJ2.EQ.0)
      AMQLU(3) = (KOBJ2.EQ.0)
      AMQLU(4) = AMQLU(3)
      AMQLU(5) = AMQLU(3)
      AMQLU(1) = (AMQLU(3).OR.KOBJ2.EQ.3)
      AMQLU(2) = AMQLU(1)
      PABSLU = KOBJ2.EQ.2

  17  CONTINUE
        IPASS1 = IPASSR
        I = I + 1

C        IF(BINARY) THEN
        IF ((BINARY) .AND. (KOBJ2.EQ.0)) THEN
 222      CONTINUE
          IF(IEND.EQ.1) GOTO 95
          READ(NL,ERR=97,END=95)
     >      IEXI,DPO,YO,TTO,ZO,PO,SO,TIMO, 
     >      DP,Y,T,Z,P,S,TIM,
     >      SIX,SIY,SIZ,SIN,SFX,SFY,SFZ,SFN,
     >      EKIN,ENERG, 
     >      ITR,IREPI,SORTI,AMQ1,AMQ2,AMQ3,AMQ4,AMQ5,RETI,DPRI,PS,
     >      BRO, IPASSR, NOELR, TDUMX,TDUM8,TDUM8,LETI

          IF(LM .NE. -1) THEN
            IF(LM .NE. NOELR) GOTO 222
          ENDIF
          
          IF(.NOT. OKKP(KP1,KP2,KP3,IPASSR,
     >                                 IEND)) THEN
            IF(IEND.EQ.1) THEN 
              IPASSR=IPASS1
              GOTO 95
            ENDIF
            GOTO 222
          ENDIF

          IF(.NOT. OKKT(KT1,KT2,ITR,
     >                              IEND)) THEN
            GOTO 222
          ENDIF

C        ELSEIF(.NOT.BINARY) THEN
        ELSE

 221      CONTINUE
          IF(IEND.EQ.1) GOTO 95

          IF  (KOBJ2.EQ.0) THEN           
            READ(NL,110,ERR=97,END=95)
     >      IEXI,DPO,YO,TTO,ZO,PO,SO,TIMO, 
     >      DP,Y,T,Z,P,S,TIM,
     >      SIX,SIY,SIZ,SIN,SFX,SFY,SFZ,SFN,
     >      EKIN,ENERG, 
     >      ITR,IREPI,SORTI,AMQ1,AMQ2,AMQ3,AMQ4,AMQ5,RETI,DPRI,PS,
     >      BRO, IPASSR, NOELR, TDUMX,TDUM8,TDUM8,LETI
            INCLUDE "FRMFAI.H"

          ELSEIF(KOBJ2.EQ.1) THEN 
C------------ Was installed for reading pion data at NuFact target
            IKAR = IKAR+1
            IF(IKAR.GT.41)  IKAR=1
C            READ(NL,*,ERR=97,END=95) Y,T,Z,P,S, DP
            IF (BINARY) THEN
              READ(NL,ERR=97,END=95) Y,T,Z,P,S, DP
            ELSE
              READ(NL,*,ERR=97,END=95) Y,T,Z,P,S, DP
            ENDIF
            TIM = 0.D0
            LETI=KAR(IKAR)
            IEXI=1
            IT = IT1
            IREPI = IT
            IPASSR =  KP1    
            BRO = BORO
            YO= 0.D0
            TTO= 0.D0
            ZO= 0.D0
            PO= 0.D0
            SO= 0.D0
            DPO= 0.D0
            TIMO= 0.D0

          ELSEIF(KOBJ2.EQ.2) THEN 
C----------- Was installed for reading e+ data provided by Rosowski/Perez. DAPNIA/SPP March 03
            IKAR = IKAR+1
            IF(IKAR.GT.41)  IKAR=1
            READ(NL,*,ERR=97,END=95) X,Y,Z,PX,PY,PZ
            BRO = BORO
            PT = SQRT(PX*PX+PY*PY+PZ*PZ)
            DP = PT/(BORO*CL9)
            T = ATAN(PY/PX)*1000.D0
            P = ATAN(PZ/SQRT(PX*PX+PY*PY))*1000.D0
            TIM = 0.D0
            LETI=KAR(IKAR)
            IEXI=1
            IT = IT1
            IREPI = IT
            IPASSR =  KP1    
            YO= 0.D0
            TTO= 0.D0
            ZO= 0.D0
            PO= 0.D0
            SO= 0.D0
            DPO= 0.D0
            TIMO= 0.D0

          ELSEIF(KOBJ2.EQ.3) THEN 
C------------ Was installed for RHS_DESIR
            IKAR = IKAR+1
            IF(IKAR.GT.41)  IKAR=1
            READ(NL,*,ERR=97,END=95) dp,Y,T,Z,P,S,time,amq1, amq2
C            TIM = 0.D0
            LETI=KAR(IKAR)
            IEXI=1
            IT = IT1
            IREPI = IT
            IPASSR =  KP1    
            BRO = BORO
            YO= 0.D0
            TTO= 0.D0
            ZO= 0.D0
            PO= 0.D0
            SO= 0.D0
            DPO= 0.D0
            TIMO= 0.D0

          ENDIF 

          IF(LM .NE. -1) THEN
            IF(LM .NE. NOELR) GOTO 221
          ENDIF
          IF(.NOT. OKKP(KP1,KP2,KP3,IPASSR,
     >                                 IEND)) THEN
            IF(IEND.EQ.1) THEN 
              IPASSR=IPASS1
              GOTO 95
            ENDIF
            GOTO 221
          ENDIF

          IF(.NOT. OKKT(KT1,KT2,ITR,
     >                              IEND)) THEN
            GOTO 221 
          ENDIF

        ENDIF

        IF(LETAG.NE.'*') THEN
          IF(LETI.NE.LETAG) IEXI=-9
        ENDIF

        IF(IEXI.LE.0) GOTO 17

        IT1 = IT1 + 1

        LET(IT1)=LETI
        IEX(IT1)=IEXI
C        FO(1,IT1)=1.D0 + DPO
        FO(1,IT1)=(1.D0 + DPO) * BRO/BORO
C        FO(1,IT1)= DPO * BRO/BORO
        FO(2,IT1)=YO
        FO(3,IT1)=TTO
        FO(4,IT1)=ZO
        FO(5,IT1)=PO
        FO(6,IT1)=SO
        FO(7,IT1)=TIMO
        IF (PABSLU) THEN
           DP0(IT1)=DP*DPFAC * BRO/BORO
C          THIS WILL NOT WORK RIGHT IF DREF.NE.0D0 AND Q IS CHANGED IN A 
C          SUBSEQUENT PARTICUL
           F(1,IT1)= DP0(IT1)/Q
        ELSE
           F(1,IT1)=(DP*DPFAC + DREF) * BRO/BORO
        ENDIF
        F(2,IT1)=  Y*YFAC  + YREF
        F(3,IT1)=  T*TFAC  + TREF
        F(4,IT1)=  Z*ZFAC  + ZREF
        F(5,IT1)=  P*PFAC  + PREF
        F(6,IT1)=  S*SFAC  + SREF
        F(7,IT1)=TIM*TIFAC + TIREF 
        IREP(IT1) = IT1
        IF (AMQLU(3)) THEN
           RET(IT1)=RETI
           DPR(IT1)=DPRI
C        IREP(IT1) = IREPI
           SORT(IT1) = SORTI
           AMQ(3,IT1) = AMQ3
           AMQ(4,IT1) = AMQ4
           AMQ(5,IT1) = AMQ5
        ENDIF
        IF (AMQLU(1)) THEN
           AMQ(1,IT1) = AMQ1
           AMQ(2,IT1) = AMQ2
        ELSE
           AMQ(1,IT1) = AAM
           AMQ(2,IT1) = Q
        ENDIF

        IF    (INITC .EQ. 1) THEN
          DO J=1,MXJ
            FO(J,IT1) = F(J,IT1)
          ENDDO
        ELSEIF(INITC .EQ. 2) THEN
          DO J=1,MXJ
            TEMP = F(J,IT1)
            F(J,IT1) = FO(J,IT1)
            FO(J,IT1) = TEMP
          ENDDO        
        ENDIF

        IF(IT1 .EQ. MXT) GOTO 169
        IF(IT1 .EQ. KT2) GOTO 169
      GOTO  17
 
 96   WRITE(TXT,FMT='(
     >'' *** SBR OBJ3 : ERROR  at  OPEN  file  "'',A,''"'')') 
     >NOMFIC(DEBSTR(NOMFIC):FINSTR(NOMFIC))
      CALL ENDJOB(TXT,-99)
 
 97   WRITE(NRES,FMT='(/,
     > '' SBR OBJ3 -> error in  reading  file '',
     >  A,'' at  event/traj #  '',I6,''/'')') NOMFIC,IT1

 95   CONTINUE
      IT1 = IT1-1
 169  CONTINUE

C      IF(KT1.GT.1) THEN
C        IMAX=IT1-1
C      ELSE
        IMAX=IT1
C      ENDIF
C-----
      CALL CNTMXT(IMAX)

      II = 0
      IS = 0
      DO 1 IT=1,IMAX
        IF(IEX(IT).LE.-1) THEN
          II = II + 1
          CALL KSTOP(ABS(IEX(IT)),IT,IEX(I),*1)
        ENDIF
        IF(LET(IT).EQ.'S') IS = IS+1
 1    CONTINUE

      IF(NRES .GT. 0) THEN
        WRITE(NRES,FMT='(/,T5,''  Reading  in  file  '',A
     >    ,/,''   ended  after  gathering '',I6,''  particles''
     >    ,''  in  requested  range :  ['',I6,'', '',I6,'']'')')
     >    NOMFIC(DEBSTR(NOMFIC):FINSTR(NOMFIC)),IMAX,KT1,KT2
        IF(IS.GT.0) WRITE(NRES,FMT='(/,T5,I6,
     >    ''  particles  are  of  secondary  type  (LET="S")'')') IS
        IF(II.GT.0) THEN 
          WRITE(NRES,
     >    FMT='(/,T5,I6,''/'',I6,'' particles  have  IEX < 0,  hence'',
     >    I6,''  only  left  to  be  ray-traced. '')') II,IMAX,IMAX-II
          WRITE(NRES,*) ' Last  pass  number  read :  ',IPASSR,'  in  ' 
     >    ,'requested range :  [',KP1,',',KP2,'], ipass-modulo=',KP3,'.'
          IF    (INITC.LT.2)  THEN
            WRITE(NRES,FMT='(/,T5,''InitC = '',I1
     >      ,'' :  The  starting  coordinates  of  this  run  are  ''
     >      ,''the  final  coordinates  read  in '',A)')
     >      INITC,NOMFIC(DEBSTR(NOMFIC):FINSTR(NOMFIC))
          ELSEIF(INITC.EQ.2)  THEN
            WRITE(NRES,FMT='(/,T5,''InitC = '',I1
     >      ,'' :  The  starting  coordinates  of  this  run  are  ''
     >      ,''the  initial  coordinates  read  in '',A)')
     >      INITC,NOMFIC(DEBSTR(NOMFIC):FINSTR(NOMFIC))
          ELSE
            STOP ' SBR obj3 : No such option InitC '
          ENDIF
        ENDIF
      ENDIF
      IDMAX = 1
      IMAXT=IMAX
 
      CLOSE(NL)
      RETURN
 999  CONTINUE
          CALL OBJERR(ABS(NRES),1,MXT,' Read error')
          CALL ENDJOB('*** Error, SBR OBJ3 -> Read error',-99)  
      RETURN
      END
