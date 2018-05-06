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
C  Upton, NY, 11973, USA
C  -------
      SUBROUTINE OBJ3(KOBJ2,BORO,
     >                           DPREF,HDPRF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     **************************************
C     CONSTITUTION DE L'OBJET INITIAL KOBJ=3
C     **************************************
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      INCLUDE "C.CHAMBR.H"     ! COMMON/CHAMBR/ LIMIT,IFORM,YLIM2,ZLIM2,SORT(MXT),FMAG,YCH,ZCH
 
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON_2.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IIP(MXL),NB,NOEL
C      PARAMETER (LNTA=132) ; CHARACTER(LNTA) TA
C      PARAMETER (MXTA=45)
      INCLUDE "C.DONT.H"     ! COMMON/DONT/ TA(MXL,MXTA)
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      CHARACTER(1) LET
      INCLUDE "C.FAISCT.H"     ! COMMON/FAISCT/ LET(MXT)
      CHARACTER(1) KAR(41)
      INCLUDE "C.KAR.H"     ! COMMON/KAR/ KAR
      INCLUDE "C.OBJET.H"     ! COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT,KZOB
      INCLUDE "C.PTICUL_2.H"     ! COMMON/PTICUL/ AAM,Q,G,TO
      INCLUDE "C.REBELO.H"   ! COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      INCLUDE "C.SPIN.H"     ! COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)
      INCLUDE "C.SYNCH.H"     ! COMMON/SYNCH/ PH(MXT), DPR(MXT), PS

      CHARACTER(LNTA) NOMFIC
      INTEGER DEBSTR,FINSTR
      LOGICAL IDLUNI
      LOGICAL BINARI, BINARY
      CHARACTER(11) FRMT
      LOGICAL OKKT3, OKKP3
      PARAMETER (KSIZ=10)
      CHARACTER(KSIZ) TDUMK
      PARAMETER (LBLSIZ=20)
      CHARACTER(LBLSIZ) TDUML

      CHARACTER(1) LETI, LETAG
      CHARACTER(130) TXT
C      CHARACTER(1) TX1
      LOGICAL OKOPN
      LOGICAL AFTREB

      SAVE NHD
      LOGICAL FIRST
 
      DATA OKOPN / .FALSE. /
      DATA AFTREB / .FALSE. /
      DATA FIRST / .TRUE. /

C----- Reset particle counter
      IF(IPASS.EQ.1) CALL CNTRST

      KT1 = NINT(A(NOEL,20)      )
      KT2 = NINT(A(NOEL,21)      )
      KT3 = NINT(A(NOEL,22))
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
        INQUIRE(FILE=NOMFIC,NUMBER=NL,OPENED=OKOPN)
        IF(.NOT. OKOPN) THEN      
          OKOPN = (IDLUNI(
     >                    NL)) 
          IF(.NOT. OKOPN) CALL ENDJOB
     >       ('Pgm obj3. Cannot get free unit for object file',-99)
          IF(NRES.GT.0) WRITE(NRES,FMT='(/,''   Opening input file  '',
     >    A,'',  in logical unit # : '',I3)') 
     >    NOMFIC(DEBSTR(NOMFIC):FINSTR(NOMFIC)),NL
          OPEN(UNIT=NL,FILE=NOMFIC,STATUS='OLD',FORM=FRMT,ERR=96, 
     >                                          IOSTAT=IOS)
            IF(IOS.NE.0) GOTO 96
            IF(NRES.GT.0) WRITE(NRES,FMT='(/,''   Reading  initial'',2X
     >      ,''conditions  in  file  '',A)') 
     >                             NOMFIC(DEBSTR(NOMFIC):FINSTR(NOMFIC))
            IF(NRES.GT.0) WRITE(NRES,FMT='(/,''   Particles  '',2X
     >      ,I9,''  to  '',I9,'',  step  '',I9)') KT1,KT2,KT3
            IF(NRES.GT.0) WRITE(NRES,FMT='(/,''   Pass #     '',2X
     >      ,I9,''  to  '',I9,'',  step  '',I9,/)') KP1,KP2,KP3
        ELSE
          REWIND(NL)
        ENDIF
C        IF(NINT(A(NOEL,11)).GE.1) THEN
Changed Jan 2018 - FM. Now managed in robjet.
        NHD = NINT(A(NOEL,12)) 
C        ELSE
C          NHD = 4
C        ENDIF
        CALL HEADER(NL,NRES,NHD,BINARY,
     >                                 IER)
C     >                              *999)
        IF(IER.GT.0) GOTO 999

      ELSE
        CALL REBEL5(
     >              AFTREB)
        IF(AFTREB) THEN
C May happen if problems are stacked, as in 
C [...]/testKEYWORDS/FIT/FITfollowedByREBELOTE/orbitsInCYCLOTRON/REBELOTE-FIT_followedByOrbits.dat
          INQUIRE(FILE=NOMFIC,NUMBER=NL,OPENED=OKOPN)
          IF(OKOPN) THEN
            REWIND(NL)
            CALL HEADER(NL,NRES,NHD,BINARY,
     >                                     IER)
C     >                                   *999)
            IF(IER.GT.0) GOTO 999
          ENDIF
        ENDIF
      ENDIF

C----- Flag for lmnt number - not used (for the moment...)
      LM = -1
C----- Traj. counter
      I = 0
      IT1 = 0
      IT2 = 0
      IKAR = 0
      IPASSR = 0
      IEND = 0

      AMQLU(3) = (BINARY.OR.KOBJ2.EQ.0)
      AMQLU(4) = AMQLU(3)
      AMQLU(5) = AMQLU(3)
      AMQLU(1) = (AMQLU(3).OR.KOBJ2.EQ.3)
      AMQLU(2) = AMQLU(1)
      PABSLU = KOBJ2.EQ.2

  17  CONTINUE
        IPASS1 = IPASSR
        I = I + 1

        IF ((BINARY) .AND. (KOBJ2.EQ.0)) THEN
C Expected to be resuming from 1 particular pass (for instance after a 
C crash), reading a zgoubi.fai style file (i.e., a file created by FAISTORE (or FAISCNL);
C KP1=KP2 is expected iin that case.

 222      CONTINUE

          IF(IEND.EQ.1) GOTO 95
          READ(NL,ERR=97,END=95)
     >      IEXI,DPO,YO,TTO,ZO,PO,SO,TIMO, 
     >      DP,Y,T,Z,P,S,TIM,
     >      SIX,SIY,SIZ,SIN,SFX,SFY,SFZ,SFN,
     >      EKIN,ENERG, 
     >      IT,IREPI,SORTI,AMQ1,AMQ2,AMQ3,AMQ4,AMQ5,PHI,DPRI,PS,
     >      BRO, IPASSR, NOELR, TDUMK,TDUML,TDUML,LETI, SRLO ,
     >      DPRF, GDPRF

            CALL SPN2(.TRUE.)

          IF(LM .NE. -1) THEN
            IF(LM .NE. NOELR) GOTO 222
          ENDIF
          
          IF(.NOT. OKKP3(KP1,KP2,KP3,IPASSR,
     >                                   IEND)) THEN
            IF(IEND.EQ.1) THEN 
              IPASSR=IPASS1
              GOTO 95
            ENDIF

            GOTO 222
          ENDIF
          
          IF(.NOT. OKKT3(KT1,KT2,KT3,IT,
     >                                  IEND)) GOTO 222

C          IF(KP1.EQ.KP2) IPASS=IPASSR
            DPREF = DPRF
            HDPRF = GDPRF
          IT2 = IT2 + 1
          IT1 = IT1 + 1

        ELSE       ! .NOT.BINARY) THEN

 221      CONTINUE
          IF(IEND.EQ.1) GOTO 95

          IF  (KOBJ2.EQ.0) THEN           
C Recommended, in case of resuming from 1 particular pass (for instance after a 
C crash), reading a zgoubi.fai style file (i.e., a file created by FAISTORE (or FAISCNL);
C KP1=KP2 is expected iin that case.
             
            READ(NL,*,ERR=97,END=95)
     >      IEXI,DPO,YO,TTO,ZO,PO,SO,TIMO, 
     >      DP,Y,T,Z,P,S,TIM, 
     >      SIX,SIY,SIZ,SIN,SFX,SFY,SFZ,SFN,
     >      EKIN,ENERG, 
     >      IT,IREPI,SORTI,AMQ1,AMQ2,AMQ3,AMQ4,AMQ5,PHI,DPRI,PS,
     >      BRO, IPASSR, NOELR,     TDUMK,
     >                              TDUML,
     >                              TDUML,     LETI, SRLO ,
     >      DPRF, GDPRF
            CALL SPN2(.TRUE.)

            IF(LM .NE. -1) THEN
              IF(LM .NE. NOELR) GOTO 221
            ENDIF

            IF(.NOT. OKKP3(KP1,KP2,KP3,IPASSR,
     >                                       IEND)) THEN
              IF(IEND.EQ.1) THEN 
                IPASSR=IPASS1
                GOTO 95
              ENDIF

              GOTO 221
            ENDIF

            IF(.NOT. OKKT3(KT1,KT2,KT3,IT,
     >                                   IEND)) GOTO 221

C            IF(KP1.EQ.KP2) IPASS=IPASSR
            DPREF = DPRF
            HDPRF = GDPRF                   
            IT2 = IT2 + 1
            IT1 = IT1 + 1

          ELSEIF(KOBJ2.EQ.1) THEN 
C------------ Was installed for reading pion data at NuFact target
            IKAR = IKAR+1
            IF(IKAR.GT.41)  IKAR=1
 171        CONTINUE
C            READ(NL,*,ERR=97,END=95) Y,T,Z,P,S, DP
            IF (BINARY) THEN
              READ(NL,ERR=97,END=95) Y,T,Z,P,S, DP
            ELSE
              READ(NL,*,ERR=97,END=95) Y,T,Z,P,S, DP
            ENDIF
            IT2 = IT2 + 1
            IF(KT3*((IT2-1)/KT3) .NE. IT2-1) GOTO 171            
            IT1 = IT1 + 1
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
            it1 = it1 + 1
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
            IT2 = IT2 + 1

          ELSEIF(KOBJ2.EQ.3) THEN 
C------------ Was installed for RHS_DESIR
            IKAR = IKAR+1
            IF(IKAR.GT.41)  IKAR=1
            READ(NL,*,ERR=97,END=95) DP,Y,T,Z,P,S,time,amq1, amq2
            it1 = it1 + 1
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
            IT2 = IT2 + 1

          ELSEIF(KOBJ2.EQ.4) THEN 
C------------ Installed for RHIC FFAG arcs
            IKAR = IKAR+1
            IF(IKAR.GT.41)  IKAR=1
 174        continue
            READ(NL,*,ERR=97,END=95) Y,T,Z,P,S, DP, leti, ien, jt2
            IT2 = JT2
            IF(IT2 .LT. KT1 .OR. IT2 .GT. KT2) GOTO 174            
            IF(IT1.GT.0) THEN 
              IF(KT3*((IT2)/KT3) .NE. IT2) GOTO 174            
            ENDIF
            IT1 = IT1 + 1
            TIM = 0.D0
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

          ELSEIF(KOBJ2.EQ.5) THEN 
C------------ Installed to read Yichao's files
            IKAR = IKAR+1
            IF(IKAR.GT.41)  IKAR=1
 175        CONTINUE
            READ(NL,*,ERR=97,END=95) Y,T,Z,P,TIM, GMA
            IF(FIRST) THEN
              AME = 0.51099892
              Q = 1.D0
              PRF = BORO*CL9*Q
            ENDIF
            PIT1 =  AME * SQRT(GMA*GMA -1.D0)
            DP = PIT1/PRF
            BTA = PIT1/SQRT(PIT1*PIT1 + AME * AME)
            S = BTA * CL * TIM * 1D2  ! cm
            TIM = TIM * 1.D6  ! mu_s
            IT2 = IT2 + 1
            IF(KT3*((IT2-1)/KT3) .NE. IT2-1) GOTO 175
            IT1 = IT1 + 1
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
          IF(.NOT. OKKP3(KP1,KP2,KP3,IPASSR,
     >                                     IEND)) THEN
            IF(IEND.EQ.1) THEN 
              IPASSR=IPASS1
              GOTO 95
            ENDIF
            GOTO 221
          ENDIF

          IF(.NOT. OKKT3(KT1,KT2,KT3,IT,
     >                                 IEND)) THEN
            
            IT1 = IT1 - 1 
            GOTO 169
          ENDIF

        ENDIF

        IF(LETAG.NE.'*') THEN
          IF(LETI.NE.LETAG) IEXI=-9
        ENDIF

        IF(IEXI.LE.0) THEN
C So to avoid problem with test on D in objets
          FO(1,IT1)= 1.D0 
          F(1,IT1) = 1.D0
C To be clean, otherwise LET(i) is undefined
          LET(IT1) = LETI
          IEX(IT1)=IEXI

          GOTO 17

        ENDIF

        LET(IT1)=LETI
        IEX(IT1)=IEXI
        FO(1,IT1)=(1.D0 + DPO) * BRO/BORO
        FO(2,IT1)=YO
        FO(3,IT1)=TTO
        FO(4,IT1)=ZO
        FO(5,IT1)=PO
        FO(6,IT1)=SO
        FO(7,IT1)=TIMO
        SI(1,IT1) = SIX
        SI(2,IT1) = SIY
        SI(3,IT1) = SIZ
        SI(4,IT1) = SIN

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
        SF(1,IT1) = SFX
        SF(2,IT1) = SFY
        SF(3,IT1) = SFZ
        SF(4,IT1) = SFN
        IREP(IT1) = IT1
        IF (AMQLU(3)) THEN
           PH(IT1)=PHI
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
C      IT1 = IT1-1

 169  CONTINUE

      IMAX=IT1

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
        WRITE(NRES,FMT='(/,T5,''  Reading in file '',A
     >    ,''has ended after gathering '',I7,'' particles''
     >    ,'' in requested range :  ['',I7,'', '',I7,'']'')')
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
 
C      CLOSE(NL)
      RETURN
 999  CONTINUE
          CALL OBJERR(ABS(NRES),1,MXT,' Read error')
          CALL ENDJOB('*** Error, SBR OBJ3 -> Read error',-99)  
      RETURN
      END
