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
      COMMON/DONT/ TA(MXL,20)
      INCLUDE "MAXCOO.H"
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),IMAX,IEX(MXT),IREP(MXT)
      CHARACTER LET
      COMMON/FAISCT/ LET(MXT)
      CHARACTER  KAR(41)
      COMMON/KAR/ KAR
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      COMMON/SYNCH/ RET(MXT), DPR(MXT),PS

      CHARACTER*80 NOMFIC
      INTEGER DEBSTR,FINSTR
      LOGICAL OKINIT, IDLUNI
      LOGICAL BINARI, BINARY
      CHARACTER*11 FRMT
      LOGICAL OKKT, OKKP
      CHARACTER TDUM8*8

      CHARACTER LETI, LETAG
      CHARACTER*130 TXT

      Q = 1.D0
      P0 = BORO*CL9*Q

C----- Reset particle counter
      IF(IPASS.EQ.1) CALL CNTRST

      KT1 = NINT(A(NOEL,20)      )
      KT2 = NINT(A(NOEL,21)      )
      KTSTP = A(NOEL,22)
      KP1 = NINT(A(NOEL,30))
      KP2 = NINT(A(NOEL,31)) 
      KPSTP = NINT(A(NOEL,32))
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
      DPREF= A(NOEL,55)
      TIREF= A(NOEL,56)
      IF(KOBJ2.EQ.0) DPREF=1.D0
      OKINIT = NINT(A(NOEL,60)) .EQ. 0
      NOMFIC=TA(NOEL,2)
      

      NOMFIC=NOMFIC(DEBSTR(NOMFIC):FINSTR(NOMFIC))
      IF(NRES.GT.0) WRITE(NRES,FMT='(/,''   Reading  initial'',2X
     >,''conditions  in  '',A)') NOMFIC

      BINARY=BINARI(NOMFIC,
     >                     IDUM)
      FRMT='FORMATTED'
      IF(BINARY) FRMT='UNFORMATTED'
      IF(IPASS .EQ. 1) THEN
        IF(IDLUNI(
     >            NL)) THEN
          OPEN(UNIT=NL,FILE=NOMFIC,STATUS='OLD',FORM=FRMT,ERR=96, 
     >                                          IOSTAT=IOS)
            IF(IOS.NE.0) GOTO 97
            CALL HEADER(NL,NRES,4,BINARY,
     >                                   *999)
        ELSE
          GOTO 96
        ENDIF
      ENDIF

C----- Flag for lmnt number - not used
      LM = -1
C----- Traj. counter
      I = 0
      IT1 = 0
      IKAR = 0
      IPASSR = 0
      IEND = 0

  17  CONTINUE
        IPASS1 = IPASSR
        I = I + 1
        IT1 = IT1 + 1

        IF(BINARY) THEN
 222      CONTINUE
          IF(IEND.EQ.1) GOTO 95
          READ(NL,ERR=97,END=95)
     >     LETI,IEXI,DPO,YO,TO,ZO,PO,SO,TIMO, 
     >     DP,Y,T,Z,P,S,TIM,EKIN, 
     >     IT,IREPI,SORTI,AMQ1,AMQ2,AMQ3,AMQ4,AMQ5,RETI,DPRI,
     >     BRO, IPASSR, TDUM8,TDUM8,TDUM8,IDUM

          IF(LM .NE. -1) THEN
            IF(LM .NE. NOELR) GOTO 222
          ENDIF
          
          IF(.NOT. OKKP(KP1,KP2,IPASSR,
     >                                 IEND)) THEN
            IF(IEND.EQ.1) THEN 
              IPASSR=IPASS1
              IT1 = IT1-1
              GOTO 95
            ENDIF
            GOTO 222
          ENDIF
          IF(.NOT. OKKT(KT1,KT2,IT1,
     >                             IEND)) GOTO 222
        ELSE
 221      CONTINUE
          IF(IEND.EQ.1) GOTO 95

          IF  (KOBJ2.EQ.0) THEN           
            READ(NL,110,ERR=97,END=95)
     >      LETI,IEXI,DPO,YO,TO,ZO,PO,SO,TIMO, 
     >      DP,Y,T,Z,P,S,TIM,EKIN, 
     >      IT,IREPI,SORTI,AMQ1,AMQ2,AMQ3,AMQ4,AMQ5,RETI,DPRI,
     >      BRO, IPASSR, TDUM8,TDUM8,TDUM8,IDUM
            INCLUDE "FRMFAI.H"

          ELSEIF(KOBJ2.EQ.1) THEN 
C------------ Was installed for reading pion data at NuFact target
            IKAR = IKAR+1
            IF(IKAR.GT.41)  IKAR=1
            READ(NL,*,ERR=97,END=95) Y,T,Z,P,S, DP
            TIM = 0.D0
            LETI=KAR(IKAR)
            IEXI=1
            IT = IT1
            IREPI = IT
            IPASSR =  KP1    
            BRO = BORO
            YO= 0.D0
            TO= 0.D0
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
            DP = PT/P0
            T = ATAN(PY/PX)*1000.D0
            P = ATAN(PZ/SQRT(PX*PX+PY*PY))*1000.D0
            TIM = 0.D0
            LETI=KAR(IKAR)
            IEXI=1
            IT = IT1
            IREPI = IT
            IPASSR =  KP1    
            YO= 0.D0
            TO= 0.D0
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
            TO= 0.D0
            ZO= 0.D0
            PO= 0.D0
            SO= 0.D0
            DPO= 0.D0
            TIMO= 0.D0

          ENDIF 

          IF(LM .NE. -1) THEN
            IF(LM .NE. NOEL) GOTO 221
          ENDIF
          IF(.NOT. OKKP(KP1,KP2,IPASSR,
     >                                 IEND)) THEN
            IF(IEND.EQ.1) THEN 
              IPASSR=IPASS1
              IT1 = IT1-1
              GOTO 95
            ENDIF
            GOTO 221
          ENDIF
          IF(.NOT. OKKT(KT1,KT2,IT1,
     >                             IEND)) GOTO 221
        ENDIF

        IF(LETAG.NE.'*') THEN
          IF(LETI.NE.LETAG) IEXI=-9
        ENDIF

        LET(IT1)=LETI
        IEX(IT1)=IEXI
C        FO(1,IT1)=1.D0 + DPO
C        FO(1,IT1)=(1.D0 + DPO) * BRO/BORO
        FO(1,IT1)= DPO * BRO/BORO
        FO(2,IT1)=YO
        FO(3,IT1)=TO
        FO(4,IT1)=ZO
        FO(5,IT1)=PO
        FO(6,IT1)=SO
        FO(7,IT1)=TIMO
        F(1,IT1)=(DP*DPFAC + DPREF) * BRO/BORO
        F(2,IT1)=  Y*YFAC  + YREF
        F(3,IT1)=  T*TFAC  + TREF
        F(4,IT1)=  Z*ZFAC  + ZREF
        F(5,IT1)=  P*PFAC  + PREF
        F(6,IT1)=  S*SFAC  + SREF
        F(7,IT1)=TIM*TIFAC + TIREF 
        RET(IT1)=RETI
        DPR(IT1)=DPRI
C        IREP(IT1) = IREPI
        IREP(IT1) = IT1
        SORT(IT1) = SORTI
        AMQ(1,IT1) = AMQ1
        AMQ(2,IT1) = AMQ2
        AMQ(3,IT1) = AMQ3
        AMQ(4,IT1) = AMQ4
        AMQ(5,IT1) = AMQ5

C TESTS COOLING NUFACT--------------------
C        AM2 = AMQ(1,IT)*AMQ(1,IT)
C        P0 = BORO*CL9*AMQ(2,IT)
C        P = P0*F(1,IT)        
C        ENRG = SQRT(P*P+AM2)-AMQ(1,IT)
C        IF(ENRG.LT.100)  iex(it) = -9
C        IF(ENRG.GT.300)  iex(it) = -9
C------------------------------------------


CC------- If a series of particle numbers have not been found :
C        IF(IT1.LT.IT) THEN 
C          DO 117 II=IT1,IT-1
C            IEX(II)=-9
C            LET(I)='*'            
C 117      CONTINUE
C          IT1=IT
C          IPASS1=IPASSR
C        ENDIF

        IF(OKINIT) THEN
           DO 116 J=1,MXJ
 116         FO(J,IT1)=F(J,IT1)
        ENDIF

        IF(IT1 .EQ. MXT) GOTO 169
        IF(IT1 .EQ. KT2) GOTO 169
C         write(*,*) ' sbr obj3  ',it,it1,ipassr
      GOTO  17
 
 96   WRITE(TXT,FMT='(
     >'' *** SBR OBJ3 : ERROR  at  OPEN  file  '',A)') NOMFIC
      IF(NRES.GT.0) WRITE(NRES,FMT='(A130)') TXT
      CALL ENDJOB(TXT,-99)
 
 97   WRITE(NRES,FMT='(/,
     > '' SBR OBJ3 -> error in  reading  file '',
     >  A,'' at  event/traj #  '',I6,''/'')') NOMFIC,IT1
 95   CONTINUE
 169  CONTINUE
C Changed for linear FFAG sudies
C      IMAX=IT

C         write(*,*) ' sbr obj3  ',it,it1,ipassr

      IMAX=IT1
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
        WRITE(NRES,FMT='(/,T5,''  Reading  in  file  '',A,
     >    ''  ended  after  gathering '',I6,''  particles'')') 
     >    NOMFIC(DEBSTR(NOMFIC):FINSTR(NOMFIC)), IMAX
        WRITE(NRES,*) '  in  requested  range :  [',KT1,',',KT2,'].'
        IF(IS.GT.0) WRITE(NRES,FMT='(/,T5,I6,
     >    ''  particles  are  of  secondary  type  (LET="S")'')') IS
        IF(II.GT.0) WRITE(NRES,
     >    FMT='(/,T5,I6,''/'',I6,'' particles  have  IEX < 0,  hence'',
     >    I6,''  only  left  to  be  ray-traced. '')') II,IMAX,IMAX-II
        WRITE(NRES,*) '  Last  pass  number  read  :  ',IPASSR, 
     >    '   in requested range  :   [',KP1,',',KP2,'].'
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
