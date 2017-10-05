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
      SUBROUTINE RCAVIT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     **************************
C     READS DATA FOR ACC. CAVITY
C     **************************
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
C      PARAMETER (ISZTA=80)
C      CHARACTER(ISZTA) TA
C      PARAMETER (MXTA=45)
      INCLUDE "C.DONT.H"     ! COMMON/DONT/ TA(MXL,MXTA)

      CHARACTER(132) TXT132
      LOGICAL STRCON, OK
      INTEGER DEBSTR, FINSTR
      PARAMETER (I9=9)
      CHARACTER(30) STRA(I9)
      LOGICAL EMPTY, ISNUM
      PARAMETER (KSIZ=10)
      CHARACTER(KSIZ) KLE

      PARAMETER (MXH=5)

C     ....IOPT -OPTION
      LINE = 1
      READ(NDAT,FMT='(A)',ERR=90,END=90) TXT132      
      READ(TXT132,*,ERR=90,END=90) A(NOEL,1)
      IOPT = NINT(A(NOEL,1))
      IF(STRCON(TXT132,'!',
     >                     IS)) TXT132 = TXT132(DEBSTR(TXT132):IS-1)

      IF(STRCON(TXT132,'PRINT',
     >                         IS) 
     >.OR.  (NINT(10.D0*A(NOEL,1)) - 10*INT(A(NOEL,1))).EQ.1) THEN
        TA(NOEL,1) = 'PRINT'        
      ELSE
        TA(NOEL,1) = ' '        
      ENDIF
      IF    (STRCON(TXT132,'CEBAF',
     >                         IS)) THEN
        IF(EMPTY(TA(NOEL,1))) THEN 
           TA(NOEL,1) = 'CEBAF'
        ELSE
           TA(NOEL,1) = 
     >     TA(NOEL,1)(debstr(TA(NOEL,1) ):finstr(TA(NOEL,1) ))//' CEBAF'
        ENDIF
      ELSEIF(STRCON(TXT132,'CornellSynch',
     >                                IS)) THEN
        IF(EMPTY(TA(NOEL,1))) THEN 
           TA(NOEL,1) = 'CornellSynch'
        ELSE
           TA(NOEL,1) = 
     >     TA(NOEL,1)(debstr(TA(NOEL,1) ):finstr(TA(NOEL,1) ))//
     >     ' CornellSynch'
        ENDIF
      ELSEIF(STRCON(TXT132,'eRHIC_RCS',
     >                                IS)) THEN
        IF(EMPTY(TA(NOEL,1))) THEN 
           TA(NOEL,1) = 'eRHIC_RCS'
        ELSE
           TA(NOEL,1) = 
     >     TA(NOEL,1)(debstr(TA(NOEL,1) ):finstr(TA(NOEL,1) ))//
     >     ' eRHIC_RCS'
        ENDIF
      ENDIF

C     ....FREQ. (Hz) or orbit length (m), H -HARMONIQUE or list {ha, hb, ...}
      LINE = LINE + 1
C      READ(NDAT,*,ERR=90,END=90) A(NOEL,10),A(NOEL,11)
      READ(NDAT,FMT='(A)',ERR=90,END=90) TXT132      
      IF(IOPT .NE. 0) then 
        IF(STRCON(TXT132,'!',
     >                     IS)) TXT132 = TXT132(DEBSTR(TXT132):IS-1)
        CALL STRGET(TXT132,I9
     >                    ,MSTR,STRA)
        IF(IOPT .EQ. 2) THEN
          IF(STRCON(TXT132,'{',
     >                        IS)) THEN
            OK = STRCON(TXT132,'}',
     >                             JS) 
            IF(.NOT. OK) GOTO 90
            TXT132 = TXT132(IS+1:JS-1)
            I = 2
            READ(TXT132,*,ERR=90,END=90) TMP
            WRITE(STRA(I),*) TMP
            DO WHILE(STRCON(TXT132,',',
     >                                 IIS)) 
              I = I + 1
              TXT132 = TXT132(IIS+1:FINSTR(TXT132))
              READ(TXT132,*,ERR=90,END=90) TMP
              WRITE(STRA(I),*) TMP
            ENDDO
            MSTR = I
            IF(MSTR-1 .GT. MXH) THEN
              WRITE(NRES,*) 
     >        'Pgm rcavit. Too many harmonics. Should be < ',MXH
              GOTO 90
            ENDIF
          ELSE
            MSTR=2       
          ENDIF
      
        ELSE
          MSTR=2       
        ENDIF

        DO I = 1, MSTR
          IF(ISNUM(STRA(I))) THEN 
            READ(STRA(I),*,ERR=90,END=90) A(NOEL,9+I)
          ELSE
            CALL ENDJOB('Pgm rcavit. Check input data, '
     >      //' non-numerical data found at line ',LINE)
          ENDIF
        ENDDO
      ENDIF
C Nb of harmonics
      NBH = MSTR-1
      IF(IOPT .EQ. 2) A(NOEL,19) = NBH

C V(Volts); PHS(rd) : dW = q*V sin( H*OMEGA*T + PHS); SR loss at pass #1 for computation of compensation (cav. 21)
C      READ(NDAT,*) A(NOEL,20),A(NOEL,21)
      LINE = LINE + 1
      READ(NDAT,FMT='(A)',ERR=90,END=90) TXT132      
      IF(IOPT .NE. 0) then 
        IF(STRCON(TXT132,'!',
     >                     IS)) TXT132 = TXT132(DEBSTR(TXT132):IS-1)
        CALL STRGET(TXT132,I9
     >                    ,MSTR,STRA)

        IF    (IOPT .EQ. 10) THEN
C        3rd data is IDMP if IOPT=10
          IF(MSTR.LE.2) THEN
            A(NOEL,22)= 2      ! Chambers matrix approximations. Default=none.
            IF(MSTR.LE.3) A(NOEL,23)= 0.D0      ! BORORef setting. Default=0.
          ENDIF
        ELSEIF(IOPT .EQ. 11) THEN
C        3rd data is energy loss at first pass, U0
          IF(MSTR.LE.2) THEN
            CALL ENDJOB(
     >      'Pgm cavite. SR loss U0[MeV] at first turn is needed. ',-99)
          ELSE
            READ(STRA(3),*) A(NOEL,22)   ! SR loss U0 at injection energy
          ENDIF
        ELSEIF(IOPT .EQ. 2) THEN
          IF(STRCON(TXT132,'{',
     >                        IS)) THEN
            TXT132 = TXT132(IS+1:FINSTR(TXT132))
            OK = STRCON(TXT132,'}',
     >                             JS) 
            IF(.NOT. OK) GOTO 90
            I = 1
            READ(TXT132(1:JS-1),*,ERR=90,END=90) TMP
            WRITE(STRA(I),*) TMP
            DO WHILE(STRCON(TXT132,',',
     >                                 IIS)) 
              I = I + 1
              TXT132 = TXT132(IIS+1:FINSTR(TXT132))
              READ(TXT132(1:JS-1),*,ERR=90,END=90) TMP
              WRITE(STRA(I),*) TMP
            ENDDO
            MSTR = I+1
            READ(TXT132(JS+1:FINSTR(TXT132)),*,ERR=90,END=90) TMP
            WRITE(STRA(MSTR),*) TMP
          ELSE
            MSTR=2       
          ENDIF
      
        ELSE
          MSTR=2       
        ENDIF

        DO I = 1, MSTR
          IF(ISNUM(STRA(I))) THEN 
            READ(STRA(I),*,ERR=90,END=90) A(NOEL,19+I)
          ELSE
            CALL ENDJOB('Pgm rcavit. Check input data, '
     >      //' non-numerical data found at line ',LINE)
          ENDIF
             
        ENDDO
      ENDIF

      RETURN

 90   CONTINUE
      CALL ZGKLEY( 
     >            KLE)
      CALL ENDJOB('*** Pgm rcavit, keyword '//KLE//' : '// 
     >'input data error, at line #',line)
      RETURN 
      END
