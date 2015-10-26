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
      SUBROUTINE ROBJET
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     *******************************************
C     READS DATA FOR OBJECT DEFINITION BY 'OBJET'
C     *******************************************
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      PARAMETER (LNTA=132) ; CHARACTER(LNTA) TA
      PARAMETER (MXTA=45)
      INCLUDE "C.DONT.H"     ! COMMON/DONT/ TA(MXL,MXTA)
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      CHARACTER(1) LET
      INCLUDE "C.FAISCT.H"     ! COMMON/FAISCT/ LET(MXT)
      INCLUDE "C.OBJET.H"     ! COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      DIMENSION IA(5)

      PARAMETER(MXJ1=MXJ-1)
      PARAMETER(I0 = 0)

C  If change MXREF, make sure 1/ to change MXD (see  line 140), 2/ that MXT is large enough
      PARAMETER(MXREF=999)

      CHARACTER(LEN=132) TXT132
      PARAMETER (MST=3)
      CHARACTER(LEN=20) STRA(MST)
      CHARACTER(LEN=80) TXT80
      LOGICAL STRCON
      INTEGER DEBSTR, FINSTR

C----- BORO
      LINE = 1
      READ(NDAT,*,ERR=99) A(NOEL,1)
C----- KOBJ - may be of the form "K.K2"
      LINE = 2
      READ(NDAT,FMT='(A)',ERR=99) TXT132
      IF(STRCON(TXT132,'!',
     >                     IS)) TXT132 = TXT132(1:IS-1)
      IF(STRCON(TXT132,'.',
     >                      IS)) THEN
        READ(TXT132(1:IS-1),*,ERR=99,END=99) K
        READ(TXT132(IS+1:FINSTR(TXT132)),*,ERR=99,END=99) K2
        A(NOEL,11) = K2
      ELSE
        K2 = 0 
        READ(TXT132,*,ERR=99,END=99) K
        A(NOEL,11) = I0
      ENDIF
      A(NOEL,10) = K
      IF(K .LT. 0) K=-K

      GOTO (1,2,3,1,5,5,7,8,9) K
      CALL ENDJOB('*** Error, SBR ROBJET -> No  such  object  KOBJ= ',K)
 
 1    CONTINUE
c      READ(NDAT,*,ERR=99) (IA(I),I=1,MXJ1)
c      DO 11 I=1,MXJ1
c 11     A(NOEL,19+I) = IA(I)
      READ(NDAT,*,ERR=99) (A(NOEL,I),I=20,25)
      READ(NDAT,*,ERR=99) (A(NOEL,I),I=30,35)
      READ(NDAT,*,ERR=99) (A(NOEL,I),I=40,45)
      RETURN
 
 2    CONTINUE
      READ(NDAT,*,ERR=99) (IA(I),I=1,2)
      A(NOEL,20) = IA(1)
      A(NOEL,21) = IA(2)
      II = 20
      DO 21 I=1,IA(1)
        II = II + 10
        READ(NDAT,*,ERR=98) (FO(J,I),J=2,MXJ1),FO(1,I),LET(I)
        IF(II .LE. 90) THEN
C--------- For allowing the use of the first 7 traj with FIT
          A(NOEL,II  ) = FO(2,I) 
          A(NOEL,II+1) = FO(3,I) 
          A(NOEL,II+2) = FO(4,I) 
          A(NOEL,II+3) = FO(5,I) 
          A(NOEL,II+4) = FO(6,I) 
          A(NOEL,II+5) = FO(1,I) 
        ENDIF
 21   CONTINUE
      READ(NDAT,*,ERR=99) (IEX(I),I=1,IA(1))
      RETURN
 
 3    CONTINUE
      CALL STRGET(TXT132,MST,
     >                     NST,STRA) 
      IF(NST .GT. MST) GOTO 90
      I = 1
      DOWHILE(STRA(I)(1:6).NE.'HEADER' .AND. I.LE.NST)
C Key is of the form 'HEADER_num', 0.le.num.le.9 
        I = I + 1
      ENDDO
      IF(I.LE.NST) THEN
        READ(STRA(I)(8:FINSTR(STRA(I))),*) NHDR
        A(NOEL,12) = NHDR
      ELSE
        A(NOEL,12) = 0.D0
      ENDIF
C----- Will read from part. #I1 to part. #I2, step
      LINE = 1
      READ(NDAT,*,ERR=99) (A(NOEL,I),I=20,22)
C----- Will read from  ipass #I1 to ipass #I2, step
      LINE = 1
      READ(NDAT,*,ERR=99) (A(NOEL,I),I=30,32)
C----- Y-,T-,Z-,P-,S-,DP-,TI-FAC, LETAG
      LINE = 1
      READ(NDAT,*,ERR=99) (A(NOEL,I),I=40,46),TA(NOEL,1)
C----- Y-,T-,Z-,P-,S-,DP-,TI-REF
      LINE = 1
      READ(NDAT,*,ERR=99) (A(NOEL,I),I=50,56)
C----- INIT flag (causes FO(j,i)=F(j,i) if A(NOEL,60)=0)
      LINE = 1
      READ(NDAT,*,ERR=99) A(NOEL,60)
C----- Name of trajectory data storage file
      LINE = 1
      READ(NDAT,100,ERR=99) TXT80
      IF(STRCON(TXT80,'!',
     >                      IS)) THEN
        TA(NOEL,2) = TXT80(1:IS-1)
      ELSE
        TA(NOEL,2) = TXT80(DEBSTR(TXT80):FINSTR(TXT80))
      ENDIF
 100  FORMAT(A80)
      RETURN
 
 5    CONTINUE
C      CALL OBJ52(K2)
      READ(NDAT,*,ERR=99) (A(NOEL,I),I=20,25)
      READ(NDAT,*,ERR=99) (A(NOEL,I),I=30,35)
      IF(K2.EQ.1) THEN
C------- Read initial beam, for possible transport by MATRIX or use by FIT
C alfy, bety,  alfz, betz,  alfx, betx, Dy, Dy', Dz, Dz'
        READ(NDAT,*,ERR=99) (A(NOEL,I),I=40,49)
      ELSEIF(K2 .GE. 2 .AND. K2 .LE. MXREF) THEN
C------- Read additional references 
        KK = 30
        DO 52 K = 2, K2
          KK = KK + 10
          IF(KK.GT.MXD-10) 
     >          CALL ENDJOB(' SBR ROBJET. MXD is too small.',-99)
          READ(NDAT,*,ERR=99) (A(NOEL,I),I=KK,KK+5)
 52       CONTINUE
      ELSEIF(K2 .GT. MXREF) THEN
        CALL ENDJOB(' SBR ROBJET. MXREF is too small.',-99)
      ENDIF
      RETURN
 
 7    CONTINUE
      RETURN
 
 8    CONTINUE
C----- IY, IZ, IX
      READ(NDAT,fmt='(a)',ERR=99) txt132
      READ(TXT132,*,ERR=99) (A(NOEL,19+I),I=1,3)
      CALL STRGET(TXT132,3,
     >                     IDUM,STRA) 
C Get number of digits after '.' in IX, IY, IZ      
      DO I = 1, 3
        IF(STRCON(STRA(I),'.'
     >                       ,IS)) THEN 
          A(NOEL,22+I) = FINSTR(STRA(I)(IS:132))-IS
        ENDIF
      ENDDO
C----- Center of ellipsoid (Y, T, Z, P, X, D
      READ(NDAT,*,ERR=99) (A(NOEL,I),I=30,35)
C----- alpha, beta, epsilon/pi, for Y, Z, X phase-spaces
      READ(NDAT,*,ERR=99) (A(NOEL,I),I=40,42)
      READ(NDAT,*,ERR=99) (A(NOEL,I),I=50,52)
      READ(NDAT,*,ERR=99) (A(NOEL,I),I=60,62)
      RETURN
 
 9    CONTINUE
C----- IY, IZ, IX = number of phase angles in Y, Z, X 
      READ(NDAT,*,ERR=99) (A(NOEL,19+I),I=1,3)
C----- Centers of ellipses (Y, T, Z, P, X, D
      READ(NDAT,*,ERR=99) (A(NOEL,I),I=30,35)
C----- alpha, beta, epsilon/pi, for Y, Z, X phase-spaces
      READ(NDAT,*,ERR=99) (A(NOEL,I),I=40,42)
      READ(NDAT,*,ERR=99) (A(NOEL,I),I=50,52)
      READ(NDAT,*,ERR=99) (A(NOEL,I),I=60,62)
      RETURN

 99   WRITE(6,*) 
     >  ' *** Execution stopped upon READ : invalid input in OBJET'
      WRITE(NRES ,*) 
     >  ' *** Execution stopped upon READ : invalid input in OBJET'
      GOTO 90
      
 98   WRITE(6,*) 
     >  ' *** Execution stopped upon READ : invalid input in OBJET',
     >  ' at particle #',I
      WRITE(NRES ,*) 
     >  ' *** Execution stopped upon READ : invalid input in OBJET',
     >  ' at particle #',I
      GOTO 90

 90   CALL ENDJOB('*** Pgm rmcobj. Input data error, at line ',line)
      RETURN
 
      END
