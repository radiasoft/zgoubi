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
      SUBROUTINE FMATCH(NLOG,NOMFIC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*(*) NOMFIC

C----- Calculate the fringe field coefficients Co-C5,
C      from a matching between, on one hand the analytical model
C        F / Fo = 1 / ( 1 + exp( P(X) ), P(X) = Co + C1*X ... + C5*X**5
C      and on the other hand numerical values read from the input 
C      file NOMFIC (Default name = fmatch.dat or last open file).
C      The matching method uses a Fletcher-Reeves alternated 
C      gradient algorithm.
C         Longitudinal coordinate should be Meters

      COMMON/CDF/ IES,IORDRE,LCHA,LIST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN

      PARAMETER(MXV=10)
      DIMENSION V(MXV), VMI(MXV), VMA(MXV)
      PARAMETER(MXC=400)
      DIMENSION CX(MXC), CY(MXC)

C      CHARACTER REP, NOMFIC*(80)
C      SAVE NOMFIC
      CHARACTER REP
      LOGICAL RFMTCH, OKOPN, EMPTY, NORM
      SAVE OKOPN
      EXTERNAL FRINGE, DFRING
      LOGICAL IDLUNI   
      SAVE NL

      DATA NORM / .TRUE. /
      DATA VMI/ 10 * -10.D0 /
      DATA VMA/ 10 * 10.D0 /

C      DATA NV / 4 /
C      DATA V/ .1D0, 6.3D0, -1.5D0, 3.6D0, -2.1D0, 1.7D0, 4 * 0.D0 /
C      DATA V/ .24D0, 1.8D0, -.6D0, .4D0, 6 * 0.D0 / 
C-------- Units are m
C      DATA XLAMB / 20.d-2 /
C      DATA XLAMB /   8.d-2 /

CC----- LHC Dipole
C      DATA NV / 6 /
C      DATA V/ .15D0, 3.87D0, -2.36D0, 2.98D0, 12.6D0, 15.D0, 4 * 0.D0 / 
C      DATA XLAMB / 11.2d-2 /
C----- RECYCLER Dipole
C      DATA NV / 4 /
C      DATA V/  0.09650,  3.76444, -0.70378,  1.31734,          6* 0.D0 / 
C      DATA XLAMB /  5.0d-2 /
C----- Fringe field Bosser/Jung
C      DATA NV / 6 /
C      DATA V/ .15D0, 3.87D0, -2.36D0, 2.98D0, 12.6D0, 15.D0, 4 * 0.D0 / 
C      DATA XLAMB / 10.d-2 /
C----- Fringe field musr bends (~ GSI quad)
      DATA NV / 6 /
      DATA V/ .015527D0, 3.875D0,-2.362D0,2.9782D0,
     >          12.6044D0,15.0257D0 , 4 * 0.D0 / 
      DATA XLAMB / 11.2d-2 /

C----- LHC Quadrupole
C      DATA NV / 6 /
C      DATA V/-.01D0, 5.46D0, .997D0, 1.57D0,-5.67D0,18.5D0, 4 * 0.D0 / 
C      DATA XLAMB / 5.6d-2 /

C----- Open data file. Normally contains two columns: CX, CY,
C       optionnally followed by  XLAMB, NV
C      NL = 20
C      CLOSE(NL)
C      OKOPN = .FALSE. 
C      OPEN(UNIT=NL, FILE=NOMFIC,STATUS='old',ERR=21)
C      OKOPN = .TRUE.

      GOTO 21
      
 20   CONTINUE      
      CALL FBGTXT
      WRITE(*,FMT='(/,'' Press RETURN for more '')') 
      READ(*,FMT='(A1)',ERR=20) REP 

 21   CONTINUE
      CALL FBGTXT
C      CALL HOMCLR

      WRITE(*,104) NV-1,EPS, XLAMB
 104  FORMAT(//,3X,60('-'),//,15X
     X,' MENU  FOR FRINGE-FIELD  MATCHING:',/
     1,/,5X,' 1  Open data file (Default = fmatch.dat)'
     2,/,5X,' 2  Degree  of  P(X) = Co + C1*X...+ C5*X5'
     X,/,5X,'                           (present = ',I1,')'
     3,/,5X,' 3  Precision   (Present = ',E10.2,')'
     4,/,5X,' 4  Lambda  (', F15.4 , ')'
     5,/,5X,' 5  Normalization to 1'
     7,/,5X,' 7  *** MATCH ***'
     8,/,5X,' 8  Print screen'
     9,/,5X,' 9  EXIT  THIS  MENU '
     2,/,5X,'12  ERASE  DISPLAY '
     >,/,3X,60('-'),//)

      WRITE(*,100)
 100  FORMAT('$  YOUR  CHOICE : ')
      READ(*,108,ERR=21) IOPT
108   FORMAT(I2)
      GOTO ( 1, 2, 3, 4, 5,21, 7, 8,98,21,21,12) IOPT  
      GOTO 21
 
 1    CONTINUE
C----- Open data file. Normally contains two columns: CX, CY,
C       optionnally followed by XLAMB, NV :
CC  FRINGE FIELD IN LHC DIPOLE
CC   -1.00000E+00  0
CC   -.800000E+00  0
CC   -.600000E+00  0
CC   -.400000E+00  0
CC        ...
CC   1.8800000E+01   6.862883E+00
CC   1.9000000E+01   6.862883E+00
CC   1.920000E+01   6.862883E+00
CC Some more data :
CC   11.2      Lambda, for LHC dipole
CC   6       P(X) of degree 5
      WRITE(*,*) 
     >  ' Give input data file name (default will be fmatch.dat) :'
      READ(*,FMT='(A)') NOMFIC
      IF(EMPTY(NOMFIC)) NOMFIC = 'fmatch.dat'
      CLOSE(NL)
      OKOPN = .FALSE.
      WRITE(*,*)
      WRITE(*,*) ' Trying to open ',NOMFIC,' ...' 
      IOS = 1
      IF (IDLUNI(NL)) OPEN(UNIT=NL, FILE=NOMFIC,STATUS='OLD',IOSTAT=IOS)
      IF(IOS .NE. 0) THEN
        WRITE(*,*) '  Cannot open ',nomfic
      ELSE
        OKOPN = .TRUE.
        WRITE(*,*) ' Ok !' 
      ENDIF
      GOTO 20

 2    CONTINUE
        WRITE(*,*) ' Degree  of  P(X) = Co + C1*X...+ C5*X5 :'
        READ(*,*,ERR=2) IDEG
        NV = IDEG + 1
        IF(NV.GT.6) GOTO 2
      GOTO 21

 3    CONTINUE
      GOTO 21

 4    CONTINUE
        WRITE(*,*) '  Give  Lambda (Presently : ',XLAMB,' m) :'
        READ(*,*,ERR=4) XLAMB
      GOTO 21

 5    CONTINUE
        WRITE(*,*) '  You want normalization to 1 ? (y/n) :'
        READ(*,*,ERR=4) REP
        IF(EMPTY(REP)) REP='y'
        NORM = REP .EQ. 'Y' .OR. REP .EQ. 'y'
      GOTO 21

 7    CONTINUE
        IF(.NOT.OKOPN) THEN
          NOMFIC = 'fmatch.dat'
          IOS = 1
          IF(IDLUNI(NL))
     >            OPEN(UNIT=NL,FILE=NOMFIC,STATUS='OLD',IOSTAT=IOS)
          IF(IOS .NE. 0) THEN
            WRITE(*,*) '  Cannot open ',nomfic
            GOTO 20
          ELSE
            OKOPN = .TRUE.
          ENDIF
        ENDIF
        IF(OKOPN) THEN
          IF ( RFMTCH(NL,MXC,NORM,
     >                       NC,CX,CY,XLAMB,NV) ) THEN

            CALL CLASS(1,NC,CX,CY,
     >                            XFB)
            WRITE(*,*) '  Estimated position of magnetic face is ',
     >                                 XFB

C----------- Normalisation coeff. XLAMB = 
C             parameter Lambda of the Zgoubi Fringe field
            CX(MXC-1) = XFB
            CX(MXC) = XLAMB

            WRITE(*,*) '  Busy,  matching...'
            CALL MATCH(FRINGE,DFRING,NV,V,VMI,VMA,NC,CX,CY,.TRUE.)
            CALL FMTCGR(FRINGE,NLOG,NV,V,NC,CX,CY,CX(MXC-1),XLAMB)


          ENDIF
        ELSE
          WRITE(*,*) ' Cannot work:'
          WRITE(*,*) '   the data file is not open. Please open it'
        ENDIF
      GOTO 20
 
 8    CONTINUE
C        CALL MENVCF
        CALL SAVPLT
      GOTO 21

 12   CONTINUE
        CALL CLSCR
      GOTO 21

 98   CLOSE(NL)
      RETURN 
      END
