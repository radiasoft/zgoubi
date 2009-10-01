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
      SUBROUTINE HMATCH(NLOG,NOMFIC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*(*) NOMFIC

C----- Calculate the coefficients Co-C5 of an Hermit series,
C      from a matching between, on one hand the analytical model
C        F / Fo = sum(Ci*exp(-x2/2/sig2)*Hi(x/sig))
C      and on the other hand numerical values read from the input 
C      file NOMFIC (Default name = hmatch.dat or last open file).
C      The matching method uses a Fletcher-Reeves alternated 
C      gradient algorithm.

      COMMON/CDF/ IES,IORDRE,LCHA,LIST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN
      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB

      PARAMETER(MXV=10)
      PARAMETER(MXC=400)
      DIMENSION V(MXV), VMI(MXV), VMA(MXV)
      DIMENSION CX(MXC), CY(MXC)

      CHARACTER REP
      LOGICAL RHMTCH, OKOPN, EMPTY
      EXTERNAL SHERMI, DSHERM
      
      DATA VMI/ 10 * -10.D0 /
      DATA VMA/ 10 * 10.D0 /

C------ Gauss profile
      DATA NV / 6 /
      DATA V/ 1.d0, 0.d0 , .5d0, 7*0.D0 /
      DATA SIG / 1.D0 /      

C----- Open data file. Normally contains two columns: CX, CY,
C       optionnally followed by  SIG, NV
      NL = 20
      CLOSE(NL)
      OKOPN = .FALSE. 
      OPEN(UNIT=NL, FILE=NOMFIC,STATUS='old',ERR=21)
      OKOPN = .TRUE.

      GOTO 21
      
 20   CONTINUE      
      CALL FBGTXT
      WRITE(*,FMT='(/,''  Press RETURN for more'')') 
      READ(*,FMT='(A1)',ERR=20) REP 

 21   CONTINUE
      CALL FBGTXT
C      CALL HOMCLR

      WRITE(*,104) NV-1,EPS, SIG
 104  FORMAT(//,3X,60('-'),//,15X
     X,' MENU  FOR FRINGE-FIELD  MATCHING:',/
     1,/,5X,' 1  Open data file (Default = hmatch.dat)'
     2,/,5X,' 2  Degree  of  the Hermit series'
     X,/,5X,'                           (present = ',I1,')'
     3,/,5X,' 3  Precision   (Present = ',E10.2,')'
     4,/,5X,' 4  Sigma of the sample distribution (', F15.4 , ')'
     7,/,5X,' 7  *** MATCH ***'
     8,/,5X,' 8  Print screen'
     9,/,5X,' 9  EXIT  THIS  MENU '
     2,/,5X,'12  ERASE  DISPLAY '
     >,/,3X,60('-'),//)

      WRITE(*,100)
 100  FORMAT('$  YOUR  CHOICE : ')
      READ(*,108,ERR=21) IOPT
108   FORMAT(I2)
      GOTO ( 1, 2, 3, 4,21,21, 7, 8,99,21,21,12) IOPT  
      GOTO 21
 
 1    CONTINUE
C----- Open data file. Normally contains two columns: CX, CY,
C       optionnally followed by SIG, NV
      CLOSE(NL)
      OKOPN = .FALSE.
      WRITE(*,*) 
     >  ' Give input data file name (default will be hmatch.dat) :'
      READ(*,FMT='(A)') NOMFIC
      IF(EMPTY(NOMFIC)) NOMFIC = 'hmatch.dat'
      WRITE(*,*)
      WRITE(*,*) ' Trying to open ',NOMFIC,' ...' 
      OPEN(UNIT=NL, FILE=NOMFIC,STATUS='OLD',IOSTAT=IOS)
      IF(IOS .NE. 0) THEN
        WRITE(*,*) '  Cannot open ',nomfic
      ELSE
        OKOPN = .TRUE.
        WRITE(*,*) ' Ok !' 
      ENDIF
      GOTO 20

 2    CONTINUE
        WRITE(*,*) ' Degree  of  H(X) = CoHo + C1*H1...+ C9*X9 :'
        READ(*,*,ERR=2) NV
        IF(NV.GT.10) GOTO 2
      GOTO 21

 3    CONTINUE
      GOTO 21

 4    CONTINUE
        WRITE(*,*) '  Give  Sigma (Presently ',SIG,') :'
        READ(*,*,ERR=4) SIG
      GOTO 21

 7    CONTINUE
        IF(OKOPN) THEN
          IF ( RHMTCH(NL,MXC,
     >                       NC,CX,CY,SIG,NV) ) THEN

C            CALL CLASS(0,NC,CX,CY,MXC,
C     >                                DUM)
            CALL CLASS(0,NC,CX,CY,
     >                            DUM)

C----------- Normalisation coeff. SIG = 
C               Sigma of the distribution
            CX(MXC-1) = 0.
            CX(MXC) = SIG

            WRITE(*,*) '  Busy,  matching...'
            CALL MATCH(SHERMI,DSHERM,NV,V,VMI,VMA,NC,CX,CY,MXV,MXC,
C     >           .TRUE.)
     >           .FALSE.)
C            CALL HMTCGR(SHERMI,NLOG,NV,V,NC,CX,CY,CX(MXC-1),
C     >            SIG,MXV,MXC)
            CALL HMTCGR(SHERMI,NLOG,NV,V,NC,CX,CY,CX(MXC-1),SIG)

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

 99   CLOSE(NL)
      RETURN 
      END
