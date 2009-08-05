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
      SUBROUTINE CLORBI(NLOG,OKECH,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL OKECH

C----- Plot the closed orbit as read from file zgoubi.averageOrbit

      COMMON/CDF/ IES,IORDRE,LCHA,LIST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN    
      COMMON/UNITS/ UNIT(6)
      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB

      PARAMETER (MXPU=10000)
      DIMENSION FCO(8,MXPU), FCO2(7,MXPU)
      DIMENSION CO1(6),CO2(6),COM(6)      

      CHARACTER TAXV(8)*11, TAXH(2)*9, TXT*80
      SAVE TAXV
      SAVE TAXH
      LOGICAL IDLUNI

      CHARACTER REP, NOMFIC*80
      SAVE NOMFIC
      SAVE KAXV, KAXH, KPU
      LOGICAL OKOPN, CHANGE
      SAVE OKOPN

      DATA TAXV / 'X_co', 'X`_co', 'Z_co', 'Z`_co', ' dist.',
     >  'dp/p', '# particles', 'PU posit.' /

      DATA TAXH / 'path  (m)', 'pick-up  ' /

      DATA NOMFIC /'none'/

      DATA KAXV, KAXH, KPU / 1, 1, -1 /

      DATA OKOPN, CHANGE / .FALSE., .FALSE. /

      GOTO 21

 20   CONTINUE      
      CALL FBGTXT
      WRITE(*,FMT='(/,'' Press RETURN for more '')') 
      READ(*,200,ERR=20) REP 
 200  FORMAT(A1)

 21   CONTINUE
      CALL FBGTXT
      CALL HOMCLR

      WRITE(*,104) NOMFIC,TAXV(KAXV),TAXH(KAXH), KPU
 104  FORMAT(5X,' MENU  -  Average orbit, closed orbit : ',/, 
     1/,5X,' 1  OPEN  FILE     -  current is ',A,
     3/,5X,' 3  To be plotted -    now,  ',A, 
     4/,5X,' 4  Horizontal axis :  path length / pick-up occurence',
     >                                                  '   - now,  ',A, 
     5/,5X,' 5  Observation pick-up #  - now (-1 is all),  ',I6, 
     7/,5X,' 7  **  PLOT  ** '             ,  
     8/,5X,' 8  Print screen',
     9/,5X,' 9  Exit  this  menu ',
     2/,5X,'12  Erase  display    ',
     >/)

      WRITE(*,100)
 100  FORMAT('$  Option  number : ')
      READ(*,108,ERR=21) IOPT
108   FORMAT(I2)
      IF(.NOT. OKOPN .AND. IOPT.EQ.7) THEN
        CALL OPNWRN(1)           
        GOTO 20
      ELSE
        GOTO ( 1, 21, 3, 4, 21, 21, 7, 8,99,21,21,12) IOPT  
        GOTO 21
      ENDIF

 1    CONTINUE    
        NOMFIC = 'zgoubi.averageOrbit'
        IF (IDLUNI(IUN)) CALL OPNMNL(IUN,
     >                                   NFCO,NOMFIC,OKOPN,CHANGE)
      GOTO 20

 3    CONTINUE
      KAXV0 = KAXV
      WRITE(6,FMT=
     >'('' To be plotted : X_co, X`_co, Z_co, Z`_co, dist.,'', 
     >''dp/p, # particles, PU posit.  (1-8) : '',T92,I1,A2)')
     > KAXV,'  '
      READ(5,FMT='(I1)',ERR=31) KAXV
      IF(KAXV .GE. 1 .AND. KAXV .LE. 8) GOTO 21
 31   KAXV = KAXV0
      GOTO 21

 4    CONTINUE
      KAXH0 = KAXH
      WRITE(*,FMT=
     > '(''  Horizontal axis :  path length / PU occur.  (1/2) : '')')
      READ(5,FMT='(I1)',ERR=41) KAXH
      IF(KAXH .GE. 1 .AND. KAXH .LE. 2) GOTO 21
 41   KAXH = KAXH0
      GOTO 21

 5    CONTINUE
      KPU0 = KPU
      WRITE(*,FMT=
     > '(''  Observation  pick-up  #  (-1 for all) : '')') 
      READ(5,FMT='(I1)',ERR=51) KPU
      IF(KPU .GE. -1) THEN 
        GOTO 21
      ELSE
        GOTO 5
      ENDIF
 51   KPU = KPU0
      GOTO 21

 7    CONTINUE
      IF(.NOT. OKOPN) THEN
        CALL OPNWRN(1)
        GOTO 20
      ENDIF

      REWIND(NFCO)
      CALL HEADER(NFCO,4,.FALSE.,*20)
      IP = 1
 71   CONTINUE
        READ(NFCO,FMT=*,END=77,ERR=77)
     >    IPU,PUPOS,(FCO(J,IP),J=1,7),IPASS
        FCO(8,IP) = PUPOS
        READ(NFCO,FMT=*,END=77,ERR=77) (FCO2(J,IP),J=1,6)

        IF(KPU.GT.0) THEN
          IF(IPU.NE.KPU) GOTO 71
        ENDIF

          DO 721 I=1,6
            FCO(I,IP) = FCO(I,IP) * UNIT(I)/ FCO(7,IP)
            FCO2(I,IP) = FCO2(I,IP) * UNIT(I)*UNIT(I)/ FCO(7,IP)
 721      CONTINUE
          FCO(8,IP) = FCO(8,IP) *  1.D-2

        IP = IP + 1

        GOTO 71

 77   CONTINUE
      NPU = IP-1
      WRITE(*,*) ' End of READ in ',NOMFIC,'.  Total # of P.U. = ',NPU

      IF( .NOT. OKECH ) THEN
        YMI = 1.D10
        YMA = -1.D10
        SIGMA = -1.D10
        DO 75 I = 1, NPU
          YM = FCO(KAXV,I) 
          SIGY = SQRT(FCO2(KAXV,I)  - YM*YM)
          IF(YM .LT. YMI) YMI = YM
          IF(YM .GT. YMA) YMA = YM
          IF(SIGY .GT. SIGMA) SIGMA=SIGY
 75     CONTINUE
        IF(KAXV.EQ.7) THEN
        ELSE
          YMI = YMI - 1.05D0 * SIGMA
          YMA = YMA + 1.05D0 * SIGMA
        ENDIF
        IF(KAXH.EQ.1) THEN 
          XMI = FCO(5,1)
          XMA = FCO(5,NPU)
        ELSE
          XMI = 1.D0
          XMA = NPU
        ENDIF
        WRITE(*,*) ' XMI-MA, YMI-MA  (m) : ',XMI,XMA,YMI,YMA

          IF(XMI .LT. XMA .AND. YMI .LT. YMA) THEN
            OKECH = .TRUE.
            CALL TXTFBG
            CALL TRAXES(XMI,XMA,YMI,YMA,-1) 
          ELSE
            WRITE(*,*) ' Scale problems ; cannot plot !!'
            GOTO 20
          ENDIF          
      ENDIF

C--------------------
      IF(KAXH.EQ.1) THEN 
        X1 = FCO(5,1)
      ELSE
        X1 = 1.D0
      ENDIF
      YM1 = FCO(KAXV,1)
      SIGY =  SQRT(FCO2(KAXV,1) - YM1*YM1)
      IF(KAXV.EQ.7) THEN
        YA1 = YM1 
        YB1 = YM1 
      ELSE
        YA1 = YM1 - SIGY
        YB1 = YM1 + SIGY
      ENDIF
      SIG1 = SIGY
      CALL LINTYP(2)
      DO 771 I = 2, NPU
        IF(KAXH.EQ.1) THEN 
          X2 = FCO(5,I)
        ELSE
          X2 = FLOAT(I)
        ENDIF
        YM2 = FCO(KAXV,I) 
        SIGY = SQRT(FCO2(KAXV,I) - YM2*YM2)
        YA2 = YM2 - SIGY
        YB2 = YM2 + SIGY
        SIG2 = SIGY
        IF(LIS .EQ. 2) CALL IMPV(NLOG,I,X1,YM1,DUM,DUM,IDUM,IDUM)
        CALL VECTPL(X1,YA1,4)
        CALL VECTPL(X2,YA2,2)
        CALL LINTYP(1)
        CALL VECTPL(X1,YM1,4)
        CALL VECTPL(X2,YM2,2)
        CALL LINTYP(2)
        CALL VECTPL(X1,YB1,4)
        CALL VECTPL(X2,YB2,2)
C        CALL VECTPL(X1,SIG1,4)
C        CALL VECTPL(X2,SIG2,2)
        X1 = X2
        YM1 = YM2
        YA1 = YA2
        YB1 = YB2
        SIG1 = SIG2
 771  CONTINUE

      CALL LOGO
      WRITE(TXT,109) TAXV(KAXV),' S (m)'
 109  FORMAT(1X,A,' v.s. ',A) 
      CALL TRTXT(120.D0,245.D0,TXT,80,0)
      WRITE(TXT,103) TAXV(KAXV),NPU
      CALL TRTXT(10.D0,21.D0,TXT,50,0)
      WRITE(TXT,*) '   '
      CALL TRTXT(.2D0,1.1D0,TXT,50,0)

      CALL FBGTXT

      WRITE(*,103) TAXV(KAXV),NPU
 103  FORMAT(1X,A,' mean orbit, ',I4,' pickups') 
      WRITE(*,101) TAXV(KAXV),TAXH(KAXH)
 101  FORMAT(A,' v.s. ',A) 
      WRITE(*,102) CO1(KAXV),CO2(KAXV),COM(KAXV)
 102  FORMAT(' Mean, Sigma, Extrem (m)',1P,3E12.4) 

      GOTO 20

 8    CONTINUE
C        CALL MENVCF
        CALL SAVPLT(*20)
      GOTO 21

 12   CONTINUE
        CALL CLSCR
        OKECH=.FALSE.
      GOTO 21

 99   CONTINUE
      RETURN

      END
