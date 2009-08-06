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
      SUBROUTINE FITGTV(NOMFIC,
     >                         FITGET)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER NOMFIC*(*)
      LOGICAL FITGET

      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
 
      CHARACTER*60 NAMFIC
      INTEGER DEBSTR,FINSTR

      LOGICAL OPN, IDLUNI
      DIMENSION AFIT(MXL,MXD), NEWVAL(MXL,MXD)
      SAVE AFIT, NMI, NMA, NEWVAL
      PARAMETER (NNEWVAL=MXL*MXD)
      DATA NEWVAL / NNEWVAL*0/

      FITGET = .FALSE.

      NAMFIC=NOMFIC(DEBSTR(NOMFIC):FINSTR(NOMFIC))

      IF(NRES .GT. 0) WRITE(NRES,100) NAMFIC
 100  FORMAT(/,10X, 
     >'Parameter values to be refreshed from FIT storage file : ',A,/)

      IF(NAMFIC.EQ.'NONE' .OR. NAMFIC.EQ.'none') GOTO 97

      INQUIRE(FILE=NAMFIC,ERR=98,OPENED=OPN,NUMBER=LUN)
      IF(OPN) CLOSE(UNIT=LUN)
        
      IF(IDLUNI(
     >          LUN)) OPEN(UNIT=LUN,FILE=NAMFIC,ERR=99)

      WRITE(NRES,101) NAMFIC
 101  FORMAT(/,10X,' Ok, opened storgae file ',A,/)

      CALL HEADER(LUN,NRES,4,.FALSE.,
     >                              *999)
      
      NMI = 999999
      NMA = -999999
      KREAD = 0
 1    CONTINUE
        READ(LUN,401,ERR=10,END=10) NUML,I,ISI,XK,XII,XI,XJ,PI
 401    FORMAT( 
     >  2X,I3,3X,I2,4X,I3,2(2X,G10.3),2X,G15.8,2(1X,G10.3))
C        WRITE(6,400) NUML,I,ISI,XK,XII,XI,XJ,PI
C 400    FORMAT(1P, 
C     >  2X,I3,3X,I2,4X,I3,2(2X,G10.3),2X,G15.8,2(1X,G10.3))
        IF(NUML.LT.NMI) NMI=NUML
        IF(NUML.GT.NMA) NMA=NUML
        AFIT(NUML,ISI) = XI
        NEWVAL(NUML,ISI) = 1
        KREAD = KREAD+1
        GOTO 1

 10   CONTINUE
      IF(KREAD.GE.1) FITGET = .TRUE.
      IF(NRES .GT. 0) WRITE(NRES,103) KREAD,NAMFIC
 103  FORMAT(/,10X, I3, 
     >  ' variables have been read from FIT storage file ',A,/)

C----- OBJECT coordinates are redefined according to FIT variables, if any
C To be completed

      RETURN

 97   CONTINUE
      RETURN 
 98   CONTINUE
      CALL ENDJOB('*** Error, SBR FITGTV -> at INQUIRE ',-99)
      RETURN 
 99   CONTINUE
      CALL ENDJOB('*** Error, SBR FITGTV -> can`t open strage file',-99)
      RETURN 
 999  CONTINUE
      CALL ENDJOB('*** Error, SBR FITGTV ->  Read error at header',-99)
      RETURN 

      ENTRY FITGT1
        DO 2 ID = 1, MXD
          IF(NEWVAL(NOEL,ID) .EQ. 1) THEN
            TEMP = A(NOEL,ID) 
            A(NOEL,ID) = AFIT(NOEL,ID)
            IF(NRES .GT. 0) WRITE(NRES,
     >      FMT='('' GETFITVAL procedure.  Former  A('',I3,
     >      '','',I3,'') = '',1P,G12.4,''   changed  to  new value '', 
     >      G12.4)') NOEL,ID,TEMP,A(NOEL,ID)
          ENDIF
 2      CONTINUE
      RETURN
      END
