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
      SUBROUTINE PLTDAT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     -------------------------
C     PLOTDATA  PURPOSES  STUFF
C     -------------------------
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
C      COMMON/DON/ A(09876,99),IQ(09876),IP(09876),NB,NOEL
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),IMAX,IEX(MXT),IREP(MXT)
      COMMON/REBELO/ NPASS,IPASS,KWRT,NNDES,STDVM

      DIMENSION FF(MXT,19,8)
      LOGICAL  IDLUNI
      SAVE SSP

      CHARACTER*10 FNAME

      DATA FNAME / 'plotda.out' /
      DATA IPASS0,SP,MAXLOC/0,0,0/

 
      LOC = A(NOEL,1)       

      IF(IPASS0 .EQ. 0) THEN
        IF(IDLUNI(
     >            LUN)) THEN
          OPEN(UNIT=LUN,FILE=FNAME,ERR=99)
        ELSE
          GOTO 99
        ENDIF
      ENDIF

      IF( IPASS .NE. IPASS0 ) THEN
C       ... NEXT LOOP
        IPASS0=IPASS
        NLOC=0
      ELSE
        SP=SP-SSP
      ENDIF
 
C     ... LOCATION NUMBER
      NLOC= NLOC+1
      IF(NLOC .GT. 20)  CALL ENDJOB(
     >       'IN SBR PLTDAT : due to  # of  locations > ',20)
      IF(IPASS .EQ. 1) THEN
        IF(NLOC .GT. MAXLOC) MAXLOC=NLOC
      ENDIF
 
C     ... SLP COUNTS HOW MANY OVER IMAX MAKE IT AT THE ACTUAL LOOP
      SLP=0
      DO 1 I=1,IMAX
        IF(IEX(I) .LT. -1) GOTO 1
        SLP=SLP+1
        FF(I,NLOC,1) = 0D0
        FF(I,NLOC,2) = F(2,I)
        FF(I,NLOC,3) = F(4,I)
        FF(I,NLOC,4) = NLOC
        FF(I,NLOC,5) = F(3,I)
        FF(I,NLOC,6) = F(5,I)
        FF(I,NLOC,7) = F(6,I)
        FF(I,NLOC,8) = F(1,I)
 1    CONTINUE
 
C     ... SSP COUNTS HOW MANY OVER IMAX MAKE IT AT THE ACTUAL LOOP
      SSP=0D0
      IF( LOC .LT. 0 ) THEN
        SLP=0
        DO 2 I= 1,IMAX
          IF( IEX(I) .LT. -1) GOTO 2
          SLP=SLP+1
          SP=SP+1
          DO 3 ILOC=1,NLOC
            WRITE(LUN,100) SP,(FF(I,ILOC,J),J=1,8),I
 100        FORMAT(9E14.6,I3)
 3        CONTINUE
 2      CONTINUE
        SSP=SLP
      ENDIF
 
      IF(NRES .GT. 0) THEN
 
        WRITE(NRES,101) LOC,NLOC,IPASS,SLP,IMAX
 101    FORMAT(5X,/,' PLOTDATA  OUTPUT',//,15X,' LOCATION  #',I4
     >  ,'        AT  POSITION  #',I2,//,15X,' LOOP  NUMBER ',I5,10X
     >  ,F4.0,'  PARTICLES  OVER  ',I3,'  MAKE  IT  TO  HERE')
 
        IF    ( LOC .GT. 0 ) THEN
          WRITE(NRES,102) SLP
 102      FORMAT(15X,F7.0,' PARTICLES  STORED  IN  ARRAY')
        ELSEIF( LOC .LT. 0 ) THEN
          WRITE(NRES,103) SP, LUN
 103      FORMAT(15X,F7.0,' PARTICLES  PRINTED  IN  UNIT=',I2)
        ENDIF
 
      ENDIF
 
      IF(IPASS .GT. 1 .AND. IPASS .EQ. NPASS+1) THEN
       IF(NLOC .EQ. MAXLOC) THEN
         MAXLOC=0
         IPASS0=0
         SP=0D0
         CLOSE(LUN)
       ENDIF
      ENDIF
 
      RETURN
 99   WRITE(6,*) ' SBR ZGOUBI : ERREUR  OPEN FICHIER ',FNAME
      RETURN

      END
