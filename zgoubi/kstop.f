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
      SUBROUTINE KSTOP(IK,II,JEX,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      INCLUDE "C.CHAMBR.H"     ! COMMON/CHAMBR/ LIMIT,IFORM,YLIM2,ZLIM2,SORT(MXT),FMAG,YCH,ZCH
 
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      INCLUDE "C.INTEG.H"     ! COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      INCLUDE 'MXFS.H'
      PARAMETER (LBLSIZ=20)
      PARAMETER (KSIZ=10)
      CHARACTER(KSIZ) FAM ; CHARACTER(LBLSIZ) LBF
      INCLUDE "C.SCALT.H"     ! COMMON/SCALT/ FAM(MXF),LBF(MXF,MLF)
      INCLUDE "C.TRAJ.H"     ! COMMON/TRAJ/ Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT

      INTEGER DEBSTR, FINSTR
      CHARACTER(42) TXT
      CHARACTER(80) TXTELT

      LOGICAL FITING
      DATA FITING / .FALSE. /

      IF(JEX.LT.-1) GOTO 99
      KEX = -IK 

      IF    (IK .EQ. 2 ) THEN
        CALL CNTNRJ
        TXT = 'too many integration steps (SBR INTEG)'
      ELSEIF( IK .EQ. 3 ) THEN
        CALL CNTNRJ
        TXT = 'going backward (SBR INTEG)'
      ELSEIF(IK .EQ. 4 ) THEN
        CALL CNTOUT
        IEX(II) = KEX
        TXT = 'hit aperture limits'
      ELSEIF(IK .EQ. 5 ) THEN
        CALL CNTNRJ
        TXT = 'Iter>itmax, or dist>1d5 in  SBR ITER'
      ELSEIF(IK .EQ. 6 ) THEN
        CALL CNTNRJ
        TXT = 'enrgy loss > partcle enrgy (SBR GASCAT)'
      ELSEIF( IK .EQ. 7 ) THEN
        CALL CNTNRJ
        TXT = 'field steps > 50% in field map' 
      ELSEIF( IK .EQ. 8 ) THEN
        CALL CNTNRJ
        TXT = 'reached field limit in optical element'
      ELSEIF(IK .EQ. 9 ) THEN
        CALL CNTNRJ
        TXT = 'particle not accounted for'
      ELSEIF(IK .EQ. 10) THEN
        CALL CNTNRJ
        TXT = 'secondary particle has decayed'
        IEX(II) = KEX
      ENDIF
 
      CALL CNTMXR(
     >            IMX)
      CALL CNTSTO(
     >            NSTOP)
      CALL ZGLMNT(
     >            TXTELT)

      IF    ( IK .EQ. 2 ) THEN
      ELSEIF( IK .EQ. 3 ) THEN
      ELSEIF( IK .EQ. 4 ) THEN
      ELSEIF( IK .EQ. 5 ) THEN
        CALL FITSTA(5,
     >                FITING)
      ELSEIF( IK .EQ. 6 ) THEN
      ELSEIF( IK .EQ. 7 ) THEN
      ELSEIF( IK .EQ. 8 ) THEN
      ELSEIF( IK .EQ. 9 ) THEN
      ELSEIF( IK .EQ. 10) THEN
      ENDIF
 
      IF(.NOT.FITING) THEN 
C        WRITE(6,100) 
C     >  'Lmnt # '//TXTELT(1:FINSTR(TXTELT))//'-> Traj. #',
C     >  II,' stopped (IK=',IK,') : ',
C     >  TXT(DEBSTR(TXT):FINSTR(TXT)),' ;  remain/launched= ',
C     >    IMX-NSTOP,'/',IMX
        CALL ZGNOEL(
     >              NUML)
        WRITE(ABS(NRES),100) 
C     >  'LMNT # '//TXTELT(1:FINSTR(TXTELT))//'  -> Traj. #',
     >  'LMNT #',NUML,'   -> Traj. #',
     >  II,'  stopped  (IK=',IK,
     >  ') : '//TXT(DEBSTR(TXT):FINSTR(TXT))//' ;  remain/launched= ',
     >    IMX-NSTOP,'/',IMX
C 100    FORMAT(A,I6,A,I4,3A,I6,A1,I6)
 100    FORMAT(A,I4,A,I4,A,I2,A,I6,A,I6)
        call flush2(abs(nres),.false.) 
      ENDIF 

      IF(NSTOP.GE.IMX) THEN
        CALL FITSTA(5,
     >                FITING)
        IF(.NOT.FITING) THEN 
          WRITE(ABS(NRES),*) ' '
          WRITE(ABS(NRES),*) 'SBR KSTOP,  IK = ', IK 
          CALL ENDJOB
     >    ('SBR KSTOP : execution stopped, all particles lost !!',-99)
        ENDIF
      ENDIF

 99   RETURN 1

      ENTRY KSTOPI(LMNTI)
      LMNT = LMNTI
      RETURN
      END
