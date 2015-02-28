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
C  Upton, NY, 11973
C  -------
      SUBROUTINE CSRI
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXSTEP.H'
      INCLUDE 'MAXTRA.H'
      INCLUDE 'CSR.H'
C      COMMON/CSR/ KTRA,KCSR,YZXB(MXSTEP,41,36),DWC(MXT)
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      COMMON/PTICUL/ AM,Q,G,TO
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM

C      CHARACTER(9)   HMS

C      DATA LUN / 20 /
C      SAVE LUN
      
      KCSR = NINT(A(NOEL,1))
      KTRA = NINT(A(NOEL,2))
C Number of macro-particles in the bunch, evenly distributed in x 
      NPART = IMAX
 
        XMI = 1.D10
        XMA = -1.D10
        DO 1 IT = 1, IMAX
          XIT = FO(6,IT)
          IF(XIT.LT.XMI) XMI = XIT
          IF(XIT.GT.XMA) XMA = XIT
 1      CONTINUE
        CALL CSRIN2(XMI,XMA-XMI)

      IF(NRES.GT.0) THEN
        IF(KCSR .EQ. 0) THEN
          WRITE(NRES,107)
 107      FORMAT(/,15X,' CSR  interaction  is  off ',/)
          GOTO 99
        ELSE
          WRITE(NRES,110) 
 110      FORMAT(/,15X,' CSR  interaction  requested. ', 
     >    ' Reference  path  is  traj. #',
     >    //,25X,' A  first  run  will  track  paths,  ',
     >        'and  store  it. ',
     >    //,25X,' A  second  run  will  compute  the  perturbation  ',
     >        'on  the  particles.' )
          WRITE(NRES,FMT='(/,15X,'' Bunch  length  :'',
     >                                   1P,G12.4,'' cm'')') XMA-XMI
        ENDIF
      ENDIF
 
      IF(AM*Q .EQ. 0.D0) THEN
        IF(NRES .GT. 0) WRITE(NRES,106)
 106    FORMAT(//,15X,' Please  give   M  &  Q  of  particles !'
     >           ,/,15X,' - use  keyword  ''PARTICUL''',/)
          GOTO 99
      ENDIF

      IF(IPASS.EQ.1) THEN
        NRBLT = 1
C------- Switch on print into zgoubi.res : 
        ANOEL2 = 1.1D0
        KWRIT = INT(ANOEL2)
C------- Switch on print to standard output :
        KWRI6=INT((ANOEL2-KWRIT)*10)
C------- Set value of macro-charge that radiates
C The model assumes NPART particles uniformely spread along the bunch, 
C e.g., OBJET, KOBJ=2, DX: -Dxmin -> Dxmax, DXstep
           QB=Q*NPART
        CALL CSRIN1(QB/DBLE(NPART),NPART)
C------- DWC is the perturbation (total variation of kinetic energy due to CSR) 
        CALL RAZ(DWC,MXT)
      ENDIF

 99   RETURN
      END
