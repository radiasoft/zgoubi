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
      SUBROUTINE FITNU(*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXV=40) 
      COMMON/CONTR/ VAT(MXV),XI(MXV)
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      COMMON/FINIT/ FINI
      COMMON/INIT/ F0(6,6),F1(6,6),EX,EZ,X0,XP0,Z0,ZP0
     >,X1,XP1,Z1,ZP1,PHIX,PHIZ,IF
      COMMON/VAR/ X(3*MXV),P(MXV)
      COMMON/VARY/ NV,IR(MXV),NC,I1(MXV),I2(MXV),V(MXV),IS(MXV),W(MXV),
     >IC(MXV),IC2(MXV),I3(MXV),XCOU(MXV),CPAR(MXV,7)
      COMMON/VAIV/IV,IBREAK,IDEB
      COMMON/VISU/INV,INVHP
      EXTERNAL FF
 
      DIMENSION Y(MXV),VI(MXV)
      SAVE MTHD
      DATA MTHD / 2 / 

      CALL FITEST(*99,IER)
      CALL FITSET
      IF(IER .NE. 0) CALL ENDJOB('End of upon FIT test procedure',-99)
      CALL FITARR(IER)
      IF(IER .NE. 1) THEN
         CALL REMPLI(0)
         CALL FBORNE
         IF(MTHD .EQ. 2) THEN
C Implemented by Scott Berg, LPSC, April 2007
           CALL MINONM(NV,X,P,VI,XI,F,FINI)
         ELSE
           CALL MINO1(FF,NV,X,P,VI,Y,XI,F,FINI)
         ENDIF
         CALL IMPAJU(F)
      ENDIF
      CALL ENDFIT(F)
      RETURN
 99   RETURN 1
      ENTRY FITNU2(MTHDI) 
      MTHD = MTHDI
      END
