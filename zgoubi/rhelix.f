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
      SUBROUTINE RHELIX(ND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER TXT*132

      READ(NDAT,*) IA
      A(NOEL,1) = IA
C----- Length (cm), Twist-pitch=distance for 2pi B rotation (cm), BO (kG), initial angle of \vecB wrt vertical axis (rad)
      READ(NDAT,*) (A(NOEL,I),I=10,13)

C IRD2= 0 or 2,4,25  for analytic or 2-,4-,5-type numerical interpolation
C mesh size= XPAS/RESOL
C      READ(NDAT,*) A(NOEL,NP+1),A(NOEL,NP+2)
      READ(NDAT,FMT='(A132)') TXT
      READ(TXT,*) AA
      IF(INT(AA).NE.0) THEN
        READ(TXT,*) A(NOEL,20),A(NOEL,21)
      ELSE
        A(NOEL,20) = AA
      ENDIF
 
C----- Step size
      ND = 30
      CALL STPSIZ(NDAT,NOEL,ND,
     >                         A)
C----- KP,XCE, YCE, ALE
      ND1=ND+10
      READ(NDAT,*) IA,(A(NOEL,I),I=ND1+1,ND1+3)
      A(NOEL,ND1)=IA

      RETURN
      END
