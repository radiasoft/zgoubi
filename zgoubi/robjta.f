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
C  Brookhaven National Laboratory                    és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE ROBJTA
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     *******************************************
C     READS DATA FOR OBJECT DEFINITION BY 'OBJET'
C     *******************************************
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      INCLUDE "C.OBJET.H"     ! COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
 
C     ....BORO
      READ(NDAT,*) A(NOEL,1)
C     ....IBODY, KOBJ
      READ(NDAT,*) A(NOEL,10),A(NOEL,11)
C     ....IMAX
      READ(NDAT,*) A(NOEL,20)
C     ....M1-M6
      READ(NDAT,*) (A(NOEL,I),I=30,35)
C     ....T1
      READ(NDAT,*) A(NOEL,40)
C     ....Yo-Do
      READ(NDAT,*) (A(NOEL,I),I=50,54)
C     ....dY-dD
      READ(NDAT,*) (A(NOEL,I),I=60,64)
C     ....XL
      READ(NDAT,*) A(NOEL,70)
C     ....IR1, IR2
      READ(NDAT,*) A(NOEL,80),A(NOEL,81)
 
      RETURN
      END
