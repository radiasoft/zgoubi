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
      SUBROUTINE RBB
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL

c          intensity=a(noel,2)
      READ(NDAT,*) (A(NOEL,I),I=1,2)
c          alfx=a(noel,10)
c          betx=a(noel,11)
c          epsnx=a(noel,12)
      READ(NDAT,*) (A(NOEL,I),I=10,12)
c          alfy=a(noel,20)
c          bety=a(noel,21)
c          epsny=a(noel,22)
      READ(NDAT,*) (A(NOEL,I),I=20,22)
c          sigz=a(noel,30)
c          dpp=a(noel,31)
      READ(NDAT,*) (A(NOEL,I),I=30,31)
c          circ=a(noel,40)
c          alfmom=a(noel,41)
      READ(NDAT,*) (A(NOEL,I),I=40,41)
c          tunex=a(noel,50)
c          tuney=a(noel,51)
c          tunez=a(noel,52)
      READ(NDAT,*) (A(NOEL,I),I=50,52)
c          ampx=a(noel,60)
c          ampy=a(noel,61)
c          ampz=a(noel,62)
      READ(NDAT,*) (A(NOEL,I),I=60,62)

      RETURN
      END
