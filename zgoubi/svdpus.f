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
      SUBROUTINE SVDPUS(NBLM, )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     -----------------------------------------------------
C     Find PUs from in A() list. Put them in PULAB
C     -----------------------------------------------------
      INCLUDE 'C.CDF.H'     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      PARAMETER (MXPUD=9,MXPU=5000)
      INCLUDE 'C.CO.H'     ! COMMON/CO/ FPU(MXPUD,MXPU),KCO,NPU,NFPU,IPU
      PARAMETER (MPULAB=5)
      INCLUDE 'C.COT.H'     ! COMMON/COT/ PULAB(MPULAB)
      INCLUDE 'MXLD.H'
      INCLUDE 'C.DON.H'     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      PARAMETER (LBLSIZ=20)
      CHARACTER(LBLSIZ) PULAB
      INCLUDE 'C.REBELO.H'   ! COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM

      LOGICAL OKPU
      parameter (I5=5)
      parameter(mxpuh =I5, mxpuv =I5)
      DATA OKPU / .FALSE. /

      read(ta(noel,6),*) npuh
      read(ta(noel,25),*) npuv

      OKPU = .FALSE.
      DO WHILE ((.NOT. OKPU) .AND. NLMP .LE. NBLM)
C Move to next PU. NLMP (1<NLMP<NBLM) is its number in the A() sequence

        NLMP = NLMP + 1
        OKPU = .false.
        do while((.NOT. OKPU) .AND. i.le.mxpuh)
          okpu = okpu .or. LABEL(NLMP,1).EQ.ta(noel,1)
        enddo                 
        do while((.NOT. OKPU) .AND. i.le.mxpuv)
          okpu = okpu .or. LABEL(NLMP,1).EQ.ta(noel,20)
        enddo                 

        IF(.NOT. OKPU) CALL ENDJOB('Pgm svdpus. Could not find any PUs.'
     >  //' Check PU names ans existence of corresponding label_1.',-99)
              
        OKPU=.FALSE.
      
        DO WHILE 
          PULAB(I), I=1,NPU)


           
      RETURN
      END
