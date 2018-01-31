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
      SUBROUTINE SVDPUS(NBLM,hpna,vpna,
     >                                 MPUL,MPU)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     -----------------------------------------------------
C     Find PUs from in A() list. Put them in PULAB
C     -----------------------------------------------------
      character(*) hpna(*), vpna(*)
      INCLUDE 'C.CDF.H'         ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      PARAMETER (MXPUD=9,MXPU=5000)
      INCLUDE 'C.CO.H'     ! COMMON/CO/ FPU(MXPUD,MXPU),KCO,NPU,NFPU,IPU
      PARAMETER (MPULAB=5)
      PARAMETER (LBLSIZ=20)
      CHARACTER(LBLSIZ) PULAB
      INCLUDE 'C.COT.H'     ! COMMON/COT/ PULAB(MPULAB)
      INCLUDE 'MXLD.H'
      INCLUDE 'C.DON.H'     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER(LBLSIZ) LABEL
      INCLUDE 'C.LABEL.H'     ! COMMON/LABEL/ LABEL(MXL,2)
      INCLUDE 'C.REBELO.H'   ! COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM

      LOGICAL OKPU
      parameter (IMON=MPULAB/2)
      parameter(mxpuh =IMON, mxpuv =IMON)
      logical deja

      DATA OKPU / .FALSE. /

c      READ(TA(NOEL,6),*) NPUH
c      READ(TA(NOEL,25),*) NPUV

      OKPU = .FALSE.
      MPUL = 0      ! Numb of PU families/labels
      MPU = 0       ! Numb of PUs
      nlm = 1    
      DO WHILE ((.NOT. OKPU) .AND. nlm .LE. NBLM)
C Move to next element in sequence. Test whether its label identifies w/ PU label.
        i = 1
        do while((.NOT. OKPU) .AND. i.le.mxpuh)
          okpu = okpu .or. (LABEL(nlm,1).EQ.hpna(i))
          i = i + 1
        enddo                 
        i = 1
        do while((.NOT. OKPU) .AND. i.le.mxpuv)
          okpu = okpu .or. (LABEL(nlm,1).EQ.vpna(i))
          i = i + 1
        enddo                 
          
        if(okpu) then
          MPU = MPU + 1
          j = 1
          deja = .false.
          dowhile (.not. deja .and. j.le. MPULAB)
            deja = deja .or. (pulab(j) .eq. LABEL(nlm,1))
            j = j+1
          enddo
          if(.not. deja ) then
            MPUL = MPUL + 1
            if(mpul .gt. mpulab) call endjob(
     >      'Pgm svdpus.  Too many PU families. Max allowed is ',MPULAB) 
            PULAB(MPUL) = LABEL(nlm,1)
          endif
          okpu = .false.
        endif
       
        nlm = nlm + 1
      enddo      

      OKPU=.FALSE.

      IF(MPUL.LE.0) CALL ENDJOB('Pgm svdpus.  None of the'
     >//' PU families listed under REBELOTE appears in the sequence. '
     >//'Check PU family names and existence of corresponding label_1.'
     >,-99)
      IF(MPU.GT.MXPU) CALL ENDJOB('Pgm svdpus. Too many PUs.'
     >//' Need re-size FPU.  Max allowed is ',MXPU)

c      write(*,*) ' svdpus  MPUL, MPU = ',MPUL,MPU
c      write(*,fmt='(a,i0,2x,a)')
c     >  (' MPUL pulab : ',i, pulab(i),i=1,MPUL)
c      read(*,*)

      RETURN
      END
