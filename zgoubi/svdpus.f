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
      SUBROUTINE SVDPUS(NBLM,HPNA,VPNA,HVPNA,
     >                                 MPUL,MPUH,MPUV,MPUHV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     -----------------------------------------------------
C     Find PUs from in A() list. Put them in PULAB
C     -----------------------------------------------------
      CHARACTER(*) HPNA(*), VPNA(*), HVPNA(*)
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
      PARAMETER (IMON=MPULAB/3)
      PARAMETER(MXPUH =IMON, MXPUV =IMON, MXPUHV =IMON)
      LOGICAL DEJA
      DIMENSION KHV(MXPU)    ! 1, 2, 3 FOR H V, HV

      DATA OKPU / .FALSE. /

      OKPU = .FALSE.
      MPUL = 0      ! Numb of PU families/labels
      MPUH = 0       ! Numb of H PUs
      MPUV = 0
      MPUHV = 0
      MPU = 0        ! total numb of PUs
      NLM = 1
      DO WHILE ((.NOT. OKPU) .AND. NLM .LE. NBLM)
C Move to next element in sequence. Test whether its label identifies w/ PU label.
        I = 1
        DO WHILE((.NOT. OKPU) .AND. I.LE.MXPUH)
          OKPU = OKPU .OR. (LABEL(NLM,1).EQ.HPNA(I))
          IF(OKPU) THEN
            MPU = MPU + 1
            KHV(MPU) = 1
            MPUH = MPUH + 1
          ENDIF
          I = I + 1
        ENDDO
        I = 1
        DO WHILE((.NOT. OKPU) .AND. I.LE.MXPUV)
          OKPU = OKPU .OR. (LABEL(NLM,1).EQ.VPNA(I))
          IF(OKPU) THEN
            MPU = MPU + 1
            KHV(MPU) = 2
            MPUV = MPUV + 1
          ENDIF
          I = I + 1
        ENDDO
        I = 1
        DO WHILE((.NOT. OKPU) .AND. I.LE.MXPUHV)
          OKPU = OKPU .OR. (LABEL(NLM,1).EQ.HVPNA(I))
          IF(OKPU) THEN
            MPU = MPU + 1
            KHV(MPU) = 3
            MPUHV = MPUHV + 1
          ENDIF
          I = I + 1
        ENDDO

        IF(OKPU) THEN
          J = 1
          DEJA = .FALSE.
          DOWHILE (.NOT. DEJA .AND. J.LE. MPULAB)
            DEJA = DEJA .OR. (PULAB(J) .EQ. LABEL(NLM,1))
            J = J+1
          ENDDO
          IF(.NOT. DEJA ) THEN
            MPUL = MPUL + 1
            IF(MPUL .GT. MPULAB) CALL ENDJOB(
     >      'Pgm svdpus.  Too many PU families. Max allowed is ',MPULAB)
            PULAB(MPUL) = LABEL(NLM,1)
          ENDIF
          OKPU = .FALSE.
        ENDIF

        NLM = NLM + 1
      ENDDO

      OKPU=.FALSE.

      IF(MPUL.LE.0) CALL ENDJOB('Pgm svdpus.  None of the'
     >//' PU families listed under REBELOTE appears in the sequence. '
     >//'Check PU family names and existence of corresponding label_1.'
     >,-99)
      IF(MPUH+MPUV+MPUHV.GT.MXPU)
     >CALL ENDJOB('Pgm svdpus. Too many PUs.'
     >//' Need re-size FPU.  Max allowed is ',MXPU)

c      write(*,*) ' svdpus  MPUL, MPU = ',MPUL,MPU
c      write(*,fmt='(a,i0,2x,a)')
c     >  (' MPUL pulab : ',i, pulab(i),i=1,MPUL)
c      read(*,*)

      KCO = 2
      NPU=MPUL
      CALL SVDPR2(MPUH,MPUV,MPUHV,KHV)
      RETURN
      END
