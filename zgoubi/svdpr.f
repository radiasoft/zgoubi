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
      SUBROUTINE SVDPR(LSVD,TXFMT,NLM,NLMC,KLE,IPASS,LABEL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER(*) KLE(*), TXFMT
      INCLUDE 'MXLD.H'
      CHARACTER(*) LABEL(MXL,2)
      PARAMETER (MXPUD=9,MXPU=5000)
      INCLUDE 'C.CO.H'     ! COMMON/CO/ FPU(MXPUD,MXPU),KCO,NPU,NFPU,IPU
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL

      DIMENSION KHV(MXPU), KHVI(MXPU)     ! 1, 2, 3 FOR H V, HV
      SAVE MPU,MPUH,MPUV,MPUHV,KHV
      DIMENSION TMP(MXPU)

      JPU = 0
      DO I = 1, MPU
        IF    (KHV(I) .EQ. 1) THEN 
          JPU = JPU + 1 
          TMP(JPU) = FPU(2,I)
        ELSEIF(KHV(I) .EQ. 2) THEN 
          JPU = JPU + 1 
          TMP(JPU) = FPU(4,I)
        ELSE
          JPU = JPU + 2
          TMP(JPU-1) = FPU(2,I)
          TMP(JPU) = FPU(4,I)
        ENDIF
      ENDDO

      IF(JPU .NE. MPU+MPUHV) CALL ENDJOB(
     >'Pgm  svdpr. Problem with numbering of PU I/O',-99) 

            WRITE(LSVD,FMT=TXFMT)
     >      (TMP(J),J=1,JPU),NLMC,IPASS-1,NLM,
     >      TRIM(KLE(IQ(NLM)))//'['//TRIM(LABEL(NLM,1))//']',A(NLM,4)
            CALL FLUSH2(LSVD,.FALSE.)

      
      RETURN

      ENTRY SVDPR2(MPUHI,MPUVI,MPUHVI,KHVI)
      MPUH = MPUHI
      MPUV = MPUVI
      MPUHV = MPUHVI
      MPU = MPUH + MPUV + MPUHV 
      KHV = KHVI
      RETURN
      END

      
