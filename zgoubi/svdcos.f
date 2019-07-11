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
      SUBROUTINE SVDCOS(NBLM,HCNA,VCNA,
     >                                 MCOL,MCOH,MCOV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     -----------------------------------------------------
C     Find COs from in A() list. Cot them in COLAB
C     -----------------------------------------------------
      CHARACTER(*) HCNA(*), VCNA(*)
      INCLUDE 'C.CDF.H'         ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      PARAMETER (MCOLAB=5)
      PARAMETER (LBLSIZ=20)
      CHARACTER(LBLSIZ) COLAB
      INCLUDE 'C.COC.H'     ! COMMON/COC/ COLAB(MCOLAB)
      INCLUDE 'MXLD.H'
      INCLUDE 'C.DON.H'     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER(LBLSIZ) LABEL
      INCLUDE 'C.LABEL.H'     ! COMMON/LABEL/ LABEL(MXL,2)
      INCLUDE 'C.REBELO.H'   ! COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM

      LOGICAL OKCO
      PARAMETER(MXCOH =5, MXCOV =5)
      LOGICAL DEJA

      DATA OKCO / .FALSE. /

      OKCO = .FALSE.
      MCOL = 0      ! Numb of c.o. families/labels
      MCOH = 0      ! Numb of H-corrs
      MCOV = 0      ! Numb of V-corrs
      NLM = 1
      DO WHILE ((.NOT. OKCO) .AND. NLM .LE. NBLM)
C Move to next element in sequence. Test whether its label identifies w/ C.O. label.

        I = 1
        DO WHILE((.NOT. OKCO) .AND. I.LE.MXCOH)
          OKCO = OKCO .OR. (LABEL(NLM,1).EQ.HCNA(I))
          IF(OKCO) MCOH = MCOH + 1
          I = I + 1
        ENDDO
        I = 1
        DO WHILE((.NOT. OKCO) .AND. I.LE.MXCOV)
          OKCO = OKCO .OR. (LABEL(NLM,1).EQ.VCNA(I))
          IF(OKCO) MCOV = MCOV + 1
          I = I + 1
        ENDDO

        IF(OKCO) THEN
          J = 1
          DEJA = .FALSE.
          DOWHILE (.NOT. DEJA .AND. J.LE. MCOLAB)
            DEJA = DEJA .OR. (COLAB(J) .EQ. LABEL(NLM,1))
            J = J+1
          ENDDO
          IF(.NOT. DEJA ) THEN
             MCOL = MCOL + 1
            IF(MCOL .GT. MCOLAB) CALL ENDJOB(
     >      'Pgm svdcos. Too many c.o. families. Must be .le. ',MCOLAB)
            COLAB(MCOL) = LABEL(NLM,1)
          ENDIF
          OKCO = .FALSE.
        ENDIF

        NLM = NLM + 1
      ENDDO

      OKCO=.FALSE.

      IF(MCOL.LE.0) CALL ENDJOB('Pgm svdcos.  None of the'
     >//' c.o. families listed under REBELOTE appears in the sequence. '
     >//'Check c.o. family names, and existence of corrspnding label_1.'
     >,-99)
C      IF(MCO.GT.MXCO) CALL ENDJOB('Pgm svdcos. Too many COs.'
C     >//' Need re-size FCO.  Max allowed is ',MXCO)

      RETURN
      END
