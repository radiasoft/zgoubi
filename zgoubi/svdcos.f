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
      SUBROUTINE SVDCOS(NBLM,hcna,vcna,
     >                                 McoL,Mcoh,mcov)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     -----------------------------------------------------
C     Find COs from in A() list. Cot them in COLAB
C     -----------------------------------------------------
      character(*) hcna(*), vcna(*)
      INCLUDE 'C.CDF.H'         ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      PARAMETER (MCOLAB=5)
      PARAMETER (LBLSIZ=20)
      CHARACTER(LBLSIZ) coLAB
      INCLUDE 'C.COC.H'     ! COMMON/COC/ coLAB(McoLAB)
      INCLUDE 'MXLD.H'
      INCLUDE 'C.DON.H'     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER(LBLSIZ) LABEL
      INCLUDE 'C.LABEL.H'     ! COMMON/LABEL/ LABEL(MXL,2)
      INCLUDE 'C.REBELO.H'   ! COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM

      LOGICAL OKCO
      parameter(mxcoh =5, mxcov =5)
      logical deja

      DATA OKCO / .FALSE. /

      OKCO = .FALSE.
      MCOL = 0      ! Numb of CO families/labels
      MCOH = 0       ! Numb of H-corrs
      MCOV = 0                  ! Numb of V-corrs
      nlm = 1    
      DO WHILE ((.NOT. OKCO) .AND. nlm .LE. NBLM)
C Move to next element in sequence. Test whether its label identifies w/ CO label.

        i = 1
        do while((.NOT. OKCO) .AND. i.le.mxcoh)
          okco = okco .or. (LABEL(nlm,1).EQ.hcna(i))
          if(okco) mcoh = mcoh + 1
          i = i + 1
        enddo                 
        i = 1
        do while((.NOT. OKCO) .AND. i.le.mxcov)
          okco = okco .or. (LABEL(nlm,1).EQ.vcna(i))
          if(okco) mcov = mcov + 1
          i = i + 1
        enddo                 
          
        if(okco) then
          j = 1
          deja = .false.
          dowhile (.not. deja .and. j.le. MCOLAB)
            deja = deja .or. (colab(j) .eq. LABEL(nlm,1))
            j = j+1
          enddo
          if(.not. deja ) then
             MCOL = MCOL + 1
            if(mcol .gt. mcolab) call endjob(
     >      'Pgm svdcos.  Too many CO families. Max allowed is ',MCOLAB) 
            COLAB(MCOL) = LABEL(nlm,1)
          endif
          okco = .false.
        endif
       
        nlm = nlm + 1
      enddo      

      OKCO=.FALSE.

      IF(MCOL.LE.0) CALL ENDJOB('Pgm svdcos.  None of the'
     >//' CO families listed under REBELOTE appears in the sequence. '
     >//'Check CO family names and existence of corresponding label_1.'
     >,-99)
C      IF(MCO.GT.MXCO) CALL ENDJOB('Pgm svdcos. Too many COs.'
C     >//' Need re-size FCO.  Max allowed is ',MXCO)

      RETURN
      END
