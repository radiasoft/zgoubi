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
      SUBROUTINE SVDOC(KLE,LABEL,
     >                           READAT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER(*) KLE(*)
      INCLUDE 'MXLD.H'
      CHARACTER(*) LABEL(MXL,2)
      LOGICAL READAT
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      PARAMETER (MXPUD=9,MXPU=5000)
      INCLUDE "C.CO.H"     ! COMMON/CO/ FPU(MXPUD,MXPU),KCO,NPU,NFPU,IPU
      PARAMETER (MCOLAB=5)
      PARAMETER (LBLSIZ=20)
      CHARACTER(LBLSIZ) COLAB
      INCLUDE 'C.COC.H'     ! COMMON/COC/ COLAB(MCOLAB)
      PARAMETER (MPULAB=5)
      CHARACTER(LBLSIZ) PULAB
      INCLUDE 'C.COT.H'     ! COMMON/COT/ PULAB(MPULAB)
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      INCLUDE "C.REBELO.H"   ! COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM

      LOGICAL OKCOR
      SAVE NBLM, OKCOR
      PARAMETER (T2KG = 10.D0)

      PARAMETER (IMON=MPULAB/3)
      PARAMETER(MXPUH =IMON, MXPUV =IMON)
      CHARACTER(LBLSIZ) HPNA(MXPUH), VPNA(MXPUV), HVPNA(MXPUV)
      CHARACTER(LBLSIZ) HPNAI(MXPUH), VPNAI(MXPUV), HVPNAI(MXPUV)
      SAVE HPNA, VPNA, HVPNA

      PARAMETER (MXCOH=5, MXCOV=5)
      CHARACTER(LBLSIZ) HCNA(MXCOH), VCNA(MXCOV)
      CHARACTER(LBLSIZ) HCNAI(MXCOH), VCNAI(MXCOV)
      SAVE HCNA, VCNA

      SAVE NLMC, NLM, KSCOR

      LOGICAL IDLUNI
C      SAVE MPU, MPUL
C      SAVE MCOH, MCOV, MCOL
      CHARACTER(70) TXFMT
      SAVE TXFMT
      SAVE LSVD

      LOGICAL FITRBL
      INTEGER(4) TODAY(3)

      SAVE NOELB

      DATA FITRBL / .FALSE. /
      DATA NOELA, NOELB / 1, MXL /

C Build SVD matrix:
C Scan all correctors. For each: 1/ change it 2/ find orbit 3/ store PUs

      CALL REBLT5(
     >            NOELA, NOELB)

      IF(IPASS.EQ.1) THEN

        KSCOR = 1                 ! H corr first
        NLMC = 0
        NLM = 0
        CALL ZGNBLM(
     >              NBLMI)
        NBLM = NBLMI

C Fill PULAB with PU names
        CALL SVDPUS(NBLM,HPNA,VPNA,HVPNA,
     >                             MPUL,MPUH,MPUV,MPUHV)
        MPU = MPUH+MPUV+MPUHV
Check corrector families
        CALL SVDCOS(NBLM,HCNA,VCNA,
     >                             MCOL,MCOH,MCOV)

        IF(IDLUNI(
     >            LSVD)) THEN
          OPEN(UNIT=LSVD,FILE='zgoubi.SVDtmp.out')
          CALL IDATE(TODAY)
          WRITE(LSVD,FMT='(''# '',2(I2.2,A1),I4.4,
     >    '' (dd-mm-yr). zgoubi.SVDtmp.out, from rebel.f/svdpr.f.'')')
     >    TODAY(1),'-',TODAY(2),'-',TODAY(3)
          WRITE(LSVD,FMT='(''# Columns: PU records 1 to '',I0
     >    ,'' (=#PUH+#PUV+2*#PUHV)'',
     >    '', ordering follows zgoubi.res optical sequence.'')')
     >    MPU+MPUHV
          WRITE(LSVD,FMT='(''# Rows: corrector 1 to '',I0,
     >    '', ordering follows zgoubi.res optical sequence.'')')
     >    MCOH+MCOV
          WRITE(LSVD,FMT='(
     >    ''# A total of '',I0,'' PUs ('',I0,''H, '',I0,''V, '',
     >    I0,''HV),'','' in '',I0,'' PU families : '',5(A,1x))') MPU,
     >    MPUH,MPUV,MPUHV,MPUL,(TRIM(PULAB(I)),I=1,MPUL)
          WRITE(LSVD,FMT='(''# and of '',I0,'' corrctrs ('',I0,
     >    '' H and '',I0,'' V),  in '',I0,'' families : ''
     >    ,5(A,1X))')
     >    MCOH+MCOV,MCOH,MCOV,MCOL,(TRIM(COLAB(I)),I=1,MCOL)
          WRITE(LSVD,FMT='(''# '',/,''# '',/,''# '')')

          WRITE(TXFMT,FMT='(A,I0,A,A)') '(1P,',
     >    MPU+MPUHV,'(E12.4,1X),','3(1X,I0),1X,A,1X,E14.6)'

        ELSE
          CALL ENDJOB('Pgm rebel. Could not open file '
     >    //'zgoubi.SVDtmp.out.',-99)
        ENDIF

        WRITE(ABS(NRES),
     >  FMT='(5X,''REBELOTE/SVD correction matrix requested:'',/,
     >  10X,''A total of '',I0,'' PUs, in '',I0,'' PU families : ''
     >  ,5(A,1x))') MPU,MPUL,(TRIM(PULAB(I)),I=1,MPUL)
        WRITE(ABS(NRES),FMT='(10X,''and of '',I0,'' corrctrs ('',I0,
     >  '' H and '',I0,'' V),  in '',I0,'' families : ''
     >  ,5(A,1X))')
     >  MCOH+MCOV,MCOH,MCOV,MCOL,(TRIM(COLAB(I)),I=1,MCOL)

      ELSE     ! IPASS.GT.1

        IF(NLMC.GE.1 .AND. NLM .GT. 0)
     >  CALL SVDPR(LSVD,TXFMT,NLM,NLMC,KLE,IPASS,LABEL)

      ENDIF

C Loop over corrector excitation, one-by-one. FIT finds the orbit each time.
      IF(KSCOR.LE.2) THEN
C There are ! 2 families of correctors at the moment

        IF(NLM .GT. 1 .AND. NLM1 .LT. NBLM) THEN
          A(NLM1,4) = ANLM1
        ENDIF

        DO WHILE ((.NOT. OKCOR) .AND. NLM .LT. NBLM)
C Move to next corrector. NLMC (1<NLMC<NBLM) is its number in the A() list
          NLM = NLM + 1
          OKCOR =
     >    (KSCOR .EQ. 1 .AND. LABEL(NLM,1).EQ.'HKIC')
     >    .OR.
     >    (KSCOR .EQ. 2 .AND. LABEL(NLM,1).EQ.'VKIC')
        ENDDO

        IF(NLM .GE. NBLM) THEN
          OKCOR=.FALSE.
          KSCOR = KSCOR + 1
          NLM = 0
        ENDIF

        IF(OKCOR) THEN
          NLMC = NLMC + 1
          OKCOR=.FALSE.
          IF(NLM .LE. NBLM) THEN
            NLM1 = NLM
            ANLM1 = A(NLM,4)
            IF    (KSCOR .EQ. 1) THEN
C              A(NLMC,4) = A(NOEL,10) / (A(NLMC,2)*1.D-2) * T2KG    ! B=(Brho==1)*kick/L
              A(NLM,4) = A(NOEL,10) / (A(NLM,2)*1.D-2) * T2KG    ! B=(Brho==1)*kick/L
            ELSEIF(KSCOR .EQ. 2) THEN
C              A(NLMC,4) =  A(NOEL,20) / (A(NLMC,2)*1.D-2) * T2KG
              A(NLM,4) =  A(NOEL,20) / (A(NLM,2)*1.D-2) * T2KG
            ELSE
              CALL ENDJOB(
     >        'SBR REBEL. NO SUCH POSSIBILITY KSCOR =',KSCOR)
            ENDIF
          ENDIF
        ENDIF

      ENDIF

C      IF(NRES .GT. 0) WRITE(NRES,100) IPASS
      WRITE(ABS(NRES),100) IPASS
 100  FORMAT(/,30X,'SVDOC. End of pass # ',
     >I0,' through the optical structure ',/)

      NRBLT = MCOH+MCOV+2
      IPASS=IPASS+1


      IF( IPASS .LE. NRBLT ) THEN

        READAT = .FALSE.
        NOEL=NOELA-1

      ELSEIF(IPASS .EQ. NRBLT+1) THEN

        READAT = .TRUE.
        KWRT = 1
        NRES = ABS(NRES)
        NOEL=NOELB

        CALL SVDINV(LSVD,MCOH+MCOV,MPU)

      ENDIF

      CALL SVDPC0

c      write(*,*) ' svdoc ',NOEL,noela-1,noelb,nrblt
c     > ,IPASS .EQ. NRBLT+1
C           read(*,*)

      RETURN

      ENTRY SVDOC2(HPNAI,VPNAI,HVPNAI,HCNAI,VCNAI)
!  HPU & HCorr name list, VPU & VCorr name list
      HPNA =       HPNAI
      VPNA =       VPNAI
      HVPNA =       HVPNAI
      HCNA =       HCNAI
      VCNA =       VCNAI
      RETURN

      END
