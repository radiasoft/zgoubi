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
      SUBROUTINE PCKUP
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     -----------------------------------------------------
C     Pick-up signal (multiturn AND multiparticle summmation)
C     at labeled elements.
C     MPULAB = max number of LABEL's. MXPU = max number of 
C     pick-ups (virtual pick-ups, positionned at indicated 
C     labeled elements!) for CO measurments.
C     -----------------------------------------------------
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      PARAMETER (MXPUD=9,MXPU=5000)
      COMMON/CO/ FPU(MXPUD,MXPU),KCO,NPU,NFPU,IPU
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      PARAMETER (LBLSIZ=10)
      CHARACTER(LBLSIZ) LABEL
      COMMON /LABEL/ LABEL(MXL,2)
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM

      DIMENSION FPUL(MXPUD,MXPU), FPUL2(MXPUD-2,MXPU), NOELPU(MXPU)
      SAVE FPUL, FPUL2, NOELPU

      PARAMETER (I4 = 4)
      DIMENSION TBT(I4,MXPU),iqpu(mxl)
      SAVE TBT, iqpu

      PARAMETER (KSIZ=10)
      CHARACTER(KSIZ)  KLOBJ

      CHARACTER(KSIZ) KLEO
      LOGICAL CHKPU

      DATA NOELPU / MXPU*0 /
      DATA CHKPU / .TRUE. /

C----- Pick-up number. Reset to 0 via ENTRY PICKP2 below, by SBR PICKUP
      IPU = IPU + 1

      IF(IPASS.EQ.1) then
        NOELPU(IPU) = NOEL
        iqpu(noel) = ipu
      endif

      IF(IPU .GT. MXPU) 
     >CALL ENDJOB('SBR PCKUP.  Too many c.o. pick-ups,  max is ',MXPU)

      NT = 0
      DO 2 I = 1, IMAX
        TBT(1,IPU) = TBT(1,IPU) + F(2,I)
        TBT(2,IPU) = TBT(2,IPU) + F(3,I)
        TBT(3,IPU) = TBT(3,IPU) + F(4,I)
        TBT(4,IPU) = TBT(4,IPU) + F(5,I)
        IF( IEX(I) .GT. 0) THEN
          NT = NT+1
          DO 1 J = 1, 7
            FJI =  F(J,I)
            IF(J.EQ.1)  FJI = FJI -1.D0
            FPUL(J,IPU) = FPUL(J,IPU) + FJI
            FPUL2(J,IPU) = FPUL2(J,IPU) + FJI*FJI
 1        CONTINUE
        ENDIF
 2    CONTINUE
      FPUL(8,IPU) = NT

C Cumulated coordinates 
      DO J = 1, 5
        FPU(J,IPU) = FPU(J,IPU) + FPUL(J,IPU)
      ENDDO
C path length
      FPU(6,IPU) = FPUL(6,IPU)
C time
      FPU(7,IPU) = FPUL(7,IPU)
C turn #
      FPU(8,IPU) = FPU(8,IPU) + NT

C------- Record pick-up position (cm)
      IF(IPASS .EQ. 1) FPU(9,IPU) = F(6,1)

      RETURN

      ENTRY PCKUP1
      DO 4 I = 1, IPU
        NT = NINT(FPUL(8,I))
        CALL ZGKLE(IQ(NOELPU(I)),
     >                           KLEO)
        WRITE(NFPU,FMT= 
     >  '(1P,2X,I4,1X,E15.7,7(1X,E12.4),I9,I7,1X,7(1X,E12.4),2X,
     >  I5,3(1X,A))')
     >  I,FPU(9,I),(FPUL(J,I),J=2,6),FPUL(1,I),FPUL(7,I),NT,IPASS,
     >  (FPUL2(J,I),J=2,7), FPUL2(1,I), 
     >  NOELPU(I),KLEO,LABEL(NOELPU(I),1),LABEL(NOELPU(I),2)
 4    CONTINUE
      RETURN

      ENTRY PCKUP2
      IPU = 0
      CALL RAZ(TBT,I4*MXPU)
      CALL RAZ(FPUL,MXPUD*MXPU)
      CALL RAZ(FPUL2,(MXPUD-2)*MXPU)
      RETURN

      ENTRY PCKUP3(NOELI)
      IPU1 = 1
      DO IL = 1, NOELI
        IF(IPU1.LE.MXPU) THEN
          IF(NOELPU(IPU1).EQ.0) THEN
            GOTO 10
          ELSE 
            IF(NOELPU(IPU1).LT.NOELI) IPU1= IPU1+1
          ENDIF
        ELSE
          CALL ENDJOB('SBR PCKUP. Problem : IPU should be .le. ',MXPU)
        ENDIF
      ENDDO
 10   CONTINUE
      IF(NOELPU(IPU1).EQ.NOELI) THEN
        IPU = IPU1        
      ELSE
        IPU = IPU1-1
      ENDIF
      CALL RAZ(FPUL,MXPUD*MXPU)
      CALL RAZ(FPUL2,(MXPUD-2)*MXPU)
      IF(NRES.GT.0) THEN
        WRITE(NRES,*) 
        WRITE(NRES,*) 
     >   ' Found ',IPU,' pick-ups prior to element #',NOELI,'.'
        IF(IPU.GT.0) WRITE(NRES,*) ' Of which last PU is # ',
     >  IPU,' (located at  element NOEL = ',NOELPU(IPU),')'
      ENDIF
      RETURN

      ENTRY PCKUP5(MPU
     >                ,YCM,ZCM)
      IF(   FPUL(8,MPU) .LE.0 )
     >  CALL ENDJOB(' SBR pckup. All particles fucked up !',-99)
      YCM = TBT(1,MPU) /    FPUL(8,MPU) 
      ZCM = TBT(3,MPU) /    FPUL(8,MPU) 
      RETURN

      ENTRY PCKUP7(NOELA,NOELB, 
     >                         IPUI,IPUF,noeli,noelf)
        CALL ENDJOB('SBR PCKUP.  Must  have  NOELA .LE. NOELB.',-99)
C Range of PUs (1st and last PU's NOEL) located between NOELA and NOELB
c      ipu1 = 0
c      DO IL = NOELA, NOELB
c          CALL ZGKLE(IQ(NOELPU(I)),
c     >                             KLEO)
cC     >  NOELPU(I),KLEO,LABEL(NOELPU(I),1),LABEL(NOELPU(I),2)
c        IF(IL.LE.MXPU) THEN
c          IF(NOELPU(IL).EQ.0) THEN
c            GOTO 17
c          ELSE 
c            IF(IPU1 .EQ. 1) noeli = NOELPU(IPU1)
c            IF(NOELPU(IPU1).LT.NOELI) IPU1= IPU1+1
c          ENDIF
c        ELSE
c          CALL ENDJOB('SBR PCKUP. Problem : IL should be .le. ',MXPU)
c        ENDIF
c      ENDDO

c 17   CONTINUE
c      CALL ENDJOB('SBR PCKUP. Problem : found a PU with NOEL=0 ',-99)   
      chkPU = .true.
      noel = noela
      do while(chkPU .and. noel.le.noelb)
        if(iqpu(noel) .gt.0) then
          chkPU = .false.
        else
          noel = noel+1
        endif
      enddo
      noeli = noel
      ipui = iqpu(noel) 
      chkPU = .true.
      noel = noelb
      do while(chkPU .and. noel .gt. noela)
        if(iqpu(noel) .gt.0) then
          chkPU = .false.
        else
          noel = noel-1
        endif
      enddo
      noelf = noel
      ipuf = iqpu(noel) 
      RETURN
      END
