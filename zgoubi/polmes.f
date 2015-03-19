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
      SUBROUTINE POLMES(SCAL,KUASEX,
     >                          BMIN,BMAX,BNORM,
     >                          XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA,NEWFIC)
      USE dynhc
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARIZ.H'
      INCLUDE "XYZHC.H"
C      COMMON//XH(MXX),YH(MXY),ZH(IZ),HC(ID,MXX,MXY,IZ,IMAP),IXMA,JYMA,KZMA
      INCLUDE "C.AIM_3.H"     ! COMMON/AIM/ ATO,AT,ATOS,RM,XI,XF,EN,EB1,EB2,EG1,EG2
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
C      INCLUDE "MAXTRA.H"
C      INCLUDE "C.CHAMBR.H"     ! COMMON/CHAMBR/ LIMIT,IFORM,YLIM2,ZLIM2,SORT(MXT),FMAG,YCH,ZCH
C 
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER(80) TA
      PARAMETER (MXTA=45)
      INCLUDE "C.DONT.H"     ! COMMON/DONT/ TA(MXL,MXTA)
      INCLUDE "C.DROITE_2.H"     ! COMMON/DROITE/ AM(9),BM(9),CM(9),IDRT
      INCLUDE "C.INTEG.H"     ! COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
C      LOGICAL ZSYM
      INCLUDE "C.TYPFLD.H"     ! COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      COMMON/ORDRES/ KORD,IRD,IDS,IDB,IDE,IDZ
 
      LOGICAL BINAR,BINARI,IDLUNI, NEWFIC
      CHARACTER(80) TITL , NOMFIC(IZ), NAMFIC
      SAVE NOMFIC, NAMFIC
 
      INTEGER DEBSTR,FINSTR
      PARAMETER (NHDF=4)
      LOGICAL STRCON
 
      DATA NOMFIC / IZ*'               '/
C FM 9 Spet. 14
C      DATA IMAP / 1 /
 
      PARAMETER (MXC = 4)
      DIMENSION AA(24+MXC-1)
      DATA AA / 27 * 0.D0 /
 
C FM 9 Spet. 14
C      CALL KSMAP(
C     >           IMAP)
 
      IF(KUASEX .EQ. 22) NDIM=2
 
C Possible SCAL change is by CAVITE
C Possible A(noel,10) change by FIT
C Apr 2013      BNORM = A(NOEL,10)*SCAL
      BNORM = A(NOEL,10)
      XNORM = A(NOEL,11)
      YNORM = A(NOEL,12)
      TITL = TA(NOEL,1)
      IXMA = A(NOEL,20)
      IF    (STRCON(TITL,'HEADER',
     >                            IS) ) THEN
        READ(TITL(IS+7:IS+7),FMT='(I1)') NHD
      ELSE
        NHD = NHDF
      ENDIF
 
      IF(IXMA.GT.MXX)
     >   CALL ENDJOB('X-dim of map is too large, max is ',MXX)
      JYMA = A(NOEL,21)
      IF(JYMA.GT.MXY )
     >   CALL ENDJOB('Y-dim of map is too large, max is ',MXY)
 
      IF(NDIM .EQ. 3) THEN
        KZMA =A(NOEL,22)
        IF(KZMA.GT.IZ )
     >   CALL ENDJOB('Z-dim of map is too large, max is ',IZ)
      ENDIF
 
      MOD = NINT(A(NOEL,22))
      MOD2 = NINT(10.D0*A(NOEL,22)) - 10*MOD
 
      IF    (NDIM.LE.2 ) THEN
        NFIC=1
        NAMFIC = TA(NOEL,2)
        NAMFIC = NAMFIC(DEBSTR(NAMFIC):FINSTR(NAMFIC))
        NEWFIC = NAMFIC .NE. NOMFIC(NFIC)
        NOMFIC(NFIC) = NAMFIC(DEBSTR(NAMFIC):FINSTR(NAMFIC))
      ELSEIF(NDIM .EQ. 3 ) THEN
        IF    (MOD .EQ. 0) THEN
C         ... 3-D map will be symmetrized wrt horizontal plane using SYMMED
          I1 = (KZMA/2) + 1
        ELSEIF(MOD .EQ. 1) THEN
C         ... No symm
          I1 = 1
          ZSYM=.FALSE.
        ENDIF
        NFIC=0
        NEWFIC = .TRUE.
        DO 129 I=I1,KZMA
          NFIC = NFIC+1
          NAMFIC = TA(NOEL,1+NFIC)
          NEWFIC = NEWFIC .AND. (NAMFIC .NE. NOMFIC(NFIC))
          NOMFIC(NFIC) = NAMFIC(DEBSTR(NAMFIC):FINSTR(NAMFIC))
 129    CONTINUE
      ENDIF
 
      CALL KSMAP4(NOMFIC,NFIC,AA(24:24+MXC-1),
     >                        NEWFIC,NBMAPS,IMAP)
 
      IF(NRES.GT.0) THEN
        IF(NEWFIC) THEN
           WRITE(NRES,209)
 209       FORMAT(/,10X,' A new field map is now used ; ',
     >     ' Name(s) of map data file(s) are : ')
C           WRITE(6   ,208) (NOMFIC(I),I=1,NFIC)
           WRITE(NRES,208) (NOMFIC(I),I=1,NFIC)
 208       FORMAT(10X,A)
        ELSE
          WRITE(NRES,210) (NOMFIC(I),I=1,NFIC)
 210      FORMAT(
     >    10X,'No  new  map  file  to  be  opened. Already  stored.',/
     >    10X,'Skip  reading  field  map  file : ',10X,A80)
        ENDIF
      ENDIF
 
      IF(KUASEX .EQ. 22 ) THEN
C------ POLARMES
C------ CARTE POLAIRE 2-D
 
        IF(NEWFIC) THEN
 
          IF(IDLUNI(
     >              LUN)) THEN
            BINAR=BINARI(NOMFIC(NFIC),IB)
            IF(BINAR) THEN
              OPEN(UNIT=LUN,FILE=NOMFIC(NFIC),FORM='UNFORMATTED'
     >        ,STATUS='OLD',ERR=96)
            ELSE
              OPEN(UNIT=LUN,FILE=NOMFIC(NFIC),STATUS='OLD',ERR=96)
           ENDIF
          ELSE
            GOTO 96
          ENDIF
 
          BMIN =  1.D10
          BMAX = -1.D10
          IRD = NINT(A(NOEL,40))
 
          CALL FMAPR(BINAR,LUN,MOD,MOD2,NHD,
     >                                 RM)
 
          DO J=1,JYMA
            DO I = 1,IXMA
              BFLD = HC(ID,I,J,1,IMAP)
              IF    (BFLD .GT. BMAX) THEN
                BMAX = BFLD
                XBMA = XH(I)
                YBMA = YH(J)
                ZBMA = 0D0
              ELSEIF(BFLD .LT. BMIN) THEN
                BMIN = BFLD
                XBMI = XH(I)
                YBMI = YH(J)
                ZBMI = 0D0
              ENDIF
            ENDDO
          ENDDO
 
          BMIN = BMIN * BNORM
          BMAX = BMAX * BNORM
          XBMA = XBMA * XNORM
          YBMA = YBMA * YNORM
          XBMI = XBMI * XNORM
          YBMI = YBMI * YNORM
 
          DO I = 1,IXMA
            XH(I)= XH(I) * XNORM
          ENDDO
          DO J=1,JYMA
            YH(J)= YH(J) * YNORM
          ENDDO
 
 
        ENDIF ! NEWFIC
 
      ENDIF
C----- KUASEX
 
      CALL CHAMK2(BNORM*SCAL)
 
      CALL MAPLI1(BMAX-BMIN)
 
      AT=XH(IXMA)-XH(1)
      ATO = 0D0
      ATOS = 0D0
C      RM=.5D0*(YH(JYMA)+YH(1))
 
      XI = XH(1)
      XF = XH(IXMA)
 
      RETURN
 
 96   WRITE(NRES,*) 'ERROR  OPEN  FILE ',NOMFIC(NFIC)
      CALL ENDJOB('ERROR  OPEN  FILE ',-99)
 
      END
