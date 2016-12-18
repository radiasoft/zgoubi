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
      SUBROUTINE TOSCAC(SCAL,NDIM,
     >                          BMIN,BMAX,BNORM,XNORM,YNORM,ZNORM,
     >               XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA,AA,NEWFIC)
      USE DYNHC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-------------------------------------------------
C     Read field map with cartesian coordinates.
C     TOSCA keyword with MOD.le.19.
C-------------------------------------------------
      LOGICAL NEWFIC
      INCLUDE 'PARIZ.H'
      INCLUDE "XYZHC.H"
      INCLUDE "C.AIM_3.H"     ! COMMON/AIM/ ATO,AT,ATOS,RM,XI,XF,EN,EB1,EB2,EG1,EG2
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
C      PARAMETER (LNTA=132) ; CHARACTER(LNTA) TA
C      PARAMETER (MXTA=45)
      INCLUDE "C.DONT.H"     ! COMMON/DONT/ TA(MXL,MXTA)
      INCLUDE "C.DROITE_2.H"     ! COMMON/DROITE/ AM(9),BM(9),CM(9),IDRT
      INCLUDE "C.INTEG.H"     ! COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      PARAMETER (LBLSIZ=20)
      CHARACTER(LBLSIZ) LABEL
      INCLUDE "C.LABEL.H"     ! COMMON/LABEL/ LABEL(MXL,2)
C      LOGICAL ZSYM
      INCLUDE "C.TYPFLD.H"     ! COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
C      INCLUDE "C.ORDRES.H"     ! COMMON/ORDRES/ KORD,IRD,IDS,IDB,IDE,IDZ
      INCLUDE 'MXFS.H'
      INCLUDE "C.SCALP.H"     ! COMMON/SCALP/ VPA(MXF,MXP),JPA(MXF,MXP)
 
      LOGICAL BINARI,IDLUNI
      LOGICAL BINAR
      LOGICAL FLIP
C      CHARACTER(80) TITL , NOMFIC(IZ), NAMFIC
      PARAMETER (MXC = 4)
      PARAMETER (NFM = IZ+MXC)
      CHARACTER(LNTA) TITL , NOMFIC(NFM), NAMFIC, NOMFI1(1)
      SAVE NOMFIC, NAMFIC
      INTEGER DEBSTR,FINSTR
      PARAMETER (NHDF=8)
 
      LOGICAL STRCON
 
      CHARACTER(20) FMTYP
      DIMENSION XXH(MXX,MMAP), YYH(MXY,MMAP), ZZH(IZ,MMAP)
      SAVE XXH, YYH, ZZH
      DIMENSION BBMI(MMAP), BBMA(MMAP), XBBMI(MMAP), YBBMI(MMAP)
      DIMENSION ZBBMI(MMAP), XBBMA(MMAP), YBBMA(MMAP), ZBBMA(MMAP)
      SAVE BBMI, BBMA, XBBMI, YBBMI, ZBBMI, XBBMA, YBBMA, ZBBMA
      DIMENSION IIXMA(MMAP), JJYMA(MMAP), KKZMA(MMAP)
      SAVE IIXMA, JJYMA, KKZMA
 
      DIMENSION DBDX(3)
 
      DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE ::
     >     HCA, HCB, HCC
C      DIMENSION HCA(ID,MXX,MXY,IZ),HCB(ID,MXX,MXY,IZ),HCC(ID,MXX,MXY,IZ)
      SAVE HCA, HCB, HCC
 
      INCLUDE 'MXSCL.H'
      DIMENSION KFM(MXSCL)
 
C      PARAMETER (MXAA2=24+MXC-1)
      PARAMETER (MXMAP=4)
      PARAMETER (MXAA2=MAX(24+MXC-1,20+10*MXMAP+9))
      DIMENSION AA(MXL,MXAA2), UU(MXAA2)

C     16/01/14 to pass the map coefficients to KSMAP4
      PARAMETER (ONE=1.D0)
      PARAMETER (MXHD=20)

C--    ERRORS
      LOGICAL ERRON
      SAVE ERRON
      PARAMETER (MXERR=MXTA)
      CHARACTER(LBLSIZ) LBL1(MXERR), LBL2(MXERR)
      CHARACTER(LBLSIZ) LBL1I, LBL2I
      SAVE LBL1, LBL2
      LOGICAL EMPTY
      PARAMETER (JPOL=6)
      DIMENSION KPOL(MXERR,JPOL), BNRM(JPOL)
      CHARACTER(2) TYPERR(MXERR,JPOL)
      CHARACTER(1) TYPAR(MXERR,JPOL),TYPDIS(MXERR,JPOL)
      CHARACTER(2) TYPERI
      CHARACTER(1) TYPAI,TYPDII
      DIMENSION ERRCEN(MXERR,JPOL),ERRSIG(MXERR,JPOL),ERRCUT(MXERR,JPOL)
      SAVE TYPERR,TYPAR,TYPDIS,ERRCEN,ERRSIG,ERRCUT
      DIMENSION DB(MXL,JPOL),DPOS(MXL,JPOL,3),TILT(MXL,JPOL,3)
      SAVE DB, DPOS, TILT
      LOGICAL OK
      LOGICAL FITING, FITFNL
      LOGICAL PRNT, PRNTI
      SAVE PRNT
      PARAMETER (KSIZ=10)
      CHARACTER(KSIZ) KLEY
      LOGICAL OKOPN
      SAVE OKOPN, LERR
      SAVE IPOL
      CHARACTER(LBLSIZ) LBL1l, LBL2l

      LOGICAL ZROBXY      
      PARAMETER (IONE = 1)

      DIMENSION NMPT(NOEL,MXC),NMPTN(MXC)
      CHARACTER(80) TXFM
      CHARACTER(1) TXT1

      DATA ZROBXY / .FALSE. /
      DATA NOMFIC / NFM*'               '/
      DATA FMTYP / ' regular' /
      DATA IALOC / 0 /
      DATA NBMAPS / 0 /
      DATA ERRON / .FALSE. /
      DATA BNRM / 6*0.D0 /
      DATA NMPTN / MXC * 1 /
      DATA NFIC / 1 /
 
      IF( .NOT.ALLOCATED( HCA ))
     >     ALLOCATE( HCA(ID,MXX,MXY,IZ), STAT = IALOC)
      IF (IALOC /= 0)
     >     CALL ENDJOB('SBR toscac Not enough memory'//
     >     ' for Malloc of HC',
     >     -99)
 
      IF( .NOT.ALLOCATED( HCB ))
     >     ALLOCATE( HCB(ID,MXX,MXY,IZ), STAT = IALOC)
      IF (IALOC /= 0)
     >     CALL ENDJOB('SBR toscac Not enough memory'//
     >     ' for Malloc of HC',
     >     -99)
 
      IF( .NOT.ALLOCATED( HCC ))
     >     ALLOCATE( HCC(ID,MXX,MXY,IZ), STAT = IALOC)
      IF (IALOC /= 0)
     >     CALL ENDJOB('SBR toscac Not enough memory'//
     >     ' for Malloc of HC',
     >     -99)
 
C Possible SCAL change is by CAVITE
C Possible A(noel,10) change by FIT
C Aug 2012      BNORM = A(NOEL,10)*SCAL
      BNORM = A(NOEL,10)
      XNORM = A(NOEL,11)
      YNORM = A(NOEL,12)
      ZNORM = A(NOEL,13)
      TITL = TA(NOEL,1)
      IF    (STRCON(TITL,'ZroBXY',
     >                            IS) ) THEN
        ZROBXY=.TRUE.
      ELSE
        ZROBXY=.FALSE.
      ENDIF
      IF    (STRCON(TITL,'HEADER',
     >                            IS) ) THEN
        READ(TITL(IS+7:IS+7),FMT='(I1)') NHD
      ELSE
        NHD = NHDF
      ENDIF
      IF(NHD .GT. MXHD) CALL ENDJOB(
     >'Pgm toscac. Field map header has too many lines. Must be .le.'
     >,MXHD)
C      IDEB = DEBSTR(TITL)
C      FLIP = TITL(IDEB:IDEB+3).EQ.'FLIP'
      FLIP = STRCON(TITL,'FLIP',
     >                          IS)
      IF(FLIP)
     >  CALL ENDJOB('SBR TOSCAC. FLIP option not implemented.',-99)
      IXMA = NINT(A(NOEL,20))
      IF(IXMA.GT.MXX)
     >   CALL ENDJOB('X-dim of map is too large,  max  is ',MXX)
      IF(NDIM .EQ. 1) THEN
        JYMA=1
        KZMA=1
      ELSE
        JYMA = NINT(A(NOEL,21))
        IF(JYMA.GT.MXY )
     >     CALL ENDJOB('Y-dim of map is too large,  max  is ',MXY)
        IF(NDIM .EQ. 3) THEN
          KZMA =NINT(A(NOEL,22))
          IF(KZMA.GT.IZ )
     >       CALL ENDJOB('Z-dim of map is too large,  max  is ',IZ)
        ELSE
          KZMA=1
        ENDIF
      ENDIF
 
      MOD = NINT(A(NOEL,23))
      MOD2 = NINT(10.D0*A(NOEL,23)) - 10*MOD
 
      IF    (NDIM.EQ.2 ) THEN
        IF    (MOD .EQ. 0) THEN
          I1=1
          I2 = 1
          NFIC=1
          NAMFIC = TA(NOEL,2)
          NOMFIC(NFIC) = NAMFIC(DEBSTR(NAMFIC):FINSTR(NAMFIC))
 
          DO LL = 24, MXAA2  ! 24+MXC-1
            UU(LL) = AA(NOEL,LL)
          ENDDO
          CALL KSMAP4(NOMFIC,NFIC,UU,
     >                          NEWFIC,NBMAPS,IMAP)
 
        ELSEIF(MOD .EQ. 3) THEN
          I1=1
          I2 = 1
          NFIC=1
          NAMFIC = TA(NOEL,2)
          NOMFIC(NFIC) = NAMFIC(DEBSTR(NAMFIC):FINSTR(NAMFIC))
 
          DO LL = 24, MXAA2     !24+MXC-1
            UU(LL) = AA(NOEL,LL)
          ENDDO
          CALL KSMAP4(NOMFIC,NFIC,UU,
     >                          NEWFIC,NBMAPS,IMAP)
 
          IF(MOD2 .EQ. 1) THEN
C TOSCA 2D map for the AGS main magnet
C dB1, dB2, dB3
            AA(NOEL,24) = A(NOEL,24)
            AA(NOEL,25) = A(NOEL,25)
            AA(NOEL,26) = 0.D0
            DUM = SCALE9(
     >                   KFM)
            DO IFM = 1, MXSCL
c            IF(KFM .GT. 0) THEN
              IF(KFM(IFM) .LE. 0) GOTO 20
                DO I = 1, JPA(KFM(IFM),MXP)
C Apply scaling to all parameters concerned
                 IF(JPA(KFM(IFM),I) .GT. MXAA2) CALL ENDJOB(
     >           'Pgm toscac. Exceeded AA size. JPA = ',JPA(KFM(IFM),I))
                 AA(NOEL,JPA(KFM(IFM),I))=
     >               AA(NOEL,JPA(KFM(IFM),I))*VPA(KFM(IFM),I)
              ENDDO
c            ENDIF
            ENDDO
 
 20         CONTINUE
 
            DBDX(1) = AA(NOEL,24)
            DBDX(2) = AA(NOEL,25)
            DBDX(3) = AA(NOEL,26)
            CALL CHAMK4(DBDX,3)
          ENDIF

        ELSEIF(MOD.EQ.15) THEN
C--------- MOD2 files are combined linearly into a single 2D map, after reading.
          IF(MOD2 .GT. MXC) CALL ENDJOB('Pgm toscac. Number of field '
     >    //'maps cannot exceed ',MXC)
          IF(MOD2.LT.1) MOD2 = 1
          I1=1
          I2 = MOD2
          DO I = 1, I2
            AA(NOEL,24-1+I) = A(NOEL,24-1+I)
          ENDDO
 
          DO LL = 24, MXAA2     !24+MXC-1
            UU(LL) = AA(NOEL,LL)
          ENDDO

          NFIC=0
          DO I=I1, I2
            NFIC = NFIC+1
            NAMFIC = TA(NOEL,1+NFIC)
            NOMFIC(NFIC) = NAMFIC(DEBSTR(NAMFIC):FINSTR(NAMFIC))
          ENDDO
 
C Store NOMFIC(1:NFIC) name series, and increment IMAP -> IMAP+1
          CALL KSMAP4(NOMFIC,NFIC,UU,
     >                          NEWFIC,NBMAPS,IMAP)

        ELSE
          CALL ENDJOB('*** Error. SBR TOSCAC. No such option MOD= ',MOD)
        ENDIF
 
      ELSEIF(NDIM .EQ. 3 ) THEN
 
        IF(MOD .EQ. 0) THEN
C--------- Several data files, normally one per XY plane
C          3-D map is  symmetrysed wrt horizontal plane
          I1 = (KZMA/2) + 1
          I2 = KZMA

        ELSEIF(MOD .EQ. 1) THEN
C--------- Several data files, normally one per XY plane
C--------  No symmetrization, map taken as is
          I1 = 1
          I2 = KZMA
          ZSYM=.FALSE.

        ELSEIF(MOD .EQ. 12) THEN
C--------- A single data file contains the all 3D volume
          I1=1
          I2 = 1

        ELSEIF(MOD .EQ. 15) THEN
C--------- MOD2 files are combined linearly into a single map after reading.
C          Each one of these files should contain the all 3D volume.
          IF(MOD2 .GT. MXC) CALL ENDJOB('Pgm toscac. Number of field '
     >    //'maps cannot exceed ',MXC)
          I1=1
          I2 = MOD2
          DO I = I1, I2
            AA(NOEL,24-1+I) = A(NOEL,24-1+I)
          ENDDO

        ELSEIF(MOD .EQ. 16) THEN
C--------- MOD2 files read for on-the-flight contribution
C          Each one of these files should contain the all 3D volume.
          IF(MOD2 .GT. MXC) CALL ENDJOB('Pgm toscac. Number of field '
     >    //'maps cannot exceed ',MXC)
          I1=1
          I2 = MOD2
          DO I = I1, I2
C Field scaling coefficients (1 per map)
            AA(NOEL,24-1+I) = A(NOEL,24-1+I)
          ENDDO
          DO I = I1, I2
C Map positioning data (1 set per map)
            DO J = 1, 5
              AA(NOEL,20+10*I +J-1) = A(NOEL,20+10*I +J-1)
            ENDDO
          ENDDO

        ELSE
          CALL ENDJOB('*** Error. SBR TOSCAC. No such option MOD= ',MOD)
        ENDIF
 
        DO LL = 24, MXAA2
          UU(LL) = AA(NOEL,LL)
        ENDDO

        NFIC=0
        DO  I = I1, I2
          NFIC = NFIC+1
          NAMFIC = TA(NOEL,1+NFIC)
          NOMFIC(NFIC) = NAMFIC(DEBSTR(NAMFIC):FINSTR(NAMFIC))
        ENDDO

C Store NOMFIC(1:NFIC) name series, and increment IMAP -> IMAP+1

        IF    (MOD .LE. 15) THEN 
          CALL KSMAP4(NOMFIC,NFIC,UU,
     >                          NEWFIC,NBMAPS,IMAP)


c              write(*,*) ' toscac MOD 15, newfic = ',newfic

        ELSEIF(MOD .EQ. 16) THEN 
          CALL KSMAP(
     >               IMP)
          NEW = IMP
          DO I = 1, NFIC
            NOMFI1(1) = NOMFIC(I)
            UU(24) = UU(24+I-1)
            CALL KSMAP4(NOMFI1, IONE ,UU,
     >                          NEWFIC,NBMAPS,IMAP)
            NMPT(NOEL,I) = IMAP
            IF(IMAP .GT. NEW) NEW = NEW +1
          ENDDO
          NEW = NEW - IMP
          DO I = 1, MXC
            NMPTN(I) = NMPT(NOEL,I)
          ENDDO

C          CALL CHAMKW(NFIC,NMPTN,UU)

        ENDIF ! MOD  .le.15, .eq.16
   
      ENDIF  ! NDIM.EQ.3

      IF    (MOD .EQ. 15) THEN
        ICOF = 24
        IFIC = I1
        DO WHILE (IFIC.LE.I2 .AND. .NOT. NEWFIC)
          IF(ICOF .GT. 27)
     >    CALL ENDJOB('SBR toscac. No such possibility ICOF=',ICOF)
          NEWFIC = NEWFIC .OR. AA(NOEL,ICOF).NE.A(NOEL,ICOF)


c              write(*,*) ' toscac test aa newfic = ',newfic
c              write(*,*) ' aa  noel,ICOF : ',AA(NOEL,ICOF),noel,ICOF
c              write(*,*) ' a ICOF  ', A(NOEL,ICOF),ICOF
c              write(*,*) ' +++++++++++++++++'
c              write(*,*) ' '
c              write(*,*) ' '
c                 read(*,*)


          ICOF = ICOF + 1
          IFIC = IFIC + 1
        ENDDO
        ICOF = 24
        DO IFIC = I1, I2
          AA(NOEL,ICOF) = A(NOEL,ICOF)
          ICOF = ICOF + 1
        ENDDO

      ELSEIF(MOD .EQ. 16) THEN
        ICOF = 24
        IFIC = I1
        DO WHILE (IFIC.LE.I2 .AND. .NOT. NEWFIC)
          IF(ICOF .GT. 27)
     >    CALL ENDJOB('SBR toscac. No such possibility ICOF=',ICOF)
          NEWFIC = NEWFIC .OR. AA(NOEL,ICOF).NE.A(NOEL,ICOF)
          ICOF = ICOF + 1
          IFIC = IFIC + 1
        ENDDO
        ICOF = 24
        DO IFIC = I1, I2
          AA(NOEL,ICOF) = A(NOEL,ICOF)
          ICOF = ICOF + 1
        ENDDO
          DO I = 1, I2
C Map positioning data (1 set per map)
            DO J = 1, 5
              AA(NOEL,20+10*I +J-1) = A(NOEL,20+10*I +J-1)
            ENDDO
          ENDDO

          CALL CHAMKW(NFIC,NMPTN,UU)

      ENDIF
 
      IF(NRES.GT.0) THEN
        WRITE(NRES,FMT='(/,5X,3(A,I0),2A,I3,A)')
     >  ,'Status of MOD.MOD2 is ', MOD,'.',MOD2,' ; NDIM = ',NDIM,' ;'
     >  ,' number of field data files used is ',NFIC,'.'
        IF    (MOD.EQ.3) THEN
          WRITE(NRES,FMT='(5X,A,I0,A)') 
     >    'Field map set stored in field array HC(IMAP = ',IMAP,')' 
          IF(MOD2.EQ.1) THEN
            WRITE(NRES,FMT='(/,5X,A,I1,2A,1P,2E15.6,/)')
     >      'MOD2 = ',MOD2,' ->  map will be perturbed using '
     >      ,'field indices db/dx, d2b/dx2 : ',dbdx(1),dbdx(2)
          ENDIF
c        ELSEIF(MOD.EQ.15) THEN
        ELSEIF(MOD.EQ.15) THEN
          WRITE(NRES,FMT='(5X,A,I0,A)') 
     >    'Field map set stored in field array HC(IMAP = ',IMAP,')' 
          IF(I2.GT.1) THEN
            WRITE(NRES,FMT='(5X,2(A,I2),2A,1P,/,15X,4E14.6)')
     >      '3-D map. MOD=',MOD,'.   Will sum up ',I2,'  field maps, ',
     >      'with field coefficient values : ',(a(noel,24+i-1),i=i1,i2)
          ELSE
            WRITE(NRES,FMT='(5X,A,I2,2A,1P,/,15X,4E14.6)')
     >      '3-D map. MOD=',MOD,'.  Single field map, ',
     >      'with field coefficient value : ',(a(noel,24+i-1),i=i1,i2)
          ENDIF
          IF(ZROBXY) THEN
            WRITE(NRES,FMT='(/,5X,A)')
     >      'ZroBXY OPTION :  found in map title. '//
     >      'BX and BY values in Z=0 plane forced to zero.'
          ELSE
            WRITE(NRES,FMT='(/,5X,A)')
     >      'ZroBXY OPTION :  absent from map title. '//
     >      'BX and BY values in Z=0 plane left as read.'
          ENDIF
        ELSEIF(MOD .EQ. 16) THEN
          IF(I2.EQ.1) THEN
            WRITE(NRES,FMT='(5X,A,I0,A)') 
     >      'Single 3-D field map, stored in field array HC(IMAP = '
     >      ,IMAP,')'
            WRITE(NRES,FMT='(5X,A,4E14.6)')
     >      'Field coefficient will apply : ',(a(noel,24+i-1),i=i1,i2)
          ELSE
            WRITE(NRES,FMT='(5X,2(A,I0),A)') 
     >      'New 3-D field map(s) stored in field array(s) HC(IMAP = ',
     >      IMP+1,'-',IMP+NEW,')' 
            WRITE(NRES,FMT='(5X,A,I2,2A,1P,/,15X,4E14.6)')
     >      'Will sum up field from ',I2,
     >      ' maps, on-the-flight at particle location. ',
     >      'Field coefficients will apply : ',(a(noel,24+i-1),i=i1,i2)
          ENDIF
          IF(ZROBXY) THEN
            WRITE(NRES,FMT='(/,5X,A)')
     >      'ZroBXY OPTION :  found in map title. '//
     >      'BX and BY values in Z=0 plane forced to zero.'
          ELSE
            WRITE(NRES,FMT='(/,5X,A)')
     >      'ZroBXY OPTION :  absent from map title. '//
     >      'BX and BY values in Z=0 plane left as read.'
          ENDIF
        ENDIF

        IF(NEWFIC) THEN
           WRITE(NRES,209)
 209       FORMAT(/,10X
     >     ,' New field map(s) opened, cartesian mesh (MOD.le.19) ; '
     >     ,/,10X,' name(s) of map data file(s) : ',/)
           WRITE(NRES,208)  (NOMFIC(I)(DEBSTR(NOMFIC(I)):
     >     FINSTR(NOMFIC(I))),I=1,NFIC)
 208       FORMAT(10X,A)
        ELSE
          WRITE(NRES,210)  (NOMFIC(I)(DEBSTR(NOMFIC(I)):
     >     FINSTR(NOMFIC(I))),I=1,NFIC)
 210      FORMAT(/,
     >    5X,'No  map  file  opened :  was  already  stored.',/
     >    5X,'Skip  reading  field  map  file : ',10X,A)
        ENDIF
        CALL FLUSH2(NRES,.FALSE.)
      ENDIF
 
      IF(NEWFIC) THEN
 
         IF    (MOD .LT. 15) THEN
           NFIC = 0
           DO 12 KZ=I1,I2
             NFIC = NFIC+1
             IF(IDLUNI(
     >                 LUN)) THEN
               BINAR=BINARI(NOMFIC(NFIC),IB)
               IF(BINAR) THEN
                 OPEN(UNIT=LUN,FILE=NOMFIC(NFIC),FORM='UNFORMATTED'
     >           ,STATUS='OLD',ERR=96)
               ELSE
                 OPEN(UNIT=LUN,FILE=NOMFIC(NFIC),STATUS='OLD',ERR=96)
               ENDIF
             ELSE
               GOTO 96
             ENDIF
 
             IF(     STRCON(NOMFIC(NFIC),'GSI',
     >                                         IS)
     >          .OR. STRCON(NOMFIC(NFIC),'gsi',
     >                                         IS)) THEN
                  NHD = 0
                  FMTYP = 'GSI'
             ELSEIF(     STRCON(NOMFIC(NFIC),'BW6',
     >                                             IS)
     >          .OR. STRCON(NOMFIC(NFIC),'bw6',
     >                                         IS)) THEN
                  NHD = 0
                  FMTYP = 'GSI'
             ENDIF
 
             IF(NRES.GT.0) THEN
               LNGTH=len(
     >         NOMFIC(NFIC)(DEBSTR(NOMFIC(NFIC)):FINSTR(NOMFIC(NFIC))))
              WRITE(NRES,FMT='(2(/,2X,A),I4,A,I0,A,/)') 
     >        ' ----',' Map file number ',NFIC,' ( of ',(I2-I1+1),') :'
               WRITE(NRES,FMT='(5X,3A)')
     >         NOMFIC(NFIC)(DEBSTR(NOMFIC(NFIC)):LNGTH)//' ',
     >         'map,  FORMAT type : ',FMTYP(DEBSTR(FMTYP):FINSTR(FMTYP))
               CALL FLUSH2(NRES,.FALSE.)
             ENDIF
 
C             IRD = NINT(A(NOEL,40))
 
C BNORM set to ONE, since sent to CHAMK below
             CALL FMAPR3(BINAR,LUN,MOD,MOD2,NHD,ZROBXY,
     >                   XNORM,YNORM,ZNORM,ONE,I1,KZ,FMTYP,
     >                                    BMIN,BMAX,
     >                                    XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)
 
 12        CONTINUE
 
C------- Store mesh coordinates
           CALL FMAPW4(IMAP,BMIN,BMAX,XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)
 
         ELSEIF(MOD .EQ. 15) THEN
 
           NFIC = 0
           DO KZ=I1,I2
             NFIC = NFIC+1
             IF(IDLUNI(
     >                 LUN)) THEN
               BINAR=BINARI(NOMFIC(NFIC),IB)
               IF(BINAR) THEN
                 OPEN(UNIT=LUN,FILE=NOMFIC(NFIC),FORM='UNFORMATTED'
     >           ,STATUS='OLD',ERR=96)
               ELSE
                 OPEN(UNIT=LUN,FILE=NOMFIC(NFIC),STATUS='OLD',ERR=96)
               ENDIF
             ELSE
               GOTO 96
             ENDIF
 
             IF(     STRCON(NOMFIC(NFIC),'GSI',
     >                                         IS)
     >          .OR. STRCON(NOMFIC(NFIC),'gsi',
     >                                         IS)) THEN
                  NHD = 0
                  FMTYP = 'GSI'
             ELSEIF(     STRCON(NOMFIC(NFIC),'BW6',
     >                                             IS)
     >          .OR. STRCON(NOMFIC(NFIC),'bw6',
     >                                         IS)) THEN
                  NHD = 0
                  FMTYP = 'GSI'
             ENDIF
 
             ICOF = 23+NFIC

             IF(NRES.GT.0) THEN
               LNGTH=len(
     >         NOMFIC(NFIC)(DEBSTR(NOMFIC(NFIC)):FINSTR(NOMFIC(NFIC))))
               WRITE(NRES,FMT='(2(/,2X,A),I4,A,I0,A,/)') 
     >         ' ----',' Map file number ',NFIC,' ( of ',(I2-I1+1),') :'
               WRITE(NRES,FMT='(5X,A,1P,E16.8,/)')
     >         NOMFIC(NFIC)(DEBSTR(NOMFIC(NFIC)):LNGTH)//' '//
     >         'map, FORMAT type : '//FMTYP(DEBSTR(FMTYP):FINSTR(FMTYP))
     >         //'.  Field multiplication factor :',AA(NOEL,ICOF)
               CALL FLUSH2(NRES,.FALSE.)
             ENDIF
 
C             IRD = NINT(A(NOEL,40))
 
             CALL FMAPW2(NFIC,AA(NOEL,ICOF))
 
             CALL FMAPR3(BINAR,LUN,MOD,MOD2,NHD,ZROBXY,
     >                   XNORM,YNORM,ZNORM,ONE,I1,KZ,FMTYP,
     >                                    BMIN,BMAX,
     >                                    XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)
 
             IF    (NFIC.EQ.1) THEN
 
               IF(MOD2.EQ.1) THEN
                 DO KKK=1,KZMA
                  DO JJJ=1,JYMA
                   DO III=1,IXMA
                     DO IID = 1, ID
                       HC(IID,III,JJJ,KKK,IMAP) =
     >                         AA(NOEL,24)*HC(IID,III,JJJ,KKK,IMAP)
                     ENDDO
                   ENDDO
                  ENDDO
                 ENDDO
 
               ELSE
 
                 DO KKK=1,KZMA
                  DO JJJ=1,JYMA
                   DO III=1,IXMA
                    DO IID = 1, ID
                      HCA(IID,III,JJJ,KKK)=AA(NOEL,24)
     >                      *HC(IID,III,JJJ,KKK,IMAP)
                    ENDDO
                   ENDDO
                  ENDDO
                 ENDDO
  
               ENDIF
               CLOSE(UNIT=LUN)
 
             ELSEIF(NFIC.EQ.2) THEN
 
               IF(MOD2.EQ.2) THEN
                 DO KKK=1,KZMA
                  DO JJJ=1,JYMA
                   DO III=1,IXMA
                     DO IID = 1, ID
                       HC(IID,III,JJJ,KKK,IMAP) = HCA(IID,III,JJJ,KKK)
     >                  +   AA(NOEL,25) * HC(IID,III,JJJ,KKK,IMAP)
                     ENDDO
                   ENDDO
                  ENDDO
                 ENDDO

               ELSE
 
                 DO KKK=1,KZMA
                  DO JJJ=1,JYMA
                   DO III=1,IXMA
                     DO IID = 1, ID
                       HCB(IID,III,JJJ,KKK) = HCA(IID,III,JJJ,KKK)
     >                  +   AA(NOEL,25) * HC(IID,III,JJJ,KKK,IMAP)
                     ENDDO
                   ENDDO
                  ENDDO
                 ENDDO
 
               ENDIF
               CLOSE(UNIT=LUN)
 
             ELSEIF(NFIC.EQ.3) THEN
 
               IF(MOD2.EQ.3) THEN
                 DO KKK=1,KZMA
                  DO JJJ=1,JYMA
                   DO III=1,IXMA
                     DO IID = 1, ID
                       HC(IID,III,JJJ,KKK,IMAP) = HCB(IID,III,JJJ,KKK)
     >                  +   AA(NOEL,26) * HC(IID,III,JJJ,KKK,IMAP)
                     ENDDO
                   ENDDO
                  ENDDO
                 ENDDO
                 CLOSE(UNIT=LUN)
 
               ELSE
                 DO KKK=1,KZMA
                  DO JJJ=1,JYMA
                   DO III=1,IXMA
                     DO IID = 1, ID
                       HCC(IID,III,JJJ,KKK) = HCB(IID,III,JJJ,KKK)
     >                  +   AA(NOEL,26) * HC(IID,III,JJJ,KKK,IMAP)
                     ENDDO
                   ENDDO
                  ENDDO
                 ENDDO
 
               endif
               CLOSE(UNIT=LUN)
 
             ELSEIF(NFIC.EQ.4) THEN
 
                 DO KKK=1,KZMA
                  DO JJJ=1,JYMA
                   DO III=1,IXMA
                     DO IID = 1, ID
                       HC(IID,III,JJJ,KKK,IMAP) = HCC(IID,III,JJJ,KKK)
     >                  +   AA(NOEL,27) * HC(IID,III,JJJ,KKK,IMAP)
                     ENDDO
                   ENDDO
                  ENDDO
                 ENDDO
                 CLOSE(UNIT=LUN)
 
             ENDIF
             CLOSE(UNIT=LUN)
 
           ENDDO

C------- Store mesh coordinates
           CALL FMAPW4(IMAP,BMIN,BMAX,XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)
 
         ELSEIF(MOD .EQ. 16) THEN
 
           NFIC = 0
           DO KZ=I1,I2
             NFIC = NFIC+1

             IF(IDLUNI(
     >                 LUN)) THEN
               BINAR=BINARI(NOMFIC(NFIC),IB)
               IF(BINAR) THEN
                 OPEN(UNIT=LUN,FILE=NOMFIC(NFIC),FORM='UNFORMATTED'
     >           ,STATUS='OLD',ERR=96)
               ELSE
                 OPEN(UNIT=LUN,FILE=NOMFIC(NFIC),STATUS='OLD',ERR=96)
               ENDIF
             ELSE
               GOTO 96
             ENDIF
 
             ICOF = 23+NFIC

             IF(NRES.GT.0) THEN
               LNGTH=len(
     >         NOMFIC(NFIC)(DEBSTR(NOMFIC(NFIC)):FINSTR(NOMFIC(NFIC))))
               WRITE(NRES,FMT='(2(/,2X,A),I4,A,I0,A,/)') 
     >         ' ----',' Map file number ',NFIC,' ( of ',(I2-I1+1),') :'
               WRITE(NRES,FMT='(5X,A,1P,E16.8,/,5x,a,4(e14.6,1x),2i4/)')
     >         NOMFIC(NFIC)(DEBSTR(NOMFIC(NFIC)):LNGTH)//' '//
     >         'map, FORMAT type : '//FMTYP(DEBSTR(FMTYP):FINSTR(FMTYP))
     >         //'.  Field multiplication factor :',AA(NOEL,23+NFIC),
     >         'Positioning (X, Y, theta, dY) : ',
     >         (AA(NOEL,JJ),JJ=20+10*KZ+1,20+10*KZ+4),
     >             20+10*KZ, 20+10*KZ+4
               CALL FLUSH2(NRES,.FALSE.)
             ENDIF
 
C             IRD = NINT(A(NOEL,40))
 
             CALL FMAPW2(NFIC,AA(NOEL,ICOF))
 
             CALL FMAPR3(BINAR,LUN,MOD,MOD2,NHD,ZROBXY,
     >                   XNORM,YNORM,ZNORM,ONE,I1,KZ,FMTYP,
     >                                    BMIN,BMAX,
     >                                    XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA) 
             CLOSE(UNIT=LUN)
 
           ENDDO

C--------- Store mesh coordinates
           CALL FMAPW4(IMAP,BMIN,BMAX,XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)

 
         ELSE
           CALL ENDJOB('SBR toscac. No such option MOD = ',MOD)
         ENDIF
 
 
      ELSE   !  .not. NEWFIC
 
C           IRD = NINT(A(NOEL,40))
 
C------- Restore mesh coordinates
           CALL FMAPR5(IMAP,
     >                   BMIN,BMAX,XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA) 
 
           IF(NRES.GT.0) THEN
             WRITE(NRES,fmt='(5X,A,I3,2A)') 
     >       'Restored mesh coordinates for field map #IMAP= ',imap,
     >       ',  name : ',
     >       NOMFIC(NFIC)(DEBSTR(NOMFIC(NFIC)):FINSTR(NOMFIC(NFIC)))
           ENDIF
 
      ENDIF ! NEWFIC
 
C----- Case erron (errors)
      IF(ERRON) THEN

        DO IRR = 1, MXERR 
          OK = (EMPTY(LBL1(IRR)) .OR. LBL1(IRR).EQ.LABEL(NOEL,1)) 
     >    .AND.(EMPTY(LBL2(IRR)) .OR. LBL2(IRR).EQ.LABEL(NOEL,2)) 

          IF(OK) THEN
            CALL REBELR(
     >                  KREB3,KDUM,KDUM)
            CALL FITSTA(5,
     >                    FITING)
            CALL FITST5(
     >                  FITFNL)         

C            IF(.NOT.FITING .AND. .NOT. FITFNL .AND. (KREB3.NE.99)) THEN
            IF(.NOT.FITING .AND. .NOT. FITFNL .AND. (KREB3.NE.99)) 
C Won't go if KREB3=99, since this is multi-turn in same lattice. 
     >        CALL TOSERR(NOEL,IRR,MXTA,BNRM, 
     >        KPOL,TYPERR,TYPAR,TYPDIS,ERRCEN,ERRSIG,ERRCUT,
     >                                             DB,DPOS,TILT)

            DO IM=1,JPOL
              IF(KPOL(IRR,IM).EQ.1) THEN
C                BNRM(IM) = BNRM(IM) + DB(NOEL,IM)
                BNRM(IM) = DB(NOEL,IM)
              ENDIF
            ENDDO

            BNORM = BNRM(1)


            IF(PRNT .AND. OKOPN) THEN
              CALL ZGKLE(IQ(NOEL), 
     >                           KLEY)
              CALL ZGIPAS( 
     >                    IPASS,NRBLT)
              IF(EMPTY(LBL1(IRR))) THEN 
                LBL1L = '.'
              ELSE
                LBL1L = LBL1(IRR)
              ENDIF
              IF(EMPTY(LBL2(IRR))) THEN 
                LBL2L = '.'
              ELSE
                LBL2L = LBL1(IRR)
              ENDIF
              WRITE(TXT1,FMT='(I1)') JPOL
              TXFM='(6(1X,I5),3(1X,A),3(1X,E16.8), 1X,E16.8,'//
     >          '6(1X,E16.8), 3(1X,A))' 
                WRITE(LERR,FMT=TXFM)
     >          IPASS,NOEL,KREB3,IRR,IPOL,KPOL(IRR,IPOL),
     >          TYPERR(IRR,IPOL), TYPAR(IRR,IPOL),TYPDIS(IRR,IPOL),
     >          ERRCEN(IRR,IPOL),ERRSIG(IRR,IPOL),ERRCUT(IRR,IPOL),
     >               DB(NOEL,IPOL),
     >          (DPOS(NOEL,IPOL,JJ),JJ=1,3),(TILT(NOEL,IPOL,JJ),JJ=1,3),
     >          KLEY,LBL1L,LBL2L
              ENDIF

          ENDIF      
        ENDDO

        IF(NRES.GT.0) THEN
          DO IRR = 1, MXERR 
            OK = (EMPTY(LBL1(IRR)) .OR. LBL1(IRR).EQ.LABEL(NOEL,1)) 
     >      .AND.(EMPTY(LBL2(IRR)) .OR. LBL2(IRR).EQ.LABEL(NOEL,2)) 
            IF(OK) THEN
              IF(.NOT.FITING .AND. .NOT.FITFNL .AND. (KREB3.NE.99)) THEN
                DO I = 1, JPOL
                  IF(KPOL(IRR,I) .EQ. 1) THEN 
                    WRITE(NRES,FMT='(/,5X,
     >              ''ERRORS are set : '')')
                    WRITE(NRES,FMT='(10X,
     >              ''- case of TOSCA key with labels : '',4A,I4,A,I4)')
     >              LBL1(IRR), ' / ',LBL2(IRR),' /  ERROR SET # ',IRR,
     >              ', # element = ',NOEL  
                    WRITE(NRES,FMT='(10X,
     >              ''- pole#, error type, A/R, G/U : '',I1,3(2X,A))')
     >              I, TYPERR(IRR,I), TYPAR(IRR,I), TYPDIS(IRR,I)
                    WRITE(NRES,FMT='(10X,A,1P,4(E14.6,2X),/)') 
     >              '- err_center, _sigma, _cutoff, '
     >              //'resulting BNORM : ',  
     >              ERRCEN(IRR,I),ERRSIG(IRR,I),ERRCUT(IRR,I),
     >              BNRM(I)
                  ENDIF
                ENDDO
              ELSE
                DO I = 1, JPOL
                  IF(KPOL(IRR,I) .EQ. 1) THEN 
                    WRITE(NRES,FMT='(/,5X,
     >              ''ERRORS found set'', 
     >              '' - maintainded unchanged so far : '')')
                    WRITE(NRES,FMT='(10X,
     >              ''- case of TOSCA key with labels : '',4A,I4,A,I4)')
     >              LBL1(IRR), ' / ',LBL2(IRR),' /  ERROR SET # ',IRR,
     >              ', # element = ',NOEL  
                    WRITE(NRES,FMT='(10X,
     >              ''- pole#, error type, A/R, G/U : '',I1,3(2X,A))')
     >              I, TYPERR(IRR,I), TYPAR(IRR,I), TYPDIS(IRR,I)
                    WRITE(NRES,FMT='(10X,A,1P,4(E14.6,2X),/)') 
     >              '- err_center, _sigma, _cutoff, '
     >              //'resulting BNORM : ',  
     >              ERRCEN(IRR,I),ERRSIG(IRR,I),ERRCUT(IRR,I),
     >              BNRM(I)
                  ENDIF
                ENDDO
              ENDIF
            ENDIF
          ENDDO
        ENDIF

      ENDIF

      CALL CHAMK2(BNORM*SCAL)
 
C Make sure this is ok with cartésien
        CALL MAPLI1(BMAX-BMIN)
C        AT=XH(IAMA)-XH(1)
C        ATO = 0.D0
C        ATOS = 0.D0
C        RM=.5D0*(YH(JRMA)+YH(1))
C        XI = XH(1)
C        XF = XH(IAMA)
 
      RETURN

 96   WRITE(ABS(NRES),*) 'Pgm toscac. Error  open  file ',
     >NOMFIC(NFIC)(DEBSTR(NOMFIC(NFIC)):FINSTR(NOMFIC(NFIC)))
      CALL ENDJOB('Leaving. ',-99)
      RETURN
 
      ENTRY TOSCA1
      DO I = 1, IZ
        NOMFIC(I) = ' '
      ENDDO
      RETURN
 
      ENTRY TOSCA2(IRRI,IPOLI,TYPERI,TYPAI,TYPDII,
     >               ERRCEI,ERRSII,ERRCUI,LBL1I,LBL2I)
      ERRON = .TRUE.
      IRR = IRRI
      IPOL = IPOLI
      KPOL(IRR,IPOL) = 1
      TYPERR(IRR,IPOL)=      TYPERI
      TYPAR(IRR,IPOL)=      TYPAI
      TYPDIS(IRR,IPOL)=      TYPDII
      ERRCEN(IRR,IPOL)=      ERRCEI
      ERRSIG(IRR,IPOL)=      ERRSII
      ERRCUT(IRR,IPOL)=      ERRCUI
      LBL1(IRR) = LBL1I
      LBL2(IRR) = LBL2I
      RETURN

      ENTRY TOSCA4
      ERRON = .FALSE.
      RETURN      

      ENTRY TOSCA8(PRNTI)
      PRNT = PRNTI
      IF(PRNT) THEN
        IF(.NOT. OKOPN) THEN
          OK = IDLUNI(
     >                 LERR)
          OPEN(UNIT=LERR,FILE='zgoubi.ERRORS.out')
          WRITE(LERR,FMT='(A,I0,A)') '# Origin of this print : toscac'
     >    //' program. This file opened with unit # ',LERR,'.'
          WRITE(LERR,FMT='(A)') '# '
          WRITE(LERR,FMT='(A)') 
     >    '# IPASS, NOEL, KREB3, IRR, IPOL, KPOL(IRR,Ipol) '
     >    //'TYPERR(IRR,Ipol), TYPAR(IRR,IPOL), TYPDIS(IRR,IPOL), '
     >    //'ERRCEN(IRR,IPOL), ERRSIG(IRR,IPOL), ERRCUT(IRR,IPOL), '
     >    //'DB(NOEL,IPOL), '
     >    //'(DPOS(NOEL,IPOL,JJ),JJ=1,3), (TILT(NOEL,IPOL,JJ),JJ=1,3), '
     >    //'KLEY, LBL1L, LBL2L'
          WRITE(LERR,FMT='(A)') '# '
          OKOPN = .TRUE.
        ENDIF
      ELSE
        IF(OKOPN) CLOSE(LERR)
      ENDIF
      RETURN      
      END
