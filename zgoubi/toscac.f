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
      SUBROUTINE TOSCAC(SCAL,NDIM,
     >                          BMIN,BMAX,BNORM,XNORM,YNORM,ZNORM,
     >               XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA,NEWFIC)
      USE dynhc
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
      CHARACTER(80) TA
      PARAMETER (MXTA=45)
      INCLUDE "C.DONT.H"     ! COMMON/DONT/ TA(MXL,MXTA)
      INCLUDE "C.DROITE_2.H"     ! COMMON/DROITE/ AM(9),BM(9),CM(9),IDRT
      INCLUDE "C.INTEG.H"     ! COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
C      LOGICAL ZSYM
      INCLUDE "C.TYPFLD.H"     ! COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      INCLUDE "C.ORDRES.H"     ! COMMON/ORDRES/ KORD,IRD,IDS,IDB,IDE,IDZ
      INCLUDE 'MXFS.H'
      INCLUDE "C.SCALP.H"     ! COMMON/SCALP/ VPA(MXF,MXP),JPA(MXF,MXP)
 
      LOGICAL BINARI,IDLUNI
      LOGICAL BINAR
      LOGICAL FLIP
      CHARACTER(80) TITL , NOMFIC(IZ), NAMFIC
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
      dimension kfm(MXSCL)
 
      DATA NOMFIC / IZ*'               '/
      DATA FMTYP / ' regular' /
 
      PARAMETER (MXC = 4)
      DIMENSION AA(MXL,24+MXC-1), UU(24+MXC-1)
c      DATA AA / MXL*(24+MXC-1)* 0.d0 /
      SAVE AA
C     16/01/14 to pass the map coefficients to KSMAP4
      PARAMETER (ONE=1.D0)
      INCLUDE 'MAPHDR.H'

      DATA IALOC / 0 /
      DATA NBMAPS/ 0 /
 
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
 
      if( .NOT.ALLOCATED( HCC ))
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
      IF    (STRCON(TITL,'HEADER',
     >                            IS) ) THEN
        READ(TITL(IS+7:IS+7),FMT='(I1)') NHD
      ELSE
        NHD = NHDF
      ENDIF
      IF(NHD .GT. MXHD) CALL ENDJOB(
     >'Pgm toscac. Field map header has too many lines. Must be .le.'
     >,MXHD)
      IDEB = DEBSTR(TITL)
C      FLIP = TITL(IDEB:IDEB+3).EQ.'FLIP'
      FLIP = STRCON(TITL,'FLIP',
     >                          IS)
      IF(FLIP)
     >  CALL ENDJOB('SBR TOSCAC. FLIP option not implemented.',-99)
      IXMA = A(NOEL,20)
      IF(IXMA.GT.MXX)
     >   CALL ENDJOB('X-dim of map is too large,  max  is ',MXX)
      IF(NDIM .EQ. 1) THEN
        JYMA=1
        KZMA=1
      ELSE
        JYMA = A(NOEL,21)
        IF(JYMA.GT.MXY )
     >     CALL ENDJOB('Y-dim of map is too large,  max  is ',MXY)
        IF(NDIM .EQ. 3) THEN
          KZMA =A(NOEL,22)
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
 
          DO LL = 24, 24+MXC-1
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
 
          DO LL = 24, 24+MXC-1
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
            CALL SCALE9(
     >                   KFM)
            do ifm = 1, MXSCL
c            IF(KFM .GT. 0) THEN
            IF(KFM(ifm) .le. 0) goto 20
              DO I = 1, JPA(KFM(IFM),MXP)
C Apply scaling to all parameters concerned
c            write(*,*) ' toscac ',
c     >        I,KFM, JPA(KFM,I), AA(NOEL,JPA(KFM,I)) , VPA(KFM,I)
c                read(*,*)
                AA(NOEL,JPA(KFM(IFM),I))=
     >                AA(NOEL,JPA(KFM(IFM),I))*VPA(KFM(IFM),I)
              ENDDO
c            ENDIF
            enddo
 
 20         continue
 
            DBDX(1) = AA(NOEL,24)
            DBDX(2) = AA(NOEL,25)
            DBDX(3) = AA(NOEL,26)
            CALL CHAMK4(DBDX,3)
          ENDIF
        ELSEIF(MOD.EQ.15) THEN
C--------- MOD2 files are combined linearly into a single 2D map, after reading.
          AA(NOEL,24) = A(NOEL,24)
          AA(NOEL,25) = A(NOEL,25)
          AA(NOEL,26) = A(NOEL,26)
          AA(NOEL,27) = A(NOEL,27)
          I1=1
          I2 = MOD2
 
          NFIC=0
          DO I=I1, I2
            NFIC = NFIC+1
            NAMFIC = TA(NOEL,1+NFIC)
            NOMFIC(NFIC) = NAMFIC(DEBSTR(NAMFIC):FINSTR(NAMFIC))

            DO LL = 24, 24+MXC-1
              UU(LL) = AA(NOEL,LL)
            ENDDO
 
            CALL KSMAP4(NOMFIC,NFIC,UU,
     >                          NEWFIC,NBMAPS,IMAP)
          ENDDO
 
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
          AA(NOEL,24) = A(NOEL,24)
          AA(NOEL,25) = A(NOEL,25)
          AA(NOEL,26) = A(NOEL,26)
          AA(NOEL,27) = A(NOEL,27)
          I1=1
          I2 = MOD2

c      call ZGNOEL(
c     >             NOEL)
c       if(noel.eq.89) then
c           write(*,*) ' toscac aa ', (a(noel,i),i=24,27)
c          endif

        ELSE
          CALL ENDJOB('*** Error. SBR TOSCAC. No such option MOD= ',MOD)
        ENDIF
 
        NFIC=0
C        NEWF = .TRUE.
        DO 129 I=I1, I2
          NFIC = NFIC+1
          NAMFIC = TA(NOEL,1+NFIC)
          NOMFIC(NFIC) = NAMFIC(DEBSTR(NAMFIC):FINSTR(NAMFIC))
 129    CONTINUE

        DO LL = 24, 24+MXC-1
          UU(LL) = AA(NOEL,LL)

c      call ZGNOEL(
c     >             NOEL)
c       if(noel.eq.89) then
c           write(*,*) ' toscac uu ',ll,uu(ll)
c          endif

        ENDDO
 
        CALL KSMAP4(NOMFIC,NFIC,UU,
     >                          NEWFIC,NBMAPS,IMAP)
 
c      call ZGNOEL(
c     >             NOEL)
c       if(noel.eq.89) then
c           write(*,*) ' toscac  nomfic ',(i1,i2,NOMFIC(iii),iii=i1,i2)
c               write(*,*) ' toscac  ',NFIC,UU,
c     >                          NEWFIC,NBMAPS,IMAP
c                 read(*,*)
c          endif
 
       IF(MOD .EQ. 15) THEN
          IFAC = 24
          IFIC = 1
          DO WHILE (IFIC.LE.I2 .AND. .NOT. NEWFIC)
            IF(IFAC .GT. 27)
     >      CALL ENDJOB('SBR toscac. No such possibility IFAC=',IFAC)
            NEWFIC = NEWFIC .OR. AA(NOEL,IFIC).NE.A(NOEL,IFAC)
c            write(*,*) newfic ,' or ', AA(NOEL,IFIC) ,' NE ',
c     >           A(NOEL,IFAC),
c     >           'IFIC',IFIC,'IFAC', IFAC
            IFAC = IFAC + 1
            IFIC = IFIC + 1
          ENDDO
          IFAC = 24
          DO IFIC = 1, I2
            AA(NOEL,IFIC) = A(NOEL,IFAC)
            IFAC = IFAC + 1
          ENDDO
        ENDIF
      ENDIF
c      write(*,*) 'after', newfic
c      read(*,*)
 
      IF(NRES.GT.0) THEN
        WRITE(NRES,FMT='(/,5X,3(A,I3,A),/,5X,A,I2,A,I2,/)')
     >  'NDIM = ',NDIM,' ;  '
     >  ,'Number of data file sets used is ',NFIC,' ;  '
     >  ,'Stored in field array # IMAP =  ',IMAP,' ;  '
     >  ,'Value of MOD.I is ', MOD,'.',MOD2
        IF    (MOD.EQ.3) THEN
          IF(MOD2.EQ.1) THEN
            WRITE(NRES,FMT='(/,5X,A,I1,2A,1P,2E15.6,/)')
     >      'MOD2 = ',MOD2,' ->  map will be perturbed using '
     >      ,'field indices db/dx, d2b/dx2 : ',dbdx(1),dbdx(2)
          ENDIF
        ELSEIF(MOD.EQ.15) THEN
          WRITE(NRES,*)
     >    ' MOD=15.   Will sum up ',I2,'  field maps.'
          WRITE(NRES,*)
     >    ' Coefficient values : ',(a(noel,24+i-1),i=i1,i2)
C     >    ' MOD=15.   Will sum up ',I2,'  3D field maps.'
        ENDIF
 
        IF(NEWFIC) THEN
           WRITE(NRES,209)
 209       FORMAT(/,10X
     >     ,' New field map(s) now used, cartesian mesh (MOD.le.19) ; '
     >     ,/,10X,' name(s) of map data file(s) : ',/)
           WRITE(NRES,208) (NOMFIC(I),I=1,NFIC)
 208       FORMAT(10X,A)
        ELSE
          WRITE(NRES,210) (NOMFIC(I),I=1,NFIC)
 210      FORMAT(
     >    10X,'No  new  map  file  to  be  opened. Already  stored.',/
     >    10X,'Skip  reading  field  map  file : ',10X,A80)
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
             LNGTH=len(
     >         NOMFIC(NFIC)(DEBSTR(NOMFIC(NFIC)):FINSTR(NOMFIC(NFIC))))
 
             IF(NRES.GT.0) THEN
               WRITE(NRES,FMT='(/,3A,/)')
     >           NOMFIC(NFIC)(DEBSTR(NOMFIC(NFIC)):LNGTH),
     >           ' map,  FORMAT type : ', FMTYP
                CALL FLUSH2(NRES,.FALSE.)
             ENDIF
 
             IRD = NINT(A(NOEL,40))
 
C BNORM set to ONE, since sent to CHAMK below
             CALL FMAPR3(BINAR,LUN,MOD,MOD2,NHD,
     >                   XNORM,YNORM,ZNORM,ONE,I1,KZ,FMTYP,
     >                                    BMIN,BMAX,
     >                                    XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)
 
 12        CONTINUE
 
C------- Store mesh coordinates
           IIXMA(IMAP) = IXMA
           DO I=1,IXMA
             XXH(I,imap) =  XH(I)
           ENDDO
           JJYMA(IMAP) = JYMA
           DO J=1,JYMA
             YYH(J,imap) =  YH(J)
           ENDDO
           KKZMA(IMAP) = KZMA
C FM Nov 2011           DO K= 2, KZMA
           DO K= 1, KZMA
             ZZH(K,imap) = ZH(K)
           ENDDO
           bBMI(imap) = BMIN
           bBMA(imap) = BMAX
           XBbMI(imap) = XBMI
           YbBMI(imap) = YBMI
           ZbBMI(imap) = ZBMI
           XbBMA(imap) = XBMA
           YBBMA(imap) = YBMA
           ZBBMA(imap) = ZBMA
 
         ELSEIF(MOD .EQ. 15) THEN
 
           NFIC = 0
           DO KZ=I1,I2
             NFIC = NFIC+1
             IF(IDLUNI(
     >                 LUN)) THEN
               BINAR=BINARI(NOMFIC(NFIC),IB)
 
c                   write(*,*) ' toscac ',i1,i2,nfic,NOMFIC(NFIC)
c                        read(*,*)
 
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
             LNGTH=len(
     >         NOMFIC(NFIC)(DEBSTR(NOMFIC(NFIC)):FINSTR(NOMFIC(NFIC))))
 
             IF(NRES.GT.0) THEN
               WRITE(NRES,FMT='(/,3A,/)')
     >           NOMFIC(NFIC)(DEBSTR(NOMFIC(NFIC)):LNGTH),
     >           ' map,  FORMAT type : ', FMTYP
                CALL FLUSH2(NRES,.FALSE.)
             ENDIF
 
             IRD = NINT(A(NOEL,40))
 
             CALL FMAPW2(NFIC,AA(NOEL,NFIC))
 
             CALL FMAPR3(BINAR,LUN,MOD,MOD2,NHD,
     >                   XNORM,YNORM,ZNORM,ONE,I1,KZ,FMTYP,
     >                                    BMIN,BMAX,
     >                                    XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)
 
             if    (nfic.eq.1) then
 
               if(mod2.eq.1) then
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
 
               else
 
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
 
c                   write(*,*) ' toscac  faca : ',faca,kzma,jyma,ixma,id
c                      read(*,*)
 
               endif
               CLOSE(UNIT=LUN)
 
             elseif(nfic.eq.2) then
 
               if(mod2.eq.2) then
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
 
c                   write(*,*) ' toscac  facb : ',facb,kzma,jyma,ixma,id
c                      read(*,*)
               else
 
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
 
               endif
               CLOSE(UNIT=LUN)
 
             elseif(nfic.eq.3) then
 
               if(mod2.eq.3) then
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
 
               else
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
 
             elseif(nfic.eq.4) then
 
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
 
             endif
             CLOSE(UNIT=LUN)
c                  write(*,*) ' sbr toscac facb ',facb
c                       read(*,*)
 
           ENDDO
c               write(*,*) ' toscac imap = ',imap
c                   read(*,*)
C------- Store mesh coordinates
           IIXMA(IMAP) = IXMA
           DO I=1,IXMA
             XXH(I,imap) =  XH(I)
           ENDDO
           JJYMA(IMAP) = JYMA
           DO J=1,JYMA
             YYH(J,imap) =  YH(J)
           ENDDO
           KKZMA(IMAP) = KZMA
C FM Nov 2011           DO K= 2, KZMA
           DO K= 1, KZMA
             ZZH(K,imap) = ZH(K)
           ENDDO
           bBMI(imap) = BMIN
           bBMA(imap) = BMAX
           XBbMI(imap) = XBMI
           YbBMI(imap) = YBMI
           ZbBMI(imap) = ZBMI
           XbBMA(imap) = XBMA
           YBBMA(imap) = YBMA
           ZBBMA(imap) = ZBMA
 
 
         ELSE
           CALL ENDJOB('SBR toscac. No such value MOD = ',MOD)
         ENDIF
 
 
      ELSE
 
           IRD = NINT(A(NOEL,40))
 
C------- Restore mesh coordinates
           IXMA = IIXMA(IMAP)
           DO I=1,IXMA
             XH(I) = XXH(I,imap)
           ENDDO
           JYMA = JJYMA(IMAP)
           DO J=1,JYMA
             YH(J) = YYH(J,imap)
           ENDDO
           KZMA = KKZMA(IMAP)
           DO K= 1, KZMA
             ZH(K) = ZZH(K,imap)
           ENDDO
           BMIN = bBMI(imap)
           BMAX = bBMA(imap)
           XBMI = XBbMI(imap)
           YBMI = YbBMI(imap)
           ZBMI = ZbBMI(imap)
           XBMA = XbBMA(imap)
           YBMA = YBBMA(imap)
           ZBMA = ZBBMA(imap)
 
           IF(NRES.GT.0) THEN
             WRITE(NRES,fmt='(2A,I3,2A)') ' SBR TOSCAC, ',
     >       ' restored mesh coordinates for field map # ',imap,
     >       ',  name : ',
     >       NOMFIC(NFIC)(DEBSTR(NOMFIC(NFIC)):FINSTR(NOMFIC(NFIC)))
           ENDIF
 
      ENDIF ! NEWFIC
 
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
 96   WRITE(ABS(NRES),*) 'Error  open  file ',NOMFIC(NFIC)
      CALL ENDJOB('Leaving... ',-99)
      RETURN
 
      ENTRY TOSCA1
      DO I = 1, IZ
        NOMFIC(I) = ' '
      ENDDO
      RETURN
 
      END
