C  ZGOUBI, a program for computing the trajectories of charged particles
C  in electric and magnetic fields
C  Copyright (C) 1988-2007  Fran�ois M�ot
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
C  Fran�ois M�ot <fmeot@bnl.gov>
C  Brookhaven National Laboratory
C  C-AD, Bldg 911
C  Upton, NY, 11973, USA
C  -------
      SUBROUTINE TOSCAP(SCAL,NDIM,
C     >                           BMIN,BMAX,BNORM,
     >                          BMIN,BMAX,BNORM,XNORM,YNORM,ZNORM,
     >                           XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA,NEWFIC)
      USE DYNHC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-------------------------------------------------
C     Read TOSCA map with cylindrical coordinates.
C     TOSCA keyword with MOD.ge.20.
C-------------------------------------------------
      INCLUDE 'PARIZ.H'
      INCLUDE "XYZHC.H"
C      COMMON//XH(MXX),YH(MXY),ZH(IZ),HC(ID,MXX,MXY,IZ),IXMA,JYMA,KZMA
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
C      LOGICAL ZSYM
      INCLUDE "C.TYPFLD.H"     ! COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      INCLUDE "C.ORDRES.H"     ! COMMON/ORDRES/ KORD,IRD,IDS,IDB,IDE,IDZ
 
      LOGICAL BINAR,BINARI,IDLUNI, NEWFIC
      CHARACTER(LNTA) TITL , NOMFIC(IZ), NAMFIC
      SAVE NOMFIC, NAMFIC
      INTEGER DEBSTR,FINSTR
      SAVE NHDF
 
      LOGICAL STRCON
 
      CHARACTER(20) FMTYP
      DIMENSION XXH(MXX,MMAP), YYH(MXY,MMAP), ZZH(IZ,MMAP)
      SAVE XXH, YYH, ZZH
      DIMENSION BBMI(MMAP), BBMA(MMAP), XBBMI(MMAP), YBBMI(MMAP)
      DIMENSION ZBBMI(MMAP), XBBMA(MMAP), YBBMA(MMAP), ZBBMA(MMAP)
      SAVE BBMI, BBMA, XBBMI, YBBMI, ZBBMI, XBBMA, YBBMA, ZBBMA
      DIMENSION IIXMA(MMAP), JJYMA(MMAP), KKZMA(MMAP)
      SAVE IIXMA, JJYMA, KKZMA
C      INCLUDE 'MAPHDR.H'
      INCLUDE 'MXHD.H'

      DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE ::
     >     HCA, HCB, HCC
      SAVE HCA, HCB, HCC

      PARAMETER (MXC = 4)
      PARAMETER (MXAA2=24+MXC-1)
      DIMENSION AA(MXL,MXAA2), UU(MXAA2)
      SAVE AA

      PARAMETER (ONE = 1.D0)

      DATA NOMFIC / IZ*'               '/ 
      DATA NHDF / 8 /
      DATA FMTYP / ' regular' /
 
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
C July 2013      BNORM = A(NOEL,10)*SCAL
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
     >'Pgm toscap. Field map header has too many lines. Must be .le.'
     >,MXHD)
      IXMA = NINT( A(NOEL,20) )
      IF(IXMA.GT.MXX)
     >   CALL ENDJOB('X-dim of map is too large, max is ',MXX)
      JYMA = NINT( A(NOEL,21) )
      IF(JYMA.GT.MXY )
     >   CALL ENDJOB('Y-dim of map is too large, max is ',MXY)
 
      KZMA =NINT( A(NOEL,22) )
      IF(KZMA.GT.IZ )
     >  CALL ENDJOB('Z-dim of map is too large, max is ',IZ)
 
      MOD = NINT(A(NOEL,23))
      MOD2 = NINT(10.D0*A(NOEL,23)) - 10*MOD

      IF    (NDIM.EQ.2 ) THEN
        IF(MOD2 .LE. 1) MOD2 = 1
C Keyword TOSCA, option cylindrical frame, with up to 4 2D maps, 
C such that total B is their linear superimposition.
        IF(MOD2 .GT. 4) CALL ENDJOB('Pgm toscap. Number of field '
     >  //'maps cannot exceed ',MXC)
        I1=1
        I2 = MOD2
        DO I = 1, I2
          AA(NOEL,24-1+I) = A(NOEL,24-1+I)
        ENDDO
 
        NFIC=0
        DO I=I1, I2
          NFIC = NFIC+1
          NAMFIC = TA(NOEL,1+NFIC)
          NOMFIC(NFIC) = NAMFIC(DEBSTR(NAMFIC):FINSTR(NAMFIC))

          DO LL = 24, MXAA2     !24+MXC-1
            UU(LL) = AA(NOEL,LL)
          ENDDO
 
          CALL KSMAP4(NOMFIC,NFIC,UU,
     >                          NEWFIC,NBMAPS,IMAP)
        ENDDO

      ELSEIF(NDIM .EQ. 3 ) THEN

        IF(MOD .EQ. 20) THEN
C--------- an option on symmetrization for creating full field map from
C          a 1-quarter map (done by  ENTRY FMAPR2).
C          Cylindrical mesh. Axis is Z
          I1 = 1
          I2 = 1
        ELSEIF(MOD .EQ. 21) THEN
C--------- another option for symmetrization by FMAPR2
          I1 = 1
          I2 = KZMA
        ELSEIF(MOD .EQ. 22) THEN
C--------- MOD2 files are combined linearly into a single map after reading.
C          Each one of these files should contain the all 3D volume.
          IF(MOD2 .GT. 4) CALL ENDJOB('Pgm toscap. Number of field '
     >    //'maps cannot exceed ',MXC)
          I1=1
          I2 = MOD2
          DO I = 1, I2
            AA(NOEL,24-1+I) = A(NOEL,24-1+I)
          ENDDO
        ELSE
          STOP ' *** Error. SBR TOSCAP. No such MOD value '
        ENDIF

        NFIC=0
        DO I=I1, I2
          NFIC = NFIC+1
          NAMFIC = TA(NOEL,1+NFIC)
          NOMFIC(NFIC) = NAMFIC(DEBSTR(NAMFIC):FINSTR(NAMFIC))
        ENDDO

        DO LL = 24, MXAA2     !24+MXC-1
          UU(LL) = AA(NOEL,LL)
        ENDDO 
        CALL KSMAP4(NOMFIC,NFIC,UU,
     >                          NEWFIC,NBMAPS,IMAP)

      ENDIF

      IF(MOD .EQ. 22) THEN
        IFAC = 24
        IFIC = 1
        DO WHILE (IFIC.LE.I2 .AND. .NOT. NEWFIC)
          IF(IFAC .GT. 27)
     >    CALL ENDJOB('SBR toscap. No such possibility IFAC=',IFAC)
          NEWFIC = NEWFIC .OR. AA(NOEL,IFIC).NE.A(NOEL,IFAC)
          IFAC = IFAC + 1
          IFIC = IFIC + 1
        ENDDO
        IFAC = 24
        DO IFIC = 1, I2
          AA(NOEL,IFIC) = A(NOEL,IFAC)
          IFAC = IFAC + 1
        ENDDO
      ENDIF

      IF(NRES.GT.0) THEN

        WRITE(NRES,FMT='(/,5X,3(A,I3,A),/,5X,A,I0,A,I0,/)')
     >  'NDIM = ',NDIM,' ;  '
     >  ,'Number of data file sets used is ',NFIC,' ;  '
     >  ,'Stored in field array # IMAP =  ',IMAP,' ;  '
     >  ,'Value of MOD.MOD2 is ', MOD,'.',MOD2
        IF(MOD.EQ.22) THEN
          WRITE(NRES,*)
     >    ' MOD=22.   Will sum up ',I2,'  field maps.'
          WRITE(NRES,*)
     >    ' Coefficient values : ',(a(noel,24+i-1),i=i1,i2)
        ENDIF
 
        IF(NEWFIC) THEN
           WRITE(NRES,209)
 209       FORMAT(/,10X
     >     ,'New field map(s) now used, polar mesh (MOD .ge. 20) ; '
     >     ,/,10X,' name(s) of map data file(s) are : ')
!           WRITE(6   ,208) (NOMFIC(I),I=1,NFIC)
           WRITE(NRES,208)  (NOMFIC(I)(DEBSTR(NOMFIC(I)):
     >     FINSTR(NOMFIC(I))),I=1,NFIC)
 208       FORMAT(10X,A)
        ELSE
          WRITE(NRES,210)  (NOMFIC(I)(DEBSTR(NOMFIC(I)):
     >     FINSTR(NOMFIC(I))),I=1,NFIC)
 210      FORMAT(
     >    10X,'No  new  map  file  to  be  opened. Already  stored.',/
     >    10X,'Skip  reading  field  map  file : ',10X,A)
        ENDIF
      ENDIF
 
      IF(NEWFIC) THEN
 
        IF    (MOD .EQ. 20 .OR. MOD .EQ. 21) THEN

          IF(IDLUNI(
     >            LUN)) THEN
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
 
          IF(NRES.GT.0) THEN
              LNGTH=len(
     >        NOMFIC(NFIC)(DEBSTR(NOMFIC(NFIC)):FINSTR(NOMFIC(NFIC))))
              WRITE(NRES,FMT='(2(/,2X,A),I4,A,I0,A,/)') 
     >        ' ----',' Map file number ',NFIC,' ( of ',(I2-I1+1),') :'
              WRITE(NRES,FMT='(5X,3A)')
     >        NOMFIC(NFIC)(DEBSTR(NOMFIC(NFIC)):LNGTH)//' ',
     >        'map,  FORMAT type : ',FMTYP(DEBSTR(FMTYP):FINSTR(FMTYP))
              CALL FLUSH2(NRES,.FALSE.)
          ENDIF

          IRD = NINT(A(NOEL,40))
 
C BNORM set to ONE, since sent to CHAMK below
          KZ = 1
          CALL FMAPR2(BINAR,LUN,MOD,MOD2,NHD,
     >                   XNORM,YNORM,ZNORM,ONE,I1,KZ,FMTYP,
     >                      BMIN,BMAX,
     >                      XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)
 
C------- Store mesh coordinates
           CALL FMAPW4(IMAP,BMIN,BMAX,XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)
 
        ELSEIF(MOD .EQ. 22) THEN
 
          NFIC = 0
          DO KZ=I1,I2
            NFIC = NFIC+1
            IF(IDLUNI(
     >                LUN)) THEN
              BINAR=BINARI(NOMFIC(NFIC),IB)
              IF(BINAR) THEN
                OPEN(UNIT=LUN,FILE=NOMFIC(NFIC),FORM='UNFORMATTED'
     >          ,STATUS='OLD',ERR=96)
              ELSE
                OPEN(UNIT=LUN,FILE=NOMFIC(NFIC),STATUS='OLD',ERR=96)
              ENDIF
            ELSE
              GOTO 96
            ENDIF
 
            IF(NRES.GT.0) THEN
              LNGTH=len(
     >        NOMFIC(NFIC)(DEBSTR(NOMFIC(NFIC)):FINSTR(NOMFIC(NFIC))))
              WRITE(NRES,FMT='(2(/,2X,A),I4,A,I0,A,/)') 
     >        ' ----',' Map file number ',NFIC,' ( of ',(I2-I1+1),') :'
              WRITE(NRES,FMT='(5X,A,1P,E16.8,/)')
     >        NOMFIC(NFIC)(DEBSTR(NOMFIC(NFIC)):LNGTH)//' '//
     >        'map,  FORMAT type : '//FMTYP(DEBSTR(FMTYP):FINSTR(FMTYP))
     >        //'.  Field multiplication factor :',AA(NOEL,23+NFIC)
              CALL FLUSH2(NRES,.FALSE.)
            ENDIF
 
            IRD = NINT(A(NOEL,40))
 
C BNORM set to ONE, since sent to CHAMK below
            CALL FMAPR2(BINAR,LUN,MOD,MOD2,NHD,
     >                   XNORM,YNORM,ZNORM,ONE,I1,KZ,FMTYP,
     >                      BMIN,BMAX,
     >                      XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)
 
            CALL FMAPW2(NFIC,AA(NOEL,NFIC))

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
 
c                   write(*,*) ' toscac  faca : ',faca,kzma,jyma,ixma,id
c                      read(*,*)
 
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
 
c                   write(*,*) ' toscac  facb : ',facb,kzma,jyma,ixma,id
c                      read(*,*)
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

        ELSE
          CALL ENDJOB('Pgm toscap. No such option MOD = ',MOD)
        ENDIF
 
      ELSE   !  .not. NEWFIC
 
C------- Restore mesh coordinates
        CALL FMAPR5(IMAP,
     >                   BMIN,BMAX,XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)
 
        IF(NRES.GT.0) THEN
          WRITE(NRES,fmt='(2A,I3,2A)') ' Pgm toscac, ',
     >    ' restored mesh coordinates for field map # ',imap,
     >    ',  name : ',
     >    NOMFIC(NFIC)(DEBSTR(NOMFIC(NFIC)):FINSTR(NOMFIC(NFIC)))
        ENDIF
 
      ENDIF
 
      CALL CHAMK2(BNORM*SCAL)
 
      CALL MAPLI1(BMAX-BMIN)
      AT=XH(IXMA)-XH(1)
      ATO = 0.D0
      ATOS = 0.D0
      RM=.5D0*(YH(JYMA)+YH(1))
      XI = XH(1)
      XF = XH(IXMA)
 
      RETURN

 96   WRITE(ABS(NRES),*) 'Pgm toscap. Error  open  file ',
     >NOMFIC(NFIC)(DEBSTR(NOMFIC(NFIC)):FINSTR(NOMFIC(NFIC)))
      CALL ENDJOB('Leaving. ',-99)
      RETURN
      END
