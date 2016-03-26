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
      SUBROUTINE EMMAC(SCAL,NDIM,
     >                          BMIN,BMAX,BNORM,XNORM,YNORM,ZNORM,
     >               XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)
      USE dynhc
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-------------------------------------------------
C     Read TOSCA map with cartesian coordinates.
C     TOSCA keyword with MOD.le.19.
C-------------------------------------------------
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
C      LOGICAL ZSYM
      INCLUDE "C.TYPFLD.H"     ! COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      INCLUDE "C.ORDRES.H"     ! COMMON/ORDRES/ KORD,IRD,IDS,IDB,IDE,IDZ
 
      LOGICAL BINARI,IDLUNI
      LOGICAL BINAR, NEWFIC
      CHARACTER(LNTA) TITL, NOMFIC(2), NAMFIC(2)
 
      SAVE NOMFIC, NAMFIC
 
      INTEGER DEBSTR,FINSTR
 
      SAVE NHDF
 
      LOGICAL STRCON
 
      DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE ::
     >     HC1, HC2, HCA, HCB
C      DIMENSION HC1(ID,MXX,MXY,IZ), HC2(ID,MXX,MXY,IZ)
C      DIMENSION HCA(ID,MXX,MXY,IZ), HCB(ID,MXX,MXY,IZ)
      SAVE HC1, HC2, HCA, HCB
 
C      DIMENSION HCU(ID,MXX,MXY,IZ)
C      SAVE HCU
 
      parameter (idmx=ID*MXX*MXY*IZ)
 
      CHARACTER(20) FMTYP
 
      DIMENSION XXH(MXX,MMAP), YYH(MXY,MMAP), ZZH(IZ,MMAP)
      SAVE XXH, YYH, ZZH
      DIMENSION BBMI(MMAP), BBMA(MMAP), XBBMI(MMAP), YBBMI(MMAP)
      DIMENSION ZBBMI(MMAP), XBBMA(MMAP), YBBMA(MMAP), ZBBMA(MMAP)
      SAVE BBMI, BBMA, XBBMI, YBBMI, ZBBMI, XBBMA, YBBMA, ZBBMA
      DIMENSION IIXMA(MMAP), JJYMA(MMAP), KKZMA(MMAP)
      SAVE IIXMA, JJYMA, KKZMA
 
C      INCLUDE 'MAPHDR.H'

      DATA NOMFIC / 2*'               '/
 
      DATA NHDF / 8 /
 
C          data hcu /  IDMX * 0.d0 /
 
      DATA FMTYP / ' regular' /
      DATA IMAP / 1 /
 
      PARAMETER (MXC = 4)
      DIMENSION AA(24+MXC-1)
      DATA AA / 27 * 0.D0 /
 
      if( .NOT.ALLOCATED( HC1 ))
     >     ALLOCATE( HC1(ID,MXX,MXY,IZ), STAT = IALOC)
      IF (IALOC /= 0)
     >     CALL ENDJOB('SBR EMMAC Not enough memory for Malloc of HC1',
     >     -99)
 
      if( .NOT.ALLOCATED( HC2 ))
     >     ALLOCATE( HC2(ID,MXX,MXY,IZ), STAT = IALOC)
      IF (IALOC /= 0)
     >     CALL ENDJOB('SBR EMMAC Not enough memory for Malloc of HC2',
     >     -99)
 
      if( .NOT.ALLOCATED( HCA ))
     >     ALLOCATE( HCA(ID,MXX,MXY,IZ), STAT = IALOC)
      IF (IALOC /= 0)
     >     CALL ENDJOB('SBR EMMAC Not enough memory for Malloc of HCA',
     >     -99)
 
      if( .NOT.ALLOCATED( HCB ))
     >     ALLOCATE( HCB(ID,MXX,MXY,IZ), STAT = IALOC)
      IF (IALOC /= 0)
     >     CALL ENDJOB('SBR EMMAC Not enough memory for Malloc of HCB',
     >     -99)
 
c      CALL KSMAP(
c     >           IMAP)
 
      BNORM = A(NOEL,10)*SCAL
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
     >'Pgm emmac. Field map header has too many lines. Must be .le.'
     >,MXHD)

      IXMA = NINT(A(NOEL,20))
      IF(IXMA.GT.MXX)
     >   CALL ENDJOB('X-dim of map is too large,  max  is ',MXX)
      JYMA = NINT(A(NOEL,21))
      IF(JYMA.GT.MXY )
     >   CALL ENDJOB('Y-dim of map is too large,  max  is ',MXY)
      IF(NDIM .EQ. 3) THEN
        KZMA =NINT(A(NOEL,22))
        IF(KZMA.GT.IZ )
     >     CALL ENDJOB('Z-dim of map is too large,  max  is ',IZ)
      ENDIF
 
      MOD = NINT(A(NOEL,23))
      MOD2 = NINT(10.D0*A(NOEL,23)) - 10*MOD
 
      AF = A(NOEL,30)
      AD = A(NOEL,31)
      DIST = A(NOEL,32)
      DIST2 = A(NOEL,33)
 
      IF    (NDIM.EQ.2 ) THEN
 
        IF(MOD .EQ. 0) THEN
C--------- Rectangular mesh
C Keyword EMMA with *two* 2D maps, one for QF, one for QD,
C a single map superimposition of both is built prior to tracking
C and used for tracking.
 
          NFIC = 2
          NAMFIC(1) = TA(NOEL,NFIC)
          CALL KSMAP4(NAMFIC(1),NFIC-1,AA(24:24+MXC-1),
     >                          NEWFIC,NBMAPS,IMAP)
          NAMFIC(2) = TA(NOEL,1+NFIC)
          NAMFIC(1) = NAMFIC(1)(DEBSTR(NAMFIC(1)):FINSTR(NAMFIC(1)))
          NAMFIC(2) = NAMFIC(2)(DEBSTR(NAMFIC(2)):FINSTR(NAMFIC(2)))
 
 
 
          NEWFIC = (NAMFIC(1) .NE. NOMFIC(1))
          NEWFIC = NEWFIC .OR. (NAMFIC(2) .NE. NOMFIC(2))
          NOMFIC(1) = NAMFIC(1)(DEBSTR(NAMFIC(1)):FINSTR(NAMFIC(1)))
          NOMFIC(2) = NAMFIC(2)(DEBSTR(NAMFIC(2)):FINSTR(NAMFIC(2)))
 
        ELSEIF(MOD .EQ. 1) THEN
C--------- Rectangular mesh
C Keyword EMMA with *two* 2D maps, one for QF, one for QD,
C the resulting single map is obtained in the following way :
C  QF_new is interpolated from QF with dr=xF,
C  QD_new is interpolated from QD with dr=xD.
C A single map superimposition of both is built prior to tracking
C and used for tracking.
          NFIC = 2
          NAMFIC(1) = TA(NOEL,NFIC)
          NAMFIC(2) = TA(NOEL,1+NFIC)
          NAMFIC(1) = NAMFIC(1)(DEBSTR(NAMFIC(1)):FINSTR(NAMFIC(1)))
          NAMFIC(2) = NAMFIC(2)(DEBSTR(NAMFIC(2)):FINSTR(NAMFIC(2)))
 
          NEWFIC = (NAMFIC(1) .NE. NOMFIC(1))
          NEWFIC = NEWFIC .OR. (NAMFIC(2) .NE. NOMFIC(2))
          NOMFIC(1) = NAMFIC(1)(DEBSTR(NAMFIC(1)):FINSTR(NAMFIC(1)))
          NOMFIC(2) = NAMFIC(2)(DEBSTR(NAMFIC(2)):FINSTR(NAMFIC(2)))
 
        ELSE
 
          CALL ENDJOB('*** Error, SBR EMMAC -> No  such  MOD= ',MOD)
 
        ENDIF
 
      ELSEIF(NDIM .EQ. 3 ) THEN
 
        CALL ENDJOB('*** Error, SBR EMMAC -> No  such  NDIM= ',NDIM)
 
      ENDIF
 
      IF(NRES.GT.0) THEN
 
        WRITE(NRES,FMT='(/,5X,2(A,I3),A,A,I3)')
     >  'NDIM = ',NDIM,' ;   value of MOD is ', MOD,' ;  '
     >  ,'Number of map files set used is ',NFIC
 
        IF(NEWFIC) THEN
           WRITE(NRES,209)
 209       FORMAT(/,10X
     >     ,' New field map(s) now used, cartesian mesh (MOD.le.19) ; '
     >     ,/,10X,' name(s) of map data file(s) : ')
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
 
        WRITE(NRES,FMT='(/,5X,A,1P,2E12.4)')
     >  ' Value of field factor coefficients AF, AD are : ', AF, AD
 
        IF(MOD .EQ. 0) THEN
          WRITE(NRES,FMT='(/,5X,A,1P,E12.4)')
     >    ' Distance between axis of quads : ', DIST
          WRITE(NRES,FMT='(/,5X,A)')
     >    'A single map is built by superimposition prior to tracking'
        ELSEIF(MOD .EQ. 1) THEN
          WRITE(NRES,FMT='(/,5X,A,1P,E14.6,A)')
     >    ' Field map 1 recomputed following radial move ',DIST,' (cm)'
          WRITE(NRES,FMT='(  5X,A,1P,E14.6,A)')
     >    ' Field map 2 recomputed following radial move ',DIST2,' (cm)'
          WRITE(NRES,FMT='(/,5X,A)')
     >    'A single map is built by superimposition prior to tracking'
        ENDIF
      ENDIF
 
      IF(NEWFIC) THEN
 
        IF(MOD .EQ. 0) THEN
C------- Rectangular mesh
C Keyword EMMA with *two* 2D maps, one for QF, one for QD,
C a single map superimposition of both is built prior to tracking
C and used for tracking.
 
             IRD = NINT(A(NOEL,50))
             KZ = 1
             I1=1
 
             NFIC = 1
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
 
             LNGTH=len(
     >         NOMFIC(NFIC)(DEBSTR(NOMFIC(NFIC)):FINSTR(NOMFIC(NFIC))))
             IF(NRES.GT.0) WRITE(NRES,FMT='(/,3A)')
     >         NOMFIC(NFIC)(DEBSTR(NOMFIC(NFIC)):LNGTH),
     >         ' map,  FORMAT type : ', FMTYP
             CALL FMAPR3(BINAR,LUN,MOD,MOD2,NHD,
     >                   XNORM,YNORM,ZNORM,BNORM,I1,KZ,FMTYP,
     >                                    BMIN,BMAX,
     >                                    XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)
 
               do jjj=1,JYMA
                 do iii=1,IXMA
                   do iid = 1, id
                     HC1(iid,iii,jjj,KZ) = HC(iid,iii,jjj,KZ,IMAP)
c                      write(86,*) iid,iii,jjj,KZ,HC(iid,iii,jjj,KZ,IMAP)
                   enddo
                 enddo
               enddo
 
             NFIC = 2
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
 
             LNGTH=len(
     >         NOMFIC(NFIC)(DEBSTR(NOMFIC(NFIC)):FINSTR(NOMFIC(NFIC))))
             IF(NRES.GT.0) WRITE(NRES,FMT='(/,3A)')
     >         NOMFIC(NFIC)(DEBSTR(NOMFIC(NFIC)):LNGTH),
     >         ' map,  FORMAT type : ', FMTYP
             CALL FMAPR3(BINAR,LUN,MOD,MOD2,NHD,
     >                   XNORM,YNORM,ZNORM,BNORM,I1,KZ,FMTYP,
     >                             BMIN,BMAX,
     >                             XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)
 
               do jjj=1,JYMA
                 do iii=1,IXMA
                   do iid = 1, id
                     HC2(iid,iii,jjj,KZ) = HC(iid,iii,jjj,KZ,IMAP)
                   enddo
                 enddo
               enddo
 
 
          do jjj=1,JYMA
            do iii=1,IXMA
              do iid = 1, id
                HC(iid,iii,jjj,KZ,IMAP) = AF * HC1(iid,iii,jjj,KZ)
     >                         + AD * HC2(iid,iii,jjj,KZ)
c                write(86,*) af,ad,iid,iii,jjj,KZ,HC(iid,iii,jjj,KZ,IMAP)
              enddo
            enddo
          enddo
 
        ELSEIF(MOD .EQ. 1) THEN
C--------- Rectangular mesh
C Keyword EMMA with *two* 2D maps, one for QF, one for QD,
C the resulting single map is obtained in the following way :
C  QF_new is interpolated from QF with dr=xF,
C  QD_new is interpolated from QD with dr=xD.
C A single map superimposition of both is built prior to tracking
C and used for tracking.
 
           IRD = NINT(A(NOEL,50))
           KZ = 1
           I1=1
 
           IF(NEWFIC) THEN
 
             NFIC = 1
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
 
             LNGTH=LEN(
     >         NOMFIC(NFIC)(DEBSTR(NOMFIC(NFIC)):FINSTR(NOMFIC(NFIC))))
             IF(NRES.GT.0) WRITE(NRES,FMT='(/,3A)')
     >         NOMFIC(NFIC)(DEBSTR(NOMFIC(NFIC)):LNGTH),
     >         ' map,  FORMAT type : ', FMTYP
             CALL FMAPR3(BINAR,LUN,MOD,MOD2,NHD,
     >                   XNORM,YNORM,ZNORM,BNORM,I1,KZ,FMTYP,
     >                                 BMIN,BMAX,
     >                                 XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)
 
             DO JJJ=1,JYMA
               DO III=1,IXMA
                 DO IID = 1, ID
                   HCA(IID,III,JJJ,KZ) = HC(IID,III,JJJ,KZ,IMAP)
                 ENDDO
               ENDDO
             ENDDO
             CLOSE(UNIT=LUN)
 
           ENDIF ! NEWFIC
 
           CALL MAPSHF(HCA,XH,YH,DIST,IXMA,JYMA,
     >                                           HC1)
 
           IF(NEWFIC) THEN
 
             NFIC = 2
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
 
             LNGTH=len(
     >         NOMFIC(NFIC)(DEBSTR(NOMFIC(NFIC)):FINSTR(NOMFIC(NFIC))))
             IF(NRES.GT.0) WRITE(NRES,FMT='(/,3A)')
     >         NOMFIC(NFIC)(DEBSTR(NOMFIC(NFIC)):LNGTH),
     >         ' map,  FORMAT type : ', FMTYP
             CALL FMAPR3(BINAR,LUN,MOD,MOD2,NHD,
     >                   XNORM,YNORM,ZNORM,BNORM,I1,KZ,FMTYP,
     >                               BMIN,BMAX,
     >                               XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)
             CLOSE(UNIT=LUN)
 
             DO JJJ=1,JYMA
               DO III=1,IXMA
                 DO IID = 1, ID
                   HCB(IID,III,JJJ,KZ) = HC(IID,III,JJJ,KZ,IMAP)
                 ENDDO
               ENDDO
             ENDDO
 
           ENDIF ! NEWFIC
 
           CALL MAPSHF(HCB,XH,YH,DIST2,IXMA,JYMA,
     >                                           HC2)
 
C Superimpose the two field maps into a single one
C          zer1 = 0.d0
          do jjj=1,JYMA
            do iii=1,IXMA
              do iid = 1, id
                HC(iid,iii,jjj,KZ,IMAP) = AF * HC1(iid,iii,jjj,KZ)
     >                         + AD * HC2(iid,iii,jjj,KZ)
c                zerA = HCA(iid,iii,jjj,KZ) - HC1(iid,iii,jjj,KZ)
c                zerB = HCB(iid,iii,jjj,KZ) - HC2(iid,iii,jjj,KZ)
              enddo
            enddo
          enddo
 
        ELSE
 
          CALL ENDJOB('*** Error, SBR EMMAC -> No  such  MOD= ',MOD)
 
        ENDIF
 
C------- Store mesh coordinates
           CALL FMAPW4(IMAP,BMIN,BMAX,XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)
 
      ELSE
C------- Restore mesh coordinates
           CALL FMAPW4(IMAP,BMIN,BMAX,XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)

           IRD = NINT(A(NOEL,40))
  
           IF(NRES.GT.0) THEN
             WRITE(NRES,fmt='(2A,I3,2A)') ' Pgm emmac, ',
     >       ' restored mesh coordinates for field map # ',imap,
     >       ',  name : ',
     >       NOMFIC(NFIC)(DEBSTR(NOMFIC(NFIC)):FINSTR(NOMFIC(NFIC)))
           ENDIF
 
      ENDIF ! NEWFIC
 
      CALL CHAMK2(BNORM*SCAL)
 
C Make sure this is ok with cartésien
        CALL MAPLI1(BMAX-BMIN)
 
      RETURN

 96   WRITE(ABS(NRES),*) 'Pgm toscap. Error  open  file ',
     >NOMFIC(NFIC)(DEBSTR(NOMFIC(NFIC)):FINSTR(NOMFIC(NFIC)))
      CALL ENDJOB('Leaving. ',-99)
      RETURN
      END
