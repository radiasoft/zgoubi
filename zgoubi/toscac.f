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
C  François Méot <meot@lpsc.in2p3.fr>
C  Service Accélerateurs
C  LPSC Grenoble
C  53 Avenue des Martyrs
C  38026 Grenoble Cedex
C  France
      SUBROUTINE TOSCAC(SCAL,NDIM, 
     >                          BMIN,BMAX,BNORM,XNORM,YNORM,ZNORM,
     >               XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-------------------------------------------------
C     Read TOSCA map with cartesian coordinates. 
C     TOSCA keyword with MOD.le.19. 
C-------------------------------------------------
      INCLUDE 'PARIZ.H'
      INCLUDE "XYZHC.H"
C      COMMON//XH(MXX),YH(MXY),ZH(IZ),HC(ID,MXX,MXY,IZ),IXMA,JYMA,KZMA
      COMMON/AIM/ ATO,AT,ATOS,RM,XI,XF,EN,EB1,EB2,EG1,EG2
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER*80 TA
      COMMON/DONT/ TA(MXL,20)
      COMMON/DROITE/ AM(9),BM(9),CM(9),IDRT
      COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      LOGICAL ZSYM
      COMMON/OPTION/ KFLD,MG,LC,ML,ZSYM
      COMMON/ORDRES/ KORD,IRD,IDS,IDB,IDE,IDZ

      LOGICAL BINARI,IDLUNI
      LOGICAL BINAR, NEWFIC
      LOGICAL FLIP
      CHARACTER TITL*80 , NOMFIC(IZ)*80, NAMFIC*80
      DATA NOMFIC / IZ*'               '/ 
      SAVE NOMFIC, NAMFIC

      INTEGER DEBSTR,FINSTR

      DATA NHDF / 8 /
      SAVE NHDF

      LOGICAL STRCON 

      CHARACTER*20 FMTYP
      DATA FMTYP / ' regular' / 

      BNORM = A(NOEL,10)*SCAL
      XNORM = A(NOEL,11)
      YNORM = A(NOEL,12)
      ZNORM = A(NOEL,13)
      TITL = TA(NOEL,1)
      IF    (STRCON(TITL,'HEADER',
     >                              IS) ) THEN
        READ(TITL(IS+7:IS+7),FMT='(I1)') NHD
             write(*,*) '  header, nhd, is ', TITL(IS:IS+7),nhd,is
      ELSE
        NHD = NHDF
      ENDIF
      IDEB = DEBSTR(TITL)
      FLIP = TITL(IDEB:IDEB+3).EQ.'FLIP'
      IXMA = A(NOEL,20)
      IF(IXMA.GT.MXX) 
     >   CALL ENDJOB('X-dim of map is too large,  max  is ',MXX)
      IF(NDIM .EQ. 1) THEN
        JYMA=1
      ELSE
        JYMA = A(NOEL,21)
        IF(JYMA.GT.MXY ) 
     >     CALL ENDJOB('Y-dim of map is too large,  max  is ',MXY)
        IF(NDIM .EQ. 3) THEN
          KZMA =A(NOEL,22)
          IF(KZMA.GT.IZ ) 
     >       CALL ENDJOB('Z-dim of map is too large,  max  is ',IZ)
        ENDIF
      ENDIF

      MOD = NINT(A(NOEL,23))
      MOD2 = NINT(10.D0*A(NOEL,23)) - 10*MOD

      IF    (NDIM.EQ.2 ) THEN
        I1=1
        I2 = 1
        NFIC=1
        NAMFIC = TA(NOEL,2)
        NAMFIC = NAMFIC(DEBSTR(NAMFIC):FINSTR(NAMFIC))
        NEWFIC = NAMFIC .NE. NOMFIC(NFIC)
!          write(*,*) 
!          write(*,*) NEWFIC 
!          write(*,*) NAMFIC 
!          write(*,*) NOMFIC(NFIC),nfic
        NOMFIC(NFIC) = NAMFIC(DEBSTR(NAMFIC):FINSTR(NAMFIC))
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
        ELSE
          STOP ' *** Error. SBR TOSCAC. No such MOD value '
        ENDIF
        NFIC=0
        NEWFIC = .TRUE.
        DO 129 I=I1, I2
          NFIC = NFIC+1
          NAMFIC = TA(NOEL,1+NFIC)
          NAMFIC = NAMFIC(DEBSTR(NAMFIC):FINSTR(NAMFIC))
          NEWFIC = NEWFIC .AND. (NAMFIC .NE. NOMFIC(NFIC))
          NOMFIC(NFIC) = NAMFIC(DEBSTR(NAMFIC):FINSTR(NAMFIC))
 129    CONTINUE
        IF(NRES.GT.0) WRITE(NRES,FMT='(/,5X,A,I3,2A,I3,/)') 
     >  'NDIM = 3 ;   Value of MOD is ', MOD,' ;  ', 
     >  'Number of field data files used is ',NFIC
      ENDIF

      IF(NRES.GT.0) THEN
        IF(NEWFIC) THEN
           WRITE(NRES,209) 
 209       FORMAT(/,10X,' A new field map is now used ; ', 
     >     ' Name(s) of map data file(s) are : ')
!           WRITE(6   ,208) (NOMFIC(I),I=1,NFIC)
           WRITE(NRES,208) (NOMFIC(I),I=1,NFIC)
 208       FORMAT(10X,A)
        ELSE
          WRITE(NRES,210) (NOMFIC(I),I=1,NFIC)
 210      FORMAT(
     >    10X,'No  new  map  file  to  be  opened. Already  stored.',/
     >    10X,'Skip  reading  field  map  file : ',10X,A80)
        ENDIF
      ENDIF 

      IF(NEWFIC) THEN

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

             IF(     STRCON(NOMFIC(NFIC),'GSI',IS)
     >          .OR. STRCON(NOMFIC(NFIC),'gsi',IS)) THEN
                  NHD = 0
                  FMTYP = 'GSI'
             ENDIF
             LNGTH=len(
     >         NOMFIC(NFIC)(DEBSTR(NOMFIC(NFIC)):FINSTR(NOMFIC(NFIC))))
             WRITE(NRES,FMT='(/,3A)') 
     >         NOMFIC(NFIC)(DEBSTR(NOMFIC(NFIC)):LNGTH), 
     >         ' map,  FORMAT type : ', FMTYP             
             IRD = NINT(A(NOEL,40))
             CALL FMAPR3(BINAR,LUN,MOD,MOD2,NHD,BNORM,I1,KZ,FMTYP,
     >                      BMIN,BMAX,
     >                      XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)

 12        CONTINUE

      ENDIF

c voir si ok avec cartésien
        CALL MAPLI1(BMAX-BMIN)
        AT=XH(IAMA)-XH(1)
        ATO = 0.D0
        ATOS = 0.D0
        RM=.5D0*(YH(JRMA)+YH(1))
        XI = XH(1)
        XF = XH(IAMA)
 
      RETURN
 96   WRITE(NRES,*) 'ERROR  OPEN  FILE ',NOMFIC(NFIC)
      CALL ENDJOB('ERROR  OPEN  FILE ',-99)
      RETURN

      ENTRY TOSCA1
      DO I = 1, IZ
        NOMFIC(I) = ' '
      ENDDO
      RETURN

      END
