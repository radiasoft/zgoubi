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
      SUBROUTINE TOSCAP(SCAL,NDIM,
C     >                           BMIN,BMAX,BNORM,
     >                          BMIN,BMAX,BNORM,XNORM,YNORM,ZNORM,
     >                           XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA,NEWFIC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-------------------------------------------------
C     Read TOSCA map with cylindrical coordinates. 
C     TOSCA keyword with MOD.ge.20. 
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
      CHARACTER(80) TA
      PARAMETER (MXTA=45)
      COMMON/DONT/ TA(MXL,MXTA)
      COMMON/DROITE/ AM(9),BM(9),CM(9),IDRT
      COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      LOGICAL ZSYM
      COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      COMMON/ORDRES/ KORD,IRD,IDS,IDB,IDE,IDZ

      LOGICAL BINAR,BINARI,IDLUNI, NEWFIC
      CHARACTER TITL*80 , NOMFIC(IZ)*80, NAMFIC*80
      SAVE NOMFIC, NAMFIC
      INTEGER DEBSTR,FINSTR
      SAVE NHDF

      LOGICAL STRCON

      CHARACTER*20 FMTYP
      DIMENSION XXH(MXX,MMAP), YYH(MXY,MMAP), ZZH(IZ,MMAP)
      SAVE XXH, YYH, ZZH
      DIMENSION BBMI(MMAP), BBMA(MMAP), XBBMI(MMAP), YBBMI(MMAP)
      DIMENSION ZBBMI(MMAP), XBBMA(MMAP), YBBMA(MMAP), ZBBMA(MMAP)
      SAVE BBMI, BBMA, XBBMI, YBBMI, ZBBMI, XBBMA, YBBMA, ZBBMA
      DIMENSION IIXMA(MMAP), JJYMA(MMAP), KKZMA(MMAP)
      SAVE IIXMA, JJYMA, KKZMA
 
      DATA NOMFIC / IZ*'               '/ 

      DATA NHDF / 8 /

      DATA FMTYP / ' regular' / 

      PARAMETER (MXC = 4)
      DIMENSION AA(24+MXC-1)
      DATA AA / 27 * 0.D0 /

      

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
      IXMA = A(NOEL,20)
      IF(IXMA.GT.MXX) 
     >   CALL ENDJOB('X-dim of map is too large, max is ',MXX)
      JYMA = A(NOEL,21)
      IF(JYMA.GT.MXY ) 
     >   CALL ENDJOB('Y-dim of map is too large, max is ',MXY)

      KZMA =A(NOEL,22)
      IF(KZMA.GT.IZ ) 
     >  CALL ENDJOB('Z-dim of map is too large, max is ',IZ)

      MOD = NINT(A(NOEL,23))
      MOD2 = NINT(10.D0*A(NOEL,23)) - 10*MOD

      IF    (NDIM.EQ.2 ) THEN
        NFIC=1
        I1 = 1
        I2 = 1
        NAMFIC = TA(NOEL,2)
        NAMFIC = NAMFIC(DEBSTR(NAMFIC):FINSTR(NAMFIC))
        NEWFIC = NAMFIC .NE. NOMFIC(NFIC)
        NOMFIC(NFIC) = NAMFIC(DEBSTR(NAMFIC):FINSTR(NAMFIC))
        CALL KSMAP4(NOMFIC,NFIC,AA(24:24+MXC-1),
     >                          NEWFIC,NBMAPS,IMAP)
      ELSEIF(NDIM .EQ. 3 ) THEN
C        IF(MOD .GE. 20) THEN
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
C--------- another option for symmetrization by FMAPR2
          I1 = 1
          I2 = 1  !!!KZMA   FM, 2009-04-10 for raccam measured 3-D field maps
        ELSE
          STOP ' *** Error. SBR TOSCAP. No such MOD value '
        ENDIF
        NFIC=0
        NEWFIC = .TRUE.
        DO 129 I=I1, I2
          NFIC = NFIC+1
          NAMFIC = TA(NOEL,1+NFIC)
          NEWFIC = NEWFIC .AND. (NAMFIC .NE. NOMFIC(NFIC))
          NOMFIC(NFIC) = NAMFIC(DEBSTR(NAMFIC):FINSTR(NAMFIC))
 129    CONTINUE
        CALL KSMAP4(NOMFIC,NFIC,AA(24:24+MXC-1),
     >                          NEWFIC,NBMAPS,IMAP)
      ENDIF
      IF(NRES.GT.0) WRITE(NRES,FMT='(/,5X,A,I1,A,I3,2A,I3,/)') 
     >'NDIM = ',NDIM,' ;   Value of MOD is ', MOD,' ;  ', 
     >'Number of data file sets used is ',NFIC

      IF(NRES.GT.0) THEN
        IF(NEWFIC) THEN
           WRITE(NRES,209) 
 209       FORMAT(/,10X  
     >     ,'New field map(s) now used, polar mesh (MOD .ge. 20) ; '
     >     ,/,10X,' name(s) of map data file(s) are : ')
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

        IF(IDLUNI(
     >            LUN)) THEN
          BINAR=BINARI(NOMFIC(NFIC),IB)
          IF(BINAR) THEN
            OPEN(UNIT=LUN,FILE=NOMFIC(NFIC),FORM='UNFORMATTED'
     >      ,STATUS='OLD',ERR=96)
          ELSE
            OPEN(UNIT=LUN,FILE=NOMFIC(NFIC),STATUS='OLD',ERR=96)
          ENDIF
        ELSE
          GOTO 96
        ENDIF

             LNGTH=len(
     >         NOMFIC(NFIC)(DEBSTR(NOMFIC(NFIC)):FINSTR(NOMFIC(NFIC))))
             WRITE(NRES,FMT='(/,3A)') 
     >         NOMFIC(NFIC)(DEBSTR(NOMFIC(NFIC)):LNGTH), 
     >         ' map,  FORMAT type : ', FMTYP             
        IRD = NINT(A(NOEL,40))
C        NHD = NHDF

C BNORM set to ONE, since sent to CHAMK below
C        CALL FMAPR2(BINAR,LUN,MOD,MOD2,NHD,BNORM,
        one = 1.d0
        kz = 1
        CALL FMAPR2(BINAR,LUN,MOD,MOD2,NHD,
     >                   XNORM,YNORM,ZNORM,ONE,I1,KZ,FMTYP,
     >                      BMIN,BMAX,
     >                      XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)

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

      ENDIF

c          write(*,*)  ' toscap BNORM SCAL ',BNORM,SCAL
      CALL CHAMK2(BNORM*SCAL)

       CALL MAPLI1(BMAX-BMIN)
       AT=XH(IXMA)-XH(1)
       ATO = 0.D0
       ATOS = 0.D0
       RM=.5D0*(YH(JYMA)+YH(1))
       XI = XH(1)
       XF = XH(IXMA)

      RETURN
 96   WRITE(ABS(NRES),*) 'ERROR  OPEN  FILE ',NOMFIC(NFIC)
      CALL ENDJOB('ERROR  OPEN  FILE ',-99)
      RETURN
      END
