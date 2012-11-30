C  ZGOUBI, a program for computing the trajectories of charged particles
C  in electric and magnetic fields
C  Copyright (C) 1988-2007  François Meot
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
C  François Meot <fmeot@bnl.gov>
C  BNL
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE MAP2D(SCAL,
     >                      BMIN,BMAX,BNORM,XNORM,YNORM,ZNORM,
     >               XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA,NEWFIC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-------------------------------------------------
C     Read TOSCA map with cartesian coordinates. 
C     TOSCA keyword with MOD.le.19. 
C-------------------------------------------------
      LOGICAL NEWFIC
      INCLUDE 'PARIZ.H'
      INCLUDE "XYZHC.H"
      COMMON/AIM/ ATO,AT,ATOS,RM,XI,XF,EN,EB1,EB2,EG1,EG2
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER*80 TA
      COMMON/DONT/ TA(MXL,40)
      COMMON/DROITE/ AM(9),BM(9),CM(9),IDRT
      COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      LOGICAL ZSYM
      COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      COMMON/ORDRES/ KORD,IRD,IDS,IDB,IDE,IDZ

      LOGICAL BINARI,IDLUNI
      LOGICAL BINAR
      LOGICAL FLIP
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

      DATA NOMFIC / IZ*' '/ 
      DATA NHDF / 8 /
      DATA FMTYP / ' regular' / 

      BNORM = A(NOEL,10)*SCAL
      XNORM = A(NOEL,11)
      YNORM = A(NOEL,12)
      ZNORM = A(NOEL,13)
      TITL = TA(NOEL,1)
      IF    (STRCON(TITL,'HEADER',
     >                              IS) ) THEN
        READ(TITL(IS+7:IS+7),FMT='(I1)') NHD
      ELSE
        NHD = NHDF
      ENDIF
      IDEB = DEBSTR(TITL)
      FLIP = TITL(IDEB:IDEB+3).EQ.'FLIP'
      IXMA = A(NOEL,20)
      IF(IXMA.GT.MXX) 
     >   CALL ENDJOB('X-dim of map is too large,  max  is ',MXX)
      JYMA = A(NOEL,21)
      IF(JYMA.GT.MXY ) 
     >   CALL ENDJOB('Y-dim of map is too large,  max  is ',MXY)

      KZMA = 1
      NFIC = 1
      NAMFIC = TA(NOEL,2)
      NOMFIC(NFIC) = NAMFIC(DEBSTR(NAMFIC):FINSTR(NAMFIC))

      CALL KSMAP4(NOMFIC,NFIC,
     >                        NEWFIC,NBMAPS,IMAP)

      IF(NRES.GT.0) THEN
C        WRITE(NRES,*) ' '
C        WRITE(NRES,*) ' map2d ', 
C     >       NAMFIC(DEBSTR(NAMFIC):FINSTR(NAMFIC)),
C     >       '       file #',nfic
        WRITE(NRES,FMT='(/,5X,2(A,I3,A),/)') 
     >  'Number of data file sets used is ',NFIC,' ;  ' 
     >  ,'Stored in field array # IMAP =  ',IMAP,' '
        IF(NEWFIC) THEN
          WRITE(NRES,209) 
 209      FORMAT(/,10X  
     >    ,' New field map now used, cartesian mesh (MOD.le.19) ; '
     >    ,/,10X,' name of map data file : ')
          WRITE(NRES,208) (NOMFIC(I),I=1,NFIC)
 208      FORMAT(10X,A)
        ELSE
          WRITE(NRES,210) (NOMFIC(I),I=1,NFIC)
 210      FORMAT(
     >    10X,'No  new  map  file  to  be  opened. Already  stored.',/
     >    10X,'Skip  reading  field  map  file : ',10X,A80)
        ENDIF
      ENDIF 

      INDEX=0
      NT = 1
      CALL PAVELW(INDEX,NT)

      IF(NEWFIC) THEN
               NFIC = 1
               IF(IDLUNI(
     >                   LUN)) THEN
                 BINAR=BINARI(NOMFIC(NFIC),IB)
                 IF(BINAR) THEN
                   OPEN(UNIT=LUN,FILE=NOMFIC(NFIC),FORM='UNFORMATTED'
     >             ,STATUS='OLD',ERR=96)
                 ELSE
                   OPEN(UNIT=LUN,FILE=NOMFIC(NFIC),STATUS='OLD',ERR=96)
                 ENDIF
               ELSE
                 GOTO 96
               ENDIF

             MOD = 0
             MOD2 = 3
             I1 = 1
             KZ = 1
             IRD = NINT(A(NOEL,40))

             CALL FMAPR3(BINAR,LUN,MOD,MOD2,NHD,
     >                   XNORM,YNORM,ZNORM,BNORM,I1,KZ,FMTYP,
     >                                    BMIN,BMAX,
     >                                    XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)

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
C FM Nov 2011   DO K= 2, KZMA
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

c           write(*,*) ' oldfic imap ', xxh(1,imap),xxh(ixma,imap),imap
c     >         ,HC(3,1,1,1,iMAP)

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

           IF(NRES.GT.0) WRITE(NRES,*) ' SBR TOSCAC, ',
     >     ' restored mesh coordinates for field map # ',imap

      ENDIF


      CALL MAPLI1(BMAX-BMIN)

      RETURN

 96   WRITE(ABS(NRES),*) 'ERROR  OPEN  FILE ',NOMFIC(NFIC)
      CALL ENDJOB('ERROR  OPEN  FILE ',-99)

      RETURN
      END 
