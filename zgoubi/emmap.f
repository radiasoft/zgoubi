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
      SUBROUTINE EMMAP(SCAL,NDIM,
     >                           BMIN,BMAX,BNORM,
     >                           XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)
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
      CHARACTER*80 TA
      COMMON/DONT/ TA(MXL,20)
      COMMON/DROITE/ AM(9),BM(9),CM(9),IDRT
      COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      LOGICAL ZSYM
      COMMON/OPTION/ KFLD,MG,LC,ML,ZSYM
      COMMON/ORDRES/ KORD,IRD,IDS,IDB,IDE,IDZ

      LOGICAL BINARI,IDLUNI
      LOGICAL BINAR, NEWFIC
      PARAMETER (MXFIC=2*40,MXFIC1 = MXFIC+1)
      CHARACTER TITL*80 , NOMFIC(MXFIC1)*80, NAMFIC(2)*80, FILEL*80
      DATA NOMFIC, NAMFIC, FILEL / MXFIC1*'   ', 2*'   ', '   '/
      SAVE NOMFIC, NAMFIC, FILEL

      INTEGER DEBSTR,FINSTR

      PARAMETER (MXDST = MXFIC/2)
      DIMENSION HCA(ID,MXX,MXY,IZ,MXDST),HCB(ID,MXX,MXY,IZ,MXDST),
     >HCC(ID,MXX,MXY,IZ,MXDST)
      SAVE HCA, HCB
      DIMENSION DISTL(MXDST)

C      LOGICAL FITING

      DATA NHDF / 8 /
      SAVE NHDF

      LOGICAL STRCON

      DATA KFIC0 / -99 /

      BNORM = A(NOEL,10)*SCAL
      XNORM = A(NOEL,11)
      YNORM = A(NOEL,12)
      ZNORM = A(NOEL,13)
      TITL = TA(NOEL,1)
      IF    (STRCON(TITL,'HEADER',
     >                            IS) ) THEN
        READ(TITL(IS+8:IS+8),FMT='(I1)') NHD
      ELSE
        NHD = NHDF
      ENDIF
      IXMA = A(NOEL,20)
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

      MOD = NINT(A(NOEL,23))
      MOD2 = NINT(100.D0*A(NOEL,23)) - 100*MOD

      AF = A(NOEL,30)
      AD = A(NOEL,31)
      DIST = A(NOEL,32)

      IF    (NDIM.EQ.2 ) THEN
        IF    (MOD .EQ. 22) THEN
C Keyword EMMA with *two* 2D maps, one for QF, one for QD, 
C such that total B is a linear superimposition of both

          NFIC = 2
          NAMFIC(1) = TA(NOEL,2)
          NAMFIC(2) = TA(NOEL,3)
          NAMFIC(1) = NAMFIC(1)(DEBSTR(NAMFIC(1)):FINSTR(NAMFIC(1)))
          NAMFIC(2) = NAMFIC(2)(DEBSTR(NAMFIC(2)):FINSTR(NAMFIC(2)))

          NEWFIC = (NAMFIC(1) .NE. NOMFIC(1))
          NEWFIC = NEWFIC .OR. (NAMFIC(2) .NE. NOMFIC(2))
          NOMFIC(1) = NAMFIC(1)
          NOMFIC(2) = NAMFIC(2)

        ELSEIF(MOD .EQ. 24) THEN
C Keyword EMMA using several pairs of 2D maps with parameter the quad axis distance DIST. 
C A pair comprises one QF & one QD map  such that total B is a linear superimposition of both. 

          NFIC = 1
          FILEL = TA(NOEL,2)
          FILEL = FILEL(DEBSTR(FILEL):FINSTR(FILEL))

          NEWFIC = (FILEL .NE. NOMFIC(MXFIC1))
          NOMFIC(MXFIC1) = FILEL

        ELSE

          CALL ENDJOB('*** Error, SBR EMMAP -> No  such  MOD= ',MOD)

        ENDIF

      ELSEIF(NDIM .EQ. 3 ) THEN

      ENDIF

      IF(NRES.GT.0) THEN

        WRITE(NRES,FMT='(/,5X,A,I3,A,I3,2A,I2,/)') 
     >  'NDIM = 3 ;   Value of MOD / MOD2 is ', MOD,' / ',MOD2,' ;  ', 
     >  'Number of field data files used is ',NFIC

        WRITE(NRES,FMT='(/,5X,A,1P,2E12.4,/)') 
     >  ' Value of multiplying coefficients AF, AD are : ', AF, AD

        WRITE(NRES,FMT='(/,5X,A,1P,E12.4,A,/)') 
     >  ' Distance between axis of quads is : ', DIST,' cm'

        IF(NEWFIC) THEN
           WRITE(NRES,209) 
 209       FORMAT(/,10X,' new field maps are now used ; ', 
     >     ' Names of files are : ')
           IF(MOD.EQ.24) THEN 
C             WRITE(6   ,208) NOMFIC(MXFIC1),' (field map file list)'
             WRITE(NRES,208) NOMFIC(MXFIC1),' (field map file list)'
           ELSE
C             WRITE(6   ,208) (NOMFIC(I),' (field map)',I=1,NFIC)
             WRITE(NRES,208) (NOMFIC(I),' (field map)',I=1,NFIC)
 208         FORMAT(10X,2A)
           ENDIF
        ELSE
          IF(MOD.EQ.24) THEN 
            WRITE(NRES,210) NOMFIC(MXFIC1)
          ELSE
            WRITE(NRES,210) (NOMFIC(I),I=1,NFIC)
 210        FORMAT(
     >      10X,'No  new  file  to  be  opened. Maps  already  stored.',/
     >      10X,'Skip  reading  file : ',10X,A80)
          ENDIF
        ENDIF
      ENDIF 


      IRD = NINT(A(NOEL,50))
      KZ = 1

      IF    (MOD .EQ. 22) THEN
C Keyword EMMA with *two* 2D maps, one for QF, one for QD, 
C such that total B is a linear superimposition of both

        IF(NEWFIC) THEN

             KFIC = 1
             IF(IDLUNI(
     >                 LUN)) THEN
               BINAR=BINARI(NOMFIC(KFIC),IB)
               IF(BINAR) THEN
                 OPEN(UNIT=LUN,FILE=NOMFIC(KFIC),FORM='UNFORMATTED'
     >           ,STATUS='OLD',ERR=96)
               ELSE
                 OPEN(UNIT=LUN,FILE=NOMFIC(KFIC),STATUS='OLD',ERR=96)
               ENDIF
             ELSE
               GOTO 96
             ENDIF
             CALL FMAPR2(BINAR,LUN,MOD,MOD2,NHD,BNORM,
     >                      BMIN,BMAX,
     >                      XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)
             DO JJJ=1,JYMA
               DO III=1,IXMA
                 DO IID = 1, ID    
                   HCA(IID,III,JJJ,KZ,1) = HC(IID,III,JJJ,KZ)
                 ENDDO
               ENDDO
             ENDDO
             CLOSE(UNIT=LUN)

             KFIC = 2
             IF(IDLUNI(
     >                 LUN)) THEN
               BINAR=BINARI(NOMFIC(KFIC),IB)
               IF(BINAR) THEN
                 OPEN(UNIT=LUN,FILE=NOMFIC(KFIC),FORM='UNFORMATTED'
     >           ,STATUS='OLD',ERR=96)
               ELSE
                 OPEN(UNIT=LUN,FILE=NOMFIC(KFIC),STATUS='OLD',ERR=96)
               ENDIF
             ELSE
               GOTO 96
             ENDIF
             CALL FMAPR2(BINAR,LUN,MOD,MOD2,NHD,BNORM,
     >                      BMIN,BMAX,
     >                      XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)
             DO JJJ=1,JYMA
               DO III=1,IXMA
                 DO IID = 1, ID    
                   HCB(IID,III,JJJ,KZ,1) = HC(IID,III,JJJ,KZ)
                 ENDDO
               ENDDO
             ENDDO
             CLOSE(UNIT=LUN)

C--------- Sum both ! 
             DO JJJ=1,JYMA
               DO III=1,IXMA
                 DO IID = 1, ID    
                   HC(IID,III,JJJ,KZ) = AF * HCA(IID,III,JJJ,KZ,1) 
     >                                + AD * HCB(IID,III,JJJ,KZ,1)
                 ENDDO
               ENDDO
             ENDDO
        ENDIF ! NEWFIC

      ELSEIF(MOD .EQ. 24) THEN
C Keyword EMMA with NFMP pairs of 2D maps (a pair is QFon and QDon
C such that total B is a linear superimposition of both). 
C NOMFIC(1) contains the names of the 2*NFMP field maps, and value of 
C parameter DSTP (quads' axis distance)

        IF(NEWFIC) THEN
     
C Open file that contains list of field map file names
             IF(IDLUNI(
     >                 LUN)) THEN
               OPEN(UNIT=LUN,FILE=NOMFIC(MXFIC1),STATUS='OLD',ERR=96)
             ELSE
               GOTO 96
             ENDIF
C Put file names into NOMFIC
             CALL GETNAM(LUN,MXDST, 
     >                             NOMFIC,NBDST)
             IF(NRES.GT.0) THEN
               WRITE(NRES,FMT='(/,5X,A,I3,A,/)') 
     >            ' EMMAP. GETNAM read  ',2*NBDST,' files : '
               WRITE(NRES,FMT='(A)') (NOMFIC(I),I=1,2*NBDST)
             ENDIF
C Put corresponding quads' axis distances into DISTL
             CALL GETDST(NOMFIC,NBDST,
     >                                DISTL)
             IF(NRES.GT.0) THEN
               WRITE(NRES,FMT='(/,5X,I4,A,80F7.2)') 
     >             NBDST,'  distances :  ',(DISTL(I),I=1,NBDST)
             ENDIF

        ENDIF  ! NEWFIC

C From the value of DIST, determines which two values in DISTL to use for interpolation
             IF    (DIST .LE. DISTL(1)) THEN
               KFIC = 1
C               WRITE(*,*) ' 1  DIST, DISTL(I): ',DIST, DISTL(1),KFIC
             ELSEIF(DIST .GE. DISTL(NBDST)) THEN
               KFIC = 2*NBDST-3
C               WRITE(*,*) ' 2  DIST, DISTL(I): ',DIST, DISTL(NBDST),KFIC
             ELSE
               DO I = 2, NBDST
                 IF(DIST .LE. DISTL(I)) THEN
                   KFIC = 2*I-3
C                   WRITE(*,*) ' 3  DIST, DISTL(I): ',DIST, DISTL(I),KFIC
                   GOTO 14
                 ENDIF
               ENDDO
             ENDIF
 14          CONTINUE

C                   WRITE(*,*) '****** DIST, AF, AD ',DIST, AF, AD 


             IF(NRES.GT.0) THEN
               WRITE(NRES,FMT='(/,5X,A,F7.2,A,I2)') 'Given that DIST=',
     >             DIST,',  then  KFIC = ',KFIC
               WRITE(NRES,FMT='(/,5X,A)') ' Files retained are :'
               WRITE(NRES,FMT='(5X,A)') (nomfic(KK),KK=KFIC,KFIC+3)
             ENDIF

C Open field map files one after the other and store field in HCA-B  tables
C             KFIC = 6           ! for test

C             CALL FITSTA(I5,
C     >                      FITING)
C                  write(*,*) ' +++++++++++ emmap fiting : ',fiting
C             IF(FITING) THEN
C               JFICA = 1
C               JFICB = 2*NBDST-1
C                  write(*,*) ' emmap jifca, jifcb 1 ', jfica, jficb 
C             ELSE
               JFICA = KFIC
               JFICB = JFICA+2
C                  write(*,*) ' emmap jifca, jifcb 2 ', jfica, jficb 
C                  write(*,*) ' kfic,kfic0 ', kfic,kfic0
C             ENDIF

           IF(KFIC.NE.KFIC0) THEN

             KFIC0 = KFIC

             JFIC = JFICA
 10          CONTINUE

               KDST = JFIC/2+1

               IF(IDLUNI(
     >                 LUN)) THEN
                 BINAR=BINARI(NOMFIC(JFIC),IB)
                 IF(BINAR) THEN
                   OPEN(UNIT=LUN,FILE=NOMFIC(JFIC),FORM='UNFORMATTED'
     >             ,STATUS='OLD',ERR=96)
                 ELSE
                   OPEN(UNIT=LUN,FILE=NOMFIC(JFIC),STATUS='OLD',ERR=96)
                 ENDIF
               ELSE
                 GOTO 96
               ENDIF
               CALL FMAPR2(BINAR,LUN,MOD,MOD2,NHD,BNORM,
     >                      BMIN,BMAX,
     >                      XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)
               IF(NRES.GT.0) 
     >           WRITE(NRES,FMT='(5X,A,2I3)') nomfic(JFIC),KDST,JFIC
               CLOSE(UNIT=LUN)
               DO JJJ=1,JYMA
                 DO III=1,IXMA
                   DO IID = 1, ID    
                     HCA(IID,III,JJJ,KZ,KDST) = HC(IID,III,JJJ,KZ)
                   ENDDO
                 ENDDO
               ENDDO

               JFIC = JFIC+1
C               write(*,*) '  emmap, nomfic 1 : ',NOMFIC(JFIC)
               IF(IDLUNI(
     >                   LUN)) THEN
                 BINAR=BINARI(NOMFIC(JFIC),IB)
                 IF(BINAR) THEN
                   OPEN(UNIT=LUN,FILE=NOMFIC(JFIC),FORM='UNFORMATTED'
     >             ,STATUS='OLD',ERR=96)
                 ELSE
                   OPEN(UNIT=LUN,FILE=NOMFIC(JFIC),STATUS='OLD',ERR=96)
                 ENDIF
               ELSE
                 GOTO 96
               ENDIF
               CALL FMAPR2(BINAR,LUN,MOD,MOD2,NHD,BNORM,
     >                        BMIN,BMAX,
     >                        XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)
               IF(NRES.GT.0) 
     >           WRITE(NRES,FMT='(5X,A,2I3)') nomfic(JFIC),KDST,JFIC
               CLOSE(UNIT=LUN)
               DO JJJ=1,JYMA
                 DO III=1,IXMA
                   DO IID = 1, ID    
                     HCB(IID,III,JJJ,KZ,KDST) = HC(IID,III,JJJ,KZ)
                   ENDDO
                 ENDDO
               ENDDO

               JFIC = JFIC + 1
C            WRITE(*,*) '+++JFIC,JFCA/B,KDST ',JFIC,JFICA,JFICB,KDST,dist

               IF(JFIC.LE.JFICB) GOTO 10
C-------------- 10 continue

           ELSE

             JFIC = JFICA
             KDST = JFIC/2+1

           ENDIF  ! KFIC  

C--------- Sum both QF and QD contributions into HCC
C           for both lower and upper DIST bounds
             KDST1 = KDST-1
             DO JJJ=1,JYMA
              DO III=1,IXMA
               DO IID = 1, ID    
                HCC(IID,III,JJJ,KZ,KDST1)=AF * HCA(IID,III,JJJ,KZ,KDST1) 
     >                                   +AD * HCB(IID,III,JJJ,KZ,KDST1)
               ENDDO
              ENDDO
             ENDDO
             DO JJJ=1,JYMA
              DO III=1,IXMA
               DO IID = 1, ID    
                HCC(IID,III,JJJ,KZ,KDST) = AF * HCA(IID,III,JJJ,KZ,KDST) 
     >                                   + AD * HCB(IID,III,JJJ,KZ,KDST)
               ENDDO
              ENDDO
             ENDDO
C--------- Interpolate the new field map that corresponds to actual DIST
             DO JJJ=1,JYMA
              DO III=1,IXMA
               DO IID = 1, ID    
                HC(IID,III,JJJ,KZ) = HCC(IID,III,JJJ,KZ,KDST1) + 
     >          (HCC(IID,III,JJJ,KZ,KDST)-HCC(IID,III,JJJ,KZ,KDST1)) / 
     >          (DISTL(KDST) - DISTL(KDST1)) * (DIST - DISTL(KDST1))
               ENDDO
              ENDDO
             ENDDO

      ELSE

        CALL ENDJOB('*** Error, SBR EMMAC -> No  such  MOD= ',MOD)

      ENDIF  ! MOD

C      ENDIF  ! NEWFIC

      CALL MAPLI1(BMAX-BMIN)
      AT=XH(IXMA)-XH(1)
      ATO = 0.D0
      ATOS = 0.D0
      RM=.5D0*(YH(JYMA)+YH(1))
      XI = XH(1)
      XF = XH(IXMA)

      RETURN
 
 96   IF(NRES.GT.0) WRITE(NRES,*) 'ERROR  OPEN  FILE ',NOMFIC(JFIC)
      CALL ENDJOB('ERROR  OPEN  FILE ',-99)

      RETURN
      END
