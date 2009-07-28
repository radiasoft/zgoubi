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
C  Fran�ois M�ot <meot@lpsc.in2p3.fr>
C  Service Acc�lerateurs
C  LPSC Grenoble
C  53 Avenue des Martyrs
C  38026 Grenoble Cedex
C  France
      SUBROUTINE EMMAC(SCAL,NDIM, 
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
      CHARACTER TITL*80 , NOMFIC(2)*80, NAMFIC(2)*80
      DATA NOMFIC / 2*'               '/ 
      SAVE NOMFIC, NAMFIC

      INTEGER DEBSTR,FINSTR

      DIMENSION HC1(ID,MXX,MXY,IZ), HC2(ID,MXX,MXY,IZ)
      DIMENSION HCA(ID,MXX,MXY,IZ), HCB(ID,MXX,MXY,IZ)
      SAVE HC1, HC2, HCA, HCB
      DIMENSION HCU(ID,MXX,MXY,IZ)
      SAVE HCU
      parameter (idmx=ID*MXX*MXY*IZ)
          data hcu /  IDMX * 0.d0 /

      BNORM = A(NOEL,10)*SCAL
      XNORM = A(NOEL,11)
      YNORM = A(NOEL,12)
      ZNORM = A(NOEL,13)
      TITL = TA(NOEL,1)
      IXMA = A(NOEL,20)
      IF(IXMA.GT.MXX) 
     >   CALL ENDJOB('X-dim of map is too large,  max  is ',MXX)
      JYMA = A(NOEL,21)
      IF(JYMA.GT.MXY ) 
     >   CALL ENDJOB('Y-dim of map is too large,  max  is ',MXY)
      IF(NDIM .EQ. 3) THEN
        KZMA =A(NOEL,22)
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
C such that resulting single map is a linear superimposition of both

          NFIC = 2
          NAMFIC(1) = TA(NOEL,NFIC)
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
C  QD_new is interpolated from QD with dr=xD, 
C  resulting single map is linear superimposition of both
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

        WRITE(NRES,FMT='(/,5X,2(A,I3),2A,I3,/)') 
     >  'NDIM = ',NDIM,' ;   Value of MOD is ', MOD,' ;  ', 
     >  'Number of field data files used is ',NFIC

        WRITE(NRES,FMT='(/,5X,A,1P,2E12.4,/)') 
     >  ' Value of field factor coefficients AF, AD are : ', AF, AD

        IF(MOD .EQ. 0) THEN
          WRITE(NRES,FMT='(/,5X,A,1P,E12.4,/)') 
     >    ' Distance between axis of quads : ', DIST
        ELSEIF(MOD .EQ. 1) THEN
          WRITE(NRES,FMT='(/,5X,A,1P,E14.6,A  )') 
     >    ' Field map 1 recomputed following radial move ',DIST,' (cm)'
          WRITE(NRES,FMT='(/,5X,A,1P,E14.6,A,/)') 
     >    ' Field map 2 recomputed following radial move ',DIST2,' (cm)'
        ENDIF

        IF(NEWFIC) THEN
           WRITE(NRES,209) 
 209       FORMAT(/,10X,' new field maps are now used ; ', 
     >     ' Names of map data files are : ')
           WRITE(6   ,208) (NOMFIC(I),I=1,NFIC)
           WRITE(NRES,208) (NOMFIC(I),I=1,NFIC)
 208       FORMAT(10X,A)
        ELSE
          WRITE(NRES,210) (NOMFIC(I),I=1,NFIC)
 210      FORMAT(
     >    10X,'No  new  map  file  to  be  opened. Already  stored.',/
     >    10X,'Skip  reading  field  map  file : ',10X,A80)
        ENDIF
      ENDIF 

C      IF(NEWFIC) THEN

        IF(MOD .EQ. 0) THEN
C------- Rectangular mesh
C Keyword EMMA with *two* 2D maps, one for QF at fixed xF, one for QD at fixed xD, 
C such that resulting single map is a linear superimposition of both

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

             CALL FMAPR3(BINAR,LUN,MOD,MOD2,BNORM,I1,KZ,
     >                      BMIN,BMAX,
     >                      XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)

               do jjj=1,JYMA
                 do iii=1,IXMA
                   do iid = 1, id    
                     HC1(iid,iii,jjj,KZ) = HC(iid,iii,jjj,KZ)
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

             CALL FMAPR3(BINAR,LUN,MOD,MOD2,BNORM,I1,KZ,
     >                      BMIN,BMAX,
     >                      XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)

               do jjj=1,JYMA
                 do iii=1,IXMA
                   do iid = 1, id    
                     HC2(iid,iii,jjj,KZ) = HC(iid,iii,jjj,KZ)
                   enddo
                 enddo
               enddo


          do jjj=1,JYMA
            do iii=1,IXMA
              do iid = 1, id    
                HC(iid,iii,jjj,KZ) = AF * HC1(iid,iii,jjj,KZ) 
     >                         + AD * HC2(iid,iii,jjj,KZ)
              enddo
            enddo
          enddo

        ELSEIF(MOD .EQ. 1) THEN
C--------- Rectangular mesh
C Keyword EMMA with *two* 2D maps, one for QF, one for QD, 
C the resulting single map is obtained in the following way : 
C  QF_new is interpolated from QF with dr=xF, 
C  QD_new is interpolated from QD with dr=xD, 
C  resulting single map is linear superimposition of both

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

             CALL FMAPR3(BINAR,LUN,MOD,MOD2,BNORM,I1,KZ,
     >                      BMIN,BMAX,
     >                      XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)

             DO JJJ=1,JYMA
               DO III=1,IXMA
                 DO IID = 1, ID    
                   HCA(IID,III,JJJ,KZ) = HC(IID,III,JJJ,KZ)
                 ENDDO
               ENDDO
             ENDDO
             CLOSE(UNIT=LUN)

            ENDIF ! NEWFIC

             CALL MAPSHF(HCA,XH,YH,DIST,IXMA,JYMA,
     >                                            HC1)

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

             CALL FMAPR3(BINAR,LUN,MOD,MOD2,BNORM,I1,KZ,
     >                      BMIN,BMAX,
     >                      XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)


             DO JJJ=1,JYMA
               DO III=1,IXMA
                 DO IID = 1, ID    
                   HCB(IID,III,JJJ,KZ) = HC(IID,III,JJJ,KZ)
                 ENDDO
               ENDDO
             ENDDO
             CLOSE(UNIT=LUN)

            ENDIF ! NEWFIC

             CALL MAPSHF(HCB,XH,YH,DIST2,IXMA,JYMA,
     >                                             HC2)

C Superimpose the two field maps into a single one
C          zer1 = 0.d0
          do jjj=1,JYMA
            do iii=1,IXMA
              do iid = 1, id    
                HC(iid,iii,jjj,KZ) = AF * HC1(iid,iii,jjj,KZ) 
     >                         + AD * HC2(iid,iii,jjj,KZ)
C                zer1 = zer1 + HC(iid,iii,jjj,KZ)
              enddo
            enddo
          enddo

C          zero = 0.d0
C          do jjj=1,JYMA
C            do iii=1,IXMA
C              do iid = 1, id    
C                zero = zero + (HCU(iid,iii,jjj,KZ) - HC(iid,iii,jjj,KZ))
C              enddo
C            enddo
C          enddo
C             write(87,*) dist,dist2,zero,zer1,' emmac '
C          do jjj=1,JYMA
C            do iii=1,IXMA
C              do iid = 1, id    
C                HCU(iid,iii,jjj,KZ) = AF * HC1(iid,iii,jjj,KZ) 
C     >                         + AD * HC2(iid,iii,jjj,KZ)
C              enddo
C            enddo
C          enddo

        ELSE

          CALL ENDJOB('*** Error, SBR EMMAC -> No  such  MOD= ',MOD)

        ENDIF

C      ENDIF

C Make sure this is ok given that field map is in cartesian frame
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
      END