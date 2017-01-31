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
C  Upton, NY, 11973, USA
C  -------
      SUBROUTINE BREVOL(SCAL,
     >                      BMIN,BMAX,BNORM,XNORM,
     >               XBMI,XBMA,NEWFIC)
      USE DYNHC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-------------------------------------------------
C     Read TOSCA map with cartesian coordinates.
C     TOSCA keyword with MOD.le.19.
C-------------------------------------------------
C      LOGICAL NEWFIC(*)
      LOGICAL NEWFIC
      INCLUDE 'PARIZ.H'
      INCLUDE "XYZHC.H"     ! COMMON// XH(MXX),YH(MXY),ZH(IZ),IXMA,JYMA,KZMA
      INCLUDE 'C.AIM_3.H'     ! COMMON/AIM/ ATO,AT,ATOS,RM,XI,XF,EN,EB1,EB2,EG1,EG2
      INCLUDE 'C.CDF.H'     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MAXTRA.H'
      INCLUDE 'C.CONST.H'     ! COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      INCLUDE 'C.DON.H'     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE 'C.DONT.H'     ! COMMON/DONT/ TA(MXL,MXTA)
      INCLUDE 'C.DROITE_2.H'     ! COMMON/DROITE/ AM(9),BM(9),CM(9),IDRT
      INCLUDE 'C.INTEG.H'     ! COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      INCLUDE 'C.TYPFLD.H'     ! COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      INCLUDE 'C.ORDRES.H'     ! COMMON/ORDRES/ KORD,IRD,IDS,IDB,IDE,IDZ
 
      LOGICAL BINARI,IDLUNI
      LOGICAL BINAR
      LOGICAL FLIP
      CHARACTER(LNTA) TITL , NOMFIC(IZ)
      SAVE NOMFIC
      INTEGER DEBSTR,FINSTR
      SAVE NHDF
 
      LOGICAL STRCON, OK, NEWF
 
C      CHARACTER(20) FMTYP
      DIMENSION XXH(MXX,MMAP), YYH(MXY,MMAP), ZZH(IZ,MMAP)
      SAVE XXH, YYH, ZZH
      DIMENSION BBMI(MMAP), BBMA(MMAP), XBBMI(MMAP), YBBMI(MMAP)
      DIMENSION ZBBMI(MMAP), XBBMA(MMAP), YBBMA(MMAP), ZBBMA(MMAP)
      SAVE BBMI, BBMA, XBBMI, YBBMI, ZBBMI, XBBMA, YBBMA, ZBBMA
 
      DIMENSION IIXMA(MMAP), JJYMA(MMAP), KKZMA(MMAP)
      SAVE IIXMA, JJYMA, KKZMA

      CHARACTER(LNTA) STRA(2) 

      DATA NOMFIC / IZ*' '/
      DATA NHDF / 8 /
 
      PARAMETER (MXC = 4)
      DIMENSION AA(24+MXC-1)
      DATA NEWF / .TRUE. /

      BNORM = A(NOEL,10)
      XNORM = A(NOEL,11)
      ZNORM = 1.D0      
      TITL = TA(NOEL,1)
      IF    (STRCON(TITL,'HEADER',
     >                              IS) ) THEN
        READ(TITL(IS+7:IS+7),FMT='(I1)') NHD
      ELSE
        NHD = NHDF
      ENDIF
      IDEB = DEBSTR(TITL)
      FLIP = TITL(IDEB:IDEB+3).EQ.'FLIP'
      IXMA = NINT(A(NOEL,20))
      IF(IXMA.GT.MXX)
     >   CALL ENDJOB('X-dim of map is too large,  max  is ',MXX)
 
      DO I = 1, IZ
        NOMFIC(I) = ' '
      ENDDO
      CALL RAZ(AA,24+MXC-1)
      KZMA = 1
      NFIC = 1
      CALL STRGET(TA(NOEL,1+NFIC),2,
     >                              IDUM,STRA)        
      NOMFIC(NFIC) = STRA(1)(DEBSTR(STRA(1)):FINSTR(STRA(1)))
      DO WHILE(IDUM.EQ.2 .AND. STRA(2) .EQ. 'SUM')
        NFIC = NFIC + 1
        CALL STRGET(TA(NOEL,1+NFIC),2,
     >                                IDUM,STRA)        
        NOMFIC(NFIC) = STRA(1)(DEBSTR(STRA(1)):FINSTR(STRA(1)))
      ENDDO

      CALL KSMAP4(NOMFIC,NFIC,AA(24:24+MXC-1),
     >                          NEWFIC,NBMAPS,IMAP)

      IF(NRES.GT.0) THEN
        WRITE(NRES,FMT='(/,5X,2(A,I3,A),/)')
     >   'Number of data files used is ',NFIC,' ;  '
     >  ,'Stored in field array # IMAP =  ',IMAP,' '

          WRITE(NRES,*) 
     >    ' Will sum up ',NFIC,'  field maps.'
 
        IF(NEWFIC) THEN
           WRITE(NRES,209)
 209       FORMAT(/,10X
     >     ,' New field map(s) now opened ; '
     >     ,' name(s) of map data file(s) : ',/)
           WRITE(NRES,208) (NOMFIC(I)(DEBSTR(NOMFIC(I)):
     >     FINSTR(NOMFIC(I))),I=1,NFIC)
 208       FORMAT(10X,A)
        ELSE
          WRITE(NRES,210) (NOMFIC(I)(DEBSTR(NOMFIC(I)):
     >     FINSTR(NOMFIC(I))),I=1,NFIC)
 210      FORMAT(
     >    10X,'No  new  map  file  to  be  opened. Already  stored.',/
     >    10X,'Skip  reading  field  map  file : ',10X,A)
        ENDIF
        CALL FLUSH2(NRES,.FALSE.)
      ENDIF
 
      IF(NEWFIC) THEN

        CALL RAZ(HC,ID*IXMA)
        JFIC = 0
        DO WHILE (JFIC .LT. NFIC)
          JFIC = JFIC+1

          OK = IDLUNI(
     >                LUN)
          IF(.NOT. OK) GOTO 96
          BINAR=BINARI(NOMFIC(JFIC),IB)
          IF(BINAR) THEN
            OPEN(UNIT=LUN,FILE=NOMFIC(JFIC),FORM='UNFORMATTED'
     >      ,STATUS='OLD',ERR=96)
          ELSE
            OPEN(UNIT=LUN,FILE=NOMFIC(JFIC),STATUS='OLD',ERR=96)
          ENDIF

          DO I = 1,IXMA
            IF(BINAR) THEN
              READ(LUN) XH(I),BMES
            ELSE
              READ(LUN,*) XH(I),BMES
            ENDIF
            IF    (BMES .GT. BMAX) THEN
              BMAX = BMES
              XBMA = XH(I)
            ELSEIF(BMES .LT. BMIN) THEN
              BMIN = BMES
              XBMI = XH(I)
            ENDIF
C--------------------- Rustine collecte e+ / Jab
C                  fac=1.D0
CCCCCCCCCCCCCCCCCCC   if(jfic.eq.2) fac = 3.5d0    ! pour obtenir le meme champ que Jab dans le 2eme soleno
C                HC(ID,I,1,1,IMAP) = HC(ID,I,1,1,IMAP) + BMES * BNORM *FAC
C--------------------------------------------------------

            HC(ID,I,1,1,IMAP) = HC(ID,I,1,1,IMAP) + BMES
            XH(I) = XH(I) * XNORM
  
          ENDDO

          CLOSE(LUN)
C--------- Will sum (superimpose) 1D field maps if 'SUM' follows map file name in zgoubi.dat. 
C          SUMAP is .T.

        ENDDO

C          Motion in this lmnt has no z-symm. 
        ZSYM=.FALSE.

C------- Store mesh coordinates
        CALL FMAPW4(IMAP,BMIN,BMAX,XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)

      ELSE
 
C------- Restore mesh coordinates
        CALL FMAPR5(IMAP,
     >                   BMIN,BMAX,XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)
        ZSYM=.FALSE.

        IF(NRES.GT.0) WRITE(NRES,*) ' Pgm brevol, ',
     >  ' restored mesh coordinates for field map # ',imap
      ENDIF

      CALL CHAMK2(BNORM*SCAL)
      
      CALL MAPLI1(BMAX-BMIN)
 
      RETURN
 
 96   WRITE(ABS(NRES),*) 'Pgm brevol. Error  open  file ',
     >NOMFIC(JFIC)(DEBSTR(NOMFIC(JFIC)):FINSTR(NOMFIC(JFIC)))
      CALL ENDJOB('Leaving. ',-99)
 
      RETURN
      END
