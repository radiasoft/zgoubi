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
      USE DYNHC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-------------------------------------------------
C     Read TOSCA map with cartesian coordinates.
C     TOSCA keyword with MOD.le.19.
C-------------------------------------------------
C      LOGICAL NEWFIC(*)
      LOGICAL NEWFIC
      INCLUDE 'PARIZ.H'
      INCLUDE 'XYZHC.H'     ! COMMON// XH(MXX),YH(MXY),ZH(IZ),IXMA,JYMA,KZMA
      INCLUDE 'C.AIM_3.H'     ! COMMON/AIM/ ATO,AT,ATOS,RM,XI,XF,EN,EB1,EB2,EG1,EG2
      INCLUDE 'C.CDF.H'     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MAXTRA.H'
      INCLUDE 'C.CONST.H'     ! COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      INCLUDE 'C.DON.H'     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
C      PARAMETER (LNTA=132) ; CHARACTER(LNTA) TA
C      PARAMETER (MXTA=45)
      INCLUDE 'C.DONT.H'     ! COMMON/DONT/ TA(MXL,MXTA)
      INCLUDE 'C.DROITE_2.H'     ! COMMON/DROITE/ AM(9),BM(9),CM(9),IDRT
      INCLUDE 'C.INTEG.H'     ! COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
C      LOGICAL ZSYM
      INCLUDE 'C.TYPFLD.H'     ! COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      INCLUDE 'C.ORDRES.H'     ! COMMON/ORDRES/ KORD,IRD,IDS,IDB,IDE,IDZ
 
      LOGICAL BINARI,IDLUNI
      LOGICAL BINAR
      LOGICAL FLIP
      CHARACTER(LNTA) TITL , NAMFIC
      character(LNTA), allocatable :: NOMFIC(:)
      SAVE NOMFIC, NAMFIC
      INTEGER DEBSTR,FINSTR
      PARAMETER (NHDF=8)
 
      LOGICAL STRCON
 
      CHARACTER(20) FMTYP
C      DIMENSION XXH(MXX,MMAP), YYH(MXY,MMAP), ZZH(IZ,MMAP)
C      SAVE XXH, YYH, ZZH
C      DIMENSION BBMI(MMAP), BBMA(MMAP), XBBMI(MMAP), YBBMI(MMAP)
C      DIMENSION ZBBMI(MMAP), XBBMA(MMAP), YBBMA(MMAP), ZBBMA(MMAP)
C      SAVE BBMI, BBMA, XBBMI, YBBMI, ZBBMI, XBBMA, YBBMA, ZBBMA
C      DIMENSION IIXMA(MMAP), JJYMA(MMAP), KKZMA(MMAP)
C      SAVE IIXMA, JJYMA, KKZMA
      

      PARAMETER (MXHD=20)
      PARAMETER (ONE=1.D0)
  
      PARAMETER (MXC = 4)
      DIMENSION AA(24+MXC-1)

      LOGICAL NEWF
      LOGICAL ZROBXY      

      DATA ZROBXY / .FALSE. /
      DATA FMTYP / ' regular' /
      DATA AA / 27 * 0.D0 /
      DATA NEWF / .TRUE. /

      call allocate_if_unallocated(NOMFIC, IZ)

C      BNORM = A(NOEL,10)*SCAL
      BNORM = A(NOEL,10)
      XNORM = A(NOEL,11)
      YNORM = A(NOEL,12)
      ZNORM = A(NOEL,13)
      IF(ZNORM .EQ. 0.D0) ZNORM = 1.D0      
      TITL = TA(NOEL,1)
      IF    (STRCON(TITL,'HEADER',
     >                              IS) ) THEN
        READ(TITL(IS+7:IS+7),FMT='(I1)') NHD
      ELSE
        NHD = NHDF
      ENDIF
      IF(NHD .GT. MXHD) CALL ENDJOB(
     >'Pgm map2d. Field map header has too many lines. Must be .le.'
     >,MXHD)
      IDEB = DEBSTR(TITL)
C      FLIP = TITL(IDEB:IDEB+3).EQ.'FLIP'
      IXMA = NINT(A(NOEL,20))

      IF(IXMA.GT.MXX)
     >   CALL ENDJOB('X-dim of map is too large,  max  is ',MXX)
      JYMA = NINT(A(NOEL,21))
      IF(JYMA.GT.MXY )
     >   CALL ENDJOB('Y-dim of map is too large,  max  is ',MXY)
 
             MOD = 0
             MOD2 = 0
C FM - Mar 2016
C             MOD2 = 3

      KZMA = 1
      NFIC = 1
      NAMFIC = TA(NOEL,2)
      NOMFIC(NFIC) = NAMFIC(DEBSTR(NAMFIC):FINSTR(NAMFIC))
 
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
                  FMTYP = 'LESB3'
             ELSE
                  FMTYP = ' regular'
             ENDIF

      CALL KSMAP4(NOMFIC,NFIC,AA(24:24+MXC-1),
     >                        NEWFIC,NBMAPS,IMAP)

      IF(NRES.GT.0) THEN

               WRITE(NRES,FMT='(/,3A,/)')
     >           (NOMFIC(I)(DEBSTR(NOMFIC(I)):
     >    FINSTR(NOMFIC(I))),I=1,NFIC), 
     >           ' map,  FORMAT type : ', FMTYP
                CALL FLUSH2(NRES,.FALSE.)

        WRITE(NRES,FMT='(/,5X,2(A,I3,A),/)')
     >  'Number of data file sets used is ',NFIC,' ;  '
     >  ,'Stored in field array # IMAP =  ',IMAP,' '

        DO I = 1, NFIC
C          IF(NEWFIC(I)) THEN
          IF(NEWFIC) THEN
           WRITE(NRES,209)NOMFIC(I)(DEBSTR(NOMFIC(I)):FINSTR(NOMFIC(I)))
 209       FORMAT(/,5X
     >     ,'New field map, cartesian mesh (MOD.le.19) ; '
     >     ,'  name of map data file : ',/,10X,A)
          ELSE
           WRITE(NRES,210)NOMFIC(I)(DEBSTR(NOMFIC(I)):FINSTR(NOMFIC(I)))
 210       FORMAT(5X,
     >     'Skip  reading  already  stored  field  map  file  ',A)
          ENDIF
        ENDDO
      ENDIF
 
      INDEX=0
      NT = 1
      CALL PAVELW(INDEX,NT)
 
      NEWF = .FALSE.
      I = 1
      DO WHILE(.NOT. NEWF .AND. I.LE.NFIC)
C        NEWF = NEWFIC(I)
        NEWF = NEWFIC
        I = I+1
      ENDDO
 
      IF(NEWF) THEN
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
 
             I1 = 1
             KZ = 1
             IRD = NINT(A(NOEL,40))
 
             CALL FMAPR3(BINAR,LUN,MOD,MOD2,NHD,ZROBXY,
     >                   XNORM,YNORM,ZNORM,ONE,I1,KZ,FMTYP,
     >                                    BMIN,BMAX,
     >                                    XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)

C------- Store mesh coordinates
        CALL FMAPW4(IMAP,BMIN,BMAX,XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)

      ELSE
 
C------- Restore mesh coordinates
        CALL FMAPR5(IMAP,
     >                   BMIN,BMAX,XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)
 
        IF(NRES.GT.0) WRITE(NRES,*) ' Pgm map2d, ',
     >  ' restored mesh coordinates for field map # ',imap
 
      ENDIF
 
      CALL CHAMK2(BNORM*SCAL)
 
      CALL MAPLI1(BMAX-BMIN)

      RETURN
 
 96   WRITE(ABS(NRES),*) 'Pgm map2d. Error  open  file ',
     >NOMFIC(I)(DEBSTR(NOMFIC(I)):FINSTR(NOMFIC(I)))
      CALL ENDJOB('Leaving. ',-99)
 
      RETURN
      contains
        subroutine allocate_if_unallocated(string_array,array_size)
          character(len=*), allocatable, intent(inout):: string_array(:)
          integer, intent(in) :: array_size
          integer istat
          if (.not. allocated(string_array)) then 
            allocate(string_array(array_size), stat=istat)
            if (istat/=0) error stop "string_array allocation failed"
          end if 
        end subroutine
      END
