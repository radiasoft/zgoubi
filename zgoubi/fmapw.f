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
      SUBROUTINE FMAPW(ACN,RFR,KART)
      USE DYNHC   ! HC(ID,MXX,MXY,IZ,IMAP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARIZ.H'     ! PARAMETERS : IZ, ID, MMAP, MXX, MXY
      INCLUDE 'XYZHC.H'     ! COMMON// XH(MXX),YH(MXY),ZH(IZ),IXMA,JYMA,KZMA
      INCLUDE 'C.CDF.H'     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'C.CONST.H'     ! COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'C.DROITE.H'     ! COMMON/DROITE/ CA(9),SA(9),CM(9),IDRT
C      LOGICAL ZSYM
      INCLUDE 'C.TYPFLD.H'     ! COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM

      DIMENSION BREAD(3)
C----- at entry FMAPR :
C      LOGICAL IDLUNI, BINARI
      LOGICAL BINAR, EMPTY
      CHARACTER(120) HDRM
      CHARACTER(20) FMTYP
      PARAMETER (MXCHAR=200)
      CHARACTER(200) TXT132
      CHARACTER(1) BE(2)
      INTEGER DEBSTR, FINSTR
      SAVE MOD, MOD2
CC 3D field maps, intermediate storage prior to summing
C      PARAMETER(MX3D=3)
C      DOUBLE PRECISION, DIMENSION(:,:,:,:,:), ALLOCATABLE :: 
C     >     HCTMP
CC      DIMENSION HCTMP(ID,MXX,MXY,IZ,MX3D)
C      SAVE HCTMP

      SAVE AA, IFIC

      DIMENSION IIXMA(MMAP),XXH(MXX,MMAP),JJYMA(MMAP),YYH(MXY,MMAP),
     >KKZMA(MMAP),ZZH(IZ,MMAP),BBMI(MMAP),BBMA(MMAP),XBBMI(MMAP),
     >YBBMI(MMAP),ZBBMI(MMAP),XBBMA(MMAP),YBBMA(MMAP),ZBBMA(MMAP)
      SAVE IIXMA,XXH,JJYMA,YYH,
     >KKZMA,ZZH,BBMI,BBMA,XBBMI,
     >YBBMI,ZBBMI,XBBMA,YBBMA,ZBBMA

      LOGICAL ZROBXY   ! Set to .TRUE. to force BX==0, BY==0, at Z=0 (some OPERA maps have residual 
                       ! non-zero, that causes particles in mid-plane to wrongly acquire Z-motion). 

      DATA BE /'B','E'/

      DATA MOD, MOD2 / 0, 0 /
      DATA IMAP / 1 /
      DATA R0 / 0.D0 /

      CALL KSMAP(
     >           IMAP) 

      IF(LF .EQ. 1) THEN
C------- Print field map in zgoubi.res
        IF(NRES.GT.0) THEN

          K = KZMA/2+1
          WRITE(NRES,
     >    FMT='(/,1X,20(''-''),/,10X,'' FIELD  MAP  (NORMALISED) :'',/)'
     $         )
          WRITE(NRES,FMT='(/,A,3(I5,1X),/)') '  IXMA, JYMA, KZMA : ',
     >    IXMA,JYMA,KZMA
          WRITE(NRES,FMT='(T3,''Y'',T16,''Z'',T29,''X'',
     >    T42,A,''Y'',T55,A,''Z'',T68,A,''X'')') 
     >    BE(KFLD),BE(KFLD),BE(KFLD)

          CALL RAZ(BREAD,3)

          DO J=1,JYMA
            DO I = 1,IXMA
              DO JD=1, ID
                BREAD(JD) = HC(JD,I,J,K,IMAP)
              ENDDO
              WRITE(NRES,FMT='(1X,1P,6(E12.4,1X),3(I5,1X))') 
     >        YH(J),ZH(K),XH(I), BREAD(2), BREAD(3), BREAD(1)
     >        ,J,K,I
            ENDDO
          ENDDO

        ENDIF
         
      ELSEIF(LF .EQ. 2) THEN
C------- Print field map in zgoubi.map

        CALL OPEN2('FMAPW',NMAP,'zgoubi.map') 

        K = KZMA/2+1
        
C Replaces the last 2 lines of header (as written by OPEN2) : 
        BACKSPACE(NMAP)
        BACKSPACE(NMAP)
        WRITE(NMAP,FMT='(1X,A11,1X,I1,1P,3G15.7)')
     >  'boundaries:',IDRT,(CA(I),SA(I),CM(I),I=1,IDRT)
        WRITE(NMAP,992) IXMA,JYMA,ACN,RFR,KART,MOD
 992    FORMAT(2(I5,1X),1P,2(E15.7,1X),1X,I1,1X,I3)

C Write the map :
        WRITE(NMAP,FMT='(A)') ' X data '
        WRITE(NMAP,991) (XH(I),I=1,IXMA)
        WRITE(NMAP,FMT='(A)') ' Y data '
        WRITE(NMAP,991) (YH(J),J=1,JYMA)
        WRITE(NMAP,FMT='(A)') ' H data '
        WRITE(NMAP,991) ((HC(ID,I,J,K,IMAP),I=1,IXMA),J=1,JYMA)
C        WRITE(NMAP,991) (XH(I),I=1,IXMA),(YH(J),J=1,JYMA)
C     >  ,((HC(ID,I,J,K,IMAP),I=1,IXMA),J=1,JYMA)
 991    FORMAT(1P,6E18.10)

        CLOSE(NMAP)

        IF(NRES.GT.0) WRITE(NRES,FMT='(/,10X,A,/)') 
     >        'Field map has been written in zgoubi.map' 

      ENDIF
      RETURN

C Called by POLMES
      ENTRY FMAPR(BINAR,LUN,MODI,MODI2,NHDI,
     >                                      RM)
      MOD = MODI
      MOD2 = MODI2
      NHD = NHDI

      CALL KSMAP(
     >           IMAP) 

      IF    (MOD2 .EQ. 0) THEN
C Ex. : Calabretta H2 field map.

        IF(BINAR) THEN

C Map data file starts with NHD-line header
C Line NHD contains IXM,JYM,ACENT,RM,KRT,MD
          DO II=1, NHD-1
            READ(LUN) HDRM
          ENDDO

          READ(LUN) IXM,JYM,ACENT,RM,KRT,MD
   
          READ(LUN) (XH(I),I=1,IXMA),(YH(J),J=1,JYMA)
     >    ,((HC(ID,I,J,1,IMAP),I=1,IXMA),J=1,JYMA)

        ELSE 

C Map data file starts with NHD-line header
C Line NHD contains IXM,JYM,ACENT,RM,KRT,MD
          DO II=1, NHD-1
            READ(LUN,FMT='(A)') HDRM
          ENDDO

          READ(LUN,*) IXM,JYM,ACENT,RM,KRT,MD

          IF(NRES.GT.0) THEN
            IF(IXM .NE. IXMA) THEN
              WRITE(6,100) 'IXMA',IXM,'IXMA',IXMA
              WRITE(NRES,100) 'IXMA',IXM,'IXMA',IXMA
 100          FORMAT(/,'SBR fmapw. WARNING - info !! ',A,' = ',I6, 
     >        ' as read in map data file .ne. ',A,' = ',I6,
     >        ' given in zgoubi.dat file.',/)
            ENDIF
            IF(JYM .NE. JYMA) THEN
              WRITE(6,100) 'JYMA','JYMA'
              WRITE(NRES,100) 'JYMA','JYMA'
            ENDIF
          ENDIF

          READ(LUN,*,END=229,ERR=98) (XH(I),I=1,IXMA),(YH(J),J=1,JYMA)
     >             ,((HC(ID,I,J,1,IMAP),I=1,IXMA),J=1,JYMA)

c          write(88,*) ' id, imap = ', id, imap
c          write(88,*)((XH(I),YH(J),HC(ID,I,J,1,IMAP),I=1,IXMA),J=1,JYMA)

 229      CONTINUE 
          GOTO 99

        ENDIF

      ELSEIF(MOD2 .EQ. 1) THEN
C Input data formatting is that of MIT Megatron field map 
C Map data file starts with NHD-line header
C Line NHD contains IXM,JYM,ACENT,RM

        IF(BINAR) THEN

          DO II=1, NHD-1
            READ(LUN) HDRM
          ENDDO
   
          READ(LUN,*) IXM,JYM,ACENT,RM

          READ(LUN) (XH(I),I=1,IXMA),(YH(J),J=1,JYMA)
     >    ,((HC(ID,I,J,1,IMAP),I=1,IXMA),J=1,JYMA)

        ELSE
 
          DO II=1, NHD-1
            READ(LUN,FMT='(A)') HDRM
          ENDDO

          READ(LUN,*) IXM,JYM,ACENT,RM

          ijma = IXMA*JYMA
          ii = 1
          i = 1
          j = 1
          do while (ii .le. ijma)
            READ(LUN,*,END=230,ERR=98) dum,dum,dum,dum,dum,
     >      HC(ID,I,J,1,IMAP),XH(I),YH(J)
c            write(88,fmt='(2i4,1p,3e17.6,3i6)') 
c     >        i,j,XH(I),YH(J),HC(ID,I,J,1,IMAP),ii,id,imap
            if(i.lt.ixma) then
              i = i+1
            else
              i = 1
              j = j+1
            endif
            ii = ii+1
          enddo

 230        CONTINUE 
          GOTO 99

        ENDIF

      ELSEIF(MOD2 .EQ. 2) THEN
C Input data formatting is that of Carol's FFAG - Sept. 2012
C Map data file starts with NHD-line header
C Line NHD contains IXM,JYM,ACENT,RM

        IF(BINAR) THEN

          DO II=1, NHD-1
            READ(LUN) HDRM
          ENDDO
   
          READ(LUN,*) IXM,JYM,ACENT,RM

          READ(LUN) (XH(I),I=1,IXMA),(YH(J),J=1,JYMA)
     >    ,((HC(ID,I,J,1,IMAP),I=1,IXMA),J=1,JYMA)

        ELSE
 
          DO II=1, NHD-1
            READ(LUN,FMT='(A)') HDRM
          ENDDO

          READ(LUN,*) IXM,JYM,ACENT,RM

          ijma = IXMA*JYMA
          ii = 1
          i = 1
          j = 1
          do while (ii .le. ijma)
            READ(LUN,*,END=231,ERR=98) dum,dum,dum,dum,dum,
     >      HC(ID,I,J,1,IMAP),XH(I),YH(J)
c            write(88,fmt='(2i4,1p,3e17.6,3i6)') 
c     >        i,j,XH(I),YH(J),HC(ID,I,J,1,IMAP),ii,id,imap
            if(i.lt.ixma) then
              i = i+1
            else
              i = 1
              j = j+1
            endif
            ii = ii+1
          enddo

 231      CONTINUE 
          GOTO 99

        ENDIF

      ELSE
        CALL ENDJOB('SBR fmapw. No such option MOD2=',IMOD)   
      ENDIF

      RETURN

C----------------------------------------------------------
C Read and interprete field maps in polar frame (MOD >= 20)
      ENTRY FMAPR2(BINAR,LUN, MODI,MODI2,NHDI,ZROBXY,
     >             XNORM,YNORM,ZNORM,BNORM,I1,KZI,FMTYP,
     >                       BMIN,BMAX,
     >                       XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)
      CALL KSMAP(
     >           IMAP) 

      MOD = MODI
      MOD2 = MODI2
      NHD = NHDI

      BMIN =  1.D10
      BMAX = -1.D10

      IF(NRES.GT.0) 
     >   WRITE(NRES,FMT='(A,I1,A)') ' HEADER  (',NHD,' lines) : '

      IF(BINAR) THEN
        IF(NHD .GE. 1) THEN
C Map data file starts with NHD-line header (NORMALLY 8)
          READ(LUN,END=97,ERR=98) HDRM
          IF(NRES.GT.0) WRITE(NRES,FMT='(5X,A120)') HDRM
          READ(HDRM,*,ERR=41,END=41) R0, DR, DTTA, DZ
          GOTO 42
 41       CONTINUE
          IF(NRES.GT.0) WRITE(NRES,*) 
     >    'Could not read R_min, DR, DTTA, DZ from line 1 of HEADER...'
 42       CONTINUE
          DO II=1, NHD-1
            READ(LUN,END=97,ERR=98) HDRM
            IF(NRES.GT.0) WRITE(NRES,FMT='(5X,A120)') HDRM
          ENDDO
          IF(NRES.GT.0) WRITE(NRES,FMT='(A)') ' '
          IF(NRES.GT.0) WRITE(NRES,fmt='(a,1p,4e14.6,/)')
     >    'R_min (cm), DR (cm), DTTA (deg), DZ (cm) : ',R0,DR,DTTA,DZ
        ELSE
C          CALL ENDJOB('Please give at least 1 line header in field map'
C     >    //' file with R_min, DR, DTTA, DZ values.',-99)
        ENDIF
      ELSE
          
        IF(NHD .GE. 1) THEN
C Map data file starts with NHD-line header (NORMALLY 8)
          READ(LUN,FMT='(A120)',END=97,ERR=98) HDRM
          IF(NRES.GT.0) WRITE(NRES,FMT='(5X,A120)') HDRM
          READ(HDRM,*,ERR=39,END=39) R0, DR, DTTA, DZ
          GOTO 40
 39       CONTINUE
          IF(NRES.GT.0) WRITE(NRES,*) 
     >    'Could not read R_min, DR, DTTA, DZ from line 1 of HEADER...'
 40       CONTINUE
          DO II=1, NHD-1
            READ(LUN,FMT='(A120)',END=97,ERR=98) HDRM
            IF(NRES.GT.0) WRITE(NRES,FMT='(5X,A120)') HDRM
          ENDDO
          IF(NRES.GT.0) WRITE(NRES,FMT='(A)') ' '
          IF(NRES.GT.0) WRITE(NRES,fmt='(a,1p,4e14.6,/)')
     >    'R_min (cm), DR (cm), DTTA (deg), DZ (cm) : ',R0,DR,DTTA,DZ
C     >                 'R_min, DR, DTTA, DZ : ',R0,DR,DTTA,DZ
        ELSE
C          CALL ENDJOB('Please give at least 1 line header in field map'
C     >    //' file with R_min, DR, DTTA, DZ values.',-99)
        ENDIF
      ENDIF

      DTTA = DTTA * RAD
             
      IF(MOD .EQ. 20 .OR. MOD .EQ. 21) THEN
C---------- Read the map. 
C Partial field map, special symmetrization applied, as follows : 
C  KEK case here (Aiba) :  position (x,y,z) and field components 
C  (Bx,By,Bz) are given in 
C  a XYZ frame centered at O=geometrical center of the magnet with  
C  Ox axis along the radius  (Ox is the symmetry axis of the triplet).  
C  Oz axis is normal to radius at O and Y is the vertical axis. 
C    The corresponding cylindrical frame is centered at O, so that the 
C  x,z position in an horizontal plane in map data relates to tta, r by  
C  x=r cos(tta),  z=r sin(tta). 
C  Coordinates x,z,y in map data file vary from fastest to slowest respectively, 
C  according to tta, r increments and to z increment. 
C  Cylindrical oordinates of the mesh nodes are computed from dR, dtheta.
 
C----------- Read the bare map data
           IRMA = JYMA
           JTMA = IXMA/2+1
           K2ZMA = KZMA/2+1

           DO 33 I=1,IRMA            
             RADIUS=R0 + DBLE(I-1) *DR 
             DO 33  K = 1,K2ZMA      
               KZC = K2ZMA-1+K
               Z = DBLE(K-1) * DZ
               DO 33 J=1,JTMA        
                 IF    (MOD .EQ. 20) THEN
                   JTC = J
                 ELSEIF(MOD .EQ. 21) THEN
                   JTC = JTMA-1+J
                 ENDIF
                 TTA =  DBLE(J-1) * DTTA
                 IF(BINAR) THEN
                   READ(LUN,END=97,ERR=98) DUM1,DUM2,DUM3,
     >                                   BREAD(1),BREAD(2),BREAD(3)
                 ELSE
                   READ(LUN,*,END=97,ERR=98) DUM1,DUM2,DUM3,
     >                                   BREAD(1),BREAD(2),BREAD(3)
C                   write(33,*) DUM1,DUM2,DUM3,BREAD(1),BREAD(2),BREAD(3)
                 ENDIF

                 BMAX0 = BMAX
                 BMAX = DMAX1(BMAX,BREAD(1),BREAD(2),BREAD(3))
                 IF(BMAX.NE.BMAX0) THEN
                   XBMA = TTA
                   YBMA = RADIUS
                   ZBMA = Z
                 ENDIF
                 BMIN0 = BMIN
                 BMIN = DMIN1(BMIN,BREAD(1),BREAD(2),BREAD(3))
                 IF(BMIN.NE.BMIN0) THEN
                   XBMI = TTA
                   YBMI = RADIUS
                   ZBMI = Z
                 ENDIF

C---------------- Horizontal component of field                 
                 BH = SQRT(BREAD(1)*BREAD(1) + BREAD(3)*BREAD(3))
                 ALP = ATAN2(BREAD(3),BREAD(1)) - TTA
C---------------- Polar components at positon (r,tta) in cyl. frame
                 BRAD = BH * COS(ALP)
                 BTTA = BH * SIN(ALP)
                 HC(1,JTC,I,KZC,IMAP) = BTTA * BNORM
                 HC(2,JTC,I,KZC,IMAP) = BRAD * BNORM
                 HC(3,JTC,I,KZC,IMAP) = BREAD(2) * BNORM
 33        CONTINUE

        BMIN = BMIN * BNORM
        BMAX = BMAX * BNORM

C------- symmetrise 3D map wrt magnet vertical symm plane
        DO 34  K=1,K2ZMA    
          KZC = K2ZMA-1+K
          IF(KZC.GT.IZ) CALL ENDJOB('Pgm FMAPW : K should be <',IZ+1)
          DO 34 I=1,IRMA      
            IF(I.GT.MXY) CALL ENDJOB('Pgm FMAPW : I should be <',MXY+1)
            DO 34 J=1,JTMA-1    
              IF    (MOD .EQ. 20) THEN
                JTC = J
                JTS = 2*JTMA-J
              ELSEIF(MOD .EQ. 21) THEN
                JTC = JTMA+J
                JTS = JTMA-J
              ENDIF
              IF(JTS.GT.MXX) CALL ENDJOB(' FMAPW, J SHOULD BE <',MXX+1)
              IF(JTC.GT.MXX) CALL ENDJOB(' FMAPW, J SHOULD BE <',MXX+1)
              HC(1,JTS,I,KZC,IMAP) = -HC(1,JTC,I,KZC,IMAP)
              HC(2,JTS,I,KZC,IMAP) =  HC(2,JTC,I,KZC,IMAP) 
              HC(3,JTS,I,KZC,IMAP) =  HC(3,JTC,I,KZC,IMAP) 
 34     CONTINUE

C------- symmetrise 3D map wrt mid-plane= bend-plane
        DO 35  K=2,K2ZMA      
          KZC = K2ZMA-1+K
          KZS = K2ZMA-K+1
          DO 35 I=1,IRMA         
            DO 35 J=1,JTMA*2-1    
              HC(1,J,I,KZS,IMAP) = -HC(1,J,I,KZC,IMAP)
              HC(2,J,I,KZS,IMAP) = -HC(2,J,I,KZC,IMAP) 
              HC(3,J,I,KZS,IMAP) = HC(3,J,I,KZC,IMAP) 
 35     CONTINUE

C------- Mesh coordinates
        DO J=1,JYMA
          YH(J) =  R0 + DBLE(J-1) *DR 
        ENDDO
        DO K= -KZMA/2,KZMA/2   
          ZH(K+KZMA/2+1) = DBLE(K) * DZ
        ENDDO
        DO I=-IXMA/2,IXMA/2
          XH(I+IXMA/2+1) =  DBLE(I) * DTTA
        ENDDO


      ELSEIF(MOD .EQ. 24) THEN
C---------- Read the map. 
C A single file contains the full 3D map. No symmetrization applied. 
C The map file gives theta, r, Z, Btheta, Br, BZ. 
C Cylindrical oordinates theta, r, Z of the mesh nodes are re-computed
C to ensure accuracy. 

C----------- Read the bare map data
           ITMA = IXMA
           JRMA = JYMA

           DO J=1,JRMA            
             RADIUS=R0 + DBLE(J-1) *DR     ! cm
             DO K = 1,KZMA      
               Z = DBLE(K-1) * DZ          ! cm
               DO I=1,ITMA        
                 TTA =  DBLE(I-1) * DTTA         ! rad
                 IF(BINAR) THEN
                   READ(LUN,END=97,ERR=98) DUM1,DUM2,DUM3,
     >                                   BREAD(1),BREAD(2),BREAD(3)
                 ELSE
                   READ(LUN,*,END=97,ERR=98) DUM1,DUM2,DUM3,
     >                                   BREAD(1),BREAD(2),BREAD(3)
C                   write(33,*) DUM1,DUM2,DUM3,BREAD(1),BREAD(2),BREAD(3)
                 ENDIF

                 BMAX0 = BMAX
                 BMAX = DMAX1(BMAX,BREAD(1),BREAD(2),BREAD(3))
                 IF(BMAX.NE.BMAX0) THEN
                   XBMA = TTA
                   YBMA = RADIUS
                   ZBMA = Z
                 ENDIF
                 BMIN0 = BMIN
                 BMIN = DMIN1(BMIN,BREAD(1),BREAD(2),BREAD(3))
                 IF(BMIN.NE.BMIN0) THEN
                   XBMI = TTA
                   YBMI = RADIUS
                   ZBMI = Z
                 ENDIF

C---------------- Polar components Btheta, Br, BZ at positon (tta,r,Z) in cyl. frame
                 IF(ZROBXY) THEN
                   IF(Z .EQ. 0.D0) THEN
                     HC(1,I,J,K,IMAP) = 0.D0
                     HC(2,I,J,K,IMAP) = 0.D0
                   ELSE
                     HC(1,I,J,K,IMAP) = BREAD(1) * BNORM
                     HC(2,I,J,K,IMAP) = BREAD(2) * BNORM
                   ENDIF
                 ELSE
                   HC(1,I,J,K,IMAP) = BREAD(1) * BNORM
                   HC(2,I,J,K,IMAP) = BREAD(2) * BNORM
                 ENDIF
                 HC(3,I,J,K,IMAP) = BREAD(3) * BNORM

               ENDDO
             ENDDO
           ENDDO

        BMIN = BMIN * BNORM
        BMAX = BMAX * BNORM

C------- Mesh coordinates
        DO I=-IXMA/2,IXMA/2
          XH(I+IXMA/2+1) =  DBLE(I) * DTTA
        ENDDO
        DO J=1,JYMA
          YH(J) =  R0 + DBLE(J-1) *DR 
        ENDDO
        DO K= -KZMA/2,KZMA/2   
          ZH(K+KZMA/2+1) = DBLE(K) * DZ
        ENDDO


      ELSEIF(MOD .EQ. 22) THEN

           IRMA = JYMA
           JTMA = IXMA
           K2ZMA = KZMA/2+1

           IRCNT = 0
           DO 133 I=1,IRMA            
             RADIUS=R0 + DBLE(I-1) *DR 
             IRCNT = IRCNT+1
             KZCNT=0       
             DO 133  K = 1,K2ZMA      
               KZC = K2ZMA-1+K
               Z = DBLE(K-1) * DZ
               KZCNT = KZCNT+1
              JTCNT=0
              DO 133 J=1,JTMA        
                   JTC = J
                 TTA =  DBLE(J-1) * DTTA
                 JTCNT = JTCNT + 1
                 IF(BINAR) THEN
                   READ(LUN,END=97,ERR=98) DUM1,ZZZ,DUM3,
     >                                   BREAD(1),BREAD(2),BREAD(3)
                 ELSE

 224                CONTINUE
                   READ(LUN,FMT='(A)') TXT132
                   IDSTR = DEBSTR(TXT132)
                   IF    (EMPTY(TXT132)) THEN
                     GOTO 224
                   ELSEIF(TXT132(IDSTR:IDSTR+1) .EQ. '%'
     >                        .OR. TXT132(IDSTR:IDSTR+1) .EQ. '#') THEN
                     GOTO 224 
                   ELSE
                     IF(FINSTR(TXT132) .EQ. MXCHAR) 
     >                  CALL ENDJOB('SBR FMAPW : # of columns in field' 
     >                  //' data file must be <',MXCHAR)
                   ENDIF

                   READ(TXT132,*,END=97,ERR=98) DUM1,ZZZ,DUM3,
     >                                   BREAD(1),BREAD(2),BREAD(3)

                 ENDIF
                 BMAX0 = BMAX
                 BMAX = DMAX1(BMAX,BREAD(1),BREAD(2),BREAD(3))
                 IF(BMAX.NE.BMAX0) THEN
                   XBMA = TTA
                   YBMA = RADIUS
                   ZBMA = Z
                 ENDIF
                 BMIN0 = BMIN
                 BMIN = DMIN1(BMIN,BREAD(1),BREAD(2),BREAD(3))
                 IF(BMIN.NE.BMIN0) THEN
                   XBMI = TTA
                   YBMI = RADIUS
                   ZBMI = Z
                 ENDIF

C---------------- Horizontal component of field                 
                 BH = SQRT(BREAD(1)*BREAD(1) + BREAD(3)*BREAD(3))
                 ALP = ATAN2(BREAD(3),BREAD(1)) - TTA
C---------------- Polar components at positon (r,tta) in cyl. frame
                 BRAD = BH * COS(ALP)
                 BTTA = BH * SIN(ALP)
C---------------- Watch the sign !! Bx and By multiplied by (-1)
                 HC(1,JTC,I,KZC,IMAP) =   -  BTTA * BNORM
                 HC(2,JTC,I,KZC,IMAP) =   -  BRAD * BNORM
                 HC(3,JTC,I,KZC,IMAP) = BREAD(2) * BNORM

C---------------- In case that non-zero Bx, By in median plane would be prohibitive
                 IF(ZZZ .EQ. 0.D0) THEN
                    HC(1,JTC,I,KZC,IMAP) = 0.D0
                    HC(2,JTC,I,KZC,IMAP) = 0.D0
                ENDIF   

C----------------  TEST RACCAM
C     >          HC(1,JTC,I,KZC,IMAP) = HC(1,JTC,I,KZC,IMAP)  *(RADIUS/348.)**(-.2)
C     >          HC(2,JTC,I,KZC,IMAP) = HC(2,JTC,I,KZC,IMAP)  *(RADIUS/348.)**(-.2)
C     >          HC(3,JTC,I,KZC,IMAP) = HC(3,JTC,I,KZC,IMAP)  *(RADIUS/348.)**(-.2)
                     
 133       CONTINUE

        BMIN = BMIN * BNORM
        BMAX = BMAX * BNORM
                   XBMA = XBMA * XNORM
                   YBMA = YBMA * YNORM
                   ZBMA = ZBMA * ZNORM
                   XBMI = XBMI * XNORM
                   YBMI = YBMI * YNORM
                   ZBMI = ZBMI * ZNORM

C------- symmetrise 3D map wrt mid-plane= bend-plane
        DO 135  K=2,K2ZMA      
          KZC = K2ZMA-1+K
          KZS = K2ZMA-K+1
          DO 135 I=1,IRMA         
            DO 135 J=1,JTMA
              HC(1,J,I,KZS,IMAP) = -HC(1,J,I,KZC,IMAP)
              HC(2,J,I,KZS,IMAP) = -HC(2,J,I,KZC,IMAP) 
              HC(3,J,I,KZS,IMAP) = HC(3,J,I,KZC,IMAP) 
 135    CONTINUE

C------- Mesh coordinates
          DO 137 J=1,JYMA
 137         YH(J) =  R0 + DBLE(J-1) *DR 
          DO K= -KZMA/2,KZMA/2   
            ZH(K+KZMA/2+1) = DBLE(K) * DZ
          ENDDO
          DO 136 I=-IXMA/2,IXMA/2
 136        XH(I+IXMA/2+1) =  DBLE(I) * DTTA

      ELSEIF(MOD .EQ. 23) THEN

           IRMA = JYMA
           JTMA = IXMA
           K2ZMA = KZMA/2+1

              JTCNT=0
              DO J=1,JTMA        
                   JTC = J
                 TTA =  DBLE(J-1) * DTTA
                 JTCNT = JTCNT + 1

             KZCNT=0       
             DO  K = 1,K2ZMA      
               KZC = K2ZMA-1+K
               Z = DBLE(K-1) * DZ
               KZCNT = KZCNT+1

           IRCNT = 0
           DO I=1,IRMA            
             RADIUS=R0 + DBLE(I-1) *DR 
             IRCNT = IRCNT+1

                 IF(BINAR) THEN
                   READ(LUN,END=97,ERR=98) DUM1,ZZZ,DUM3,
     >                                   BREAD(1),BREAD(2),BREAD(3)
                 ELSE

 223                CONTINUE
                   READ(LUN,FMT='(A)') TXT132
                   IDSTR = DEBSTR(TXT132)
                   IF    (EMPTY(TXT132)) THEN
                     GOTO 223
                   ELSEIF(TXT132(IDSTR:IDSTR+1) .EQ. '%'
     >                        .OR. TXT132(IDSTR:IDSTR+1) .EQ. '#') THEN
                     GOTO 223 
                   ELSE
                     IF(FINSTR(TXT132) .EQ. MXCHAR) 
     >                  CALL ENDJOB('SBR FMAPW : # of columns in field' 
     >                  //' data file must be <',MXCHAR)
                   ENDIF

                   READ(TXT132,*,END=97,ERR=98) DUM1,ZZZ,DUM3,
     >                                   BREAD(1),BREAD(2),BREAD(3)

                 ENDIF
                 BMAX0 = BMAX
                 BMAX = DMAX1(BMAX,BREAD(1),BREAD(2),BREAD(3))
                 IF(BMAX.NE.BMAX0) THEN
                   XBMA = TTA
                   YBMA = RADIUS
                   ZBMA = Z
                 ENDIF
                 BMIN0 = BMIN
                 BMIN = DMIN1(BMIN,BREAD(1),BREAD(2),BREAD(3))
                 IF(BMIN.NE.BMIN0) THEN
                   XBMI = TTA
                   YBMI = RADIUS
                   ZBMI = Z
                 ENDIF

C---------------- Horizontal component of field                 
                 BH = SQRT(BREAD(1)*BREAD(1) + BREAD(3)*BREAD(3))
                 ALP = ATAN2(BREAD(3),BREAD(1)) - TTA
C---------------- Polar components at positon (r,tta) in cyl. frame
                 BRAD = BH * COS(ALP)
                 BTTA = BH * SIN(ALP)
C---------------- Watch the sign !! Bx and By multiplied by (-1)
                 HC(1,JTC,I,KZC,IMAP) =   -  BTTA * BNORM       !!!!!!!!!!!! sign TBC
C                       write(*,*) '   !!!!!!!!!!!! sign TBC '
                 HC(2,JTC,I,KZC,IMAP) =   -  BRAD * BNORM       !!!!!!!!!!!! sign TBC
                 HC(3,JTC,I,KZC,IMAP) = BREAD(2) * BNORM

C---------------- In case that non-zero Bx, By in median plane would be prohibitive
                 IF(ZZZ .EQ. 0.D0) THEN
                    HC(1,JTC,I,KZC,IMAP) = 0.D0
                    HC(2,JTC,I,KZC,IMAP) = 0.D0
                ENDIF   

              ENDDO
            ENDDO
          ENDDO

        BMIN = BMIN * BNORM
        BMAX = BMAX * BNORM
                   XBMA = XBMA * XNORM
                   YBMA = YBMA * YNORM
                   ZBMA = ZBMA * ZNORM
                   XBMI = XBMI * XNORM
                   YBMI = YBMI * YNORM
                   ZBMI = ZBMI * ZNORM

C------- symmetrise 3D map wrt mid-plane= bend-plane
        DO K=2,K2ZMA      
          KZC = K2ZMA-1+K
          KZS = K2ZMA-K+1
          DO I=1,IRMA         
            DO J=1,JTMA
              HC(1,J,I,KZS,IMAP) = -HC(1,J,I,KZC,IMAP)
              HC(2,J,I,KZS,IMAP) = -HC(2,J,I,KZC,IMAP) 
              HC(3,J,I,KZS,IMAP) = HC(3,J,I,KZC,IMAP) 
            ENDDO
          ENDDO
        ENDDO

C------- Mesh coordinates
        DO J=1,JYMA
           YH(J) =  R0 + DBLE(J-1) *DR 
        ENDDO
        DO K= -KZMA/2,KZMA/2   
          ZH(K+KZMA/2+1) = DBLE(K) * DZ
          ENDDO
        DO I=-IXMA/2,IXMA/2
          XH(I+IXMA/2+1) =  DBLE(I) * DTTA
        ENDDO

      ELSEIF(MOD .EQ. 25) THEN

           IRMA = JYMA
           JTMA = IXMA
           KZC = 1

           IRCNT = 0
           DO I=1,IRMA            
             RADIUS=R0 + DBLE(I-1) *DR 
             IRCNT = IRCNT+1
              JTCNT=0
              DO J=1,JTMA        
                 JTC = J
                 TTA =  DBLE(J-1) * DTTA
                 JTCNT = JTCNT + 1
                 IF(BINAR) THEN
                   READ(LUN,END=97,ERR=98) DUM1,DUM2,zzz,
     >                                   BREAD(1),BREAD(2),BREAD(3)
                 ELSE

 227                CONTINUE
                   READ(LUN,FMT='(A)') TXT132
                   IDSTR = DEBSTR(TXT132)
                   IF    (EMPTY(TXT132)) THEN
                     GOTO 227
                   ELSEIF(TXT132(IDSTR:IDSTR+1) .EQ. '%'
     >                        .OR. TXT132(IDSTR:IDSTR+1) .EQ. '#') THEN
                     GOTO 227 
                   ELSE
                     IF(FINSTR(TXT132) .EQ. MXCHAR) 
     >                  CALL ENDJOB('SBR FMAPW : # of columns in field' 
     >                  //' data file must be <',MXCHAR)
                   ENDIF

                   READ(TXT132,*,END=97,ERR=98) dum1,dum2,zzz,
     >                                   BREAD(1),BREAD(2),BREAD(3)

                 ENDIF
                 BMAX0 = BMAX
                 BMAX = DMAX1(BMAX,BREAD(1),BREAD(2),BREAD(3))
                 IF(BMAX.NE.BMAX0) THEN
                   XBMA = TTA
                   YBMA = RADIUS
                   ZBMA = Z
                 ENDIF
                 BMIN0 = BMIN
                 BMIN = DMIN1(BMIN,BREAD(1),BREAD(2),BREAD(3))
                 IF(BMIN.NE.BMIN0) THEN
                   XBMI = TTA
                   YBMI = RADIUS
                   ZBMI = Z
                 ENDIF

                 HC(3,JTC,I,KZC,IMAP) = BREAD(3) * BNORM

               ENDDO
           ENDDO

        BMIN = BMIN * BNORM
        BMAX = BMAX * BNORM
                   XBMA = XBMA * XNORM
                   YBMA = YBMA * YNORM
                   ZBMA = ZBMA * ZNORM
                   XBMI = XBMI * XNORM
                   YBMI = YBMI * YNORM
                   ZBMI = ZBMI * ZNORM

C------- Mesh coordinates
          DO J=1,JYMA
            YH(J) =  R0 + DBLE(J-1) *DR 
          ENDDO
          DO I=-IXMA/2,IXMA/2
            XH(I+IXMA/2+1) =  DBLE(I) * DTTA
          ENDDO

      ELSE
        CALL ENDJOB('SBR FMAPW/FMAPR2, no such option MOD = ',MOD)
      ENDIF 

      RETURN

C-------------------------------------------------------------
C Read and interprete field maps in cartesian frame (MOD < 20)
      ENTRY FMAPR3(BINAR,LUN,MODI,MODI2,NHDI,ZROBXY,
     >             XNORM,YNORM,ZNORM,BNORM,I1,KZI,FMTYP,
     >                       BMIN,BMAX,
     >                       XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)
      
C      IF( .NOT.ALLOCATED( HCTMP )) 
C     >     ALLOCATE( HCTMP(ID,MXX,MXY,IZ,MX3D), STAT = IALOC)
C      IF (IALOC .NE. 0) 
C     >     CALL ENDJOB('SBR FMAPW Not enough memory'//
C     >     ' for Malloc of HCTMP',
C     >     -99)

Check current number of field map, IMAP in HC(*,*,*,*,IMAP)
      CALL KSMAP(
     >           IMAP) 
      MOD = MODI
      MOD2 = MODI2
      NHD = NHDI
      KZ = KZI
      BMIN =  1.D10
      BMAX = -1.D10
      ZBMA = 0.D0
      ZBMI = 0.D0
C Map data file starts with NHD-line header
      IF(NRES.GT.0) 
     >  WRITE(NRES,FMT='(A,I1,A)') ' HEADER  (',NHD,' lines) : '

      IF(BINAR) THEN

        IF(NHD .GE. 1) THEN
C Map data file starts with NHD-line header
          READ(LUN,END=97,ERR=98) HDRM
          IF(NRES.GT.0) WRITE(NRES,FMT='(5X,A120)') HDRM
          IF(MOD .EQ. 15) THEN
          ELSE
            READ(HDRM,*,ERR=51,END=51) R0, DR, DTTA, DZ
          ENDIF
          GOTO 52
 51       CONTINUE
          IF(NRES.GT.0) WRITE(NRES,*) 
     >    'Could not read R_min, DR, DTTA, DZ from line 1 of HEADER...'
 52       CONTINUE
          IF(NHD .GE. 2) THEN
            READ(LUN) HDRM
            IF(NRES.GT.0) WRITE(NRES,FMT='(5X,A120)') HDRM
            IF(NHD .GE. 3) THEN
              DO II=1, NHD-2
                READ(LUN) HDRM
                IF(NRES.GT.0) WRITE(NRES,FMT='(5X,A120)') HDRM
              ENDDO
              IF(NRES.GT.0) THEN
                WRITE(NRES,FMT='(A,1P,4E17.8,/)') 
     >          'R_min (cm), DR (cm), DTTA (deg), DZ (cm) : ',
     >           R0,DR,DTTA,DZ
C     >          ' R_min, DR, DTTA, DZ : ',R0, DR, DTTA, DZ
                CALL FLUSH2(NRES,.FALSE.)
              ENDIF
            ENDIF
          ENDIF
        ELSE
C          CALL ENDJOB('Please give at least 1 line header in field map'
C     >    //' file with R_min, DR, DTTA, DZ values.',-99)
        ENDIF
      ELSE

        IF(NHD .GE. 1) THEN

          READ(LUN,FMT='(A120)') HDRM
          IF(NRES.GT.0) WRITE(NRES,FMT='(5X,A120)') HDRM
          READ(HDRM,*,END=49,ERR=49) R0, DR, DX, DZ
          GOTO 50
 49       CONTINUE
          IF(NRES.GT.0) WRITE(NRES,*) 
     >    'Could not read R_min, DR, DX, DZ from line 1 of HEADER...'
 50       CONTINUE
          IF(NHD .GE. 2) THEN
            READ(LUN,FMT='(A120)') HDRM
            IF(NRES.GT.0) WRITE(NRES,FMT='(5X,A120)') HDRM
            IF(NHD .GE. 3) THEN
              DO II=1, NHD-2
                READ(LUN,FMT='(A120)') HDRM
                IF(NRES.GT.0) WRITE(NRES,FMT='(5X,A120)') HDRM
              ENDDO
              IF(NRES.GT.0) THEN
                WRITE(NRES,FMT='(A,1P,4E17.8,/)') 
     >          'R_min (cm), DR (cm), DTTA (deg), DZ (cm) : ',
     >          R0,DR,DX,DZ
C                WRITE(NRES,*) ' R_min, DR, DX, DZ : ',R0,DR,DX,DZ
                CALL FLUSH2(NRES,.FALSE.)
              ENDIF
            ENDIF
          ENDIF
        ELSE
C          CALL ENDJOB('Please give at least 1 line header in field map'
C     >    //' file with  R_min, DR, DX, DZ values.',-99)
        ENDIF

      ENDIF

      IF(MOD .EQ. 12) THEN
C Read 3D field map of (part or total of) magnet, contained in a single file

        IF(MOD2 .EQ. 0) THEN    
C------- Default : upper half of magnet, symmetrise 3D map wrt (X,Y)=mid-plane=bend-plane
           JTCNT=0
           IRCNT = 0
           KZCNT=0       

           DO 121 I=1,IXMA            
             IRCNT = IRCNT+1
             DO 121  K = 1,KZMA      
               KZC = KZMA-1+K
               KZCNT = KZCNT+1
               DO 121 J=1,JYMA        
                 JTC = J
                 JTCNT = JTCNT + 1

                 IF(BINAR) THEN
                   READ(LUN)YH(J),ZH(K),XH(I),BREAD(2),BREAD(3),BREAD(1)
                 ELSE
                   READ(LUN,FMT='(6E20.2)') YH(J),ZH(K),XH(I), 
     >                                    BREAD(2),BREAD(3),BREAD(1)
                 ENDIF
                 BMAX0 = BMAX
                 BMAX = DMAX1(BMAX,BREAD(1),BREAD(2),BREAD(3))
                 IF(BMAX.NE.BMAX0) THEN
                   XBMA = XH(I)
                   YBMA = YH(J)
                   ZBMA = ZH(K)
                 ENDIF
                 BMIN0 = BMIN
                 BMIN = DMIN1(BMIN,BREAD(1),BREAD(2),BREAD(3))
                 IF(BMIN.NE.BMIN0) THEN
                   XBMI = XH(I)
                   YBMI = YH(J)
                   ZBMI = ZH(K)
                 ENDIF

                 HC(1,I,JTC,KZC,IMAP) = BREAD(1) * BNORM
                 HC(2,I,JTC,KZC,IMAP) = BREAD(2) * BNORM
                 HC(3,I,JTC,KZC,IMAP) = BREAD(3) * BNORM

 121          CONTINUE

          BMIN = BMIN * BNORM
          BMAX = BMAX * BNORM

C------- symmetrise 3D map wrt (X,Y) plane
          DO 122  K=2,KZMA      
            KZC = KZMA-1+K
            KZS = KZMA-K+1
            DO 122 I=1,IXMA         
              DO 122 J=1,JYMA
                HC(1,J,I,KZS,IMAP) = -HC(1,J,I,KZC,IMAP)
                HC(2,J,I,KZS,IMAP) = -HC(2,J,I,KZC,IMAP) 
                HC(3,J,I,KZS,IMAP) = HC(3,J,I,KZC,IMAP) 
 122       CONTINUE

C------- Mesh coordinates
          DY = YH(KZMA*IXMA) - YH(1)
          DO 127 J=2,JYMA
 127          YH(J) =  YH(1) + DY
          DZ = ZH(IXMA) - ZH(1)
          DO 128 K= 2, 2*KZMA-1
 128         ZH(K) = ZH(1) + DZ

        ELSEIF(MOD2 .EQ. 1) THEN
C Full map of magnet. Differs from MOD2=0 by absence of symmetrization
C Used for instance for AGS helical snake maps
           JTCNT=0
           IRCNT = 0
           KZCNT=0       
           DO J=1,JYMA        
             JTC = J
             JTCNT = JTCNT + 1
             DO  K = 1,KZMA      
               KZC = K
               KZCNT = KZCNT+1
               DO I=1,IXMA            
                 IRCNT = IRCNT+1

                 IF(BINAR) THEN
                   READ(LUN)YH(J),ZH(K),XH(I),BREAD(2),BREAD(3),BREAD(1)
                 ELSE
                   READ(LUN,*) YH(J),ZH(K),XH(I), 
     >                                    BREAD(2),BREAD(3),BREAD(1)
                 ENDIF
                 BMAX0 = BMAX
                 BMAX = DMAX1(BMAX,BREAD(1),BREAD(2),BREAD(3))
                 IF(BMAX.NE.BMAX0) THEN
                   XBMA = XH(I)
                   YBMA = YH(J)
                   ZBMA = ZH(K)
                 ENDIF
                 BMIN0 = BMIN
                 BMIN = DMIN1(BMIN,BREAD(1),BREAD(2),BREAD(3))
                 IF(BMIN.NE.BMIN0) THEN
                   XBMI = XH(I)
                   YBMI = YH(J)
                   ZBMI = ZH(K)
                 ENDIF

                 HC(1,I,JTC,KZC,IMAP) = BREAD(1) * BNORM
                 HC(2,I,JTC,KZC,IMAP) = BREAD(2) * BNORM
                 HC(3,I,JTC,KZC,IMAP) = BREAD(3) * BNORM

               ENDDO
             ENDDO
           ENDDO

           BMIN = BMIN * BNORM
           BMAX = BMAX * BNORM

C------- Mesh coordinates
           DX = (XH(2) - XH(1))*Xnorm
           XH(1) = XH(1)*Xnorm
           DO J=2,IXMA
             XH(J) =  XH(J-1) + DX
           ENDDO
           DY = (YH(2) - YH(1))*ynorm
           YH(1) = YH(1)*ynorm
           DO J=2,JYMA
             YH(J) =  YH(J-1) + DY
           ENDDO
           IF(IZ .EQ. 1) THEN 
             IZ1 = 1
           ELSE
             IZ1 = 2
           ENDIF
           DZ = (ZH(IZ1) - ZH(1))*Znorm
           ZH(1) = ZH(1)*Znorm
           DO K= 2, KZMA
             ZH(K) = ZH(K-1) + DZ
           ENDDO

        ELSEIF(MOD2 .EQ. 2) THEN
C Upper right bottom 1/8th of magnet. Symmetrization wrt (i) (X,Y) plane, (ii) (Y,Z) plane, 
C (iii) (X,Z) plane. Used for D5 IPM magnet in AGS.
           jtcnt=0
           ircnt = 0
           kzcnt=0       
           DO J = JYMA/2+1, JYMA
             jtcnt = jtcnt + 1
             DO  K = KZMA/2+1, KZMA
               kzcnt = kzcnt+1
               DO I = IXMA/2+1, IXMA
                 ircnt = ircnt+1

                 IF(BINAR) THEN
                   READ(LUN)YH(J),ZH(K),XH(I),BREAD(2),BREAD(3),BREAD(1)
                 ELSE
                   READ(LUN,*) YH(J),ZH(K),XH(I), 
     >                                    BREAD(2),BREAD(3),BREAD(1)
                 ENDIF

                 BMAX0 = BMAX
                 BMAX = DMAX1(BMAX,BREAD(1),BREAD(2),BREAD(3))
                 IF(BMAX.NE.BMAX0) THEN
                   XBMA = XH(I)
                   YBMA = YH(J)
                   ZBMA = ZH(K)
                 ENDIF
                 BMIN0 = BMIN
                 BMIN = DMIN1(BMIN,BREAD(1),BREAD(2),BREAD(3))
                 IF(BMIN.NE.BMIN0) THEN
                   XBMI = XH(I)
                   YBMI = YH(J)
                   ZBMI = ZH(K)
                 ENDIF

                 HC(1,I,J,K,IMAP) = BREAD(1) * BNORM
                 HC(2,I,J,K,IMAP) = BREAD(2) * BNORM
                 HC(3,I,J,K,IMAP) = BREAD(3) * BNORM

               ENDDO
             ENDDO
           ENDDO

          BMIN = BMIN * BNORM
          BMAX = BMAX * BNORM

C--------- symmetrise 3D map wrt (Y,Z) plane= longitudinal plane
          DO I = 1, IXMA/2   
            IXC = IXMA/2+1+I
            IXS = IXMA/2+1-I
            DO  K = KZMA/2+1, KZMA
              DO J = JYMA/2+1, JYMA
                HC(1,IXS,J,K,IMAP) = - HC(1,IXC,J,K,IMAP)
                HC(2,IXS,J,K,IMAP) =   HC(2,IXC,J,K,IMAP) 
                HC(3,IXS,J,K,IMAP) =   HC(3,IXC,J,K,IMAP) 
              ENDDO
            ENDDO
          ENDDO
C--------- symmetrise 3D map wrt (X,Z) plane= longitudinal plane
          DO J = 1, JYMA/2
            JYC = JYMA/2+1+J
            JYS = JYMA/2+1-J
            DO  K = KZMA/2+1, KZMA
              DO I = 1, IXMA
                HC(1,I,JYS,K,IMAP) =   HC(1,I,JYC,K,IMAP)
                HC(2,I,JYS,K,IMAP) = - HC(2,I,JYC,K,IMAP) 
                HC(3,I,JYS,K,IMAP) =   HC(3,I,JYC,K,IMAP) 
              ENDDO
            ENDDO
          ENDDO
C--------- symmetrise 3D map wrt (X,Y) plane= longitudinal plane
          DO K = 1, KZMA/2
            KZC = KZMA/2+1+K
            KZS = KZMA/2+1-K
            DO  J = 1, JYMA     
              DO I = 1, IXMA 
                HC(1,I,J,KZS,IMAP) = - HC(1,I,J,KZC,IMAP)
                HC(2,I,J,KZS,IMAP) = - HC(2,I,J,KZC,IMAP) 
                HC(3,I,J,KZS,IMAP) =   HC(3,I,J,KZC,IMAP) 
              ENDDO
            ENDDO
          ENDDO

C--------- Mesh coordinates
           DX = (XH(IXMA) - XH(IXMA-1))*XNORM
           XH(1) = -XH(IXMA)*XNORM
           DO I=2, IXMA
             XH(I) = XH(I-1) + DX
           ENDDO
           DY = (YH(JYMA) - YH(JYMA-1))*YNORM
           YH(1) = -YH(JYMA)*YNORM
           DO J=2, JYMA
             YH(J) = YH(J-1) + DY
           ENDDO
           DZ = (ZH(KZMA) - ZH(KZMA-1))*ZNORM
           ZH(1) = -ZH(KZMA)*ZNORM
           DO K= 2, KZMA
             ZH(K) = ZH(K-1) + DZ
           ENDDO

        ENDIF ! MOD2=0, 1, 2

      ELSEIF(MOD .EQ. 0 .OR. MOD .EQ. 1) THEN
C------- EMMAC, GSI spectro for instance. 
C------- MAP2D (2D map)
C        Keyword EMMA with two 2D maps 
C        MOD=0 : # of files is NF=1+|IZ/2|, from z=0 to z_max, symmetrizing wrt median plane
C        MOD=1 : # of files is NF= IZ, from -z_max to +z_max, no symmetrizing

           I = KZ

             DO 10 J=1,JYMA
               DO  10  K = 1,IXMA

                 IF( BINAR ) THEN
                   READ(LUN)YH(J),ZH(I),XH(K),BREAD(2),BREAD(3),BREAD(1)

                 ELSE

 14                CONTINUE
                   READ(LUN,FMT='(A)') TXT132

                   IDSTR = DEBSTR(TXT132)
                   IF    (EMPTY(TXT132)) THEN
                     GOTO 14
                   ELSEIF(TXT132(IDSTR:IDSTR+1) .EQ. '%'
     >                        .OR. TXT132(IDSTR:IDSTR+1) .EQ. '#') THEN
                     GOTO 14
                   ELSE
                     IF(FINSTR(TXT132) .EQ. MXCHAR) 
     >                  CALL ENDJOB('SBR FMAPW : # of columns in field' 
     >                  //' data file must be <',MXCHAR)
                   ENDIF

                   IF    (MOD2 .EQ. 1) THEN 
                     READ(TXT132,FMT='(1X,6E11.2)',ERR=96)  
     >                                        YH(J),ZH(I),XH(K), 
     >                                        BREAD(2),BREAD(3),BREAD(1)
                   ELSEIF(MOD2 .EQ. 2) THEN 
                     READ(TXT132,FMT='(1X,6E12.2)',ERR=96) 
     >                                        YH(J),ZH(I),XH(K), 
     >                                        BREAD(2),BREAD(3),BREAD(1)
                   ELSEIF(MOD2 .EQ. 3) THEN 
                     READ(TXT132,*,ERR=96) YH(J),ZH(I),XH(K), 
     >                                        BREAD(2),BREAD(3),BREAD(1)
                   ELSE

C-------------------- Default MOD2
                     IF(FMTYP .EQ. 'GSI') THEN

                       READ(TXT132,FMT='(1X,6E11.2)',ERR=96) 
     >                                      YH(J),ZH(I),XH(K),  
     >                                      BREAD(2),BREAD(3),BREAD(1)

                     ELSEIF(FMTYP .EQ. 'LESB3') THEN
                       READ(TXT132,FMT='(1X,6E11.2)',ERR=96) 
     >                                      YH(J),ZH(I),XH(K),  
     >                                      BREAD(2),BREAD(3),BREAD(1)

                     ELSE
                       READ(TXT132,*,ERR=96) YH(J),ZH(I),XH(K),  
CCCCCCCCCCCCC                       READ(TXT132,FMT='(1X,6E11.2)') YH(J),ZH(I),XH(K),  
     >                                      BREAD(2),BREAD(3),BREAD(1)
                     ENDIF
                   ENDIF

CC----  Manip VAMOS GANIL, oct. 2001 
C                     YH(J) = YH(J) * 1.D2
C                     zH(i) = zH(i) * 1.D2
C                     xH(k) = xH(k) * 1.D2

                 ENDIF
                 BMAX0 = BMAX
                 BMAX = DMAX1(BMAX,BREAD(1),BREAD(2),BREAD(3))
                 IF(BMAX.NE.BMAX0) THEN
                   XBMA = XH(K)
                   YBMA = YH(J)
                   ZBMA = ZH(I)
                 ENDIF
                 BMIN0 = BMIN
                 BMIN = DMIN1(BMIN,BREAD(1),BREAD(2),BREAD(3))
                 IF(BMIN.NE.BMIN0) THEN
                   XBMI = XH(K)
                   YBMI = YH(J)
                   ZBMI = ZH(I)
                 ENDIF
C FM 11/03 
                 IF(ID .EQ. 1) THEN
                   HC(ID,K,J,I,IMAP) =  BREAD(3) * BNORM
                 ELSE
                   DO 187 LHC=1,ID
 187                  HC(LHC,K,J,I,IMAP) = BREAD(LHC) * BNORM
                 ENDIF

                 YH(J) = YH(J) * YNORM
                 ZH(I) = ZH(I) * ZNORM
                 XH(K) = XH(K) * XNORM
                    
                 IF(I1 .GT. 1) THEN
C------------------- Symmetrize 3D map wrt XY plane (IZ>1)
                   IF( I .GT. I1 ) THEN
                     ZH(2*I1-I) = -ZH(I)
                     FAC = -1.D0
                     HC(1,K,J,2*I1-I,IMAP) =   -BREAD(1)* BNORM
                     II=2
                     HC(II,K,J,2*I1-I,IMAP) =   -BREAD(2)* BNORM
                     II=3
                     HC(II,K,J,2*I1-I,IMAP) =   BREAD(3)* BNORM
                   ENDIF
                 ENDIF
 10          CONTINUE
             CLOSE(UNIT=LUN)
C 12        CONTINUE
 
             BMIN = BMIN * BNORM
             BMAX = BMAX * BNORM
                   XBMA = XBMA * XNORM
                   YBMA = YBMA * YNORM
                   ZBMA = ZBMA * ZNORM
                   XBMI = XBMI * XNORM
                   YBMI = YBMI * YNORM
                   ZBMI = ZBMI * ZNORM

      ELSEIF(MOD .EQ. 3) THEN
C------- AGS magnet maps (2D map, half-magnet, symmetrized wrt YZ plane - at 45inches in A and C type)

           K = 1
           
           READ(LUN,*,ERR=96) JYMA
           IF(JYMA.GT.MXY) CALL 
     >       ENDJOB(' SBR fmapw. In PARIZ.H need MXY .ge. ',JYMA)
           READ(LUN,*,ERR=96) (YH(J),J=1,JYMA)
           DO J=1,JYMA
             YH(J) = YH(J) * YNORM
           ENDDO

           READ(LUN,*,ERR=96) IXMA2
           IXMA = 2*IXMA2 - 1
           IF(IXMA.GT.MXX) CALL 
     >       ENDJOB(' SBR fmapw. In PARIZ.H need MXX .ge. ',IXMA)
           READ(LUN,*,ERR=96) (XH(I),I=1,IXMA2)
           DO I=1,IXMA2
             XH(I) = XH(I) * XNORM
           ENDDO
           DX = XH(2) - XH(1)
           DO I=IXMA2+1,IXMA
             XH(I) = XH(IXMA2) + DBLE(I-IXMA2)*DX
           ENDDO

           READ(LUN,*,ERR=96) NTOT  

           DO I = 1, IXMA2
             READ(LUN,FMT='(A132)',ERR=96) TXT132
             READ(LUN,*,ERR=96) (HC(ID,I,J,K,IMAP), J=1,JYMA)
             DO J = 1, JYMA
               BREAD(3) = HC(ID,I,J,K,IMAP)
               BMAX0 = BMAX
               BMAX = DMAX1(BMAX,BREAD(3))
               IF(BMAX.NE.BMAX0) THEN
                 XBMA = XH(I)
                 YBMA = YH(J)
               ENDIF
               BMIN0 = BMIN
               BMIN = DMIN1(BMIN,BREAD(3))
               IF(BMIN.NE.BMIN0) THEN
                 XBMI = XH(I)
                 YBMI = YH(J)
               ENDIF
               TEMP =  BREAD(3) * BNORM
               HC(ID,       I,J,K,IMAP) = TEMP
C-------------- Set the other X-half of magnet
               HC(ID,IXMA+1-I,J,K,IMAP) = TEMP
             ENDDO
           ENDDO

             CLOSE(UNIT=LUN)
C 12        CONTINUE
 
             BMIN = BMIN * BNORM
             BMAX = BMAX * BNORM
                   XBMA = XBMA * XNORM
                   YBMA = YBMA * YNORM
                   ZBMA = ZBMA * ZNORM
                   XBMI = XBMI * XNORM
                   YBMI = YBMI * YNORM
                   ZBMI = ZBMI * ZNORM

      ELSEIF(MOD .EQ. 15) THEN
C Full 3-D map of magnet as in MOD=12. However can sum up several maps. 
C Used for instance for AGS cold snake = helix map + solenoid map
           JTCNT=0
           IRCNT = 0
           KZCNT=0       

           DO J=1,JYMA        
             JTC = J
             JTCNT = JTCNT + 1
             DO  K = 1,KZMA      
               KZC = K
               KZCNT = KZCNT+1
               DO I=1,IXMA            
                 IRCNT = IRCNT+1

                 IF(BINAR) THEN
                   READ(LUN,ERR=96)
     >             YH(J),ZH(K),XH(I),BREAD(2),BREAD(3),BREAD(1)
                 ELSE
                   READ(LUN,*,ERR=96) YH(J),ZH(K),XH(I), 
     >                                    BREAD(2),BREAD(3),BREAD(1)
                 ENDIF
                 BMAX0 = BMAX
                 BMAX = DMAX1(BMAX,BREAD(1),BREAD(2),BREAD(3))
                 IF(BMAX.NE.BMAX0) THEN
                   XBMA = XH(I)
                   YBMA = YH(J)
                   ZBMA = ZH(K)
                 ENDIF
                 BMIN0 = BMIN
                 BMIN = DMIN1(BMIN,BREAD(1),BREAD(2),BREAD(3))
                 IF(BMIN.NE.BMIN0) THEN
                   XBMI = XH(I)
                   YBMI = YH(J)
                   ZBMI = ZH(K)
                 ENDIF


                 IF(ZROBXY) THEN
                   IF(ZH(K) .EQ. 0.D0) THEN
                     HC(1,I,JTC,KZC,IMAP) = 0.D0
                     HC(2,I,JTC,KZC,IMAP) = 0.D0
                   ELSE
                     HC(1,I,JTC,KZC,IMAP) = BREAD(1) * BNORM
                     HC(2,I,JTC,KZC,IMAP) = BREAD(2) * BNORM
                   ENDIF
                 ELSE
                   HC(1,I,JTC,KZC,IMAP) = BREAD(1) * BNORM
                   HC(2,I,JTC,KZC,IMAP) = BREAD(2) * BNORM
                 ENDIF
                 HC(3,I,JTC,KZC,IMAP) = BREAD(3) * BNORM

               ENDDO
             ENDDO
           ENDDO

           BMIN = BMIN * BNORM
           BMAX = BMAX * BNORM

C------- Mesh coordinates
           DX = (XH(2) - XH(1))*XNORM
           XH(1) = XH(1)*XNORM
           DO J=2,IXMA
             XH(J) =  XH(J-1) + DX
           ENDDO
           DY = (YH(2) - YH(1))*YNORM
           YH(1) = YH(1)*YNORM
           DO J=2,JYMA
             YH(J) =  YH(J-1) + DY
           ENDDO
           IF(IZ .EQ. 1) THEN 
             IZ1 = 1
           ELSE
             IZ1 = 2
           ENDIF
           DZ = (ZH(IZ1) - ZH(1))*ZNORM
           ZH(1) = ZH(1)*ZNORM
           DO K= 2, KZMA
             ZH(K) = ZH(K-1) + DZ
           ENDDO

      ELSEIF(MOD .EQ. 16) THEN
C     Full 3-D map of magnet. Can sum up field as read from several maps at particle
C     location, on-the-flight (maps are used as provided, NOT summed up into a single one). 
           JTCNT=0
           IRCNT = 0
           KZCNT=0       

           DO J=1,JYMA        
             JTC = J
             JTCNT = JTCNT + 1
             DO  K = 1,KZMA      
               KZC = K
               KZCNT = KZCNT+1
               DO I=1,IXMA            
                 IRCNT = IRCNT+1

                 IF(BINAR) THEN
                   READ(LUN,ERR=96)
     >             YH(J),ZH(K),XH(I),BREAD(2),BREAD(3),BREAD(1)
                 ELSE
                   READ(LUN,*,ERR=96) YH(J),ZH(K),XH(I), 
     >                                    BREAD(2),BREAD(3),BREAD(1)
                 ENDIF
                 BMAX0 = BMAX
                 BMAX = DMAX1(BMAX,BREAD(1),BREAD(2),BREAD(3))
                 IF(BMAX.NE.BMAX0) THEN
                   XBMA = XH(I)
                   YBMA = YH(J)
                   ZBMA = ZH(K)
                 ENDIF
                 BMIN0 = BMIN
                 BMIN = DMIN1(BMIN,BREAD(1),BREAD(2),BREAD(3))
                 IF(BMIN.NE.BMIN0) THEN
                   XBMI = XH(I)
                   YBMI = YH(J)
                   ZBMI = ZH(K)
                 ENDIF


                 IF(ZROBXY) THEN
                   IF(ZH(K) .EQ. 0.D0) THEN
                     HC(1,I,JTC,KZC,IMAP) = 0.D0
                     HC(2,I,JTC,KZC,IMAP) = 0.D0
                   ELSE
                     HC(1,I,JTC,KZC,IMAP) = BREAD(1) * BNORM
                     HC(2,I,JTC,KZC,IMAP) = BREAD(2) * BNORM
                   ENDIF
                 ELSE
                   HC(1,I,JTC,KZC,IMAP) = BREAD(1) * BNORM
                   HC(2,I,JTC,KZC,IMAP) = BREAD(2) * BNORM
                 ENDIF

C               write(*,*) ' fmapw ',I,JTC,KZC,IMAP,BNORM,BREAD(3)

                 HC(3,I,JTC,KZC,IMAP) = BREAD(3) * BNORM

               ENDDO
             ENDDO
           ENDDO

           BMIN = BMIN * BNORM
           BMAX = BMAX * BNORM

C------- Mesh coordinates
           DX = (XH(2) - XH(1))*XNORM
           XH(1) = XH(1)*XNORM
           DO J=2,IXMA
             XH(J) =  XH(J-1) + DX
           ENDDO
           DY = (YH(2) - YH(1))*YNORM
           YH(1) = YH(1)*YNORM
           DO J=2,JYMA
             YH(J) =  YH(J-1) + DY
           ENDDO
           IF(IZ .EQ. 1) THEN 
             IZ1 = 1
           ELSE
             IZ1 = 2
           ENDIF
           DZ = (ZH(IZ1) - ZH(1))*ZNORM
           ZH(1) = ZH(1)*ZNORM
           DO K= 2, KZMA
             ZH(K) = ZH(K-1) + DZ
           ENDDO

C FM. Comented May 2016. Seems redundat with summation in toscac, TBC...
C           IF    (IFIC .EQ. 1) THEN
C             DO I = 1, IXMA
C               DO J = 1, JYMA
C                 DO K = 1, KZMA
C                   HCTMP(1,I,J,K,IFIC) = HC(1,I,J,K,IMAP) 
C                   HCTMP(2,I,J,K,IFIC) = HC(2,I,J,K,IMAP) 
C                   HCTMP(3,I,J,K,IFIC) = HC(3,I,J,K,IMAP) 
C                 ENDDO
C               ENDDO
C             ENDDO
C           ELSE
C             IF    (IFIC.LT.MOD2) THEN
C               DO I = 1, IXMA
C                 DO J = 1, JYMA
C                   DO K = 1, KZMA
C                     HCTMP(1,I,J,K,IFIC) = HC(1,I,J,K,IMAP) 
C     >                 + HCTMP(1,I,J,K,IFIC) 
C                     HCTMP(2,I,J,K,IFIC) = HC(2,I,J,K,IMAP) 
C     >                 + HCTMP(2,I,J,K,IFIC) 
C                     HCTMP(3,I,J,K,IFIC) = HC(3,I,J,K,IMAP) 
C     >                 + HCTMP(3,I,J,K,IFIC) 
C                   ENDDO
C                 ENDDO
C               ENDDO
C             ELSEIF(IFIC .EQ. MOD2) THEN
C               DO I = 1, IXMA
C                 DO J = 1, JYMA
C                   DO K = 1, KZMA
C                     HC(1,I,J,K,IMAP) = HC(1,I,J,K,IMAP) 
C     >               + HCTMP(1,I,J,K,IFIC) 
C                     HC(2,I,J,K,IMAP) = HC(2,I,J,K,IMAP) 
C     >               + HCTMP(2,I,J,K,IFIC) 
C                     HC(3,I,J,K,IMAP) = HC(3,I,J,K,IMAP) 
C     >               + HCTMP(3,I,J,K,IFIC) 
C                   ENDDO
C                 ENDDO
C               ENDDO
C             ELSE
C               CALL ENDJOB(
C     >         'Pgm fmapw. No such possibility IFIC=',IFIC)
C             ENDIF
C           ENDIF

      ENDIF  

      IF(NRES.GT.0) THEN
        WRITE(NRES,*)' Sbr fmapw/fmapr3 : completed job of reading map.'
        CALL FLUSH2(NRES,.FALSE.)
      ENDIF
      RETURN

 96   CONTINUE
      CALL ENDJOB('Pgm fmapw. Problem reading into field map.'
     >//' Check field map data formatting',-99)
      RETURN

 97   CONTINUE
      IF(NRES.GT.0) THEN
        WRITE(NRES,*) ' WARNING !  Field map may be uncomplete'
        WRITE(NRES,FMT='(A,3(1X,I0))') '   End of file encountered '
     >  //' at jtcnt-1, ircnt,  kzcnt : ', jtcnt-1, ircnt,kzcnt
        WRITE(NRES,*) ' Check IXMAX and IRMAX ?'
        WRITE(6,*) ' WARNING !  Field map may be uncomplete'
        WRITE(6,FMT='(A,3(1X,I0))') '   End of file encountered '
     >  //' at jtcnt-1, ircnt,  kzcnt : ', jtcnt-1, ircnt,kzcnt
        WRITE(6,*) ' Check IXMAX and IRMAX ?'
      ENDIF
      RETURN

 98   CONTINUE
      IF(NRES.GT.0) THEN
        WRITE(NRES,*) ' WARNING !  Field map may be uncomplete'
        WRITE(NRES,FMT='(A,3(1X,I0))') '   Error encountered '
     >  //' at jtcnt-1, ircnt,  kzcnt : ', jtcnt-1, ircnt,kzcnt
        WRITE(NRES,*) '   Check map data file'
        WRITE(6,*) ' WARNING !  Field map may be uncomplete'
        WRITE(6,FMT='(A,3(1X,I0))') '   Error encountered '
     >  //' at jtcnt-1, ircnt,  kzcnt : ', jtcnt-1, ircnt,kzcnt
        WRITE(6,*) '   Check map data file'
      ENDIF
      RETURN

 99   CONTINUE
      RETURN

      ENTRY FMAPW2(IFICI,AAI)
      IFIC = IFICI
      AA = AAI
      RETURN

      ENTRY FMAPW4(IMAPI,BMIN,BMAX,XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)
           IMAP = IMAPI
           IIXMA(IMAP) = IXMA
           DO I=1,IXMA
             XXH(I,IMAP) =  XH(I)
           ENDDO
           JJYMA(IMAP) = JYMA
           DO J=1,JYMA
             YYH(J,IMAP) =  YH(J)
           ENDDO
           KKZMA(IMAP) = KZMA
           DO K= 1, KZMA
             ZZH(K,IMAP) = ZH(K)
           ENDDO
           BBMI(IMAP)  = BMIN
           BBMA(IMAP)  = BMAX
           XBBMI(IMAP) = XBMI
           YBBMI(IMAP) = YBMI
           ZBBMI(IMAP) = ZBMI
           XBBMA(IMAP) = XBMA
           YBBMA(IMAP) = YBMA
           ZBBMA(IMAP) = ZBMA
      RETURN

      ENTRY FMAPR5(IMAPI,
     >                  BMIN,BMAX,XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)
           IMAP = IMAPI
           IXMA = IIXMA(IMAP)
           DO I=1,IXMA
             XH(I) = XXH(I,IMAP)
           ENDDO
           JYMA = JJYMA(IMAP)
           DO J=1,JYMA
             YH(J) = YYH(J,imap)
           ENDDO
           KZMA = KKZMA(IMAP)
           DO K= 1, KZMA
             ZH(K) = ZZH(K,imap)
           ENDDO
           BMIN = BBMI(IMAP)
           BMAX = BBMA(IMAP)
           XBMI = XBBMI(IMAP)
           YBMI = YBBMI(IMAP)
           ZBMI = ZBBMI(IMAP)
           XBMA = XBBMA(IMAP)
           YBMA = YBBMA(IMAP)
           ZBMA = ZBBMA(IMAP)
      RETURN
      END
