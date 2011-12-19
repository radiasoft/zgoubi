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
C  Brookhaven National Laboratory               és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE FMAPW(ACN,RFR,KART)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARIZ.H'
      INCLUDE "XYZHC.H"
C      COMMON//XH(MXX),YH(MXY),ZH(IZ),HC(ID,MXX,MXY,IZ,IMAP),IXMA,JYMA,KZMA
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      COMMON/DROITE/ CA(9),SA(9),CM(9),IDRT
      LOGICAL ZSYM
      COMMON/OPTION/ KFLD,MG,LC,ML,ZSYM

      DIMENSION BREAD(3)
C----- at entry FMAPR :
C      LOGICAL IDLUNI, BINARI
      LOGICAL BINAR, EMPTY
      CHARACTER*120 TITL
      CHARACTER*20 FMTYP
      PARAMETER (MXCHAR=200)
      CHARACTER TXT132*200
      CHARACTER BE(2)
      INTEGER DEBSTR, FINSTR
      SAVE MOD, MOD2
      
      DATA BE /'B','E'/

      DATA MOD, MOD2 / 0, 0 /
      DATA IMAP / 1 /

      CALL KSMAP(
     >           IMAP) 

      IF(LF .EQ. 1) THEN
C------- Print field map in zgoubi.res
        IF(NRES.GT.0) THEN

          K = KZMA/2+1
          WRITE(NRES,
     >    FMT='(/,1X,20(''-''),/,10X,'' FIELD  MAP  (NORMALISED) :'',/)'
     $         )
          WRITE(NRES,FMT='(/,A,3I2,/)') '  IXMA,JYMA,KZMA ',IXMA,JYMA,K
          WRITE(NRES,FMT='(
     >    5X,''Y'',4X,''Z'',4X,''X'',4X,A,''y'',4X,A,''z'',4X,A,''x'')') 
     >    BE(KFLD),BE(KFLD),BE(KFLD)

          CALL RAZ(BREAD,3)

             DO J=1,JYMA
               DO I = 1,IXMA
                 DO JD=1, ID
                   BREAD(JD) = HC(JD,I,J,K,IMAP)
                 ENDDO
                   WRITE(NRES,FMT='(1X,1P,6G11.2)') YH(J),ZH(K),XH(I),
     >                                  BREAD(2), BREAD(3), BREAD(1)
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
 992    FORMAT(2I3,1P,2E15.7,1X,I1,I3)

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
      ENTRY FMAPR(BINAR,LUN,
     >                      RM)

      CALL KSMAP(
     >           IMAP) 
      MOD = 0 
      MOD2 = 0 

      IF(BINAR) THEN

C Map data file starts with 4-line header
C        DO 21 II=1, 4
        DO II=1, 3
          READ(LUN) TITL
        ENDDO
        READ(LUN) IXM,JYM,ACENT,RM,KRT,MD
  
        READ(LUN) (XH(I),I=1,IXMA),(YH(J),J=1,JYMA)
     >    ,((HC(ID,I,J,1,IMAP),I=1,IXMA),J=1,JYMA)

      ELSE
C Map data file starts with 4-line header
C        DO 22 II=1, 4
        DO II=1, 3
          READ(LUN,FMT='(A)') TITL
        ENDDO
        READ(LUN,*) IXM,JYM,ACENT,RM,KRT,MD

        IF(NRES.GT.0) THEN
          IF(IXM .NE. IXMA) THEN
            WRITE(6,100) 'IXM','IXMA'
            WRITE(NRES,100) 'IXM','IXMA'
 100        FORMAT('WARNING - info !! ',A, 
     >       ' from map file .ne. ',A,' from .dat file.')
          ENDIF
          IF(JYM .NE. JYMA) THEN
            WRITE(6,100) 'JYM','JYMA'
            WRITE(NRES,100) 'JYM','JYMA'
          ENDIF
        ENDIF

        READ(LUN,*,END=229,ERR=98) (XH(I),I=1,IXMA),(YH(J),J=1,JYMA)
     >             ,((HC(ID,I,J,1,IMAP),I=1,IXMA),J=1,JYMA)

 229    CONTINUE 
        GOTO 97

      ENDIF
      RETURN

C----------------------------------------------------------
C Read and interprete field maps in polar frame (MOD >= 20)
      ENTRY FMAPR2(BINAR,LUN, MODI,MODI2,NHDI, BNORM,
     >                       BMIN,BMAX,
     >                       XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)

      CALL KSMAP(
     >           IMAP) 

      IF(NRES.GT.0) WRITE(NRES,*)'IMAP : ',IMAP
C      IF(IMAP.LT.0) CALL ENDJOB('SBR FMAPW, IMAP has to be .ge.0',-99)
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
          READ(LUN,END=97,ERR=98) TITL
          READ(TITL,*,ERR=41,END=41) R0, DR, DTTA, DZ
          IF(NRES.GT.0) WRITE(NRES,*)'R0, DR, DTTA, DZ : ',R0,DR,DTTA,DZ
          GOTO 42
 41       CONTINUE
          IF(NRES.GT.0) WRITE(NRES,*) 
     >        'Could not read R0, DR, DTTA, DZ from line 1 of HEADER...'
 42       CONTINUE
          DO II=1, NHD-1
            READ(LUN,END=97,ERR=98) TITL
            IF(NRES.GT.0) WRITE(NRES,FMT='(5X,A120)') TITL
          ENDDO
          IF(NRES.GT.0) WRITE(NRES,FMT='(A)') ' '
        ENDIF
      ELSE
        IF(NHD .GE. 1) THEN
C Map data file starts with NHD-line header (NORMALLY 8)
          READ(LUN,FMT='(A120)',END=97,ERR=98) TITL
          READ(TITL,*,ERR=39,END=39) R0, DR, DTTA, DZ
          IF(NRES.GT.0) WRITE(NRES,*)'R0, DR, DTTA, DZ : ',R0,DR,DTTA,DZ
          GOTO 40
 39       CONTINUE
          IF(NRES.GT.0) WRITE(NRES,*) 
     >        'Could not read R0, DR, DTTA, DZ from line 1 of HEADER...'
 40       CONTINUE
          DO II=1, NHD-1
            READ(LUN,FMT='(A120)',END=97,ERR=98) TITL
            IF(NRES.GT.0) WRITE(NRES,FMT='(5X,A120)') TITL
          ENDDO
          IF(NRES.GT.0) WRITE(NRES,FMT='(A)') ' '
        ENDIF
      ENDIF

      DTTA = DTTA * RAD

      IF(MOD .NE. 22 .AND. MOD .NE. 24) THEN
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

C----------- Read the bare map data
           IRMA = JYMA
           JTMA = IXMA/2+1
           KKZMA = KZMA/2+1

           DO 33 I=1,IRMA            
             RADIUS=R0 + DBLE(I-1) *DR 
C             WRITE(*,*)'fmapw R0, radius : ',R0,radius
             DO 33  K = 1,KKZMA      
               KZC = KKZMA-1+K
               Z = DBLE(K-1) * DZ
               DO 33 J=1,JTMA        
                 IF    (MOD.EQ.20) THEN
                   JTC = J
                 ELSEIF(MOD.EQ.21) THEN
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
        DO 34  K=1,KKZMA    
          KZC = KKZMA-1+K
          if(kzc.gt.iz) call endjob(' FMAPW, K should be <',IZ+1)
          DO 34 I=1,IRMA      
            if(i.gt.mxy) call endjob(' FMAPW, I should be <',mxy+1)
            DO 34 J=1,JTMA-1    
              IF    (MOD.EQ.20) THEN
                JTC = J
                JTS = 2*JTMA-J
              ELSEIF(MOD.EQ.21) THEN
                JTC = JTMA+J
                JTS = JTMA-J
              ENDIF
              if(jts.gt.mxx) call endjob(' FMAPW, J should be <',mxx+1)
              if(jtc.gt.mxx) call endjob(' FMAPW, J should be <',mxx+1)
              HC(1,JTS,I,KZC,IMAP) = -HC(1,JTC,I,KZC,IMAP)
              HC(2,JTS,I,KZC,IMAP) =  HC(2,JTC,I,KZC,IMAP) 
              HC(3,JTS,I,KZC,IMAP) =  HC(3,JTC,I,KZC,IMAP) 
 34     CONTINUE

C------- symmetrise 3D map wrt mid-plane= bend-plane
        DO 35  K=2,KKZMA      
          KZC = KKZMA-1+K
          KZS = KKZMA-K+1
          DO 35 I=1,IRMA         
            DO 35 J=1,JTMA*2-1    
              HC(1,J,I,KZS,IMAP) = -HC(1,J,I,KZC,IMAP)
              HC(2,J,I,KZS,IMAP) = -HC(2,J,I,KZC,IMAP) 
              HC(3,J,I,KZS,IMAP) = HC(3,J,I,KZC,IMAP) 
C              WRITE(*,*)'fmapw j,i,kzs,kzc : ',j,i,kzs,kzc 
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


      ELSEIF(MOD .EQ. 22 .OR. MOD .EQ. 24) THEN

           IRMA = JYMA
           JTMA = IXMA
           KKZMA = KZMA/2+1

           ircnt = 0
           DO 133 I=1,IRMA            
             RADIUS=R0 + DBLE(I-1) *DR 
             ircnt = ircnt+1
             kzcnt=0       
             DO 133  K = 1,KKZMA      
               KZC = KKZMA-1+K
               Z = DBLE(K-1) * DZ
               kzcnt = kzcnt+1
              jtcnt=0
              DO 133 J=1,JTMA        
                   JTC = J
                 TTA =  DBLE(J-1) * DTTA
                 jtcnt = jtcnt + 1
                 IF(BINAR) THEN
                   READ(LUN,END=97,ERR=98) DUM1,ZZZ,DUM3,
     >                                   BREAD(1),BREAD(2),BREAD(3)
                 ELSE
                   READ(LUN,*,END=97,ERR=98) DUM1,ZZZ,DUM3,
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
                 IF(ZZZ.EQ.0.D0) THEN
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

C------- symmetrise 3D map wrt mid-plane= bend-plane
        DO 135  K=2,KKZMA      
          KZC = KKZMA-1+K
          KZS = KKZMA-K+1
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

      ENDIF  ! MOD
      RETURN

C-------------------------------------------------------------
C Read and interprete field maps in cartesian frame (MOD < 20)
      ENTRY FMAPR3(BINAR,LUN,MODI,MODI2,NHDI,
     >             XNORM,YNORM,ZNORM,BNORM,I1,KZI,FMTYP,
     >                       BMIN,BMAX,
     >                       XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)
      
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
C Map data file starts with NHD-line header (NORMALLY 8)
          READ(LUN,END=97,ERR=98) TITL
          READ(TITL,*,ERR=51,END=51) R0, DR, DTTA, DZ
          GOTO 52
 51       CONTINUE
          IF(NRES.GT.0) WRITE(NRES,*) 
     >        'Could not read R0, DR, DTTA, DZ from line 1 of HEADER...'
 52       CONTINUE
          IF(NHD .GE. 2) THEN
            READ(LUN) TITL
            IF(NRES.GT.0) WRITE(NRES,FMT='(5X,A120)') TITL
            IF(NHD .GE. 3) THEN
              DO II=1, NHD-2
                READ(LUN) TITL
                IF(NRES.GT.0) WRITE(NRES,FMT='(5X,A120)') TITL
              enddo
              IF(NRES.GT.0) THEN
                WRITE(NRES,*) ' R0, DR, DX, DZ : ',R0,DR,DX,DZ
                CALL FLUSH2(NRES,.FALSE.)
              ENDIF
            ENDIF
          ENDIF
        ENDIF

      ELSE

        IF(NHD .GE. 1) THEN
          READ(LUN,FMT='(A120)') TITL
          IF(NRES.GT.0) WRITE(NRES,FMT='(5X,A120)') TITL
          READ(TITL,*,END=49,ERR=49) R0, DR, DX, DZ
          GOTO 50
 49       CONTINUE
          IF(NRES.GT.0) WRITE(NRES,*) 
     >        'Could not read R0, DR, DX, DZ from line 1 of HEADER...'
 50       CONTINUE
          IF(NHD .GE. 2) THEN
            READ(LUN,FMT='(A120)') TITL
            IF(NRES.GT.0) WRITE(NRES,FMT='(5X,A120)') TITL
            IF(NHD .GE. 3) THEN
              DO II=1, NHD-2
                READ(LUN,FMT='(A120)') TITL
                IF(NRES.GT.0) WRITE(NRES,FMT='(5X,A120)') TITL
              enddo
              IF(NRES.GT.0) THEN
                WRITE(NRES,*) ' R0, DR, DX, DZ : ',R0,DR,DX,DZ
                CALL FLUSH2(NRES,.FALSE.)
              ENDIF
            ENDIF
          ENDIF
        ENDIF

      ENDIF

      IF(MOD .EQ. 12) THEN
C Read 3D field map contained in a single file

        IF(MOD2 .EQ. 0) THEN         ! Default
           jtcnt=0
           ircnt = 0
           kzcnt=0       

           DO 121 I=1,IXMA            
             ircnt = ircnt+1
             DO 121  K = 1,KZMA      
               KZC = KZMA-1+K
               kzcnt = kzcnt+1
               DO 121 J=1,JYMA        
                 JTC = J
                 jtcnt = jtcnt + 1

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

C------- symmetrise 3D map wrt mid-plane= bend-plane
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
C Differs from MOD2=0 by just READ format and absence of symmetrization
C Used for AGS snakes maps
           jtcnt=0
           ircnt = 0
           kzcnt=0       
C           write(*,*) 'sbr fmapr3 ixma, kzma, jyma: ',ixma,kzma,jyma 
           DO J=1,JYMA        
             JTC = J
             jtcnt = jtcnt + 1
             DO  K = 1,KZMA      
               kzc = k
               kzcnt = kzcnt+1
               DO I=1,IXMA            
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
           IF(IZ.EQ.1) THEN 
             IZ1 = 1
           ELSE
             IZ1 = 2
           ENDIF
           DZ = (ZH(IZ1) - ZH(1))*Znorm
           ZH(1) = ZH(1)*Znorm
           DO K= 2, KZMA
             ZH(K) = ZH(K-1) + DZ
           ENDDO


        ENDIF ! MOD2=0, 1

      ELSEIF(MOD.EQ.0 .OR. MOD.EQ.1) THEN
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
                     IF(FINSTR(TXT132).EQ.MXCHAR) 
     >                  CALL ENDJOB('SBR FMAPW : # of columns in field' 
     >                  //' data file must be <',MXCHAR)
                   ENDIF

                   IF    (MOD2.EQ.1) THEN 
                     READ(TXT132,FMT='(1X,6E11.2)') YH(J),ZH(I),XH(K), 
     >                                        BREAD(2),BREAD(3),BREAD(1)
                   ELSEIF(MOD2.EQ.2) THEN 
                     READ(TXT132,FMT='(1X,6E12.2)') YH(J),ZH(I),XH(K), 
     >                                        BREAD(2),BREAD(3),BREAD(1)
                   ELSEIF(MOD2.EQ.3) THEN 
                     READ(TXT132,*) YH(J),ZH(I),XH(K), 
     >                                        BREAD(2),BREAD(3),BREAD(1)
                   ELSE
C-------------------- Default MOD2
                     IF(FMTYP.EQ.'GSI') THEN
                       READ(TXT132,FMT='(1X,6E11.2)') YH(J),ZH(I),XH(K),  
     >                                      BREAD(2),BREAD(3),BREAD(1)
                     ELSEIF(FMTYP.EQ.'LESB3') THEN
                       READ(TXT132,FMT='(1X,6E11.2)') YH(J),ZH(I),XH(K),  
     >                                      BREAD(2),BREAD(3),BREAD(1)
                     ELSE
                       READ(TXT132,*) YH(J),ZH(I),XH(K),  
C                       READ(TXT132,FMT='(1X,6E11.2)') YH(J),ZH(I),XH(K),  
     >                                      BREAD(2),BREAD(3),BREAD(1)
                     ENDIF
                   ENDIF
c                   write(89,fmt='(1p,6e12.4,6i8)') YH(J),ZH(I),XH(K), 
c     >                BREAD(2),BREAD(3),BREAD(1),j,i,k,kz,jyma,ixma
c        write(89,*) j,i,k,jyma,KZ,ixma,imap,' J,I,k,jyma,KZ,ixma,imap'

CC----  Manip VAMOS ganil, oct 2001 
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
                 IF(ID.EQ.1) THEN
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

      ELSEIF(MOD.EQ.3) THEN
C------- AGS magnet maps (2D map, half-magnet, 2b symmetrized wrt YZ plane at 45inches)

           K = 1
           
           READ(LUN,*) JYMA
           IF(JYMA.GT.MXY) CALL 
     >       ENDJOB(' SBR fmapw. In PARIZ.H need MXY .ge. ',JYMA)
           READ(LUN,*) (YH(J),J=1,JYMA)
           DO J=1,JYMA
             YH(J) = YH(J) * YNORM
           ENDDO

           READ(LUN,*) IXMA2
           IXMA = 2*IXMA2 - 1
           IF(IXMA.GT.MXX) CALL 
     >       ENDJOB(' SBR fmapw. In PARIZ.H need MXX .ge. ',IXMA)
           READ(LUN,*) (XH(I),I=1,IXMA2)
           DO I=1,IXMA2
             XH(I) = XH(I) * XNORM
C                write(*,*) ' fmapw i, xh : ', i,xh(i)
           ENDDO
           DX = XH(2) - XH(1)
           DO I=IXMA2+1,IXMA
             XH(I) = XH(IXMA2) + DBLE(I-IXMA2)*DX
C                write(*,*) ' fmapw i, xh : ', i,xh(i)
           ENDDO

C           DO I=1,IXMA
C                write(*,*) ' fmapw i, xh(i) : ',i, xh(i)  
C           ENDDO

           READ(LUN,*) NTOT  

           DO I = 1, IXMA2
             READ(LUN,FMT='(A132)') TXT132
C             IF(NRES.GT.0) WRITE(NRES,*) 'Read from field map : ',TXT132
             READ(LUN,*) (HC(ID,I,J,K,IMAP), J=1,JYMA)
C             write(*,*) (HC(ID,I,J,K,IMAP), J=1,JYMA)
C                 stop
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
      ENDIF  

      IF(NRES.GT.0) THEN
        WRITE(NRES,*)' SBR FMAPW/FMAPR3 : completed job of reading map.'
        CALL FLUSH2(NRES,.FALSE.)
      ENDIF
      RETURN

 97   CONTINUE
      IF(NRES.GT.0) THEN
        WRITE(NRES,*) ' WARNING !  Field map may be uncomplete'
        WRITE(NRES,*) '   End of file encountered before',
     >  ' jtcnt-1, ircnt,  kzcnt : ', jtcnt-1, ircnt,kzcnt
        WRITE(NRES,*) ' Check IXMAX and IRMAX ?'
        WRITE(6,*) ' WARNING !  Field map may be uncomplete'
        WRITE(6,*) '   End of file encountered before',
     >  ' jtcnt-1, ircnt,  kzcnt : ', jtcnt-1, ircnt,kzcnt
        WRITE(6,*) ' Check IXMAX and IRMAX ?'
      ENDIF
      RETURN

 98   CONTINUE
      IF(NRES.GT.0) THEN
        WRITE(NRES,*) ' WARNING !  Field map may be uncomplete'
        WRITE(NRES,*) '   Error encountered during map reading'
        WRITE(NRES,*) '   Check map data file'
        WRITE(6,*) ' WARNING !  Field map may be uncomplete'
        WRITE(6,*) '   Error encountered during map reading'
        WRITE(6,*) '   Check map data file'
      ENDIF
      RETURN
      END
