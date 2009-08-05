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
      SUBROUTINE FMAPW(ACN,RFR,KART)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARIZ.H'
      INCLUDE "XYZHC.H"
C      COMMON//XH(MXX),YH(MXY),ZH(IZ),HC(ID,MXX,MXY,IZ),IXMA,JYMA,KZMA
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      COMMON/DROITE/ CA(9),SA(9),CM(9),IDRT
      LOGICAL ZSYM
      COMMON/OPTION/ KFLD,MG,LC,ML,ZSYM

      DIMENSION BREAD(3)
C----- at entry FMAPR :
C      LOGICAL IDLUNI, BINARI
      LOGICAL BINAR
      CHARACTER*120 TITL
      CHARACTER*20 FMTYP
      CHARACTER BE(2)
      SAVE MOD, MOD2

      DATA BE /'B','E'/

      DATA MOD, MOD2 / 0, 0 /

      IF(LF .EQ. 1) THEN
C------- Print field map in zgoubi.res
        IF(NRES.GT.0) THEN

          K = KZMA/2+1
          WRITE(NRES,
     >    FMT='(/,1X,20(1H-),/,10X,'' FIELD  MAP  (NORMALISED) :'',/)')
          WRITE(NRES,FMT='(/,A,3I2,/)') '  IXMA,JYMA,KZMA ',IXMA,JYMA,K
          WRITE(NRES,FMT='(
     >    5X,''Y'',4X,''Z'',4X,''X'',4X,A,''y'',4X,A,''z'',4X,A,''x'')') 
     >    BE(KFLD),BE(KFLD),BE(KFLD)

          CALL RAZ(BREAD,3)

             DO 14 J=1,JYMA
              DO  14  I = 1,IXMA
               DO 15 JD=1, ID
 15              BREAD(JD) = HC(JD,I,J,K)
                 WRITE(NRES,FMT='(1X,1P,6G11.2)') YH(J),ZH(K),XH(I),
     >                                  BREAD(2), BREAD(3), BREAD(1)
 14           CONTINUE
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
        WRITE(NMAP,991) ((HC(ID,I,J,K),I=1,IXMA),J=1,JYMA)
C        WRITE(NMAP,991) (XH(I),I=1,IXMA),(YH(J),J=1,JYMA)
C     >  ,((HC(ID,I,J,K),I=1,IXMA),J=1,JYMA)
 991    FORMAT(1P,6E18.10)

        CLOSE(NMAP)

        IF(NRES.GT.0) WRITE(NRES,FMT='(/,10X,A,/)') 
     >        'Field map has been written in zgoubi.map' 

      ENDIF
      RETURN

      ENTRY FMAPR(BINAR,LUN,
     >                      RM)

      MOD = 0 
      MOD2 = 0 

      IF(BINAR) THEN

C Map data file starts with 4-line header
C        DO 21 II=1, 4
        DO 21 II=1, 3
 21       READ(LUN) TITL
        READ(LUN) IXM,JYM,ACENT,RM,KRT,MD
  
        READ(LUN) (XH(I),I=1,IXMA),(YH(J),J=1,JYMA)
     >    ,((HC(ID,I,J,1),I=1,IXMA),J=1,JYMA)

      ELSE
C Map data file starts with 4-line header
C        DO 22 II=1, 4
        DO 22 II=1, 3
 22       READ(LUN,FMT='(A)') TITL
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
     >             ,((HC(ID,I,J,1),I=1,IXMA),J=1,JYMA)

 229    continue 
        goto 97

      ENDIF
      RETURN

C----------------------------------------------------------
C Read and interprete field maps in polar frame (MOD >= 20)
      ENTRY FMAPR2(BINAR,LUN, MODI,MODI2,NHDI, BNORM,
     >                       BMIN,BMAX,
     >                       XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)
      MOD = MODI
      MOD2 = MODI2
        NHD = NHDI

      BMIN =  1.D10
      BMAX = -1.D10

      IF(NRES.GT.0) 
     >   WRITE(NRES,FMT='(A,I1,A)') ' HEADER  (',NHD,' lines) : '

      IF(BINAR) THEN
C Map data file starts with 1-line header
        READ(LUN) R0, DR, DTTA, DZ    !    430  1. .1 1.  (cm cm deg cm)
        IF(NRES.GT.0) WRITE(NRES,*) '  R0, DR, DTTA, DZ :',R0,DR,DTTA,DZ
      ELSE
C Map data file starts with NHD-line header (NORMALLY 8)
        IF(NHD .GT. 0) THEN
          READ(LUN,FMT='(A120)',END=97,ERR=98) TITL
          READ(TITL,*) R0, DR, DTTA, DZ
          IF(NRES.GT.0) WRITE(NRES,*)'R0, DR, DTTA, DZ : ',R0,DR,DTTA,DZ
          DO 39 II=1, NHD-1
            READ(LUN,FMT='(A120)',END=97,ERR=98) TITL
 39         IF(NRES.GT.0) WRITE(NRES,FMT='(5X,A120)') TITL
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
C                   write(33) DUM1,DUM2,DUM3,BREAD(1),BREAD(2),BREAD(3)
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
                 HC(1,JTC,I,KZC) = BTTA * BNORM
                 HC(2,JTC,I,KZC) = BRAD * BNORM
                 HC(3,JTC,I,KZC) = BREAD(2) * BNORM
 33        CONTINUE

        BMIN = BMIN * BNORM
        BMAX = BMAX * BNORM

C------- symmetrise 3D map wrt magnet vertical symm plane
        DO 34  K=1,KKZMA    
          KZC = KKZMA-1+K
          DO 34 I=1,IRMA      
            DO 34 J=1,JTMA-1    
              IF    (MOD.EQ.20) THEN
                JTC = J
                JTS = 2*JTMA-J
              ELSEIF(MOD.EQ.21) THEN
                JTC = JTMA+J
                JTS = JTMA-J
              ENDIF
C              HC(1,JTS,I,KZC) = HC(1,JTC,I,KZC)
              HC(1,JTS,I,KZC) = -HC(1,JTC,I,KZC)
              HC(2,JTS,I,KZC) = HC(2,JTC,I,KZC) 
              HC(3,JTS,I,KZC) = HC(3,JTC,I,KZC) 
 34     CONTINUE

C------- symmetrise 3D map wrt mid-plane= bend-plane
        DO 35  K=2,KKZMA      
          KZC = KKZMA-1+K
          KZS = KKZMA-K+1
          DO 35 I=1,IRMA         
            DO 35 J=1,JTMA*2-1    
              HC(1,J,I,KZS) = -HC(1,J,I,KZC)
              HC(2,J,I,KZS) = -HC(2,J,I,KZC) 
              HC(3,J,I,KZS) = HC(3,J,I,KZC) 
 35     CONTINUE

C------- Mesh coordinates
          DO 37 J=1,JYMA
 37          YH(J) =  R0 + DBLE(J-1) *DR 
          DO 38 K= -KZMA/2,KZMA/2   
 38         ZH(K+KZMA/2+1) = DBLE(K) * DZ
          DO 36 I=-IXMA/2,IXMA/2
 36         XH(I+IXMA/2+1) =  DBLE(I) * DTTA

      ELSEIF(MOD .EQ. 22 .OR. MOD .EQ. 24) THEN

           IRMA = JYMA
           JTMA = IXMA
           KKZMA = KZMA/2+1

           jtcnt=0
           ircnt = 0
           kzcnt=0       
C             write(*,*) 'irma,kkzma,jtma : ', irma,kkzma,jtma
           DO 133 I=1,IRMA            
             RADIUS=R0 + DBLE(I-1) *DR 
             ircnt = ircnt+1
             DO 133  K = 1,KKZMA      
               KZC = KKZMA-1+K
               Z = DBLE(K-1) * DZ
               kzcnt = kzcnt+1
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
C                 write(33,*) DUM1,zzz,DUM3,BREAD(1),BREAD(2),BREAD(3)

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
c                    write(88,*) sqrt(dum1*dum1 + dum3*dum3 ), BREAD(2), 
c     >                 tta*180/(4.*atan(1.)), alp*180/(4.*atan(1.)), 
c     >                         ' fmapr'
C---------------- Polar components at positon (r,tta) in cyl. frame
                 BRAD = BH * COS(ALP)
                 BTTA = BH * SIN(ALP)
C---------------- Watch the sign !! Bx and By multiplied by (-1)
                 HC(1,JTC,I,KZC) =   -  BTTA * BNORM
                 HC(2,JTC,I,KZC) =   -  BRAD * BNORM
                 HC(3,JTC,I,KZC) = BREAD(2) * BNORM

C---------------- In case that non-zero Bx, By in median plane would be prohibitive
                 IF(ZZZ.EQ.0.D0) THEN
                    HC(1,JTC,I,KZC) = 0.D0
                    HC(2,JTC,I,KZC) = 0.D0
C                 ELSE
C                    WRITE(*,*) ' FMAPR2  Z, HC ',ZZZ,HC(3,JTC,I,KZC),kkzma,z
                 ENDIF   

C----------------  TEST RACCAM
C     >          HC(1,JTC,I,KZC) = HC(1,JTC,I,KZC)  *(RADIUS/348.)**(-.2)
C     >          HC(2,JTC,I,KZC) = HC(2,JTC,I,KZC)  *(RADIUS/348.)**(-.2)
C     >          HC(3,JTC,I,KZC) = HC(3,JTC,I,KZC)  *(RADIUS/348.)**(-.2)
C                    write(*,*) ' HC_1-3 multiplied by (r/r_0)^dK'
                     
 133       CONTINUE

        BMIN = BMIN * BNORM
        BMAX = BMAX * BNORM

C------- symmetrise 3D map wrt mid-plane= bend-plane
        DO 135  K=2,KKZMA      
          KZC = KKZMA-1+K
          KZS = KKZMA-K+1
          DO 135 I=1,IRMA         
            DO 135 J=1,JTMA
              HC(1,J,I,KZS) = -HC(1,J,I,KZC)
              HC(2,J,I,KZS) = -HC(2,J,I,KZC) 
              HC(3,J,I,KZS) = HC(3,J,I,KZC) 
 135    CONTINUE

C------- Mesh coordinates
          DO 137 J=1,JYMA
 137         YH(J) =  R0 + DBLE(J-1) *DR 
          DO 138 K= -KZMA/2,KZMA/2   
 138        ZH(K+KZMA/2+1) = DBLE(K) * DZ
          DO 136 I=-IXMA/2,IXMA/2
 136        XH(I+IXMA/2+1) =  DBLE(I) * DTTA

C        write(89,fmt='(a,3i4)') ' itta, jr, kz =',ixma,jyma,kzma 
C        do  kk=1,kzMA
C          do jj=1,jyma 
C            do  i=1,ixMA
C             write(89,fmt='(6E14.6,3i4)') 
C     >        HC(1,i,jj,kkzma),HC(2,i,jj,kkzma),
C     >        HC(3,i,jj,kkzma)
C     >        ,xh(i),yh(jj),zh(kk),i,jj,kk
C            enddo
C          enddo
C          write(89,*) ' ' 
C        enddo

      ENDIF  ! MOD
      RETURN

C-------------------------------------------------------------
C Read and interprete field maps in cartesian frame (MOD < 20)
      ENTRY FMAPR3(BINAR,LUN,MODI,MODI2,NHDI,BNORM,I1,KZI,FMTYP,
     >                       BMIN,BMAX,
     >                       XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)
      
      MOD = MODI
      MOD2 = MODI2
        NHD = NHDI
      KZ = KZI
      BMIN =  1.D10
      BMAX = -1.D10

C Map data file starts with NHD-line header
      IF(NRES.GT.0) 
     >  WRITE(NRES,FMT='(A,I1,A)') ' HEADER  (',NHD,' lines) : '

      IF(NHD .GE. 1) THEN
        READ(LUN,FMT='(A120)') TITL
        IF(NRES.GT.0) WRITE(NRES,FMT='(5X,A120)') TITL
        READ(TITL,*,END=49,ERR=49) R0, DR, DX, DZ
        GOTO 50
 49     CONTINUE
        IF(NRES.GT.0) WRITE(NRES,*) 
     >         ' Error upon readinng R0, DR, DX, DZ'
 50     CONTINUE
        IF(NHD .GE. 2) THEN
          READ(LUN,FMT='(A120)') TITL
          IF(NRES.GT.0) WRITE(NRES,FMT='(5X,A120)') TITL
          IF(NHD .GE. 2) THEN
            DO 32 II=1, NHD-2
              READ(LUN,FMT='(A120)') TITL
              IF(NRES.GT.0) WRITE(NRES,FMT='(5X,A120)') TITL
 32         CONTINUE
            IF(NRES.GT.0) THEN
              WRITE(NRES,*) ' R0, DR, DX, DZ : ',R0,DR,DX,DZ
            ENDIF
          ENDIF
        ENDIF
      ENDIF
 
      IF(MOD .EQ. 12) THEN

           jtcnt=0
           ircnt = 0
           kzcnt=0       
C             write(*,*) 'iXma, jYma, kzma : ', iXma,jYma,kzma

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
C                   READ(LUN,*) YH(J),ZH(K),XH(I), 
     >                                    BREAD(2),BREAD(3),BREAD(1)
C                   write(*,*) YH(J),ZH(K),XH(I), 
C     >                       BREAD(2),BREAD(3),BREAD(1),i,j,k
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

                 HC(1,I,JTC,KZC) = BREAD(1) * BNORM
                 HC(2,I,JTC,KZC) = BREAD(2) * BNORM
                 HC(3,I,JTC,KZC) = BREAD(3) * BNORM

 121          CONTINUE

        BMIN = BMIN * BNORM
        BMAX = BMAX * BNORM


C------- symmetrise 3D map wrt mid-plane= bend-plane
        DO 122  K=2,KZMA      
          KZC = KZMA-1+K
          KZS = KZMA-K+1
          DO 122 I=1,IXMA         
            DO 122 J=1,JYMA
              HC(1,J,I,KZS) = -HC(1,J,I,KZC)
              HC(2,J,I,KZS) = -HC(2,J,I,KZC) 
              HC(3,J,I,KZS) = HC(3,J,I,KZC) 
 122       CONTINUE

C------- Mesh coordinates
          DY = YH(KZMA*IXMA) - YH(1)
          DO 127 J=2,JYMA
 127          YH(J) =  YH(1) + DY
          DZ = ZH(IXMA) - ZH(1)
          DO 128 K= 2, 2*KZMA-1
 128         ZH(K) = ZH(1) + DZ
C            write(*,*) 'sbr fmapw ',xh(1),xh(2),yh(1),yh(2),zh(1),zh(2)

      ELSEIF(MOD.EQ.0 .OR. MOD.EQ.1) THEN
C------- GSI spectro for instance. 
C        Keyword EMMA with two 2D maps 
C        MOD=0 : # of files is NF=1+|IZ/2|, from z=0 to z_max, symmetrizing wrt median plane
C        MOD=1 : # of files is NF= IZ, from -z_max to +z_max, no symmetrizing

           I = KZ

c                 write(*,*) 
c           write(*,*) ' kz, jyma, ixma, lun',kz, jyma, ixma,bnorm,lun
c                 write(*,*) 

             DO 10 J=1,JYMA
               DO  10  K = 1,IXMA
                 IF( BINAR ) THEN
                   READ(LUN)YH(J),ZH(I),XH(K),BREAD(2),BREAD(3),BREAD(1)
                 ELSE

C----  Manip VAMOS ganil, oct 2001
C                   READ(LUN,FMT='(1X,6G12.2)') YH(J),ZH(I),XH(K), 

                   IF    (MOD2.EQ.1) THEN 
                     READ(LUN,FMT='(1X,6E11.2)') YH(J),ZH(I),XH(K), 
     >                                      BREAD(2),BREAD(3),BREAD(1)
                   ELSE
C-------------------- Default MOD, also case MOD=0
                     IF(FMTYP.EQ.'GSI') THEN
                       READ(LUN,fmt='(1x,6e11.2)') YH(J),ZH(I),XH(K),  
     >                                      BREAD(2),BREAD(3),BREAD(1)
                     ELSEIF(FMTYP.EQ.'LESB3') THEN
                       READ(LUN,fmt='(1x,6e11.2)') YH(J),ZH(I),XH(K),  
     >                                      BREAD(2),BREAD(3),BREAD(1)
                     ELSE
                       READ(LUN,*) YH(J),ZH(I),XH(K),  
     >                                      BREAD(2),BREAD(3),BREAD(1)
                     ENDIF
                   ENDIF
C                   write(89,fmt='(1p,6e12.4,3i4)') YH(J),ZH(I),XH(K), 
C     >                BREAD(2),BREAD(3),BREAD(1),j,i,k

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
                   HC(ID,K,J,I) =  BREAD(3) * BNORM
                 ELSE
                   DO 187 LHC=1,ID
 187                  HC(LHC,K,J,I) = BREAD(LHC) * BNORM
                 ENDIF

                 IF(I1 .GT. 1) THEN
C------------------- Symmetrize 3D map wrt XY plane (IZ>1)
                   IF( I .GT. I1 ) THEN
                     ZH(2*I1-I) = -ZH(I)
                     FAC = -1.D0
C                     DO 188 LHC=1,ID
C                       IF(LHC .EQ. 3) FAC=1.D0
C 188                   HC(LHC,K,J,2*I1-I) =  FAC * BREAD(LHC)* BNORM
                     HC(1,K,J,2*I1-I) =   -BREAD(1)* BNORM
                     II=2
                     HC(II,K,J,2*I1-I) =   -BREAD(2)* BNORM
                     II=3
                     HC(II,K,J,2*I1-I) =   BREAD(3)* BNORM
                   ENDIF
                 ENDIF
 10          CONTINUE
             CLOSE(UNIT=LUN)
C 12        CONTINUE

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
