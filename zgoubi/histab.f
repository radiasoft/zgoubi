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
C  USA
C  -------
      SUBROUTINE HISTAB(IC1,IC2,NLIN,KAR,MODE,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     -------------------------------------------------------------
C     Y EST UN HISTOGRAMME DE VALEURS
C     ON TRACE Y SUR LISTING, DANS LES LIMITES UTILISATEUR
C       IC1, IC2( <120 ) EN CANAL , XMIN, XMAX EN PHYSIQ.
C     LA FENETRE 'ECRAN' EST IC1-IC2 (NUMEROS DE COLONNES) EN
C       HORIZONTAL, ET NLIN (NOMBRE DE LIGNES LISTING) EN VERTICAL
C     MODE=0 (PAS DE NORMALISTION) OU MODE=1 (NORMALISATION /NCMAX)
C     -------------------------------------------------------------
 
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      PARAMETER (JH=24,KH=5)
      COMMON/HISTO/ ICTOT(JH,KH),MOYC(JH,KH) ,CMOY(JH,KH),JMAX(JH,KH)
      COMMON/HISTOG/NC(JH,120,KH),JJ,NH,XMI(JH,KH),XMO(JH,KH),XMA(JH,KH)
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
 
      DIMENSION ISTO(130)
      CHARACTER(1) KARL(130),KARO(10),KAR,BLANC
      DATA KARO / '1','2','3','4','5','6','7','8','9','0'/
      DATA BLANC/ ' '/
C     ** numero de la 1-er colonne sur  le listing
      DATA IC0/ 10 /
 
      IC10=IC0+IC1
      IC20=IC0+IC2
 
C     *** COMPTAGE, MOYENNE
      CMOY(JJ,NH)=0D0
      MOYC(JJ,NH)=0
      ICTOT(JJ,NH)=0
      DO 6 IC=IC1,IC2
        ICTOT(JJ,NH)=ICTOT(JJ,NH) + NC(JJ,IC,NH)
 6      CMOY(JJ,NH) = CMOY(JJ,NH) + IC*NC(JJ,IC,NH)
 
      IF(ICTOT(JJ,NH) .GT. 0) THEN
        MOYC(JJ,NH) = NINT(CMOY(JJ,NH)/ICTOT(JJ,NH))
 
        DO 5 IC=IC1,IC2
 5        ISTO(IC+IC0)=NC(JJ,IC,NH)
 
        IF(MODE .EQ. 1) THEN
C         ** NORMALISE NC (DANS LA FENETRE IC1-IC2)
          YMAX=0D0
          DO 4 IC=IC1,IC2
            IF(YMAX .LT. NC(JJ,IC,NH)) YMAX=NC(JJ,IC,NH)
 4        CONTINUE
          IF(YMAX .NE. 0.D0) THEN
            DO 3  IC=IC10,IC20
 3            ISTO(IC) = .9D0 * ISTO(IC) * NLIN/YMAX
          ENDIF
        ENDIF
 
C       ** TRACE L'HISTO
        DO 1 LINE=NLIN,1,-1
          DO 2 IC=IC10,IC20
            IF(ISTO(IC) .GE. LINE) THEN
              IF( (LINE/10)*10 .EQ. LINE ) THEN
                KARL(IC)='0'
              ELSE
                KARL(IC)=KAR
              ENDIF
            ELSE
              KARL(IC)=' '
            ENDIF
 2        CONTINUE
          WRITE(NRES,100) LINE, (BLANC,IC=6,IC10-1)
     >    , (KARL(JC),JC=IC10,IC20)
 100      FORMAT(I5,126A1)
 1      CONTINUE
 
        WRITE(NRES,103)
     >  (BLANC,II=1,IC10-1),(KARO(I-(I/10)*10+1),I=IC10-1,IC20-1)
 103    FORMAT(/,131A1)
        WRITE(NRES,107) (BLANC,II=1,(IC10/10)*10-1)
     >  ,( ( BLANC,II=(I-(I/10)*10),8 )
     >  ,KARO((I)/10) , I=(IC10/10)*10,(IC20/10-1 )*10,10)
 107    FORMAT(1X,131A1)
 
        WRITE(NRES,104) ICTOT(JJ,NH),JMAX(JJ,NH)
     >  ,MOYC(JJ,NH),NC(JJ,MOYC(JJ,NH),NH)
 104    FORMAT(
     >  //,15X,' TOTAL  COMPTAGE                 : ',I7,'  SUR',I7
     >  ,/,15X,' NUMERO   DU  CANAL  MOYEN       : ',I7
     >  ,/,15X,' COMPTAGE  AU   "      "         : ',I7 )
 
      ELSEIF(ICTOT(JJ,NH) .EQ. 0) THEN
 
        WRITE(NRES,108)
 108    FORMAT(////,15X,' TOTAL  COMPTAGE                 : 0',///)
        RETURN 1
 
      ENDIF
 
      RETURN
      END
