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
      SUBROUTINE CHAMB
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ****************************************
C     CHAMBRE DE LIMITES TRANSVERSALES YLIM2 ET ZLIM2.
C     ( ATTENTION: EN POLAIRES  Y EST UN RAYON...)
C        LE COMPTAGE EST EFFECTUE A CHAQUE PAS D'INTEGRATION, ENTRE
C     LE 1-ER ET LE 2-EME APPEL A CHAMBR. CETTE CHAMBR A DONC UNE
C     ETENDUE longitudinale NON NULLE, EGALE A LA DIFFERENCE 
C            D'ABSCISSES ENTRE LES
C     2 POINTS DE LA STRUCTURE OU ONT LIEU LES APPELS a CHAMBR.
C        LE TEST DE SORTIE A LIEU DANS LE SPGM COFIN.
C        LE 1-ER APPEL A CHAMBR AVEC LIMIT=1 DECLENCHE LE COMPTEUR,
C     LES APPELS SUIVANTS AVEC LIMIT=1 PERMETTENT DE MODIFIER LES LI-
C     MITES DE LA CHAMBRE SANS ARRETER LE COMPTAGE. UN APPEL AVEC LIMITE
C     =2 LISTE LE BILAN DU COMPTAGE ET RAZ LE COMPTEUR.
C     ***************************************
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      COMMON/CHAMBR/ LIMIT,IFORM,YLIM2,ZLIM2,SORT(MXT),FMAG,BMAX
     > ,YC,ZC
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      COMMON/REBELO/ NPASS,IPASS,KWRT,NNDES,STDVM
      COMMON/UNITS/ UNIT(MXJ)
 
      DIMENSION T(MXT)
      SAVE N1

C     .... LIMIT = 0 : ELEMENT INACTIF
C          LIMIT=1 : DEFINIT LA GEOMETRIE DE LA CHAMBRE, ET
C     DECLENCHE (OU POURSUIT) LE COMPTAGE.
C          LIMIT=2 : LIST DES TRAJECTOIRES COMPTEES HORS LIMITES
C     DE LA CHAMBRE, RAZ DU COMPTEUR, ET ARRET DU COMPTAGE.
C     IFORM = 1 : CHAMBRE RECTANGULAIRE
C     IFORM = 2 : CHAMBRE ELLIPTIQUE
C
      LIMIT = A(NOEL,1)
      IFORM = A(NOEL,10)
      YL = A(NOEL,11)
      ZL = A(NOEL,12)
      YC = A(NOEL,13)
      ZC = A(NOEL,14)
 
      YLIM2=YL*YL
      ZLIM2=ZL*ZL
 
 
        IF( LIMIT .EQ. 1 ) THEN
          CALL CNTOUR(
     >                N1)
          IF(NRES .GT. 0) THEN
           IF    (IFORM .EQ. 1) THEN
             WRITE(NRES,101)
 101         FORMAT(
     >       /,20X,' Start  of  rectangular  chamber  with  size :')
             WRITE(NRES,103) YL*UNIT(1),ZL*UNIT(3)
     >                      ,YC*UNIT(1),ZC*UNIT(3)
 103         FORMAT(25X,' YL = +/-',1P,G12.4,' m'
     >            ,5X,' ZL = +/-',   G12.4,' m'
     >       ,/,20X,' centered  on :',/
     >           ,25X,' YC =    ',   G12.4,' m'
     >            ,5X,' ZC =    ',   G12.4,' m')
             WRITE(NRES,107)  (YC-YL)*UNIT(1),(YC+YL)*UNIT(1), 
     >                        (ZC-ZL)*UNIT(3),(ZC+ZL)*UNIT(3)
 107         FORMAT(/,20X,
     >       ' => Max.  accepted  coordinates  Y,  Z  such  that :',1P,
     >       /,44X,G12.4,' < Y < ',G12.4,'  m',
     >       /,44X,G12.4,' < Z < ',G12.4,'  m')
           ELSEIF(IFORM .EQ. 2) THEN
             WRITE(NRES,102)
 102         FORMAT(/,20X,
     >               ' Start  of  elliptical  chamber  with  axes :')
             WRITE(NRES,103) YL*UNIT(1),ZL*UNIT(3),
     >                       YC*UNIT(1),ZC*UNIT(3)
C             WRITE(NRES,108) YL*UNIT(1),ZL*UNIT(3),
C     >                       YC*UNIT(1),ZC*UNIT(3)
             WRITE(NRES,108) YC*UNIT(1),YL*UNIT(1)
     >                      ,ZC*UNIT(3),ZL*UNIT(3)
 108         FORMAT(/,20X,
     >       ' => Max.  accepted  coordinates  Y,  Z  such  that :',1P,
     >       /,25X,'((Y-',G10.3,')/',G10.3,')^2 + ((Z-',G10.3,')/',
     >       G10.3,')^2 < 1')
           ENDIF
          ENDIF
 
        ELSEIF( LIMIT .EQ. 2 ) THEN
 
         CALL CNTOUR(
     >               NOUT)
         CALL CNTMXR(
     >               IMX) 
         CALL CNTSTO(
     >               NTOT) 
         IF(NRES .GT. 0) THEN
          WRITE(NRES,106)
 106      FORMAT(/,15X,' End  of  CHAMBRE  option ')
          WRITE(NRES,109) NOUT-N1, NOUT, IMX-NTOT, IMX, 
     >                            (100.D0*(IMX-NTOT))/IMX
 109      FORMAT(/,T20,' Number  of  particles  out  of  acceptance  '
     >    ,/,T25,' -  in  the  sense  of  CHAMBR  since  last  call : ',
     >                                                                I6
     >    ,/,T25,' -  by  collimations,  from  the  beginning : ',I6,
     >    /,T25,' Overall  survival  :  ',I9,' / ',I9,'  (',G10.3,'%)')

          IF(IPASS .EQ. 1) THEN
 
C------------ TRI DES VALEURS P0MAX(T0)
            IF(IFORM .EQ. 10) THEN
              WRITE(NRES,105)
105           FORMAT(//,10X,10('+')
     >        ,' Max.  P0  values  accepted  as  a  function  of  T0 :'
     >        ,10('+'),/)
              DO 6 K=1,MXT
6               T(K)=1.D10
              KTRA=1
              DO 3 ITRA=1,IMAX
                IF( IEX(ITRA) .EQ. -4) THEN
                  DO 5 K=1,KTRA
                     IF(FO(3,ITRA) .EQ. T(K)) GOTO 3
5                 CONTINUE
                  T(KTRA)=FO(3,ITRA)
                  P=0D0
                  DO 4 JTRA=1,IMAX
                    IF(IEX(JTRA) .NE. -4) THEN
C                     ** PARTICULE NON STOPPEES PAR CHAMBRE
                      IF(FO(3,JTRA) .EQ. T(KTRA)) THEN
                        IF(ABS(FO(5,JTRA)) .GT. ABS(P)) P=FO(5,JTRA)
                      ENDIF
                    ENDIF
 4                CONTINUE
                  WRITE(NRES,104) T(KTRA)*UNIT(2) , P*UNIT(4)
 104              FORMAT(20X,'  Max-P0 ( T0 =',1P,G11.3,' )   = ',
     >                                                 G11.3,'  rd ')
                  KTRA=KTRA+1
                ENDIF
 3            CONTINUE
            ENDIF
C
          ENDIF
         ENDIF
        ELSEIF(LIMIT .EQ. 0) THEN
         IF(NRES .GT. 0) THEN
          WRITE(NRES,100)
 100      FORMAT(/,20X,' +++++++  ELEMENT  CHAMBR  INACTIF  +++++++',/)
          RETURN
         ENDIF
        ENDIF
 
      RETURN
      END
