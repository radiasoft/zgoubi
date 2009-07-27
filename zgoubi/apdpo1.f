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
C  A P D P O 1 / RUBRIQUE P , P038                         MINI BBL MATH
C  A P D P O 1 / RUBRIQUE P , P038                         MINI BBL MATH
C  FORTRAN
      SUBROUTINE APDPO1(A,EPS,N,FP,I)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C ----------------------------------------------------------------------
C  CALCUL DES PROBABILITES DE LA LOI DE POISSON
C
C --------------------------------------------------- APDPO1 -----------
      DIMENSION FP(N)
C
      N0=N/2
      K0=INT(A)
      AK=AINT(A)
      IF(A-15)10,11,11
C ***CAS OU A EST INFERIEUR A 15
 10   P0=EXP(-A)
      FP(1)=P0
      IF(K0.EQ.0)GO TO 50
C **CALCUL DES FP DE FP(1) A FP(K0+1) ET CUMUL DE LA SOMME
      B=0D0
      SIGMA=FP(1)
      FPP=0D0
      DO 5 J=1,K0
        B=B+1.D0
        FP(J+1)=FP(J)*A/B
  5     SIGMA=SIGMA+FP(J+1)
C **PROB MAX EN P0=FP(K0+1) ET EPSP= LIMITE ABSOLUE IMPOSEE A FP .
      P0=FP(K0+1)
      EPSP=P0*EPS
      I=0
      JL=K0+1-N0
C **CUMUL EVENTUEL POUR LES PREMIERS INDICES,SELON EPS ET LA DIMENSION
C RESERVEE POUR LE TABLEAU .
  6   IF((FP(I+2).GE.EPSP).AND.(JL.LE.1))GO TO 7
      FPP=FPP+FP(I+2)
      JL=JL-1
      I=I+1
      GO TO 6
C **RENVOI AU CALCUL DES FP POUR I>K0+1,APRES DECALAGE EVENTUEL .
  7   JP=K0+1-I
      IF(I.EQ.0)GO TO 40
      FP(1)=FP(1)+FPP
      DO 8 JJ=2,JP
  8     FP(JJ)=FP(I+JJ)
      GO TO 40
C *CAS OU K0=0
 50   JP=1
      EPSP=EPS*P0
      I=0
      SIGMA=P0
      GO TO 40
C ***CAS OU A > OU = 15
C **CALCUL DE P0 PAR LA FORMULE DE STIRLING
 11   B=1.D0
      AKA=AK-A
      F=-AKA*AKA/(2D0*AK)
      DF=F
 13   B=B+1D0
      DF=DF*AKA*B/(AK*(B+1D0))
      F=F+DF
      IF((DF.LE.-1.D-8).OR.(DF.GE.1.D-8))GO TO 13
      G=2.506628D0*SQRT(AK)
      G=G*(1D0+(1D0+(1D0-(1D0+8.558153D-3/AK)*.7722222D0/AK)*
     $     4.166667D-2/AK)*8.333333D-2/AK)
      P0=EXP(F)/G
C **IND=1,NMIN=N0-1 SI N0<K0+1,SINON IND=2 ET NMIN=K0-1
      IND=1
      NMIN=N0-1
      IF(N0.GE.(K0+1))IND=2
      NMIN=NMIN+(IND-1)*(K0-1-NMIN)
      SIGMA=P0
      FP(N0)=P0
      EPSP=EPS*P0
C **CALCUL DES PROB. POUR I<N0.ARRET A FP(N0-I)<EPSP ET RENVOI A CUMUL
      DO 22 I=1,NMIN
        FP(N0-I)=FP(N0-I+1)*(K0-I+1)/A
        SIGMA=SIGMA+FP(N0-I)
 22     IF(FP(N0-I).LT.EPSP)GO TO 30
C **SORTIE NORMALE DE LA BOUCLE DE CALCUL AIGUILLAGE SELON N0:K0+1
      GO TO(20,21),IND
C *ARRET DU CALCUL DE FP A (N-1)<EPSP.INITIALISATION DU CUMUL .
 30   JLIM=K0-I+1
      IND=2
      JP=I+1
      FPP=FP(N0-I)
      FP(1)=FPP
      GO TO 31
C ***CAS N0> OU =(K0+1).DECALAGE NECESSAIRE APRES CALCUL DE FP(1)
C  RENVOI AU CALCUL POUR I>K0+1
 21   FP(1)=FP(N0-NMIN)/A
      JP=K0+1
 26   DO 23 JJ=2,JP
 23     FP(JJ)=FP(N0-JP+JJ)
      GO TO 40
C ***CAS N0<K0+1.INITIALISATION DU CUMUL
 20   JLIM=K0-N0+2
      FPP=FP(1)
      JP=NMIN+1
C ***CUMUL PROPREMENT DIT
 31   DO 24 JJ=1,JLIM
        FPP=FPP*(JLIM-JJ)/A
        SIGMA=SIGMA+FPP
        FP(1)=FP(1)+FPP
C AIGUILLAGE EN FIN DE CUMULSELON DECALAGE NECESSAIRE AVANT LE
C  PASSAGE AUX INDICES CROISSANTS .
 24     IF(FPP.LT.1.D-8)GO TO (40,26),IND
C *CALCUL DES FP POUR I>K0+1
 40   JLIM=N-JP-1
      DO 41 JJ=1,JLIM
        FP(JP+JJ)=FP(JJ+JP-1)*A/(K0+JJ)
        SIGMA=SIGMA+FP(JP+JJ)
        IF(FP(JJ+JP).LT.EPSP)GO TO 42
 41     CONTINUE
C *RETOUR APRES CALCUL DE N VALEURS DE SOMME 1.0 .
      FP(N)=1D0-SIGMA
      I=K0+1-JP
      RETURN
C *RETOUR AVEC UN TABLEAU DE MOINS DE N VALEURS(ON A TOUJOURS S=1.0)
 42   FP(JJ+JP)=1D0-SIGMA+FP(JJ+JP)
      N=JJ+JP
      I=K0+1-JP
      RETURN
      END
