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
      LOGICAL FUNCTION ITSENS(FONC,N,X,XMIN,XMAX,P,V,F0,D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(*),XMIN(*),XMAX(*),P(*),V(*)
C----------------------------------------------------------
C                                                     
C   Soit F={Fj} 1<=j<=M / Fj=Fj(Xi) 1<=i<=N             
C                                                         
C   On cherche {Xi} / {Fj(Xi)}={Cj}                     
C              
C   Soit G=(Somme j(Pj*(Fj(Xi)-Cj)**2))=G(Fj)  
C                                              
C   On se ramene a la recherche du minimum de G
C                                              
C----------------------------------------------------------
C                                                          
C          Recherche du minimum de G par deplacement       
C                                                          
C----------------------------------------------------------
C                                                          
C  1) Determination du sens de deplacement                 
C                                                          
C     a) U={Xj} 1<=j<=N ; G0=G(U) ; k=0                    
C                                                          
C     b) i=1                                               
C                                                          
C     c) U'={Xj} /Xj = Xj +/- Pj * d(i-j), XjMin <= Xj <= XjMax ; G1=G(U')
C                                                              
C     d) si G1 < G0 => Vi = Pi ; G0=G1                         
C                                                              
C        sinon         Vi = 0  ; Xi = Xi - Pi                  
C                                                              
C        i=i+1 si i <=N => reprendre en 1c)                    
C                                                              
C        sinon  1e)                                            
C                                                              
C     e) k=k+1 si k < 4 => reprendre en 1f)                    
C                                                              
C        sinon arret "recherche de chemin impossible"          
C                                                              
C     f) si V = {0} => Pj = Pj/D 1<=j<N ; reprendre en 1b)
C                                                         
C        sinon  2)                                        
C                                                         
C  2) Avance dans la direction trouvee                    
C                                                         
C     a) U{Xj = Xj + Vj} / 1<=j<=N , XjMin <= Xj <= XjMax ; G1=G(U)
C                                                                  
C        si G1 < epsilon arret "minimum trouve"                    
C                                                                  
C        sinon si G1 < G0 => G0=G1 reprendre en 2a)                
C                                                                  
C              sinon Xj= Xj - Vj 1<=j<=N ; reprendre en 1)         
C                                                                  
C------------------------------------------------------------------
C                                                                  
C     SUBROUTINE ITMINO(FONC,N,X,)                              
C                                                                
C     FONC   :FONCTION A TRAITER                                 
C                                                                
C     N      :NOMBRES DE  VARIABLES                              
C     X      :VALEURS DES VARIABLES   ATTEINTES Xi    1<=i<=N    
C     XMIN   :MINIMUM DES VARIABLES             XMINi 1<=i<=N    
C     XMAX   :MAXIMUM DES VARIABLES             XMAXi 1<=i<=N    
C                                                                
C                                                                
C----------------------------------------------------------------
      LOGICAL OK
      PARAMETER (NITER=3)
      EXTERNAL FONC
 
      NI=1
      OK=.FALSE.
1     CONTINUE
      DO 2 I=1,N
         V(I)=0.0
         IF(P(I) .NE. 0) THEN
            X1=X(I)
            X2=P(I)
            X(I)=X1+X2
            IF(X(I).GT.XMAX(I)) THEN
               X(I)=XMAX(I)
               X2=0.0
            ENDIF
            CALL CPTFCT(FONC,F1)
            IF(F1.GE.F0) THEN
               X2=-P(I)
               X(I)=X1+X2
               IF(X(I).LT.XMIN(I)) THEN
                  X(I)=XMIN(I)
                  X2=0.0
               ENDIF
               CALL CPTFCT(FONC,F1)
               IF(F1.GE.F0) THEN
                  X(I)=X1
               ELSE
                  OK=.TRUE.
                  V(I)=X2
                  F0=F1
               ENDIF
            ELSE
               OK=.TRUE.
               V(I)=X2
               F0=F1
            ENDIF
         ENDIF
2     CONTINUE
      IF(.NOT.OK) THEN
         NI=NI+1
         IF(NI.LE.NITER) THEN
            DO 3 I=1,N
               P(I)=P(I)/D
3           CONTINUE
            GOTO 1
         ENDIF
      ENDIF
      ITSENS=OK
      RETURN
      END
