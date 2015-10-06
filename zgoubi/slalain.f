C  A L A I N / ALGEBRE LINEAIRE RUBRIQUE L,L301A          MINI BBL MATH         
C  FORTRAN                                                                      
      SUBROUTINE ALAIN ( ID,NZ,Z,N0,B,IER )                                     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C       ------ INVERSION D'UNE MATRICE REELLE Z D ORDRE NZ ------               
C  ID    1ERE DIMENSION DU BLOC BIDIMENSIONNE Z AVEC ID >= NZ                   
C  NZ    ORDRE DE LA MATRICE Z                                                  
C  Z     MATRICE A INVERSER, REMPLACEE PAR SON INVERSE Z**-1                    
C  N0,B  BLOCS DE TRAVAIL DE DIMENSION SUPERIEURE OU EGALE A NZ                 
C  IER   ENTIER PERMETTANT D ARRETER OU NON L EXECUTION                         
C        DU PROGRAMME APPELANT EN CAS DE MATRICE SINGULIERE                     
C                                                                               
C  LES ELEMENTS DE L'INVERSE INFERIEURS A 1.E-16 SONT ANNULES.                  
C  LES CALCULS SUR Z SONT FAITS AVEC 2 INDICES                                  
C  SOUS-PROGRAMME APPELE : ABS                                                  
C  ----------------------------------------- A L A I N -----------------        
      DIMENSION Z(ID,1),N0(1),B(1)                                              
      N=NZ                                                                      
      IF(N.EQ.1) GO TO 170                                                      
C                                                                               
      EPS=1.E-16                                                                
      DO 10 I=1,N                                                               
  10  N0(I)=I                                                                   
      N1=N-1                                                                    
      NP1 = N + 1                                                               
      DO 120 L=1,N                                                              
      I = NP1 - L                                                               
      J0=1                                                                      
      B1=ABS(Z(1,1))                                                            
C                                                                               
      IF(L-N) 20,40,160                                                         
   20 DO 30 J=2,I                                                               
      AZ = ABS( Z(1,J) )                                                        
      IF ( B1.GE.AZ ) GO TO 30                                                  
      B1 = AZ                                                                   
      J0=J                                                                      
   30 CONTINUE                                                                  
C                                                                               
   40 IF(B1) 160,150,50                                                         
C                                                                               
   50 IF(J0-1)160,80,60                                                         
   60 I=L+J0-1                                                                  
      J=N0(L)                                                                   
      N0(L)=N0(I)                                                               
      N0(I)=J                                                                   
      DO 70 I=1,N                                                               
      B1=Z(I,1)                                                                 
      Z(I,1)=Z(I,J0)                                                            
   70 Z(I,J0)=B1                                                                
C                                                                               
   80 DO 90 J=1,N1                                                              
   90 B(J)=Z(1,J+1)/Z(1,1)                                                      
      B(N)=1.0/Z(1,1)                                                           
      DO 110 I=1,N1                                                             
      IP1 = I + 1                                                               
      DO 100 J=1,N1                                                             
      Z(I,J) = Z(IP1,J+1) - Z(IP1,1) * B(J)                                     
      IF ( ABS( Z(I,J) ).LT.EPS ) Z(I,J) = 0.                                   
  100 CONTINUE                                                                  
      Z(I,N) = -Z(IP1,1) * B(N)                                                 
      IF ( ABS( Z(I,N) ).LT.EPS ) Z(I,N) = 0.                                   
  110 CONTINUE                                                                  
C                                                                               
      DO 120 J=1,N                                                              
  120 Z(N,J)=B(J)                                                               
C                                                                               
      DO 140 L=1,N1                                                             
      DO 140 I=L,N                                                              
      IF ( N0(I).NE.L ) GO TO 140                                               
      DO 130 J = 1,N                                                            
      B1=Z(I,J)                                                                 
      Z(I,J)=Z(L ,J)                                                            
  130 Z(L ,J)=B1                                                                
C                                                                               
      J=N0(I)                                                                   
      N0(I)=N0(L)                                                               
      N0(L)=J                                                                   
  140 CONTINUE                                                                  
      RETURN                                                                    
C                                                                               
  150 PRINT 1001                                                                
 1001 FORMAT (//' ALAIN * MATRICE SINGULIERE ****')                             
      IF ( IER.NE.-1 ) STOP                                                     
      IER = 1                                                                   
  160 RETURN                                                                    
C                                                                               
  170 IF(Z(1,1).EQ.0.) GO TO 150                                                
      Z(1,1)=1./Z(1,1)                                                          
      RETURN                                                                    
      END                                                                       
