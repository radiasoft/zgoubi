C  D L G A U / ALGEBRE LINEAIRE RUBRIQUE L,L201A          MINI BBL MATH         
C  FORTRAN    
      SUBROUTINE DLGAU ( ID,NA,M,A,K,IER )                
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C             
C         RESOLUTION DE SYSTEMES LINEAIRES A COEFFICIENTS REELS                 
C             AVEC POSSIBILITE DE PLUSIEURS SECONDS MEMBRES                     
C                        PAR L'ALGORITHME DE GAUSS        
C             AVEC RECHERCHE DE PIVOT MAXIMUM ( PAR LIGNE )                     
C         -----------------------------------------------------                 
C   ID 1ERE DIMENSION DU BLOC A  AVEC ID >= NA            
C   NA ORDRE DU SYSTEME LINEAIRE    
C   M  NOMBRE DE SECONDS MEMBRES    
C   A  MATRICE ET SECONDS MEMBRES ( REMPLACES PAR LES SOLUTIONS )               
C   K  BLOC DE TRAVAIL DE DIMENSION NA ( K : VECTEUR DE PERMUTATION )           
C             
C   SI MATRICE A REGULIERE : IER INCHANGE                 
C   SI IER = -1  ET SI MATRICE A SINGULIERE : IER = 1 ET RETOUR                 
C   SI IER.NE.-1 ET SI MATRICE A SINGULIERE : STOP        
C   SI IER = -1  ET SI MATRICE A QUASI SINGULIERE : IER = 2 ET RETOUR           
C   SI IER.NE.-1 ET SI MATRICE A QUASI SINGULIERE : STOP  
C  ------------------------------------------- D L G A U ---------------        
      DIMENSION A(ID,1),K(1)        
      DATA EPS /1.D-15/             
      N=NA    
      NDEB=N+1
      NM=N+M  
C             
      DO 10 I=1,N                   
   10 K(I)=I  
C             
      DO 120 I=1,N                  
      AMXX= DABS(A(I,I))            
      JMAX=I  
      I1=I+1  
C             
      IF(I.GE.N) GO TO 30           
C             
      DO 20 J=I1,N                  
      AR= DABS(A(I,J))              
      IF(AMXX.GT.AR) GO TO 20       
      AMXX=AR 
      JMAX=J  
   20 CONTINUE
C             
   30 IF ( AMXX.GT.0.D0 ) GO TO 40  
      PRINT 1001                    
 1001 FORMAT ( //' DLGAU * MATRICE SINGULIERE ****')      
      IF ( IER.NE.-1 ) STOP         
      IER=1   
      RETURN  
   40 IF ( JMAX.LE.I ) GO TO 60     
C             
      DO 50 I2=1,N                  
      AUX=A(I2,I)                   
      A(I2,I)=A(I2,JMAX)            
      A(I2,JMAX)=AUX                
   50 CONTINUE
C             
   60 IF ( I.LE.1 ) GO TO 80        
C             
      S=0.D0  
      T=0.D0  
      IN=I-1  
      DO 70 IT=1,IN                 
      P=A(IT,I)*A(I,IT)             
      S=S+P   
   70 T=T+ DABS(P)                  
      ERA=   EPS*(T+ DABS(A(I,I)-S))
C             
      IF ( AMXX.GT.ERA ) GO TO 80   
      PRINT 1002                    
 1002 FORMAT ( //' DLGAU * MATRICE QUASI SINGULIERE ****')
      IF ( IER.NE.-1 ) STOP         
      IER=2   
      RETURN  
C             
   80 DO 90 J2=I1,NM                
   90 A(I,J2)=A(I,J2)/A(I,I)        
C             
      IF ( I.GE.N ) GO TO 110       
C             
      DO 100 I3=I1,N                
      DO 100 J3=I1,NM               
  100 A(I3,J3)=A(I3,J3)-A(I3,I)*A(I,J3)                   
C             
  110 IF(JMAX.LE.I) GO TO 120       
C             
      NAB=K(JMAX)                   
      K(JMAX)=K(I)                  
      K(I)=NAB
C             
  120 CONTINUE
      IF ( N.EQ.1 ) RETURN          
C             
      DO 150 KC=NDEB,NM             
      J=N     
  130 I=J-1   
  140 A(I,KC)=A(I,KC)-A(J,KC)*A(I,J)
      I=I-1   
      IF ( I.GT.0 ) GO TO 140       
      J=J-1   
      IF ( J.GT.1 ) GO TO 130       
  150 CONTINUE
C             
      DO 180 I=1,N                  
  160 J=K(I)  
      IF ( J.LE.I ) GO TO 180       
      K(I)=K(J)                     
      K(J)=J  
C             
      DO  170 MP=NDEB,NM            
      AUX=A(J,MP)                   
      A(J,MP)=A(I,MP)               
  170 A(I,MP)=AUX                   
      GO TO 160                     
C             
  180 CONTINUE
      RETURN  
      END       
