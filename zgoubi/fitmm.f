      SUBROUTINE FITMM(Y,T,Z,P,SAR,DP,TAR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "MAXCOO.H"
      INCLUDE "MXLD.H"
      DIMENSION FMI(MXJ,MXL), FMA(MXJ,MXL)

      DIMENSION FMIO(MXJ,MXL), FMAO(MXJ,MXL)
      SAVE FMIO, FMAO

      CALL ZGNOEL(
     >            NOEL)

      FMA(1,NOEL) = AMAX1(DP,FMA(1,NOEL))
      FMA(2,NOEL) = AMAX1(Y,FMA(2,NOEL))
      FMA(3,NOEL) = AMAX1(T,FMA(3,NOEL))
      FMA(4,NOEL) = AMAX1(Z,FMA(4,NOEL))
      FMA(5,NOEL) = AMAX1(P,FMA(5,NOEL))
      FMA(6,NOEL) = AMAX1(SAR,FMA(6,NOEL))
      FMA(7,NOEL) = AMAX1(TAR,FMA(7,NOEL))
      FMI(1,NOEL) = AMIN1(DP,FMI(1,NOEL))
      FMI(2,NOEL) = AMIN1(Y,FMI(2,NOEL))
      FMI(3,NOEL) = AMIN1(T,FMI(3,NOEL))
      FMI(4,NOEL) = AMIN1(Z,FMI(4,NOEL))
      FMI(5,NOEL) = AMIN1(P,FMI(5,NOEL))
      FMI(6,NOEL) = AMIN1(SAR,FMI(6,NOEL))
      FMI(7,NOEL) = AMIN1(TAR,FMI(7,NOEL))
          
      DO 1 JJ = 1, MXJ
          FMIO(JJ,NOEL) = FMI(JJ,NOEL) 
          FMAO(JJ,NOEL) = FMA(JJ,NOEL)
C         write(*,*) jj,noel,FMIO(JJ,NOEL),FMaO(JJ,NOEL)
 1    ENDDO

C        write(*,*) ' fitmm y,t,noel ',y,t,noel
C        write(*,*) ' fitmm y,t,noel ',FMAO(2,NOEL),FMAO(3,NOEL),noel

      RETURN

      ENTRY FITMM1(JI,LI,MIMA,
     >                        VAL)
      IF(MIMA.EQ.1) THEN
        VAL = FMIO(JI,LI)
      ELSEIF(MIMA.EQ.2) THEN
        VAL = FMAO(JI,LI)
      ENDIF
C        write(*,*) ' fitmm ji,li,mima,val ',ji,li,mima,val
C         stop
      RETURN

      ENTRY FITMM2
      DO 11 LL = 1, MXL
        DO 11 JJ = 1, MXJ
          FMI(JJ,LL) = 1.D99
          FMA(JJ,LL) = -1.D99
 11   ENDDO
      RETURN

      END      
