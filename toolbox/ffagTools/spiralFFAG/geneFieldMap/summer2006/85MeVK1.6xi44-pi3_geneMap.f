      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NA=100, NR=60)
      DIMENSION ATAB(NA), RTAB(NR), BZ(NA,NR)
      DIMENSION XB(NA,NR), YB(NA,NR)

      OPEN(UNIT=10,FILE='FIELDMAP',FORM='FORMATTED')

C     DATA DEB
      PI = 4.D0 * ATAN(1.D0)
      RAD2D=180.D0/PI
      D2RAD = 1.D0/RAD2D

C----------------------- Data hypothesis
      B0 = 15.D0        ! KG
      AK = 1.6D0
      XI=44.D0  * D2RAD
      R0 = 2.D0
      RMIN = R0 - .85D0 
      RMAX = 3.  !!R0 + 0.05D0 
      PF=0.454
      NCELL = 6 
      PRIOD = 2.D0 * PI / NCELL
      TTF = PRIOD * PF    ! sector angle of magnet
C total angle of field map
      AT =  1.5D0 * PRIOD  

c Positionning of field region inside map
      DAT = AT/5.D0

      G=0.06D0           ! gap size
      C0=0.1455D0        ! fringe coefficients
      C1=2.267D0
      C2=-.6395D0
      C3=1.1558D0
      C4=0.D0
      C5=0.D0    
C----------------------------------------

      XACC=0.0001D0     ! accuracy for newton zero method

      RO1=R0
      B1=1.D0/TAN(XI)
      RO2=R0
      B2=B1

C Positionning of spiral EFB #1 (equiv. OMEGA+ in zgoubi)
      TTA1= TTF  /2.D0
C Positionning of spiral EFB #2 (equiv. OMEGA- in zgoubi)
      TTA2=  -TTF  /2.D0

        chrfe =  -((AT - PRIOD )/2.D0 + DAT)
        chrfs =    (AT - PRIOD )/2.D0 - DAT
        write(6,*) ' RE ALE RS ALS =  0. ',chrfe,' 0. ',chrfs
      AMIN = -AT/2.D0 - DAT       !-AT/2. - DAT
      AMAX =  AT/2.D0 - DAT       ! AT/2. - DAT

      WRITE(6,*) ' Amin, Amax : ', amin,amax
      WRITE(6,*) ' Total sector angle : ',AT*RAD2D,' deg.'
      WRITE(6,*) ' Extnt of magnetic sector: ',(TTA1-TTA2)*RAD2D,' deg.'
      SETVAL = FFSPI1(TTA1,TTA2)

      ASTP =(AMAX-AMIN)/FLOAT(NA-1)
      RSTP = (RMAX-RMIN)/FLOAT(NR-1)

      WRITE(6,*) ' ASTP, RSTP : ',ASTP*D2RAD,RSTP,AMIN*RAD2D,AMAX*RAD2D
      CALL CPU_TIME(TIMSEC)

      IRSTP = 1
      IASTP = 1    
      DO 1 IR = 1, NR, IRSTP
        RTAB(IR) = RMIN + (IR-1.D0)*RSTP
 1    CONTINUE
      DO 2 IA = 1, NA, IASTP
        ATAB(IA) = AMIN + (IA-1.D0)*ASTP
 2    CONTINUE
      DO 3 IR = 1, NR, IRSTP
        RB = RTAB(IR) 
        DO 3 IA = 1, NA, IASTP
          XB(IA,IR) = RB * COS(ATAB(IA))
          YB(IA,IR) = RB * SIN(ATAB(IA))
 3    CONTINUE

      DO 4 IR = 1, NR, IRSTP

        WRITE(88,*) 
        WRITE(11,*) 

        DO 4 IA = 1, NA, IASTP

          BZ(IA,IR)=FFSPIF(ATAB(IA),RTAB(IR),XB(IA,IR),YB(IA,IR),B0,
     >       R0,AK,AMIN,AMAX,XACC,RO1,B1,RO2,B2,G,C0,C1,C2,C3,C4,C5,
     >                                                  DE,FE,DS,FS)
        WRITE(11,FMT='(1P, 3G12.4, 2I4,6G12.4)') ATAB(IA),RTAB(IR),
     >    BZ(IA,IR), IA,IR,DE,FE,DS,FS,XB(IA,IR),YB(IA,IR)
        
 4    CONTINUE

      TEMP = TIMSEC
      CALL CPU_TIME(TIMSEC)
      WRITE(   6,*) '  CPU time, total :  ',  TIMSEC-TEMP


        WRITE(10,*)  ' carte de champ aimant spiral'
        WRITE(10,*)  ' carte de champ aimant spiral'
        WRITE(10,*)  ' gggggggggggggggg   '
            ACN = 0.D0
            KART = 1
      WRITE(10,992) NA,NR,ACN,(RMAX+RMIN)/2.D0*100.,KART
 992  FORMAT(2I3,1P,2G15.7,1X,I1)
      WRITE(10,*) (ATAB(I),I=1,NA),(RTAB(J)*100.,J=1,NR)
     >  ,((BZ(I,J),I=1,NA),J=1,NR)

      STOP
      END

      FUNCTION FFSPIF(A,R,XB,YB,B0,R0,AK,AMIN,AMAX,XACC,
     >                      RO1,B1,RO2,B2,G,C0,C1,C2,C3,C4,C5,
     >                                                DE,FE,DS,FS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (PLIM=20.D0)
      SAVE TTA1,TTA2

C Rustine carte pi/2
         IF(R.GT.1.05D0*R0)  THEN
           FFSPIF = 0.D0
           goto 88
         ENDIF

C  DISTANCE TO FIRST EFB, AND FIELD FACTOR
         DE=-DSTEFB(XB,YB,RO1,B1,AMIN,AMAX,XACC,TTA1,
     >                                             YN)
         IF(YB.GT.YN) DE = -DE
         GE=DE/G         
C             write(*,*) ge
         P = (C0+GE*(C1+GE*(C2+GE*(C3+GE*(C4+GE*C5)))))
            IF    (P .GE.  PLIM) THEN
              FE = 0.D0
            ELSEIF(P .LE. -PLIM) THEN
              FE = 1.D0
            ELSE
              FE = 1.D0/(1.D0+EXP(P))
            ENDIF

C  DISTANCE TO FIRST EFB, AND FIELD FACTOR
         DS=DSTEFB(XB,YB,RO2,B2,AMIN,AMAX,XACC,TTA2,
     >                                             YN)
         IF(YB.GT.YN) DS = -DS
         GS=DS/G
C             write(*,*) gs
         P = (C0+GS*(C1+GS*(C2+GS*(C3+GS*(C4+GS*C5)))))
            IF    (P .GE.  PLIM) THEN
              FS = 0.D0
            ELSEIF(P .LE. -PLIM) THEN
              FS = 1.D0
            ELSE
              FS = 1.D0/(1.D0+EXP(P))
            ENDIF
            
         FMIN=1.D-60
         IF (FE.LT.FMIN) FE=0.D0
         IF (FS.LT.FMIN) FS=0.D0
       
         FFSPIF = B0 * (R/R0)**AK*FE*FS 

 88      WRITE(88,FMT='(1P, 6G14.6)') A*R, R, FFSPIF

      RETURN
      
      ENTRY FFSPI1(TTA1I,TTA2I)
      TTA1=TTA1I
      TTA2=TTA2I
      RETURN

      END

      FUNCTION DSTEFB(XB,YB,D,E,AMIN,AMAX,XACC,TTARF,
     >                                             YN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      TZ=RTNEWT(XB,YB,D,E,AMIN,AMAX,XACC,TTARF)
C      write(*,*) ' FCTN DSTEFB   TZ : ',tz
      XN = D*EXP(E*TZ)*COS(TZ+TTARF)
      YN = D*EXP(E*TZ)*SIN(TZ+TTARF)
      IF(( (XB - XN)**2 + (YB - YN)**2 ) .LT. 0.D0)
     >   STOP ' *** FCTN DSTEFB     Negative distance !!  '
      DSTEFB=SQRT( (XB - XN)**2 + (YB - YN)**2 )
      RETURN
      END
      
      FUNCTION RTNEWT(XB,YB,D,E,X1,X2,XACC,TTARF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (JMAX=10000)
      RTNEWT=.5D0*(X1+X2)
      DO 11 J=1,JMAX
         XTNEWT = RTNEWT
         CALL FUNCD(RTNEWT,XB,YB,D,E,F,DF,TTARF)
         DX=F/DF
C         write(*,*) ' FCTN RTNEWT   RTNEWT, DX, F, DF : ',RTNEWT,DX,F,DF
         RTNEWT=RTNEWT-DX
C         IF((X1-RTNEWT)*(RTNEWT-X2).LT.0.D0)PAUSE'jumped outof brackets'
         IF((X1-RTNEWT)*(RTNEWT-X2).LT.0.D0) THEN
           RTNEWT=XTNEWT
           RETURN
         ENDIF
         IF(ABS(DX).LT.XACC) RETURN 
 11      CONTINUE
      WRITE(6,*) ' ** WARNING :   rtnewt exceeding maximum iterations'
      RETURN
      END

      SUBROUTINE FUNCD(X,XB,YB,D,E,FN,DF,TTARF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      FN=FCTN(XB,YB,D,E,X,TTARF)
      DF=DFCTN(XB,YB,D,E,X,TTARF)
      RETURN
      END
      
      FUNCTION FCTN(XB,YB,R0,B,T,TTARF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE EBT
      EBT = EXP(B*T)
      COST=COS(TTARF+T)
      SINT=SIN(TTARF+T)
      FCTN= 
     - (XB - EBT*R0*COST)*(B*EBT*R0*COST - EBT*R0*SINT) + 
     - (YB - EBT*R0*SINT)*(EBT*R0*COST + B*EBT*R0*SINT)
C        WRITE(*,*) ' (FCTN) T, FCTN = ',T, FCTN,B,XB,YB,E,TTARF,R0
      RETURN
      ENTRY DFCTN(XB,YB,R0,B,T,TTARF)
       DFCTN= 
     - (B*EBT*R0*COST - EBT*R0*SINT)*(-(B*EBT*R0*COST) + EBT*R0*SINT) + 
     - (XB - EBT*R0*COST)*(-(EBT*R0*COST) + B**2*EBT*R0*COST- 
     - 2*B*EBT*R0*SINT) + (-(EBT*R0*COST) - B*EBT*R0*SINT)*
     - (EBT*R0*COST + B*EBT*R0*SINT) + (YB - EBT*R0*SINT)*
     - (2*B*EBT*R0*COST - EBT*R0*SINT + B**2*EBT*R0*SINT)
      RETURN
      END

     
