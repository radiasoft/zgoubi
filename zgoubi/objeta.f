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
C  Upton, NY, 11973, USA
C  -------
      SUBROUTINE OBJETA
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     **********************************************
C     GENERE DES LEPTONS DE DECROISSANCE DU ETA
C     To , Po , et D  SONT ALEATOIRES ; Y0 = Z0 = 0.
C     EXAMPLE :
C                    BEAM    TARGET
C       MASSES (GEV): P   +    D   ->  3He  + ETA   &  ETA-> MU+  +  MU-
C                    AM1     ,AM2,     AM3,   AM4       ,    AM5  ,  AM6
C                  .93828  1.87563   2.80892 .5488,        .10566  .10566
C     **********************************************
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      CHARACTER(1) LET
      INCLUDE "C.FAISCT.H"     ! COMMON/FAISCT/ LET(MXT)
      CHARACTER(1) KAR(41)
      INCLUDE "C.KAR.H"     ! COMMON/KAR/ KAR
      INCLUDE "C.OBJET.H"     ! COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      INCLUDE "C.REBELO.H"   ! COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
 
      DIMENSION B(3),B4(3),PX1(5),PX2(5),P3(5),P4(5)
      DIMENSION FP(6,MXT)
 
      DIMENSION AMAX(5),CENTRE(5)
      CHARACTER  KTIR(2)*8 , BODY(3)*2

      SAVE STIR

      DATA KTIR /'Uniform', 'Gaussian'/
      DATA BODY /'M3' , 'M5' , 'M6' /

      AMQLU(1) = .FALSE.
      AMQLU(2) = .FALSE.
      AMQLU(3) = .FALSE.
      AMQLU(4) = .FALSE.
      AMQLU(5) = .FALSE.
      PABSLU = .FALSE.

C  .... MAGNETIC  RIGIDITY (KG*CM), MASSE (MeV/c/2)
C  .... Magnetic rigidity (Kg*cm), Mass (MeV/c/2)
      BORO = A(NOEL,1)
 
      IF(NRES.GT.0) WRITE(NRES,103) BORO
 103  FORMAT(25X,' MAGNETIC  RIGIDITY =',F15.3,' kG*cm')
C 103  FORMAT('1',15X,' RIGIDITE  MAGNETIQUE =',F15.3,' KG*CM')

      CALL REBELR(KREB3,KREB31,KDUM)
      IF(KREB3 .EQ. 99) THEN
C       ... SET TO 99 IN SBR REBELOTE
        IF(NRES.GT.0) WRITE(NRES,133)
 133    FORMAT(//,15X,
     >  'FINAL  COORDINATES  TAKEN  AS  INITIAL  COORDINATES')
        RETURN
      ENDIF
 
 
C     ... GENERE M3(IBODY=1) OU M5( PAR DECROISSANCE,IBODY=2 )
C         OU M6 ( = PARTNERS OF M5 PRVIOUSLY GENERATED; IBODY=3)
C         KOBJ=LOI DE Y,Z = 1(UNIF) OU 2(GAUSS)
      IBODY = A(NOEL,10)
      IF    (IBODY.LT.10) THEN
C       ... ANGULAR LIMITS T AND P =ZGOUBI ANGLES
        KAXE=1
      ELSE
C       ... ANGULAR LIMITS T AND P =POLAR
        IBODY=IBODY-10
        KAXE=2
      ENDIF
      KOBJ = A(NOEL,11)
      IMAX = A(NOEL,20)
      IDMAX=1
      IMAXT=IMAX/IDMAX
C     ... REST MASSES  (GEV/C2)
      AM1 = A(NOEL,30)
      AM2 = A(NOEL,31)
      AM3 = A(NOEL,32)
      AM4 = A(NOEL,33)
      AM5 = A(NOEL,34)
      AM6 = A(NOEL,35)
C     ... KINETIC-E (GEV)  INCIDENT  OF  BEAM AM1 (.896) :
      T1 = A(NOEL,40)
 
C   ...LECTURE CENTRE ET FRONTIERES DU TIRAGE
      CENTRE(2) = A(NOEL,50)
      CENTRE(3) = A(NOEL,51)
      CENTRE(4) = A(NOEL,52)
      CENTRE(5) = A(NOEL,53)
      CENTRE(1) = A(NOEL,54)
C     ** FRONTIERES  +/-   :
      AMAX(2) = A(NOEL,60)
      AMAX(3) = A(NOEL,61)
      AMAX(4) = A(NOEL,62)
      AMAX(5) = A(NOEL,63)
      AMAX(1) = A(NOEL,64)
 
      TMIN=CENTRE(3)-AMAX(3)
      TMAX=CENTRE(3)+AMAX(3)
      PMIN=CENTRE(5)-AMAX(5)
      PMAX=CENTRE(5)+AMAX(5)
      DMIN=CENTRE(1)-AMAX(1)
      DMAX=CENTRE(1)+AMAX(1)
 
C     ** DEMI-LONGUEUR DE LA CIBLE
      XL = A(NOEL,70)
      IF(IPASS .EQ. 1) THEN
        IRAND = A(NOEL,80)
        IRAND2= A(NOEL,81)
      ENDIF
 
      IF(NRES.GT.0) THEN
        WRITE(NRES,102) IMAX,BODY(IBODY),AM1,AM2,AM3,AM4,AM4,AM5,AM6
 102    FORMAT(15X,' GENERATION  DE ',I10,'  PARTICULES  DU  GENRE   '
     >  ,A,//,15X,' SELON  LES  REACTIONS   M1 + M2 --> M3 + M4 :'
     >  , /,20X,1P,G12.4,' + ',G12.4,'  -->  ',G12.4,' + ',G12.4
     >  ,   /,15X,' ET  M4  -->  M5 + M6 :'
     >  , /,20X,G12.4,'  -->  ',G12.4,' + ',G12.4,/)
        IF    (IBODY .NE. 3) THEN
          WRITE(NRES,104) T1
 104      FORMAT(15X,' M1  kinetic  energy :',1P,G12.4,' GeV',/)
 
          WRITE(NRES,100) (CENTRE(J),J=2,5),CENTRE(1),KTIR(KOBJ)
 100      FORMAT(15X,' CENTRE  DU  TIRAGE   : '
     >    ,/,11X,' YO, TO, ZO, PO,  BR/BORO : ',T50,1P,5G12.4
     >    ,/,11X,'( ',A,'  EN  Y  ET  Z )')
 
          WRITE(NRES,109) (AMAX(J),J=2,5),AMAX(1)
 109      FORMAT(15X,' FRONTIERES  DU  TIRAGE  (  +/- ) :'
     >    ,/,11X,' DY, DT, DZ, DP, DBR/BORO : ',T50,1P,5G12.4)
 
          WRITE(NRES,110) XL
 110      FORMAT(15X,' FRONTIERE  LONGITUDINALE  :'
     >    ,/,11X,' XL   =    +/-',F8.4,'  CM',/)
        ELSEIF(IBODY .EQ. 3) THEN
          WRITE(NRES,200)
 200      FORMAT(15X,' OBJET  CONSTITUE  DES  PARTNERS  ( M6 )',/)
        ENDIF
 
      ENDIF
 
      GOTO(1,1,2) IBODY
 
 1    CONTINUE
C     ++++ GENERE AM3(IBODY=1)
C          OU AM5 ( = LEPTON DE DECROISSANCE DE AM4; IBODY=2 )
 
      E1=AM1+T1
      ECM2=AM1**2+AM2**2+2.D0*AM2*E1
      ECM=SQRT(ECM2)
      E3CM=(ECM2+AM3**2-AM4**2)/2.D0/ECM
      E4CM=ECM-E3CM
      AK=SQRT(E3CM**2-AM3**2)
      PP1=SQRT(T1*(T1+2.D0*AM1))
      B(1)=PP1/(E1+AM2)
      B(2)=0D0
      B(3)=0D0
C
C     *** CONSTITUTION DU FAISCEAU
C
      IF(IPASS .EQ. 1) STIR=0.D0
      IKAR = 1
      I = 1
 
 11   CONTINUE
        STIR = STIR+1
 
        CALL DESBIN(AM4,AM5,AM6,PX1,PX2)
        CALL GENETA(AK,E3CM,E4CM,P3,P4,B4)
 
        IF(IBODY .EQ. 1) THEN
 
C         ... M3 IN LAB
          CALL BOOST(B,P3,3)
 
          T=GANG(1,P3,KAXE)
          IF( T .LT. TMIN .OR. T .GT. TMAX ) GOTO 11
 
          P=GANG(2,P3,KAXE)
          IF( P .LT. PMIN .OR. P .GT. PMAX ) GOTO 11
 
          DPP = P3(5)/(BORO*CL*1.D-12)
          IF( DPP .LT. DMIN .OR. DPP .GT. DMAX ) GOTO 11
 
        ELSEIF(IBODY .EQ. 2) THEN
 
C         ... M5
 
C         PX1 = 4-VECTEUR QdM DE M5 : PX1(1,2,3)=PX,Y,Z ,
C         PX1(4)=Etot=T+Masse ,  PX1(5)=P=SQRT(PX**2+PY**2+PZ**2) .
          CALL BOOST(B4,PX1,3)
          CALL BOOST(B,PX1,3)
 
          T=GANG(1,PX1,KAXE)
          IF( T .LT. TMIN .OR. T .GT. TMAX ) GOTO 11
 
          P=GANG(2,PX1,KAXE)
          IF( P .LT. PMIN .OR. P .GT. PMAX ) GOTO 11
 
          DPP = PX1(5)/(BORO*CL*1.D-12)
          IF( DPP .LT. DMIN .OR. DPP .GT. DMAX ) GOTO 11
        ENDIF
 
        FO(1,I) = DPP
        FO(3,I) = T
        FO(5,I) = P
        IF    (KOBJ .EQ. 1) THEN
C         ... TIRAGE Y , Z  Uniform
          DO 14 J=2,4,2
            FO(J,I)= 2.D0*(RNDM()-.5D0)*AMAX(J)
 14       CONTINUE
        ELSEIF(KOBJ .EQ. 2) THEN
C         ... TIRAGE Y , Z  Gaussian
          DO 15 J=2,4,2
            FO(J,I)=APHERF(CENTRE(J),AMAX(J))
 15       CONTINUE
        ENDIF
        IF ( XL .EQ. 0.D0) THEN
          FO(6,I) = 0D0
        ELSEIF ( XL .NE. 0.D0) THEN
C         ... CIBLE LONGUE : ON PROJETE LES COORDONNEES
C         DANS LE PLAN ORTHOGONAL A L'AXE LONGITUDINAL
C         AU CENTRE DE LA CIBLE
          AL = 2.D0*XL*( RNDM() - .5D0 )
          FO(2,I)=FO(2,I) - AL*TAN(FO(3,I)*1.D-3)
          FO(4,I)=FO(4,I) - AL*TAN(FO(5,I)*1.D-3)/COS(FO(3,I)*1.D-3)
C          FO(6,I)= - AL/(COS(FO(3,I)*1.D-3)*COS(FO(5,I)*1.D-3))
          FO(6,I)= AL/(COS(FO(3,I)*1.D-3)*COS(FO(5,I)*1.D-3))
        ENDIF
 
        IF(IBODY .EQ. 2) THEN
 
C         ... STORE M6 = PARTNER OF M5, FOR FURTHER USE OF IBODY=3
C        .... PX2 = 4-VECTEUR QdM DE M6
 
          CALL BOOST(B4,PX2,3)
          CALL BOOST(B,PX2,3)
 
          FP(3,I)=GANG(1,PX2,KAXE)
          FP(5,I)=GANG(2,PX2,KAXE)
          FP(1,I) = PX2(5)/(BORO*CL*1.D-12)
          FP(2,I)=FO(2,I)
          FP(4,I)=FO(4,I)
          FP(6,I)=FO(6,I)
        ENDIF
 
        IREP(I) = I
        IEX(I) = 1
        IF(IKAR.GT.41)  IKAR=1
        LET(I)=KAR(IKAR)
        IKAR = IKAR+1
        DO 12 J=1,6
          F(J,I)=FO(J,I)
 12     CONTINUE
 
        I = I+1
      IF(I .LE. IMAX) GOTO 11
 
 
      IF(NRES.GT.0) WRITE(NRES,105)STIR,IMAX*IPASS
     X ,IMAX*IPASS/STIR*1.D2
 105  FORMAT(/,10X,' Total number of sortings  = ',F15.0
     X      ,/,10X,' for ',I5,'  particles in  Window '
     X      ,1X,'( ',F7.2,' % )',/)
 
      GOTO 99
 
 2    CONTINUE
C     ... OBJECT = STORED M6 ( = PARTNERS OF NOT STOPPED M5 )
 
      II = 0
      DO 21 I=1,IMAX
C       ... IEX<-1 <=> PARTICULE STOPPEE
        IF(IEX(I) .LT. -1) GOTO 21
        II = II+1
 
C       .... SYM MIROIR <=> ROTATION DE 180 DEG / AXE X
        IF(KAXE .EQ. 1) THEN
          FP(3,I) = -FP(3,I)
        ELSE
          FP(5,I) = ASIN( SIN(.001D0*FP(5,I)) )*1.D3
        ENDIF
 
        DO  22 J=1,6
          FO(J,II) = FP(J,I)
          F(J,II) =  FP(J,I)
 22     CONTINUE
        IEX(II) = IEX(I)
        IREP(II) = II
        LET(II) = LET(I)
 21   CONTINUE
      IF(IPASS .EQ. 1) STIRP=0D0
      STIRP = STIRP+II
      IMAX0 = IMAX
      IMAX = II
      IMAXT=IMAX/IDMAX
 
      IF(NRES.GT.0) THEN
        WRITE(NRES,105) STIR , IMAX0*IPASS ,IMAX0*IPASS/STIR*1.D2
        WRITE(NRES,205) STIRP, STIRP/STIR*1.D2
 205    FORMAT(/,10X,' NOMBRE  TOTAL  DE  PARTNERS  = ',F15.0
     X    ,1X,'( ',F7.2,' % )',/)
      ENDIF
 
 99   CONTINUE
      IF (IBODY.EQ.1) THEN
         AM = AM3
      ELSE IF (IBODY.EQ.2) THEN
         AM = AM5
      ELSE IF (IBODY.EQ.3) THEN
         AM = AM6
      ELSE
         CALL ENDJOB('OBJETA: INVALID VALUE FOR IBODY: ',IBODY)
         AM = AAM
      ENDIF
      DO 993 I=1,IMAX
         AMQ(1,I) = AM
 993     AMQ(2,I) = Q
      AMQLU(1) = .TRUE.
      
      IF(IPASS.EQ.1) CALL CNTMXW(IMAX)
      RETURN
      END
