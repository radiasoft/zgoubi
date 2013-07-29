C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
C From b_zgoubi.fai, will : 
C  - Get ellipse parameters for given particle, as matched over NPTMAX 
C    points evenly spread over turn #kpa to kpa+ksmpl-1
C  - Get particle coordinates contained in a elementary surface  
C    (epsilon+/-d_epsilon)*(phi+/-d_phi), and their spin coordinates, 
C    from turn #kpa to kpa + ksmpl -1, and store that in stroboscopy.out. 

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL OKOPN, CHANGE
C----- PLOT SPECTRUM     
      COMMON/CDF/ IES,IORDRE,LCHA,LIST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN
      DIMENSION BORNE(3)

      LOGICAL IDLUNI

      CHARACTER*80 NOMFIC
      logical exs
      character HV*2

      INCLUDE 'MAXNPT.H'

      data HV / '  ' /
      data borne / 0.05, 0.d0, 0.00872664625997  /
c Fourier transf. ksmpl turns strarting from kpa
      data kpa, ksmpl / 1, 99999999 /
      data nt / -1 /

      pi = 4.d0*atan(1.d0)

      write(*,*) ' '
      write(*,*) '----------------------------------------------------'
      write(*,*) 'Now running pgm stroboscopy... '
      write(*,*) '----------------------------------------------------'
      write(*,*) ' '
 
C Range of turns to be considered may be specified using stroboscopy.in
      INQUIRE(FILE='stroboscopy.in',exist=EXS)

        IF (IDLUNI(NLU)) THEN
          open(unit=nlu,file='stroboscopy.in')
        ELSE
          stop 'Pgm stroboscopy :   No idle unit number ! '
        ENDIF

      if(.NOT.exs) then
        write(*,*)'WARNING : File stroboscopy.in does not exist'
        write(*,*)'Pgm stroboscopy creates one from default values'

        write(nlu,fmt='(i6,t60,a)')  nt,' !  particl #'
        write(nlu,fmt='(2(i9,1x),t60,a)')  kpa, ksmpl
     >  ,' ! kpa, ksmpl : stroboscop ksmpl turns starting from kpa'
        write(nlu,fmt='(3(f8.4,1x),t60,a)') (borne(i),i=1,3)
     >  ,' ! dS loc and size : depsilon/epsilon, phi, dphi'
      endif

      rewind(nlu)

        read(nlu,*,err=11,end=11) NT                 ! particle #
        read(nlu,*,err=11,end=11) kpa, ksmpl         ! turn #kpa to kpa+ksmpl-1
        read(nlu,*,err=11,end=11) (borne(i),i=1,3)   ! surface element in phase-space : d_epsln, phi, d_phi

      close(nlu)
        write(*,*) ' Read following data from stroboscopy.in :'
 
      kpb = kpa + ksmpl -1
      kpc = ksmpl/nptmax
      if(kpc.eq.0) kpc = 1
      write(*,*) ' Particle # : ','   all ' 
      write(*,*) 
     >' Turn # range :  ',ksmpl,' turns,  from ',kpa,' to ',kpb
      write(*,*) ' Boundaries of surface lmnt (d_epsilon/phi/d_phi) : '
     >             ,(borne(i),i=1,3)
      call READC2B(KPa,KPb,KPc)

      call INIGR(
     >           NLOG, LM, NOMFIC)
      nl = nfai
      okopn = .false.
      change = .true.
C Get parameters of particle's lips, as matched over NPTMAX points evenly spread over turn #kpa to kpa+ksmpl-1
      call LPSPRM(NLOG,NL,LM,OKOPN,CHANGE,HV,NT,kpa,kpb
     >                              ,YM,YPM,YMX,YPMX,U,A,B)
C Get coordinates of particle and its spin within +/-dU/U/2 at U, and +/-d-phi/2 at phi, and store that in stroboscope.out

        write(*,*)  ' YM,YPM,YMX,YPMX,U,A,B : '
     >                              ,YM,YPM,YMX,YPMX,U,A,B


      dUU = 0.2d0
      phi = 0.d0
      dphi = 2.d0*pi / 720.d0
C      call STROB(borne,YM,YPM,YMX,YPMX,U,A,B,dUU,phi,dphi)
      
      stop ' Pgm "stroboscopy" ended correctly '

 11   continue
      stop 'Error during read from stroboscopy.in'
      end

      SUBROUTINE LPSPRM(NLOG,NL,LM,OKOPN,CHANGE,HV,NT,kpa,kpb
     >                              ,YM,YPM,YMX,YPMX,U,A,B)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL OKOPN, CHANGE
      character*(*) HV
C----- PLOT SPECTRUM     
      COMMON/CDF/ IES,IORDRE,LCHA,LIST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN
      INCLUDE 'MAXNPT.H'
      COMMON/TRACKM/COOR(NPTMAX,9),NPTS,NPTR

      DIMENSION YM(3), YPM(3), U(3), A(3), B(3), YNU(3)
      DIMENSION YMX(6), YPMX(6)
 
      LOGICAL OKECH
      CHARACTER REP, NOMFIC*80
      LOGICAL BINARY, BINARF
      CHARACTER HVL(3)*12

      LOGICAL OPN
      DATA OPN / .FALSE. /
      LOGICAL IDLUNI, OKKT5

      INCLUDE 'FILFAI.H'

      SAVE NPASS

      DATA HVL / 'Horizontal', 'Vertical', 'Longitudinal' /
    
      IF(.NOT.OKOPN) 
     > CALL OPNDEF(NFAI,FILFAI,NL,
     >                            NOMFIC,OKOPN) 

      NPTR=NPTS

      OKECH = .FALSE.

      IF(NT.EQ.-1) THEN
        CALL READC6B(1,NPTS)
      ELSE
        CALL READC6B(NT,NT)        
      ENDIF

 6    CONTINUE
        IF(.NOT. OPN) THEN
         IF(NT.EQ.-1) THEN
          IF (IDLUNI(IUN)) THEN
           call system('cat stroboscopy.out >>stroboscopy.out_old')
           call system('rm -f stroboscopy.out')
           OPEN(UNIT=IUN,FILE='stroboscopy.out',ERR=699)
           OPN = .TRUE.
           WRITE(IUN,*) '//////////////////////////////////////////'
           WRITE(IUN,*) 
     >     '% Stroboscopy considers turn#',kpa,' to turn# ',kpb
           WRITE(IUN,*) 
     >     '% XN, XPN, X, XP, Y, YP, L, DP, IEX, KT '
          ELSE
            GOTO 698
          ENDIF
         ENDIF
        ENDIF

          IF(NT.EQ.-1) THEN
            KPR = 2
            KT = 1
            CALL READC6B(KT,KT)
          ELSE
            KPR = 1
          ENDIF
         
C--------- Coordinate reading/storing loop
 62       CONTINUE
          IF(CHANGE) THEN
            NPTR = NPTMAX
            NPTS=NPTR
            I1 = 1
            CALL STORCO(NL,LM,I1  ,BINARY, 
     >                                     NPASS)
            CHANGE=.FALSE.
            IF(NPTR.GT.0) THEN
              IF(NPTS.GT. NPTR) NPTS=NPTR
            ELSE
              NPTR = NPTS
              IF(NT.EQ.-1) KT = 1 
             GOTO 69
            ENDIF
          ENDIF
          IF(NPTR .GT. 0) THEN
            CALL LPSFIT(NLOG,KPR,LM,
     >                              YM,YPM,YMX,YPMX,U,A,B,*60,*60)
 60         CONTINUE

          ENDIF

          IF(NT.EQ.-1) THEN
            KT = KT+1
            CALL READC6B(KT,KT)
            CALL READC5(
     >                  KT1,KT2)
            CHANGE = .TRUE.
            IF(KT.LE.KT2) GOTO 62             !--------- Coordinate reading/storing loop
          ENDIF

 69       CONTINUE
          IF(OPN) THEN
            CALL FLUSH2(IUN,.FALSE.)
            CLOSE(IUN)
            CHANGE = .TRUE.
          ENDIF
      GOTO 99


 698  WRITE(6,*) ' *** Problem : No idle unit for stroboscopy.out '
      stop ' *** Problem : No idle unit for stroboscopy.out '
c      GOTO 99
 699  WRITE(6,*) ' *** Problem at OPEN stroboscopy.out '
      stop ' *** Problem : No idle unit for stroboscopy.out '
c      GOTO 99

 99   RETURN
      END

      SUBROUTINE LPSFIT(NLOG,KPR,LM,
     >                              YM,YPM,YMX,YPMX,U,A,B,*,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION YM(*), YPM(*), U(*), A(*), B(*)
      DIMENSION YMX(*), YPMX(*)
      INCLUDE 'MAXNPT.H'
      COMMON/TRACKM/COOR(NPTMAX,9),NPTS,NPTR

      DIMENSION G(3)
      CHARACTER TXT(3)*5, REP*1

      INCLUDE 'MAXCOO.H'

      DIMENSION UM(3), UMI(3), UMA(3), SMEAR(3)
      CHARACTER*1 KLET

      PARAMETER ( PI=3.1415926536, SQ2 = 1.414213562)
C      PARAMETER ( PI=4.D0*ATAN(1.D0), SQ2 = SQRT(2.D0) )

      DATA TXT/ 'HORI.', 'VERT.', 'LONG.'/

      DO 7 J=1,MXJ-1,2
        YMX(J) = 1.D10
        YMX(J+1) = -1.D10
        YPMX(J) = 1.D10
        YPMX(J+1) = -1.D10
 7    CONTINUE

      DO 2 J=1,MXJ-1,2
C------- JJ = 1 , 2 or 3  for  Y-T, Z-P or T-P(time-momentum) planes
        J1 = J
        J2 = J+1
        JJ= J2 / 2
        XM=0.D0
        XPM=0.D0
        SNPT = 0.D0
        DO 21 I=1,NPTS
            X = COOR(I,J1)
            XP = COOR(I,J2)
            XM = XM + X
            XPM = XPM + XP
c              write(*,fmt='(a,2g12.4,3i4)') 'blibli ',x,xp,i,npts,jj
C----------- Min-max of the distribution
            IF(YMX(J1) .GT. X) YMX(J1) = X
            IF(YMX(J2) .LT. X) YMX(J2) = X
            IF(YPMX(J1) .GT. XP) YPMX(J1) = XP
            IF(YPMX(J2) .LT. XP) YPMX(J2) = XP
            SNPT = SNPT + 1
 21     CONTINUE
        XM = XM/SNPT
        XPM = XPM/SNPT
        YM(JJ)=XM
        YPM(JJ)=XPM
 2    CONTINUE

      DO 25 J=1,MXJ-1,2
C------- JJ = 1 , 2 or 3  for  Y-T, Z-P or T-P(time-momentum) planes
        J1 = J
        J2 = J+1
        JJ= (J+1) / 2
        X2=0.D0
        XP2=0.D0
        XXP=0.D0
        SNPT = 0.D0
        DO 26 I=1,NPTS
            X = COOR(I,J1)
            XP = COOR(I,J2)
            X2 = X2 + (X-YM(JJ))**2
            XP2 = XP2 + (XP-YPM(JJ))**2
            XXP = XXP + (X-YM(JJ))*(XP-YPM(JJ))
            SNPT = SNPT + 1
c              write(*,fmt='(a,3g12.4,i4)') 'blublu ',x,x2,ym(jj),jj
 26     CONTINUE
        X2  = X2/SNPT
        XP2 = XP2/SNPT
        XXP = XXP/SNPT

C G. Leleux : surface de l'ellipse S=4.pi.sqrt(DELTA)
C Soit d11=X2/sqrt(DELTA), d12=XXP/sqrt(DELTA), d22=XP2/sqrt(DELTA), alors 
C  d22.x^2-2.d12.x.x'+d11.x'^2=S/pi=4sqrt(DELTA), ce qui permet d'ecrire 
C   gamma=d22=XP2/sqrt(DELTA), -alpha=d12=XXP/sqrt(DELTA), beta=d11=X2/sqrt(DELTA). 
C En outre, par definition des dij, 
C     2.sigma_x=sqrt(d11.S/pi),  2.sigma_x'=sqrt(d22.S/pi). 
C En outre, frontiere : 
C          <x^2>_frontiere=2.(sigma_x)^2,    <x'^2>_frontiere=2.(sigma_x')^2

        SQ = SQRT(X2*XP2-XXP*XXP) 
        IF(SQ .GT. 0.D0) THEN
          B(JJ)=  X2/SQ
C Error  FM 03/02
C          A(JJ)=  XXP/SQ
          A(JJ)=  -XXP/SQ
          G(JJ)=  XP2/SQ
        ENDIF
C------- Courant invariant at 1 sigma is U=4.sqrt(DELTA)=Eps/pi (consistant with zgoubi !!) :
C Eps=ellipse surface
C        U(JJ) = 4.D0*SQ
        U(JJ) = SQ

c        write(*,fmt='(a,3g12.4,2i6)') 'blabla ',x2,xp2,u(jj),npts,jj

 25   CONTINUE

      IF(KPR .EQ. 0) RETURN
     
C----- SMEAR
      DO 3 J=1,MXJ-1,2
C------- JJ = 1 , 2 or 3  for  Y-T, Z-P or T-P(time-momentum) planes
        J1 = J
        J2 = J+1
        JJ= (J+1) / 2
        UMA(JJ) = -1.D10
        UMI(JJ) = 1.D10
        UM(JJ)=0.D0
        U2M=0.D0
        SNPT = 0.D0
        DO 31 I=1,NPTS
C--------- Normalized coordinates (*Beta), for phase-space point I:
            X = ( COOR(I,J1) - YM(JJ) )
            XP = ( COOR(I,J2) - YPM(JJ) )
            XN = X 
            XPN = ( A(JJ) * XN + B(JJ) * XP )
C----------- Courant invariant Epsilon/pi at phase-space point I:
            UI = ( XN * XN + XPN * XPN )/ B(JJ)

            IF(UI .GT. UMA(JJ)) UMA(JJ) = UI
            IF(UI .LT. UMI(JJ)) UMI(JJ) = UI

            UM(JJ) = UM(JJ) + UI
            U2M = U2M + UI * UI
            SNPT = SNPT + 1
 31     CONTINUE

        UM(JJ) = UM(JJ)/SNPT
        U2M = U2M/SNPT
        SMEAR(JJ) = SQRT( U2M - UM(JJ) * UM(JJ) )

 3    CONTINUE

C----- Twiss parameters and emittance
      CALL READC5(NT1,NT2)
      CALL READC9(KKEX,KLET)
      WRITE(6,FMT='(/,'' Particle # '',I6,'' to '',I6,''    @ lmnt # '',
     >I5,''. '',I6,'' points,  IEX='',I2)') NT1,NT2, LM, NINT(SNPT),KKEX
      WRITE(6,100) 
 100  FORMAT(/,10X,'Ellipse parameters and center  (MKSA) :',//,T12,
     >'BETA',T24,'ALPHA',T36,'GAMMA',T48,'EPS/PI',T63,'CENTER')

      DO 6 JJ=1,3
        I=2*JJ-1
C        IF(U(JJ).LE. 1.D-30) THEN
        IF(U(JJ).EQ. 0.D0) THEN
          WRITE(6,121) TXT(JJ),U(JJ)
 121      FORMAT(A,T20,'*** UNDETERMINED ***',T49,1P,G12.4)
        ELSE
          WRITE(6,120) TXT(JJ),B(JJ),A(JJ),G(JJ),U(JJ),YM(JJ),YPM(JJ)
 120      FORMAT(A,T10,1P,G12.4,T22,G12.4,T34,G12.4, 
     >    T46,G12.4,T58,2G12.4)
        ENDIF
 6    CONTINUE

      WRITE(6,FMT='(/,A)') ' Frontier values : '
      DO 66 JJ=1,3
        I=2*JJ-1
        FAC = 1.D0
        IF(U(JJ).EQ. 0.D0) THEN
          WRITE(6,121) TXT(JJ),U(JJ)
        ELSE
          WRITE(6,120) TXT(JJ),B(JJ),A(JJ),G(JJ),U(JJ),YM(JJ),YPM(JJ)
        ENDIF
 66   CONTINUE

C----- Smear
      WRITE(6,130) 
 130  FORMAT(/,10X,'Smear :',/,T17,
     >'<E/pi>',T26,'Sigma(E/pi)/<E/pi>',T48,'Min(E/pi)',T60,'Max(E/pi)')

      DO 5 JJ=1,3
        I=2*JJ-1
C        IF(U(JJ).LE. 1.D-30) THEN
        IF(U(JJ).EQ. 0.D0) THEN
          WRITE(6,131) TXT(JJ),U(JJ)
 131      FORMAT(A,T13,1P,G12.4,'      *** UNDETERMINED ***')
        ELSE
          WRITE(6,132) TXT(JJ),U(JJ),SMEAR(JJ)/U(JJ),UMI(JJ),UMA(JJ)
 132      FORMAT(A,T13,1P,G12.4,T28,G12.4,T45,G12.4,T57,G12.4)
        ENDIF
 5    CONTINUE

      IF(KPR.EQ.1) THEN   
 20     WRITE(6,*)
        WRITE(6,*) '  PRINT IN tunes.log (Y/N)?' 
        READ(*,FMT='(A1)',ERR=20) REP        
        IF(REP .NE. 'N' .AND. REP .NE. 'n') REP = 'y'
      ELSEIF(KPR.EQ.2) THEN
        REP = 'y'
      ENDIF

      IF(REP.EQ. 'y') THEN

C----- Twiss parameters and emittance
        WRITE(6,*) '  Ellipse parameter calculated will'
     >  ,' be printed in tunes.log'

        CALL READC9(KKEX,KLET)
        WRITE(NLOG,FMT='(/,'' Particle # '',I6,'' to '',I6,
     >  ''  @ lmnt # '',I5,''. '',I6,'' points,  IEX='',I2)') 
     >  NT1,NT2,LM,NINT(SNPT),KKEX
        WRITE(NLOG,100) 

        DO 4 JJ=1,3
          IF(U(JJ).LE. 1.D-15) THEN
            WRITE(NLOG,121) TXT(JJ),U(JJ) 
C           WRITE(6,121) TXT(JJ),U(JJ) 
          ELSE
            WRITE(NLOG,120) 
     >         TXT(JJ),B(JJ),A(JJ),G(JJ),U(JJ),YM(JJ),YPM(JJ)
          ENDIF
 4      CONTINUE

C------- Smear
        WRITE(NLOG,130) 
        DO 51 JJ=1,3
          I=2*JJ-1
          WRITE(NLOG,131) TXT(JJ),U(JJ),SMEAR(JJ),UMI(JJ),UMA(JJ)
 51     CONTINUE

        RETURN 1
      ENDIF      

      RETURN 2
      END
      SUBROUTINE BLOCK
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

C----- NUMERO DES UNITES LOGIQUES D'ENTREES-SORTIE
      COMMON/CDF/ 
     >      IES,IORDRE,LCHA,LIST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN
 
C----- CONSTANTES
      COMMON/CONST/ CL,PI,DPI,RAD,DEG,QE,AH

      PARAMETER (MRD=9)
      COMMON/DROITE/ AM(MRD),BM(MRD),CM(MRD),IDRT

      COMMON/EFBS/ AFB(MRD), BFB(MRD), CFB(MRD), IFB

      CHARACTER  KAR(41)
      COMMON/KAR/ KAR
 
      COMMON/OBJET/ FO(6,1),KOBJ,IDMAX,IMAXT
 
      LOGICAL ZSYM
      COMMON/OPTION/ KORD,KFLD,MG,LC,ML,ZSYM
 
      COMMON/PTICUL/ AAM,Q,G,TO
 
      COMMON/REBELO/ NRBLT,IPASS,KWRI,NNDES,STDVM
 
      CHARACTER FAM*8,KLEY*10
      PARAMETER (MXF=30,MXC=10) 
      COMMON/SCAL/
     >  SCL(MXF,MXC),TIM(MXF,MXC),FAM(MXF),KTI(MXF),KSCL,KLEY
 
C----- CONVERSION COORD. (CM,MRD) -> (M,RD)
      PARAMETER (MXJ=7)
      COMMON/UNITS/ UNIT(MXJ-1)

      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB

      DATA NDAT,NRES,NPLT,NFAI,NMAP,NSPN/3,4,1,2,8,9/
 
      DATA CL, PI, DPI, RAD, DEG, QE, AH /
     >2.99792458D8 , 3.141592653589D0, 6.283185307178D0,
     > .01745329252D0 , 57.29577951D0, 1.602176487D-19, 6.626075D-34 /
 
      DATA IDRT / 0 /

      DATA IFB / 0 /

      DATA (KAR(I),I=1,41) /
     > 'O','A','B','C','D','E','F','G','H','I','J','K','L','M','N'
     >,'P','Q','R','T','U','V','W','X','Y','Z','2','3','4','5','6'
     >,'7','8','(',')','+','-','/','=','"','0','*'/
 
      DATA KOBJ /0/
 
      DATA MG,LC,ML,ZSYM/ 1,2,3,.TRUE./
 
      DATA Q / 1.60217733D-19  /
 
      DATA NRBLT,IPASS/0, 1/
 
      DATA (FAM(I),I=1,MXF)/
     > 'AIMANT' , 'QUADRUPO', 'SEXTUPOL', 'QUADISEX' , 'SEXQUAD'
     >,'TOSCA3D', 'OCTUPOLE', 'DECAPOLE', 'DODECAPO'
     >, 'TOSCA' , 'MULTIPOL' , 'DIPOLE'
     >, 'BEND'    , 'SOLENOID' , 'CAVITE'
     >,'POISSON', 14*' ' /
 
C                    Y     T     Z       P     X,S   dp/p
      DATA UNIT / .01D0,.001D0,.01D0, .001D0, .01D0, 1.D0 /

C      DATA KX,KY,IAX,LIS,NB /6, 2, 1, 1, 100 /
      DATA KX,KY,IAX,LIS,NB /2, 3, 1, 1, 100 /

      RETURN

      ENTRY UNITR(KXI,KYI,
     >                    UXO,UYO)
      IF(KXI.EQ.1) THEN
        UXO = UNIT(6)
      ELSEIF(KXI.LE.MXJ) THEN
        UXO = UNIT(KXI-1)
      ELSE
        UXO = 1.D0
      ENDIF
      IF(KYI.EQ.1) THEN
        UYO = UNIT(6)
      ELSEIF(KYI.LE.MXJ) THEN
        UYO = UNIT(KYI-1)
      ELSE
        UYO = 1.D0
      ENDIF
      RETURN
      END
      SUBROUTINE READCO(NL,LM,
     >                        KART,LET,YZXB,NDX,*,*,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ----------------------------------------------------
C     Look for and read coordinates, etc. of particle # NT
C     ----------------------------------------------------
      CHARACTER*1 LET
      INCLUDE 'MXVAR.H'
      DIMENSION YZXB(MXVAR),NDX(5)

      INCLUDE 'MXLD.H'
      COMMON/LABCO/ ORIG(MXL,6) 
      COMMON/CDF/ IES,IORDRE,LCHA,LIST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN
c      COMMON/LUN/ NDAT,NRES,NPLT,NFAI,NMAP,NSPN
      INCLUDE 'MAXNPT.H'          
      COMMON/TRACKM/COOR(NPTMAX,9),NPTS,NPTR
      PARAMETER (MXJ=7)
      COMMON/UNITS/ UNIT(MXJ-1) 
      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB

      PARAMETER (MXS=4)
      DIMENSION FO(MXJ),F(MXJ),SI(MXS),SF(MXS)
      PARAMETER (MXT=10)
      DIMENSION SX(MXT),SY(MXT),SZ(MXT)
      
      CHARACTER*8 LBL1, LBL2
      CHARACTER KLEY*10

      LOGICAL BINARY,BINAR,OKKP,OKKT,OKKL

      CHARACTER*1 KLET, KLETO, KLETI, TX1

      SAVE KP1, KP2, KP3, BINARY
      SAVE KL1, KL2
      SAVE KT1, KT2
      SAVE KKEX, KLET

      DATA KP1, KP2, KP3 / 1, NPTMAX, 1 /
      DATA KT1, KT2 / 1, 1 /
      DATA KL1, KL2 / 1, 999999 /
      DATA KKEX, KLET / 1, '*' / 

      IF(NL .EQ. NSPN) THEN
      ELSE
        IF(NL .EQ. NFAI) THEN
C--------- read in zgoubi.fai type storage file

          IMAX = 0
          IF(BINARY) THEN
 222        CONTINUE
            READ(NL,ERR=222,END=10) 
     >      KEX,(FO(J),J=1,7),
     >      (F(J),J=1,7), 
     >      (SI(J),J=1,4),(SF(J),J=1,4),
     >      ENEKI, ENERG, 
     >      IT, IREP, SORT, AMQ1,AMQ2,AMQ3,AMQ4,AMQ5, RET, DPR, PS,
     >      BORO, IPASS,NOEL, KLEY,LBL1,LBL2,LET
C            IF(LM .NE. -1) THEN
C              IF(LM .NE. NOEL) GOTO 222
C            ENDIF


            IF(.NOT. OKKT(KT1,KT2,IT,KEX,LET,
     >                             IEND)) GOTO 222
            IF(.NOT. OKKP(KP1,KP2,KP3,IPASS,
     >                                IEND)) GOTO 222
            IF(.NOT. OKKL(KL1,KL2,NOEL,
     >                                IEND)) GOTO 222
            IF(IEND.EQ.1) GOTO 91

          ELSE
 21         READ(NL,110,ERR=21,END=10)
     >      KEX,(FO(J),J=1,7),
     >      (F(J),J=1,7), 
     >      (SI(J),J=1,4),(SF(J),J=1,4),
     >      ENEKI, ENERG, 
     >      IT, IREP, SORT, AMQ1,AMQ2,AMQ3,AMQ4,AMQ5, RET, DPR, PS,
     >      BORO, IPASS,NOEL, 
     >      TX1,KLEY,TX1,TX1,LBL1,TX1,TX1,LBL2,TX1,TX1,LET,TX1
C     >      KLEY,LBL1,LBL2,LET
            INCLUDE "FRMFAI.H"
CCCCCCCCCCC           if(it.eq.1) yref = f(2)
C            IF(LM .NE. -1) THEN
C              IF(LM .NE. NOEL) GOTO 21
C            ENDIF

            IF(.NOT. OKKT(KT1,KT2,IT,KEX,LET,
     >                             IEND)) GOTO 21

            IF(.NOT. OKKP(KP1,KP2,KP3,IPASS,
     >                                IEND)) GOTO 21

            IF(.NOT. OKKL(KL1,KL2,NOEL,
     >                           IEND)) GOTO 21

            IF(IEND.EQ.1) GOTO 91

          ENDIF

        ELSEIF(NL .EQ. NPLT) THEN
C--------- read in zgoubi.plt type storage file

          IMAX = 0
          IF(BINARY) THEN
 232         CONTINUE
            READ(NL,ERR=232,END=10) 
     >      KEX,(FO(J),J=1,7),
     >      (F(J),J=1,7), DS, 
     >      KART, IT, IREP, SORT, XX, BX, BY, BZ, RET, DPR, 
     >      EX, EY, EZ, BORO, IPASS,NOEL, KLEY,LBL1,LBL2,LET
C            IF(LM .NE. -1) THEN
C              IF(LM .NE. NOEL) GOTO 232
C            ENDIF

            IF(.NOT. OKKT(KT1,KT2,IT,KEX,LET,
     >                             IEND)) GOTO 232
            IF(.NOT. OKKP(KP1,KP2,KP3,IPASS,
     >                                IEND)) GOTO 232

            IF(.NOT. OKKL(KL1,KL2,NOEL,
     >                           IEND)) GOTO 232

            IF(IEND.EQ.1) GOTO 91

          ELSE
 31         READ(NL,100,ERR=31,END=10)
     >      KEX,(FO(J),J=1,MXJ),
     >      (F(J),J=1,MXJ), DS,
     >      KART, IT, IREP, SORT, XX, BX, BY, BZ, RET, DPR,
     >      EX, EY, EZ, BORO, IPASS, NOEL, 
     7      TX1,KLEY,TX1,TX1,LBL1,TX1,TX1,LBL2,TX1,
     8                                            TX1,LET,TX1
C     >      KLEY,LBL1,LBL2,LET
            INCLUDE "FRMPLT.H"
CCCCCCCCCCC           if(it.eq.1) yref = f(2)

C            IF(LM .NE. -1) THEN
C              IF(LM .NE. NOEL) GOTO 31
C            ENDIF

            IF(.NOT. OKKT(KT1,KT2,IT,KEX,LET,
     >                             IEND)) GOTO 31

            IF(.NOT. OKKP(KP1,KP2,KP3,IPASS,
     >                                IEND)) GOTO 31

            IF(.NOT. OKKL(KL1,KL2,NOEL,
     >                           IEND)) GOTO 31

            IF(IEND.EQ.1) GOTO 91

          ENDIF
        ENDIF        

        IF(KP1.GE.0.AND.IPASS.GT.KP2) GOTO 91

C------- dp/p
        J = 1
        JU = 6
        YZXB(J)   =   F(J)   * UNIT(JU)     
C        YZXB(J)   =   1.D0 + F(J)  
        YZXB(J+10) =  FO(J)   * UNIT(JU)        ! dp/p_initial

C------- Y, T, Z, P, S, Time
        DO 20 J=2,MXJ
          JU = J-1
          YZXB(J)   =  F(J)   * UNIT(JU)     
CCCCCCC          if(j.eq.2) YZXB(J) = (f(j)-yref) * UNIT(JU)
          YZXB(J+10) = FO(J)  * UNIT(JU) 
 20     CONTINUE

C           write(66,*) ' sbr readco it ipass, f7 :',it,ipass,yzxb(7)

C------- KART=1 : Cartesian coordinates, X is current x-coordinate
C        KART=2 : Cylindrical coordinates, X is current angle
        YZXB(8) = XX

        IF(KART .EQ. 1)  YZXB(8) = YZXB(8) * UNIT(5)
C         step size :
        YZXB(9) = DS       * UNIT(5)
C         r = sqrt(y^2+z^2) :
        YZXB(10) = SQRT(YZXB(2)*YZXB(2) + YZXB(4)*YZXB(4))
        YZXB(18) = RET
C------- (p_ps)/ps
        YZXB(19) = DPR            
C-------- momentum
C        YZXB(19) = BORO * (1.D0+F(1))*0.299792458D0   
        YZXB(20) = ENERG
C         convert B from kG to T
        YZXB(30) = BX      * .1D0
        YZXB(31) = BY      * .1D0
        YZXB(32) = BZ      * .1D0
        YZXB(33) = SQRT(BY*BY +  BZ*BZ) * .1D0

        YZXB(34) = EX
        YZXB(35) = EY     !!!/YZXB(2)
        YZXB(36) = EZ 
        YZXB(37) = SQRT(EY*EY +  EZ*EZ)

        PI2 = 2.D0*ATAN(1.D0)
        YINL = F(2)* UNIT(1) !!!!!!* AMAG  
        ZINL = F(4)* UNIT(3) !!!!!!* AMAG
C FM, Dec. 05       XINL = XX* UNIT(5) - ORIG(NOEL,5)
        XINL = XX* UNIT(5) + ORIG(NOEL,5)
        YZXB(44) = ZINL 
        PHI = ORIG(NOEL,6)
        CT = COS(ORIG(NOEL,4)+PHI) 
        ST = SIN(ORIG(NOEL,4)+PHI)
        YZXB(48) = ( XINL*CT - YINL*ST) + ORIG(NOEL,1)
        YZXB(42) = ( XINL*ST + YINL*CT) + ORIG(NOEL,2)
      ENDIF

C      Location about where particle was lost
      YZXB(38) = SORT * 1.D-2
      YZXB(39) = IPASS 
      YZXB(57) = NOEL

      NDX(1)=KEX
      NDX(2)=IT
      NDX(3)=IREP
      NDX(4)=IMAX
      NDX(5)=NOEL
      RETURN

 91   RETURN 3

C------------------ Pass #, KP1 to KP2
      ENTRY READC1(KP1O,KP2O)
C------- Read pass #,  KP1 to KP2
      KP1O=KP1
      KP2O=KP2
      RETURN
C--
      ENTRY READC2(LN)
C------- Write pass #,  KP1 to KP2
 12   WRITE(6,FMT='(''  Option status is now : KP1='',I6,
     >   '',   KP2='',I6)') KP1, KP2
        WRITE(6,FMT='(''     Available options are : '',
     >   /,10X,'' KP1, KP2 > 0 : will plot within range [KP1,KP2]'', 
     >   /,10X,'' KP1=-1, KP2 > 0 : will plot every KP2 other pass '')')
        WRITE(6,FMT='(/,
     >        '' Enter desired values KP1, KP2  ( 0 0 to exit ) : '')')
        READ(LN,*,ERR=12) KP1NEW, KP2NEW
      GOTO 11
C--
      ENTRY READC2B(KP1W,KP2W,KP3W)
C------- Pass KP1 to pass KP2, step KP3
        KP1NEW=KP1W
        KP2NEW=KP2W
        KP3NEW=KP3W
 11     CONTINUE
        IF(KP1NEW.NE.0) KP1=KP1NEW
        IF(KP2NEW.NE.0) KP2=KP2NEW
        IF(KP3NEW.NE.0) KP3=KP3NEW
      RETURN
C----------------------------

C------------------ Element #, KL1 to KL2
      ENTRY READC3(KL1O,KL2O)
C------- Read lmnt #,  KL1 to KL2
      KL1O=KL1
      KL2O=KL2
      RETURN
C--
      ENTRY READC4(LN)
        KLA = KL1
        KLB = KL2
C------- Write lmnt #,  KL1 to KL2
 30     WRITE(6,FMT='('' Observation now at elements  KL1='',I6,
     >   ''to   KL2='',I6)') KL1, KL2
        WRITE(6,FMT='(''     Available options for elment are : '',
     >   /,10X,'' 0 < KL1 < KL2  : will plot within range [KL1,KL2]'', 
     >   /,10X,'' KL1=-1, KL2 > 0 : will plot every KL2 other pass '')')
        WRITE(6,FMT='(/,
     >        '' Enter desired values KL1, KL2  : '')')
        READ(LN,fmt='(2I3)',ERR=33) KL1NEW, KL2NEW
      GOTO 32
 33     KL1NEW = KLA
        KL2NEW = KLB
      GOTO 32
C--
      ENTRY READC4B(KL1W,KL2W)
C------- Lmnt  KL1 to lmnt KL2
        KL1NEW=KL1W
        KL2NEW=KL2W
 32     CONTINUE
        IF(KL1NEW.NE.0) KL1=KL1NEW
        IF(KL2NEW.NE.0) KL2=KL2NEW
      RETURN
C----------------------------

C------------------ Traj # KT1 to traj # KT2
      ENTRY READC5(
     >             KT1O,KT2O)
C------- Read traj.  KT1 to KT2
        KT1O=KT1
        KT2O=KT2
      RETURN
C--
      ENTRY READC6B(KT1W,KT2W)
C------- Traj KT1 to traj KT2
        KT1NEW=KT1W
        KT2NEW=KT2W
 51     CONTINUE
        IF(KT1NEW.NE.0) KT1=KT1NEW
        IF(KT2NEW.NE.0) KT2=KT2NEW
      RETURN
C----------------------------

      ENTRY READC7(BINAR)
      BINAR=BINARY
      RETURN

      ENTRY READC8(BINAR)
      BINARY=BINAR
      RETURN

      ENTRY READC9(KKEXO,KLETO)
        KKEXO=KKEX
        KLETO=KLET
      RETURN

      ENTRY READCX(KKEXI,KLETI)
        KKEX=KKEXI
        KLET=KLETI
      RETURN      

 10   CONTINUE
      KT2 = IT
      RETURN 1

 99   RETURN 2
      END
      SUBROUTINE OPNDEF(LU2O,DEFN2O,LUO,
     >                                  FNO,OKOPN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*(*) DEFN2O, FNO
      LOGICAL OKOPN

      LOGICAL EXS, OPN, BINARY 
      CHARACTER*11 FRMT

      IF( LU2O .EQ. -1) THEN
C------- Looks for a free LUO starting from #10
        OKOPN = .FALSE.
        LU2O = 9
 1      CONTINUE  
        LU2O = LU2O+1
        IF(LU2O.EQ.100) GOTO 96
        INQUIRE(UNIT=LU2O,EXIST=EXS,OPENED=OPN,IOSTAT=IOS)
        IF( IOS .EQ. 0 .AND. OPN ) GOTO 1
      ELSEIF(LU2O .EQ. LUO) THEN
C------- Check whether LUO is already open
        INQUIRE(UNIT=LUO,EXIST=EXS,OPENED=OPN,NAME=FNO,IOSTAT=IOS)
        IF(IOS .EQ. 0) THEN
          OKOPN = OPN
        ELSE
          OKOPN = .FALSE.
        ENDIF
      ELSE
        IF(LUO.GT.0) CLOSE(LUO)
        OKOPN = .FALSE.
      ENDIF

      IF(.NOT. OKOPN) THEN
C--------- Check existence of DEFN2O
        INQUIRE(FILE=DEFN2O,EXIST=EXS,IOSTAT=IOS)

        IF(IOS .EQ. 0) THEN

          BINARY=DEFN2O(1:2).EQ.'B_' .OR. DEFN2O(1:2).EQ.'b_'
          IF(BINARY) THEN 
            FRMT='UNFORMATTED'
          ELSE
            FRMT='FORMATTED'
          ENDIF
          IF(EXS) THEN
            OPEN(UNIT=LU2O,FILE=DEFN2O,STATUS='OLD',ERR=99,IOSTAT=IOS,
     >           FORM=FRMT)
            IF(IOS.NE.0) GOTO 97
            I4=4
            IPRNT = 0
            CALL HEADER(LU2O,I4,IPRNT,BINARY,*99)
          ELSE
            OPEN(UNIT=LU2O,FILE=DEFN2O,STATUS='NEW',ERR=99,IOSTAT=IOS,
     >           FORM=FRMT)
            IF(IOS.NE.0) GOTO 97
          ENDIF

          FNO = DEFN2O
          LUO = LU2O
          OKOPN = .TRUE.
        ELSE

          WRITE(6,*) ' Exec error occured in Subroutine OPNDEF : '
          WRITE(6,*) ' IOS NON-ZERO '

        ENDIF
        CALL READC8(BINARY)
      ENDIF

      IF(OKOPN) WRITE(6,FMT='(2A)') 'Opened file is ',DEFN2O

      RETURN

 96   WRITE(6,*) ' Exec error occured in Subroutine OPNDEF : ' 
      WRITE(6,*) ' Logical unit # exceeds 100 ; CANNOT OPEN ',DEFN2O
      RETURN

 97   WRITE(6,*) ' Error occured while in Subroutine OPNDEF : '
      WRITE(6,*) '      IOS NON-ZERO ; CANNOT OPEN ',DEFN2O
 99   CONTINUE
      RETURN

C 98   CONTINUE
C      CLOSE(LU2O)
C      RETURN

      END
      SUBROUTINE INIGR(
     >                 NLOG, LM, NOMFIC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*(*) NOMFIC

      LOGICAL OKECH, OKVAR, OKBIN
      COMMON/ECHL/OKECH, OKVAR, OKBIN

      INCLUDE 'MXVAR.H'
      CHARACTER KVAR(MXVAR)*10, KPOL(2)*9, KDIM(MXVAR)*10
      COMMON/INPVR/ KVAR, KPOL, KDIM

      INCLUDE 'MAXNPT.H'
      COMMON/TRACKM/COOR(NPTMAX,9),NPTS,NPTR

      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB

      CHARACTER * 9   DMY
      CHARACTER*80 TXT
      CHARACTER LOGOT*18, TEMP*80

      DATA  OKECH ,OKVAR, OKBIN / .FALSE., .TRUE., .FALSE.  /

      DATA KVAR/
     >' dp/p ','   Y  ','   T  ','   Z  ','   P  ','   S  ',' Time ',
     >'   X  ',' Step ','   r  ',
     >'dp/p|o','  Yo  ','  To  ','  Zo  ','  Po  ','  So  ',' Time ',
     >'Phase ',' dp/p ','KinEnr',
     >'  Sx  ','  Sy  ','  Sz  ',' <S>  ',
     >' <Sx> ',' <Sy> ',' <Sz> ','COUNT ','      ',
     >'  Bx  ','  By  ','  Bz  ','  Br  ',
     >'  Ex  ','  Ey  ','  Ez  ','  Er  ',
     >' S_out',' Pass#'  ,2*'      ',
     >' Y_Lab','      ',' Z_Lab','      ','      ','      ',' X_Lab',
     >8*' ',
     >' lmnt#' ,
     >13*' '
     >/
      DATA KPOL/ 'CARTESIAN' , 'CYLINDR.' /

      DATA NPTS / NPTMAX /

      DATA KDIM/
     >'       ','  (m)  ',' (rad) ','  (m)  ',' (rad) ','  (m)  ',
     >'(mu_s) ','  (m)  ','  (m)  ','  (m)  '                  ,
     >'       ','  (m)  ',' (rad) ','  (m)  ',' (rad) ','  (m)  ',
     >'       ',' (rad) ','       ',' (MeV) ',9*'       ',
     > 4*'  (T)  ', 4*'(eV/m) ' ,
     >'  (m)  ',3*' ',
     >'  (m)  ','      ','  (m)  ','      ','      ','      ','  (m)  ',
     >22*'   '/


      DATA KARSIZ / 4 /
      SAVE KARSIZ 

C----- Number of the lmnt concerned by the plot (-1 for all)
      LM = -1

C----- tunes log unit
      NLOG = 30
C----- Input data file name
C      Normally, default is zgoubi.fai, .plt, .spn, .map...
      NOMFIC = 'b_zgoubi.fai'

      RETURN
      END

      SUBROUTINE STORCO(NL,LM,KPS,BINARY,
     >                                   NPASS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ---------------------------------------------------
C     Read coordinates from zgoubi output file, and store  
C     ---------------------------------------------------
      LOGICAL BINARY
      INCLUDE 'MAXNPT.H'          
      COMMON/TRACKM/COOR(NPTMAX,9),NPTS,NPTR

      CHARACTER LET 

      INCLUDE 'MXVAR.H'
      DIMENSION YZXB(MXVAR),NDX(5)

c      logical first(mxt)
c      save first
c      data first / mxt * .true. /

      CALL RAZ(COOR,NPTMAX*9)

      CALL REWIN2(NL,*96)
      WRITE(6,*) '  REWIND-ing...' 
      WRITE(6,*) '  READING  AND  STORING  COORDINATES...'

      NOC=0
      NRBLT = -1 
C----- BOUCLE SUR READ FICHIER NL 
 44   CONTINUE
        CALL READCO(NL,LM,
     >                    KART,LET,YZXB,NDX,*10,*44,*17)

C----- NDX: 1->KEX, 2->IT, 3->IREP, 4->IMAX

        NOC=NOC+1
        IF(NINT(YZXB(39)) .GE. NRBLT+1) NRBLT = NINT(YZXB(39)) -1
        CALL FILCOO(KPS,NOC,YZXB,NDX)
        IF(NOC.EQ. NPTR) GOTO 10

      GOTO 44             
C     ----------------------------------

 99   CONTINUE
      WRITE(6,*) ' *** Coordinates  storage  stopped: error during',
     > ' read of event # ',NOC+1
      GOTO 11

 10   CONTINUE
      WRITE(6,*) ' END OF FILE  encountered, read ',noc,' points'
      GOTO 11

 17   CONTINUE
      WRITE(6,*) ' Required # of passes has been read, ',noc,' points'
      GOTO 11

 11   CONTINUE
      NPASS = NRBLT + 1
      NPTR=NOC
      CALL READC5(KT1,KT2)
      write(*,*) ' Pgm stroboscopy, trjctries kt1:kt2 : ',kt1,':',kt2
      IF(KT1 .EQ. -1 .OR. KT2 .GT. KT1) THEN
        WRITE(6,*) '  Analysis of particles from a set '
        IF(KPS.EQ. 0) WRITE(6,*) '  Initial  phase-space'
        IF(KPS.EQ. 1) WRITE(6,*) '  Final  phase-space'
        WRITE(6,*) '  # of turns  in the structure :',NPASS
      ELSEIF(KT1 .EQ. KT2) THEN
        IF(NPASS.EQ.0) THEN
        ELSE
          WRITE(6,*) '  A single  particle  analized          :',KT1
          WRITE(6,*) '  # of turns in the structure   :',NPASS
        ENDIF
      ENDIF

      WRITE(6,*) ' ',NOC,' points have been stored'

 96   RETURN                  
      END
      FUNCTION BINARF(LUN)
      LOGICAL BINARF
      CHARACTER*7 QUID
      INQUIRE(UNIT=LUN,UNFORMATTED=QUID)
      BINARF=QUID.EQ.'YES'
      RETURN
      END
      FUNCTION OKKL(KL1,KL2,NOEL,
     >                           IEND)
      LOGICAL OKKL

      INCLUDE "OKKL.H"

      RETURN
      END
      FUNCTION OKKP(KP1,KP2,KP3,IPASS,
     >                            IEND)
      LOGICAL OKKP

      INCLUDE "OKKP.H"

      RETURN
      END
      FUNCTION OKKT(KT1,KT2,IT,KEX,LET,
     >                                 IEND)
      LOGICAL OKKT
      CHARACTER*1 LET,KLETO
      LOGICAL OKKT5

      INCLUDE 'MAXNPT.H'          
      LOGICAL LOST(NPTMAX)
      SAVE LOST
      DATA LOST / NPTMAX* .FALSE. / 

      INCLUDE "OKKT.H"

      CALL READC9(KEXO,KLETO)

C------ Plot or store only as long as KEX is correct 
      IF(KEXO.NE.99) OKKT = OKKT .AND. (KEX.EQ.KEXO)
      IF(KLETO.NE.'*') THEN
        IF(KLETO.EQ.'S') THEN
C------ Only secondary particles (from decay process using MCDESINT) can be plotted
          OKKT = OKKT .AND. (LET.EQ.'S')
        ELSEIF(KLETO.EQ.'P') THEN
C------ Only parent particles (if decay process using MCDESINT) can be plotted
          OKKT = OKKT .AND. (LET.NE.'S')
        ENDIF
      ENDIF

      IF(.NOT.LOST(IT)) LOST(IT) = KEX.LE.0
      RETURN

      ENTRY OKKT5(KT)      
      OKKT5 = .NOT. LOST(KT)
      RETURN
      END
      SUBROUTINE FLUSH2(IUNIT,BINARY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      LOGICAL BINARY
      CHARACTER*80 TXT80
      BACKSPACE(IUNIT)
      IF(.NOT.BINARY) THEN
        READ(IUNIT,FMT='(A80)') TXT80
      ELSE
        READ(IUNIT) TXT80
      ENDIF
      RETURN
      END
      SUBROUTINE REWIN2(NL,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL BINAR
      REWIND(NL)
C------- Swallow the header (4 lines)
      CALL READC7(BINAR)
      I4 = 4
      IPRNT = 0
      CALL HEADER(NL,I4,IPRNT,BINAR,*99)
      RETURN
 99   RETURN 1
      END
      FUNCTION IDLUNI(LN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL IDLUNI

      LOGICAL OPN

      I = 20
 1    CONTINUE
        INQUIRE(UNIT=I,ERR=99,IOSTAT=IOS,OPENED=OPN)
        I = I+1
        IF(I .EQ. 100) GOTO 99
        IF(OPN) GOTO 1
        IF(IOS .GT. 0) GOTO 1
      
      LN = I-1
      IDLUNI = .TRUE.
      RETURN

 99   CONTINUE
      LN = 0
      IDLUNI = .FALSE.
      RETURN
      END
      FUNCTION EMPTY(STR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL EMPTY
      CHARACTER*(*) STR
      INTEGER FINSTR
      EMPTY = FINSTR(STR) .EQ. 0
      RETURN
      END
      SUBROUTINE RAZ(TAB,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION TAB(1)
      DO 1 I=1,N
 1      TAB(I) = 0.D0
      RETURN
      END
      SUBROUTINE FILCOO(KPS,NOC,YZXB,NDX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION YZXB(*), NDX(*)
      INCLUDE 'MAXNPT.H'
      PARAMETER (NTR=NPTMAX*9)
      COMMON/TRACKM/COOR(NPTMAX,9),NPTS,NPTR

      DIMENSION KXYL(6)
      SAVE KXYL

C      DATA KXYL / 2,3,4,5,7,1 /
C                 Y  T  Z  P  time energ
C      DATA KXYL / 2, 3, 4, 5, 7,   20 /
C                 Y  T  Z  P  time dp/p
      DATA KXYL / 2, 3, 4, 5, 7,   1 /
C      DATA KXYL / 2,3,4,5,18,19 /

C----------- Current coordinates
          II = 0
C----------- Initial coordinates
          IF(KPS.EQ. 0) II = 10

          COOR(NOC,1)=YZXB(KXYL(1)+II )
          COOR(NOC,2)=YZXB(KXYL(2)+II )
          COOR(NOC,3)=YZXB(KXYL(3)+II )
          COOR(NOC,4)=YZXB(KXYL(4)+II )
          COOR(NOC,5)=YZXB(KXYL(5)+II )
          COOR(NOC,6)=YZXB(KXYL(6) )
      RETURN
      END
      SUBROUTINE HEADER(NL,N,IPRNT,BINARY,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL BINARY
      CHARACTER*80 TXT80
      WRITE(6,FMT='(/,A,I2,A)') ' Now reading  ',N,'-line  file  header'
      IF(.NOT.BINARY) THEN
        READ(NL,FMT='(A80)',ERR=99,END=99) TXT80
        if(IPRNT.eq.1) 
     >        WRITE(6,FMT='(A)') TXT80
        READ(NL,FMT='(A80)',ERR=99,END=99) TXT80
        if(IPRNT.eq.1) 
     >        WRITE(6,FMT='(A)') TXT80
      ELSE
        READ(NL,ERR=99,END=89) TXT80
        if(IPRNT.eq.1) 
     >        WRITE(6,FMT='(A)') TXT80
        READ(NL,ERR=99,END=89) TXT80
        if(IPRNT.eq.1) 
     >        WRITE(6,FMT='(A)') TXT80
      ENDIF
      IF(.NOT.BINARY) THEN
        DO 1 I=3, N
           READ(NL,FMT='(A)',ERR=99,END=89) TXT80
           if(IPRNT.eq.1) 
     >           WRITE(6,FMT='(A)') TXT80
 1      CONTINUE
      ELSE
        DO 2 I=3, N
           READ(NL,          ERR=99,END=89) TXT80
           if(IPRNT.eq.1) 
     >           WRITE(6,FMT='(A)') TXT80
 2      CONTINUE
      ENDIF
      RETURN
 89   CONTINUE
      WRITE(6,*) 'END of file reached while reading data file header'
      RETURN 1
 99   CONTINUE
      WRITE(6,*) '*** READ-error occured while reading data file header'
      WRITE(6,*) '        ... Empty file ?'
      RETURN 1
      END
      FUNCTION FINSTR(STR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER FINSTR
      CHARACTER * (*) STR
C     -----------------------------------
C     Renvoie dans FINSTR le rang du
C     dernier caractere non-blanc de STR.
C     Renvoie 0 si STR est vide ou blanc.
C     -----------------------------------

      FINSTR=LEN(STR)+1
1     CONTINUE
         FINSTR=FINSTR-1
         IF(FINSTR.EQ. 0) RETURN
         IF (STR(FINSTR:FINSTR).EQ. ' ') GOTO 1
      RETURN
      END
