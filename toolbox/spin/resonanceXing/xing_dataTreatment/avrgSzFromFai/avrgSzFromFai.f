C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL OKOPN, CHANGE
C----- PLOT SPECTRUM     
      COMMON/CDF/ IES,IORDRE,LCHA,LIST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN

      parameter(mst=4)
      CHARACTER NOMFIC*80, txt132*132, stra(mst)*20
      logical strcon
      integer debstr
      character*2 HV
      data HV / ' ' /

      write(6,*) ' '
      write(6,*) '----------------------------------------------------'
      write(6,*) 'NOW RUNNING PGM AVRGSZFROMFAI...'
      write(6,*) '----------------------------------------------------'
      write(6,*) ' '

C Get number of turns, from zgoubi.res
      open(unit=nres,file='zgoubi.res')
 2    continue
        read(nres,fmt='(a)',end=10) txt132
        txt132 = txt132(debstr(txt132):132)

        if(strcon(txt132,'''REBELOTE''',10,
     >                                     IS)) then
            read(nres,*) npass
            write(6,*) ' PGM avrgSzFromFai.f, got from zgoubi.res npass'
     >          ,' = ',npass
            goto 10
        endif
        goto 2
 10   continue

      call INIGR(
     >           NLOG, LM, NOMFIC)
      NL = NFAI
      okopn = .false.
      change = .true.

      write(6,*) '-------'
      write(6,*) ' Pgm avrgSzFromFai. Now computing initial polar.'
      ipass1 = 1
      ipass2 = npass/20         
      write(6,*) ' using passes ',ipass1,' to ',ipass2 
      write(6,*)
      call averag(NLOG,NL,LM,OKOPN,CHANGE,HV,ipass1,ipass2)

      write(6,*) '-------'
      write(6,*) ' Pgm avrgSzFromFai. Now computing final polar.'
      ipass1 = 1+npass*19/20
      ipass2 = npass
      write(6,*) ' using passes ',ipass1,' to ',ipass2 
      write(6,*)
      call averag(NLOG,NL,LM,OKOPN,CHANGE,HV,ipass1,ipass2)
 
      stop
      end
      SUBROUTINE averag(NLOG,NL,LM,OKOPN,CHANGE,HV,ip1,ip2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL OKOPN, CHANGE
      character*(*) HV
C----- PLOT SPECTRUM     
      COMMON/CDF/ IES,IORDRE,LCHA,LIST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN
      PARAMETER (NCANAL=2500)
      COMMON/SPEDF/BORNE(6),SPEC(NCANAL,3),PMAX(3),NC0(3)
      PARAMETER (NTRMAX=99999,MXT=99999)
      COMMON/TRACKM/COOR(NTRMAX,9),NPTS,NPTR

      DIMENSION SM(4),SM2(4),SPK(8)
 
      LOGICAL OKECH
      CHARACTER REP, NOMFIC*80
      LOGICAL BINARY, BINARF
      CHARACTER HVL(3)*12

      SAVE NT
      DATA NT / -1 /

      LOGICAL OPN
      DATA OPN / .FALSE. /
      save opn, iun

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
        CALL READC6B(-NT,NPTS)
      ELSE
        CALL READC6B(NT,NT)        
      ENDIF

 6    CONTINUE
      IF(.NOT. OPN) THEN
          IF (IDLUNI(IUN)) THEN
            call system('cat avrgSzFromFai.out >>avrgSzFromFai.out_old')
            OPEN(UNIT=IUN,FILE='avrgSzFromFai.out',ERR=699)
            OPN = .TRUE.
            WRITE(6,*) 'Now opened file avrgSzFromFai.out',iun
            WRITE(IUN,*) 
     >      '%Spin components, average and peak values at start and end'
          ELSE
            GOTO 698
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
            NPTR = NTRMAX
            NPTS=NPTR
            I1 = 1
            CALL STORCO(NL,LM,I1  ,BINARY,ip1,ip2, 
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
            CALL AVRGS(IUN,LM,ip1,ip2,
     >                              SM,SM2,SPK)
               write(6,*) ' ///////////// '
               write(6,*) ' ///////////// '
          ENDIF

          IF(NT.EQ.-1) THEN
            KT = KT+1
            CALL READC5(
     >                  KT1,KT2)
            CALL READC6B(KT,KT)
            CHANGE = .TRUE.
            IF(KT.LE.KT2) GOTO 62             !--------- Coordinate reading/storing loop
          ENDIF

 69       CONTINUE
            CALL FLUSH2(IUN,.FALSE.)
            CHANGE = .TRUE.

      GOTO 99


 698  WRITE(6,*) ' *** Problem : No idle unit for avrgSzFromFai.out '
      stop
C      GOTO 99
 699  WRITE(6,*) ' *** Problem at OPEN avrgSzFromFai.out '
      stop
C      GOTO 99

 99   RETURN
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
 
      PARAMETER (MXT=99999)
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
     >                        KART,LET,YZXB,NDX,*,*)
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
      PARAMETER (NTRMAX=99999,MXT=99999)          
      COMMON/TRACKM/COOR(NTRMAX,9),NPTS,NPTR
      PARAMETER (MXJ=7)
      COMMON/UNITS/ UNIT(MXJ-1) 
      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB

      PARAMETER (MXS=4)
      DIMENSION FO(MXJ),F(MXJ),SI(MXS),SF(MXS)
      DIMENSION SX(MXT),SY(MXT),SZ(MXT)
      
      CHARACTER*8 LBL1, LBL2
      CHARACTER KLEY*10

      LOGICAL BINARY,BINAR,OKKP,OKKT,OKKL

      CHARACTER*1 KLET, KLETO, KLETI, TX1

      SAVE KP1, KP2, KP3, BINARY
      SAVE KL1, KL2
      SAVE KT1, KT2
      SAVE KKEX, KLET

      DATA KP1, KP2, KP3 / 1, 99999, 1 /
      DATA KT1, KT2 / 1, MXT /
      DATA KL1, KL2 / 1, 99999 /
      DATA KKEX, KLET / 1, '*' / 


      IF(NL .EQ. NSPN) THEN
      ELSE
        IF(NL .EQ. NFAI) THEN
C--------- read in zgoubi.fai type storage file

          IMAX = 0
          IF(BINARY) THEN
 222        CONTINUE
            READ(NL,ERR=99,END=10) 
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
 21         READ(NL,110,ERR=99,END=10)
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
            READ(NL,ERR=99,END=10) 
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
 31         READ(NL,100,ERR=99,END=10)
     >      KEX,(FO(J),J=1,MXJ),
     >      (F(J),J=1,MXJ), DS,
     >      KART, IT, IREP, SORT, XX, BX, BY, BZ, RET, DPR,
     >      EX, EY, EZ, BORO, IPASS,NOEL, 
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

        IF(KP1.GE.0.AND.IPASS.GT.KP2) RETURN 1

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
        YZXB(21) = SF(1)
        YZXB(22) = SF(2)    
        YZXB(23) = SF(3)    
        YZXB(24) = SF(4)
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

 91   RETURN 1

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
      ENTRY READC2B(KP1W,KP2W)
C------- Pass KP1 to pass KP2
        KP1NEW=KP1W
        KP2NEW=KP2W
 11     CONTINUE
        IF(KP1NEW.NE.0) KP1=KP1NEW
        IF(KP2NEW.NE.0) KP2=KP2NEW
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

C------------------ Traj #, KT1 to traj KT2
      ENTRY READC5(KT1O,KT2O)
C------- Read traj.,  KT1 to KT2
        KT1O=KT1
        KT2O=KT2
      RETURN
C--
      ENTRY READC6(LN)
C------- Specify traj. #,  KT1 to KT2
 50     WRITE(6,FMT='(''  Option status is now : KT1='',I6,
     >   '',   KT2='',I6)') KT1,KT2
        WRITE(6,FMT='(''     Available options are : '',
     >  /,9X,'' KT1, KT2 > 0  : will treat range [KT1,KT2]'', 
     >  /,9X,'' KT1=-1, KT2 > 0 '',
     >                    '' : will treat every KT2 other particle '')') 
        WRITE(6,FMT='(/,
     >        '' Enter desired values KT1, KT2  ( 0 0 to exit ) : '')')
        READ(LN,*,ERR=50) KT1NEW, KT2NEW
        IF(KT2NEW.GT.MXT) THEN
          WRITE(6,FMT='(9X,'' N2  cannot exceed '',I6)') MXT
          GOTO 50
        ENDIF
        IF(KT1NEW.GT.0) THEN
          IF(KT2NEW.LT.KT1NEW) GOTO 50
        ELSEIF(KT1NEW.NE.-1) THEN
          GOTO 50
        ENDIF
      GOTO 51
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
C            I4=4
C            IPRNT = 0
C            CALL HEADER(LU2O,I4,IPRNT,BINARY,*99)
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
C        IF(OKOPN) WRITE(6,FMT='(2A)') 'Opened file is ',DEFN2O
C      ELSE
C         WRITE(6,*) ' Already opened file ',DEFN2O
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

      PARAMETER (NCANAL=2500)
      COMMON/SPEDF/BORNE(6),SPEC(NCANAL,3),PMAX(3),NC0(3)

      PARAMETER (NTRMAX=99999,MXT=99999)
      COMMON/TRACKM/COOR(NTRMAX,9),NPTS,NPTR

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

C      DATA BORNE/ .01D0, .99D0, .01D0, .99D0, .001D0, .999D0 /
      DATA BORNE/ .0D0, .5D0, .0D0, .5D0, .001D0, .999D0 /
c      DATA NC0/ 2000, 2000, 2000 /
      DATA NC0/ 1000, 1000, 1000 /

      DATA NPTS / NTRMAX /

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

C-----  log unit
      NLOG = 30
C----- Input data file name
C      Normally, default is zgoubi.fai, .plt, .spn, .map...
      NOMFIC = 'b_zgoubi.fai'

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

      PARAMETER (NTRMAX=99999,MXT=99999)          
      LOGICAL LOST(NTRMAX)
      SAVE LOST
      DATA LOST / NTRMAX* .FALSE. / 

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
c         write(89,*) 'oeiuvc',okkt, okkt5,it
c         pause
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
      PARAMETER (NTRMAX=99999,MXT=99999)
      PARAMETER (NTR=NTRMAX*9)
      COMMON/TRACKM/COOR(NTRMAX,9),NPTS,NPTR

      DIMENSION KXYL(6)
      SAVE KXYL

C      DATA KXYL / 2,3,4,5,7,1 /
C                 Y  T  Z  P  time energ
C      DATA KXYL / 2, 3, 4, 5, 7,   20 /
C      DATA KXYL / 2,3,4,5,6,1 /
C      DATA KXYL / 2,3,4,5,18,19 /
      DATA KXYL / 21, 22, 23, 24, 23, 23 /

C----------- Current coordinates
          II = 0
C----------- Initial coordinates
          IF(KPS.EQ. 0) II = 10

          COOR(NOC,1)=YZXB(KXYL(1)+II )
          COOR(NOC,2)=YZXB(KXYL(2)+II )
          COOR(NOC,3)=YZXB(KXYL(3)+II )
          COOR(NOC,4)=YZXB(KXYL(4)+II )
C          COOR(NOC,5)=YZXB(KXYL(5)+II )
C          COOR(NOC,6)=YZXB(KXYL(6) )
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
      FUNCTION STRCON(STR,STRIN,NCHAR,
     >                                IS)
      implicit double precision (a-h,o-z)
      LOGICAL STRCON
      CHARACTER STR*(*), STRIN*(*)
C     ------------------------------------------------------------------------
C     .TRUE. if the string STR contains the string STRIN with NCHAR characters
C     at least once.
C     IS = position of first occurence of STRIN in STR
C     ------------------------------------------------------------------------

      INTEGER DEBSTR,FINSTR

      II = 0
      DO 1 I = DEBSTR(STR), FINSTR(STR)
        II = II+1
        IF( STR(I:I+NCHAR-1) .EQ. STRIN ) THEN
          IS = II
          STRCON = .TRUE.
          RETURN
        ENDIF
 1    CONTINUE
      STRCON = .FALSE.
      RETURN
      END
      FUNCTION DEBSTR(STRING)
      implicit double precision (a-h,o-z)
      INTEGER DEBSTR
      CHARACTER * (*) STRING

C     --------------------------------------
C     RENVOIE DANS DEBSTR LE RANG DU
C     1-ER CHARACTER NON BLANC DE STRING,
C     OU BIEN 0 SI STRING EST VIDE ou BLANC.
C     --------------------------------------

      DEBSTR=0
      LENGTH=LEN(STRING)
C      LENGTH=LEN(STRING)+1
1     CONTINUE
        DEBSTR=DEBSTR+1
C        IF(DEBSTR .EQ. LENGTH) RETURN
C        IF (STRING(DEBSTR:DEBSTR) .EQ. ' ') GOTO 1
        IF (STRING(DEBSTR:DEBSTR) .EQ. ' ') THEN
          IF(DEBSTR .EQ. LENGTH) THEN
            DEBSTR = 0
            RETURN
          ELSE
            GOTO 1
          ENDIF
        ENDIF

      RETURN
      END
      FUNCTION FINSTR(STRING)
      implicit double precision (a-h,o-z)
      INTEGER FINSTR
      CHARACTER * (*) STRING
C     --------------------------------------
C     RENVOIE DANS FINSTR LE RANG DU
C     DERNIER CHARACTER NON BLANC DE STRING,
C     OU BIEN 0 SI STRING EST VIDE ou BLANC.
C     --------------------------------------

      FINSTR=LEN(STRING)+1
1     CONTINUE
        FINSTR=FINSTR-1
        IF(FINSTR .EQ. 0) RETURN
        IF (STRING(FINSTR:FINSTR) .EQ. ' ') GOTO 1

      RETURN
      END
      SUBROUTINE IMP
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NTRMAX=99999,MXT=99999)
      dimension fo(7,mxt)
      save fo
      INCLUDE 'MXVAR.H'
      DIMENSION YZXB(MXVAR)

      RETURN

      entry speim2(yzxb,it)
      fo(2,it) = yzxb(12)
      fo(3,it) = yzxb(13)
      fo(4,it) = yzxb(14)
      fo(5,it) = yzxb(15)
      return

      RETURN
      END
      SUBROUTINE STORCO(NL,LM,KPS,BINARY,IP1,IP2,
     >                                           NPASS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ---------------------------------------------------
C     Read coordinates from zgoubi output file, and store  
C     ---------------------------------------------------
      LOGICAL BINARY
      PARAMETER (NTRMAX=99999,MXT=99999)          
      COMMON/TRACKM/COOR(NTRMAX,9),NPTS,NPTR

      CHARACTER LET 

      INCLUDE 'MXVAR.H'
      DIMENSION YZXB(MXVAR),NDX(5)

      logical first(mxt)
      save first
      data first / mxt * .true. /

      CALL RAZ(COOR,NTRMAX*9)

      CALL REWIN2(NL,*96)
      WRITE(6,*) '  READING  AND  STORING  COORDINATES...'

      NOC=0
      NRBLT = -1 
C----- BOUCLE SUR READ FICHIER NL 
 44   CONTINUE
        CALL READCO(NL,LM,
     >                    KART,LET,YZXB,NDX,*10,*44)

C----- NDX: 1->KEX, 2->IT, 3->IREP, 4->IMAX

        NOC=NOC+1

        IF(NINT(YZXB(39)) .lt. IP1) goto 44
        IF(NINT(YZXB(39)) .gt. IP2) goto 10

        IF(NINT(YZXB(39)) .GE. NRBLT+1) NRBLT = NINT(YZXB(39)) -1
        CALL FILCOO(KPS,NOC,YZXB,NDX)
        IF(NOC.EQ. NPTR) GOTO 10

        if(first(ndx(2))) then
          call speim2(yzxb,ndx(2))
          first(ndx(2)) = .false. 
        endif

      GOTO 44             
C     ----------------------------------

 10   CONTINUE
      WRITE(6,*) ' READ  OK; END  OF  FILE  ENCOUNTERED'
      NPASS = NRBLT + 1
      NPTR=NOC
      CALL READC5(KT1,KT2)
      write(6,*) ' Pgm avrSzFromFai, particles kt1:kt2 : ',kt1,':',kt2
      WRITE(6,*) ' ',NOC,' points have been stored'

 96   RETURN                  
      END
      SUBROUTINE avrgs(IUN,LM,ip1,ip2,
     >                              SM,SM2,SPK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION SM(*), SM2(*), SPK(*)
      PARAMETER (NTRMAX=99999,MXT=99999)
      COMMON/TRACKM/COOR(NTRMAX,9),NPTS,NPTR

      DIMENSION G(3)
      CHARACTER TXT(3)*5, REP*1

      INCLUDE 'MAXCOO.H'

      DIMENSION UM(3), UMI(3), UMA(3), SMEAR(3)
      CHARACTER*1 KLET

C      PARAMETER ( PI=3.1415926536, SQ2 = 1.414213562)
      PARAMETER ( PI=4.D0*ATAN(1.D0), SQ2 = SQRT(2.D0) )

      character TXTS(4)*3
      data TXTS / ' SX', ' SY', ' SZ', '|S|' / 

      DATA TXT/ 'HORI.', 'VERT.', 'LONG.'/

      DO 7 J=1,4
        SPK(2*J-1) = 1.D10
        SPK(2*J) = -1.D10
 7    CONTINUE

      DO 2 J=1,4
C------- JJ = 1, 2, 3, 4  for  Sx, SY, SZ, |S|
        XM=0.D0
        SNPT = 0.D0
        DO 21 I=ip1, ip2
            X = COOR(I,J)
            XM = XM + X
C----------- Min-max of the distribution
            IF(SPK(2*J-1) .GT. X) SPK(2*J-1) = X
            IF(SPK(2*J) .LT. X) SPK(2*J) = X
            SNPT = SNPT + 1
 21     CONTINUE
        SM(J) = XM/SNPT
 2    CONTINUE

      DO 25 J=1,4
C------- JJ = 1, 2, 3, 4  for  Sx, SY, SZ, |S|
        X2=0.D0
        SNPT = 0.D0
        DO 26 I=ip1, ip2
            X = COOR(I,J)
            X2 = X2 + (X-SM(J))**2
            SNPT = SNPT + 1
 26     CONTINUE
        SM2(J) = X2/SNPT
 25   CONTINUE

C----- Twiss parameters and emittance
      CALL READC5(NT1,NT2)
      CALL READC9(KKEX,KLET)
      WRITE(6,FMT='(/,'' Particle # '',I6,'' to '',I6,''    @ lmnt # '',
     >I5,''. '',I6,'' points, from pass# '',I6,'' to pass# '',I6,
     >'',  IEX='',I2)') NT1,NT2, LM, NINT(SNPT),ip1,ip2,KKEX
      write(iun,*) ip1, ip2, '        ip1 - ip2, analysis turn# range '
      DO 6 J=1,4
          WRITE(iun,123) 
     >     SM(J),SM2(J),SPK(2*J-1),SPK(2*J),(SPK(2*J-1)+SPK(2*J))/2.d0
     >     ,TXTS(J),'_M,_M2,_MI,_MA (_MI+_MA)/2 '
          WRITE(6,123) 
     >     SM(J),SM2(J),SPK(2*J-1),SPK(2*J),(SPK(2*J-1)+SPK(2*J))/2.d0
     >     ,TXTS(J),'_M,_M2,_MI,_MA (_MI+_MA)/2 '
 123      FORMAT(1P,5(1X,G12.4),'  ',2A)
 6    CONTINUE
        write(6,*) '**************************'
        write(6,*) '**************************'
        write(6,*) '**************************'
        write(6,*) '**************************'
      RETURN 
      END
