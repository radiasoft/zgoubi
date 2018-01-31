C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CDF/ 
     >      IES,IORDRE,LCHA,LIST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IIP(MXL),NB,NOEL
      CHARACTER(80) TA
      COMMON/DONT/ TA(MXL,40)
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      CHARACTER LET
      COMMON/FAISCT/ LET(MXT)
C----- KAR: tagging letter ( 'S'  is reserved for tagging secondary particles 
C            as resulting from decay (keyword 'MCDESINT')
      CHARACTER KAR(41)
      COMMON/KAR/ KAR
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      COMMON/PTICUL/ AAM,Q,G,TO
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      COMMON/SYNCH/ RET(MXT), DPR(MXT),PS
C----- CONVERSION DES COORD. (CM,MRD) -> (M,RD)
      COMMON/UNITS/ UNIT(MXJ)
 
      PARAMETER(MXJ1=MXJ-1)
      DIMENSION DE(5,MXT),IDE(5),JDE(5),P(MXJ)
      EQUIVALENCE (IDE(2),IYMAX),(IDE(3),ITMAX),(IDE(4),IZMAX),
     > (IDE(5),IPMAX),(IDE(1),IMAXD)
      EQUIVALENCE (JDE(2),IY   ),(JDE(3),IT   ),(JDE(4),IZ   ),
     > (JDE(5),IP   ),(JDE(1),ID)
 
      DIMENSION REF(MXJ)
      DIMENSION istp(mxt)
      parameter (mxstp=10000)
      DIMENSION numel(mxstp)

      LOGICAL IDLUNI

      logical exs, ok, gttext

      CHARACTER(200) TXT200,TXT
      CHARACTER(800) TXT8
      parameter (mstp=1000)
      dimension ddr(mxt,mstp)
      dimension yyr(mxt,mstp),ttr(mxt,mstp),zzr(mxt,mstp)
      dimension ppr(mxt,mstp),ssr(mxt,mstp),timr(mxt,mstp)
      dimension ddi(mxt)
      dimension yyi(mxt),tti(mxt),zzi(mxt),ppi(mxt),ssi(mxt),timi(mxt)

      INTEGER DEBSTR, FINSTR
      
      data istp / mxt*0 /

      noel = 1

      write(*,*) ' '
      write(*,*) '----------------------------------------------------'
      write(*,*) 'Now running pgm betaFromPlt... '
      write(*,*) '----------------------------------------------------'
      write(*,*) ' '

      KPa=1 ; KPb=1 ; KPc=1
      call READC2B(KPa,KPb,KPc)
c      if(oksav) call spsav2(kpa,kpb,ksmpl)

      ok = idluni(
     >            nres)
      open(unit=NRES,file='zgoubi.res')

      OK = GTTEXT(6,nres,'''OBJET''',
     >                                txt)
      read(NRES,*) Brho
      read(NRES,*) xobj
      a(noel,10) = xobj
      if(nint(10*xobj) .ne. 51) then
        write(*,*) 
        write(*,*) ' Sorry, I need a zgoubi.res with kobj = 5.1'
        write(*,*)
        stop
      endif
      KOBJ = INT(A(NOEL,10))
      KOBJ2 = NINT(10*A(NOEL,10)) - 10*KOBJ
      write(*,*) 'Brho, Kobj, Kogbj2 = ',Brho,kobj,kobj2 

      read(NRES,*) (a(noel,i),i=20,25)
c      write(*,*) ' (a(noel,i),i=20,25) ',(a(noel,i),i=20,25)
      read(NRES,*) (a(noel,i),i=30,35)
c      write(*,*) ' (a(noel,i),i=30,35) ',(a(noel,i),i=30,35)
      read(NRES,*) (a(noel,i),i=40,49)
c      write(*,*) ' (a(noel,i),i=40,49) ',(a(noel,i),i=40,49)

      close(nres)

C Re-create the object
C with nbref*11 traj. for 1st order matrix calculation
      CALL RAZ(FO,MXJ*MXT)
      CALL OBJ52(KOBJ2)
      nbref = kobj2
      CALL OBJ5
      
      imax = kobj2*11
      write(*,*) 'Kobj, Kobj2 : ',kobj, kobj2
      write(*,*) 'IMAX = ', imax

      write(*,*)  ' '
      write(*,*)  ' Initial object (d,y,t,z,p,s) : '
      write(*,fmt='(6(f12.4,2x))') ((fo(i,j),i=1,6),j=1,imax)

      ok = idluni(
     >             lunout)
      open(unit=lunout,file='betaFromPlt.out')
      WRITE(LUNout,fmt='(a)') '# F0(1,2),  F0(1,1),  F0(3,4),  F0(3,3) '
     >//',  F0(5,6),  F0(5,5)'
     >//',  F0(1,6),  F0(2,6),  F0(3,6),  F0(4,6)'
     >//',  PHY/(2.D0*PI),  PHZ/(2.D0*PI),  SCUM,  NOEL'
     >//',  R(1,6),  R(2,6),  R(3,6),  R(4,6)'
      
c      ok = idluni(
c     >             lunout2)
c      open(unit=lunout2,file='betaFromPlt.out2')
      ok = idluni(
     >             lunin)
      open(unit=lunin,file='zgoubi.plt')

      write(*,*) ' :  Header '
      read(lunin,fmt='(a)') txt200
      write(*,*) txt200(debstr(txt200):80)
      read(lunin,fmt='(a)') txt200
      write(*,*) txt200(debstr(txt200):80)
      read(lunin,fmt='(a)') txt200
      write(*,*) txt200(debstr(txt200):80)
      read(lunin,fmt='(a)') txt200
      write(*,*) txt200(debstr(txt200):80)

      itra = 0
 1    continue
          read(lunin,fmt='(a)',err=10,end=10) txt8 
          read(txt8,*,err=10,end=10) 
     >    KEX,  dbroi, aaa, bbb, ccc, ddd, eee, fff, 
     >          dbro, ppp, qqq, rrr, sss, ttt, uuu,
     >         bta, ds, kart, it
        
        if(it.ne.itra) then
          itra = it
          ddi(itra)=dbroi
          yyi(itra)=aaa
          tti(itra)=bbb
          zzi(itra)=ccc
          ppi(itra)=ddd
          ssi(itra)=eee
          timi(itra)=eee
        endif

        istp(itra) = istp(itra) +1
        read(txt8(701:705),*) numel(istp(itra))

        ddr(itra,istp(itra))=dbro
        yyr(itra,istp(itra))=ppp
        ttr(itra,istp(itra))=qqq
        zzr(itra,istp(itra))=rrr
        ppr(itra,istp(itra))=sss
        ssr(itra,istp(itra))=ttt
        timr(itra,istp(itra))=uuu

c        write(lunout2,fmt='(1p,13(e12.4,1x),2(i9,1x))') 
c     >  dbro,tti(itra),zzi(itra),ppi(itra),
c     >                                  ssi(itra),timi(itra),
c     >  dbro,yyr(itra,istp(itra)),ttr(itra,istp(itra)),
c     >    zzr(itra,istp(itra)),ppr(itra,istp(itra)),
c     >    ssr(itra,istp(itra)),timr(itra,istp(itra)),
c     >  itra,istp(itra)
  
      goto 1

 10   continue
      
c      do i = 1, itra
c        write(*,*) ' IT, istp :',i,istp(i)
c        write(*,*) yi(i),yr(i,1),yr(i,istp(i))
c      enddo

      KOPTIP=1
      do ist = 1, istp(1)
        do i = 1, 11

          fo(1,i) = 1.d0+ddi(i) 
          fo(2,i) = yyi(i) 
          fo(3,i) = tti(i) 
          fo(4,i) = zzi(i) 
          fo(5,i) = ppi(i) 
          fo(6,i) = ssi(i) 
          fo(7,i) = tti(i) 
        
          f(1,i) = 1.d0+ddr(i,ist) 
          f(2,i) = yyr(i,ist) 
          f(3,i) = ttr(i,ist) 
          f(4,i) = zzr(i,ist) 
          f(5,i) = ppr(i,ist) 
          f(6,i) = ssr(i,ist) 
          f(7,i) = ttr(i,ist) 
        enddo

        call scums(f(6,1))
C        call opticc(lunout,noel,KOPTIP)
        call opticc(lunout,numel(ist),KOPTIP)

      enddo

      close(lunout)
c      close(lunout2)
      close(lunin)
      stop

 11   continue
      stop 'Error during read from betaFromPlt.In.'
      end

      SUBROUTINE BLOCK
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

C----- NUMERO DES UNITES LOGIQUES D'ENTREES-SORTIE
      COMMON/CDF/ 
     >      IES,IORDRE,LCHA,LIST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN
 
C----- CONSTANTES
      COMMON/CONST/ CL,PI,DPI,RAD,DEG,QE,AH, CM2M

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
      COMMON/UNITS/ UNIT(MXJ)

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
 
C                    Y     T     Z       P     X,S   dp/p  time
      DATA UNIT / .01D0,.001D0,.01D0, .001D0, .01D0, 1.D0, 1.d0 /

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
      COMMON/UNITS/ UNIT(MXJ) 
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
C--------- read in zgoubi.plt type storage file

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
      SUBROUTINE STORCO(NL,LM,KPS,BINARY,
     >                                   NPASS,energ)
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
      LOGICAL IDLUNI

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

        if(noc.eq.0) energ = yzxb(20)

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
      IF (IDLUNI(ITMP)) open(unit=itmp,file='temp_ipmx')
      write(itmp,*) -1, noc,' max # of turns in .plt file'
      close(itmp)
      GOTO 11

 17   CONTINUE
      WRITE(6,*) ' Required # of passes has been read, ',noc,' points'
      GOTO 11

 11   CONTINUE
      NPASS = NRBLT + 1
      NPTR=NOC
      CALL READC5(KT1,KT2)
C      write(*,*) ' Pgm betaFromPlt, trjctries kt1:kt2 : ',kt1,':',kt2
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
      if(noc.eq.0) STOP ' No more points to analize '

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
      DATA KXYL / 2, 3, 4, 5, 7,   20 /
C                 Y  T  Z  P  time dp/p
C      DATA KXYL / 2, 3, 4, 5, 7,   1 /
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
      FUNCTION GTTEXT(NRES,LUNR,TXT,
     >                              TXTIN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL GTTEXT
      CHARACTER(*) TXT
      CHARACTER(*) TXTIN
      LOGICAL STRCON
      INTEGER DEBSTR, FINSTR

      READ(LUNR,FMT='(A)',ERR=99,END=98)  TXTIN

      DOWHILE(.NOT.STRCON(TXTIN,TXT(DEBSTR(TXT):FINSTR(TXT)),
     >                                                       IS))
        READ(LUNR,FMT='(A)',ERR=99,END=98)  TXTIN
      ENDDO

      GTTEXT = .TRUE.
      GOTO 10

 99   CONTINUE
      IF(NRES.GT.0) THEN
        WRITE(NRES,*) ' Pgm zgoubi/gttext - ERR upon read.' 
        WRITE(NRES,*) ' Text was : ',TXT(DEBSTR(TXT):FINSTR(TXT))
      ENDIF
      GTTEXT = .FALSE.
      GOTO 10

 98   CONTINUE
      IF(NRES.GT.0) THEN
        WRITE(NRES,*) ' Pgm zgoubi/gttext - EOF upon read.' 
        WRITE(NRES,*) ' Text was : ',TXT(DEBSTR(TXT):FINSTR(TXT))
      ENDIF
      GTTEXT = .FALSE.
      GOTO 10

 10   CONTINUE
      RETURN
      END
      SUBROUTINE REFER(IO,IORD,IFOC,IT1,IT2,IT3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     **********************************************
C     SETS THE REFERENCE FRAME
C     **********************************************
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      COMMON/CONST2/ ZERO, UN
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
 
      SAVE XI,YI,ALE,PATHL
      DATA PATHL / -9999.D0 / 
      GOTO (1,2) IO
 
 1    CONTINUE
C----- GOES TO NEW REFERENCE FRAME + computes new coordinates there
 
      IF(NRES .GT. 0) THEN
        WRITE(NRES,FMT='(/,1P,''  Reference, absolute '',
     >  ''(part #'',I4,'')  : '',
     >   1P,6E13.5,1X,E13.5)') IT1,-1.D0+F(1,IT1),(F(J,IT1),J=2,7)
      ENDIF 

      IF    (IFOC .EQ. 0) THEN
        XI=0.D0
        YI=F(2,IT1)
        ALE=F(3,IT1)*.001D0
      ELSEIF(IFOC .EQ. 1) THEN
C------- RECHERCHE DES COORDONNEES DU POINT DE FOCALISATION
        IF    (IORD .EQ. 1) THEN
           CALL FOCAL1(IT1,IT2,IT3,XI,YI)
        ELSEIF(IORD .EQ. 2) THEN
           CALL FOCAL1(IT1,IT2,IT3,XI,YI)
        ENDIF
C------- LE SYSTEME DE REFERENCE POUR LE CALCUL DES COEFFICIENTS
C        DE TRANSFERT S'APPUIE SUR LA DIRECTION DE LA TRAJECTOIRE #1
        ALE=F(3,IT1)*.001D0
      ENDIF
      IF(IFOC.LE.1) THEN
        DO 8 I=1,IMAX
           CALL INITRA(I)
           CALL CHAREF(.FALSE.,XI,YI,ALE)
           CALL MAJTRA(I)
 8      CONTINUE
        IF(NRES .GT. 0) WRITE(NRES,100) XI,YI,ALE*DEG,ALE
 100    FORMAT(/,10X,' Frame for MATRIX calculation moved by :'
     >        ,/,10X,'  XC =',F9.3,' cm , YC =',F9.3,' cm ,   A ='
     >        ,F9.5,' deg  ( =',F9.6,' rad )',/)
      ENDIF
      PATHL = F(6,IT1)
      IF(NRES .GT. 0) THEN
        WRITE(NRES,FMT='(/,1P,''  Reference particle '',
     >  ''(#'',I4,''), path length :'',G16.8,'' cm'', 
     >  ''  relative momentum : '',G14.6)') IT1, F(6,IT1), F(1,IT1)
      ENDIF 
      RETURN
 
 2    CONTINUE
C----- COMES BACK TO OLD FRAME + OLD COORDINATES
 
      IF(IFOC.LE.1) THEN
         DO 81 I=1,IMAX
            CALL INITRA(I)
            CALL CHAREF(.FALSE.,ZERO,ZERO,-ALE)
            CALL CHAREF(.FALSE.,-XI,-YI,ZERO)
            CALL MAJTRA(I)
81       CONTINUE
      ENDIF
 
      RETURN

      ENTRY REFER1(
     >             PATHLO)
      PATHLO = PATHL
      RETURN
      
      END
      SUBROUTINE OBJ5
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     **************************************
C     CONSTITUTION DE L'OBJET INITIAL KOBJ=5
C     **************************************
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IIP(MXL),NB,NOEL
      CHARACTER*80 TA
      COMMON/DONT/ TA(MXL,40)
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      CHARACTER LET
      COMMON/FAISCT/ LET(MXT)
C----- KAR: LETTRES AFFECTEES AUX TRAJECTOIRES ( 'S'  EST RESERVEE
C      POUR ETIQUETER LES PARTICULES SECONDAIRES -OPTION 'MCDESINT')
      CHARACTER  KAR(41)
      COMMON/KAR/ KAR
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
 
      PARAMETER(MXJ1=MXJ-1)

      DIMENSION IDE(5),JDE(5),P(MXJ)
      EQUIVALENCE (IDE(2),IYMAX),(IDE(3),ITMAX),(IDE(4),IZMAX),
     > (IDE(5),IPMAX),(IDE(1),IMAXD)
      EQUIVALENCE (JDE(2),IY   ),(JDE(3),IT   ),(JDE(4),IZ   ),
     > (JDE(5),IP   ),(JDE(1),ID)
 
      PARAMETER(MXREF=99)
      DIMENSION REF(MXJ,MXREF)
      DIMENSION FI(6,6)

      SAVE NBREF
      DATA NBREF / 1 /
 
      IMAX=11 * NBREF
      IDMAX=1
      IMAXT=IMAX/IDMAX
      P(2) = A(NOEL,20)
      P(3) = A(NOEL,21)
      P(4) = A(NOEL,22)
      P(5) = A(NOEL,23)
      P(6) = A(NOEL,24)
      P(1) = A(NOEL,25)
           
      IREF = 0

 1    CONTINUE
        IREF = IREF + 1
        IREF1 = IREF-1
        K = 30 + 10 * IREF1
        DO 52 J = 2,MXJ1
          REF(J,IREF) = A(NOEL,K)
          K = K + 1
 52     CONTINUE
        REF(1,IREF) = A(NOEL,K)

        I = 11 * IREF1 
        DO 53 J=2,5
          I=I+2
          DX = P(J)
          FO(J,I   ) = DX
          FO(J,I+1 ) = - DX
 53     CONTINUE
        FO(1,11*IREF-1) = P(1)
        FO(1,11*IREF) = - P(1)   

        IKAR = 1
        DO 51 I=11*IREF1+1,11*IREF
          IEX(I) = 1
          IREP(I) = I
          LET(I) = KAR(IKAR)
          IKAR = IKAR+1
          IF(IKAR .GT. 41) IKAR = 1
          DO 51 J = 1, 6
            FO(J,I) = FO(J,I) + REF(J,IREF)
            F(J,I) = FO(J,I)
 51     CONTINUE

      IF(IREF.LT.NBREF) GOTO 1

C----- Alpha_y, beta_y, *_z, *_d
      FI(1,1) = A(NOEL,41)  
      IF(FI(1,1) .EQ. 0.D0) FI(1,1) = 1.D0    
      FI(2,1) = A(NOEL,40)      
      FI(1,2) = FI(2,1)
      FI(2,2) = (1.D0+FI(2,1)*FI(2,1))/FI(1,1)
      FI(3,3) = A(NOEL,43)  
      IF(FI(3,3) .EQ. 0.D0) FI(3,3) = 1.D0    
      FI(4,3) = A(NOEL,42)      
      FI(3,4) = FI(4,3)
      FI(4,4) = (1.D0+FI(4,3)*FI(4,3))/FI(3,3)
      FI(5,5) = A(NOEL,45)  
      IF(FI(5,5) .EQ. 0.D0) FI(5,5) = 1.D0    
      FI(6,5) = A(NOEL,44)      
      FI(5,6) = FI(6,5)
      FI(6,6) = (1.D0+FI(6,5)*FI(6,5))/FI(5,5)        
C Dy, Dy', Dz, Dz'
      FI(1,6) = A(NOEL,46)      
      FI(6,1) = FI(1,6)
      FI(2,6) = A(NOEL,47)      
      FI(6,2) = FI(2,6)
      FI(3,6) = A(NOEL,48)      
      FI(6,3) = FI(3,6)
      FI(4,6) = A(NOEL,49)      
      FI(6,4) = FI(4,6)
      CALL BEAMA1(FI)

      IF(NRES.GT.0) THEN
        WRITE(NRES,100) KOBJ,IMAX
  100   FORMAT(/,41X,'CALCUL  DES  TRAJECTOIRES',//,30X,'OBJET  (',I1,
     >  ')  FORME  DE ',I6,' POINTS ',//)
        WRITE(NRES,FMT='(/,T33,''Y (cm)'',T48,''T (mrd)'',T62,
     >  ''Z (cm)'',T76,''P (mrd)'',T90,''S (cm)'',T103,'' dp/p '')')
       WRITE(NRES,FMT='(14X,'' Sampling : '',T30, 
     >  5(4X,G10.2),4X,G12.4)') (P(J), J=2,6), P(1)
        IREF = 1
        WRITE(NRES,FMT='(2X,  ''Reference trajectory #'',I1,'' : '',T30,
     >                         5(4X,G10.2),4X,G12.4)') 
     >  IREF,(REF(J,IREF), J=2,6), REF(1,IREF)
        DO 20 IREF=2, NBREF
          WRITE(NRES,FMT='(2X,''                     #'',I1,'' : '',T30,
     >                           5(4X,G10.2),4X,G12.4)') 
     >    IREF,(REF(J,IREF), J=2,6), REF(1,IREF)
 20     CONTINUE
      ENDIF

      RETURN

      ENTRY OBJ51(
     >            NBREFO)
      NBREFO = NBREF
      RETURN

      ENTRY OBJ52(KOBJ2)
      IF(KOBJ2 .EQ. 0) THEN
        NBREF = 1
      ELSEIF(KOBJ2 .EQ. 1) THEN
        NBREF = 1
      ELSEIF(KOBJ2 .GE. 2 .AND. KOBJ2 .LE.MXREF) THEN
        NBREF = KOBJ2
      ELSE
        CALL ENDJOB('OBJ5, wrong value KOBJ2 in OBJET. Max is ',MXREF)
      ENDIF 
      RETURN
      END
      SUBROUTINE MAT1(R,T,IT1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(6,*) , T(6,6,*)
C     -------------------------------------------------
C     OPTION  IORD = 1 :
C       MATRICES ORDRE 1 ET 2, SERIES DE TAYLOR ORDRE 3
C     -------------------------------------------------
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
 
      CALL RAZ(R,6*6)
      R(5,5) = 1.D0
      R(6,6) = 1.D0
      CALL RAZ(T,5*6*6)
      S1 = 2.D0 * F(6,IT1)
 
C..............................................
C             Y     T     Z     P     L     D
C
C        Y   R11   R12   R13   R14   R15   R16
C        T   R21   R22   R23   R24   R25   R26
C        Z   R31   R32   R33   R34   R35   R36
C        P   R41   R42   R43   R44   R45   R46
C        L   R51   R52   R53   R54   R55   R56
C        D   R61   R62   R63   R64   R65   R66
C..............................................

      I10 = IT1+9
      I11 = IT1+10
      DP = ( FO(1,I10) - FO(1,I11) ) / (.5D0*( FO(1,I10) + FO(1,I11) ) )
      DP2 = DP*DP
      DO 11 J=2,5
        R(J-1,6)  = (F(J,I10) - F(J,I11)) /DP
        DO 11 I=1,4
          I2 = 2*I +IT1-1
          I3 = I2+1 
          UO = FO(I+1,I2)-FO(I+1,I3)
          R(J-1,I)  = (F(J,I2) - F(J,I3)) / UO
          IF(J .EQ. 5) THEN
            R(5,I)  = ( F(6,I2) - F(6,I3) ) / UO
          ENDIF
C          write(*,*) j,I2,I3, F(J,I2), F(J,I3)
 11   CONTINUE
      R(5,6)  = ( F(6,I10) - F(6,I11) ) / DP
 
      IF(IMAX.NE.13) RETURN
C     ... Compute Ri5, i=1,5. 
      IF(FO(6,13)-FO(6,12) .NE. 0.D0) THEN
        DL=FO(6,13)-FO(6,12)
        R(1,5)  = ( F(2,13) - F(2,12) ) / DL
        R(2,5)  = ( F(3,13) - F(3,12) ) / DL
        R(3,5)  = ( F(4,13) - F(4,12) ) / DL
        R(4,5)  = ( F(5,13) - F(5,12) ) / DL
        R(5,5)  = ( F(6,13) - F(6,12) ) / DL
      ENDIF
      RETURN
      END
      SUBROUTINE MKSA(IORD,R,T,T3,T4)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(6,*) , T(6,6,*)
      DIMENSION  T3(5,*), T4(5,*)

      INCLUDE "MAXCOO.H"
      COMMON/UNITS/ UNIT(MXJ)
 
      DO 1 IA=1,6
        DO 1 IB=1,6
          R(IA,IB) = R(IA,IB)*UNIT(IA)/UNIT(IB)
C          DO 1 IC=1,IB
          DO 1 IC=1,6
            T(IA,IC,IB) = T(IA,IC,IB)*UNIT(IA)/UNIT(IC)/UNIT(IB)
 1    CONTINUE
 
      IF(IORD .EQ. 1) RETURN
 
      DO 2 IA=1,5
       DO 2 IB=1,6
        T3(IA,IB) = T3(IA,IB)*UNIT(IA)/UNIT(IB)**3
        T4(IA,IB) = T4(IA,IB)*UNIT(IA)/UNIT(IB)**4
 2    CONTINUE
 
      RETURN
      END
      SUBROUTINE MATIMP(R) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(6,*) , T(6,6,*)
      DIMENSION  T3(5,*)

      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU

      LOGICAL KWRI, KWRMAT, IDLUNI
      CHARACTER FNAME*17
      LOGICAL EXS, OPN
      SAVE KWRMAT

      DATA FNAME / 'zgoubi.MATRIX.out' /
  
      DETY=R(1,1)*R(2,2)-R(1,2)*R(2,1)
      DETZ=R(3,3)*R(4,4)-R(3,4)*R(4,3)
      RIJ = R(2,2)
      IF(RIJ .EQ. 0.D0) RIJ = 1.D-10
      SFH = - R(1,2)/RIJ
      RIJ = R(4,4)
      IF(RIJ .EQ. 0.D0) RIJ = 1.D-10
      SFZ = - R(3,4)/RIJ
 
      I=1
      IF(NRES.GT.0) THEN
        WRITE(NRES,103) I
 103    FORMAT(//,18X,'TRANSFER  MATRIX  ORDRE',I3,'  (MKSA units)',/)
        WRITE(NRES,104) (( R(IA,IB) , IB=1,6) , IA=1,6)
 104    FORMAT(6X,1P,6G16.6)
        WRITE(NRES,112) DETY-1.D0,DETZ-1.D0
112     FORMAT(/,10X,'DetY-1 = ',F18.10,',',4X,'DetZ-1 = ',F18.10)
        WRITE(NRES,FMT='(/,10X,''R12=0 at '',G12.4,'' m, '',7X, 
     >                       ''R34=0 at '',G12.4,'' m'')') SFH,SFZ
      ENDIF

      CALL SYMPL(R)

      IF(KWRMAT) THEN
        IF(IDLUNI(
     >            LNWRT)) THEN
          INQUIRE(FILE=FNAME,EXIST=EXS,OPENED=OPN,IOSTAT=I)
          IF(OPN) THEN
            CLOSE(LNWRT)
          ENDIF
          IF(EXS) THEN
            OPEN(UNIT=LNWRT, FILE=FNAME, status='OLD',ERR=96)
            CALL GO2END(LNWRT)
          ELSE
            OPEN(UNIT=LNWRT, FILE=FNAME, status='NEW',ERR=96)
            WRITE(LNWRT,*) '% R11 R12 R13 R14 R21 R22 R23 ... R43 R44'
            WRITE(LNWRT,*) '%  '
          ENDIF
c          write(lnwrt,*)'%  transport coefficients',
c     >              ' ((R(IA,IB),IB=1,4),IA=1,4)'
C This will stack results from stacked jobs, or will stack with earlier results
        ELSE
          GOTO 95
        ENDIF
        WRITE(LNWRT,FMT='(1P,16(1X,E12.4))') ((R(IA,IB),IB=1,4),IA=1,4)
        CLOSE(LNWRT)
        KWRMAT = .FALSE.
      ENDIF

      RETURN

      ENTRY MATIM2(R,T,T3)

C MODIFIED, FM, 04/97
C       ** CHANGE MATRICE TRIANGULAIRE EN CARREE SYMMETRIQUE/DIAG
        DO 11 IA=1,6
          DO 11 IB=1,6
            IC1=IB+1
            DO 11 IC=IC1,6
              T(IA,IC,IB)=T(IA,IB,IC)
 11     CONTINUE

      I=2
      IF(NRES.GT.0) THEN
        WRITE(NRES,103) I
        DO 16 IA=1,6
          IF(IA.GT.1) WRITE(NRES,107)
 107      FORMAT(/)
          DO 16 IB=1,6
C          WRITE(NRES,108) ( IA,IC,IB, T(IA,IC,IB)  , IC=1,IB )
C MODIFIED, FM, 04/97
            WRITE(NRES,108) ( IA,IC,IB, T(IA,IC,IB)  , IC=1,6 )
 108        FORMAT( 6(I4,I2,I1,1P,G11.3) )
 16     CONTINUE
      ENDIF

      CALL SYMPL2(R,T)
 
      IF(NRES.GT.0) THEN
        WRITE(NRES,123) T3(1,1),T3(1,2),T3(1,3),T3(1,4)
 123    FORMAT(//,15X,'COEFFICIENTS  D''ORDRE  SUPERIEUR  ( MKSA ):'
     >  ,//,10X,' Y/Y3   ',5X,1P,G14.5
     >  , /,10X,' Y/T3   ',5X,   G14.5
     >  , /,10X,' Y/Z3   ',5X,   G14.5
     >  , /,10X,' Y/P3   ',5X,   G14.5,/)
        WRITE(NRES,124) T3(2,1),T3(2,2),T3(2,3),T3(2,4)
 124    FORMAT(
     >     10X,' T/Y3   ',5X,1P,G14.5
     >  ,/,10X,' T/T3   ',5X,   G14.5
     >  ,/,10X,' T/Z3   ',5X,   G14.5
     >  ,/,10X,' T/P3   ',5X,   G14.5,/)
        WRITE(NRES,125) T3(3,1),T3(3,2),T3(3,3),T3(3,4)
 125    FORMAT(
     >     10X,' Z/Y3   ',5X,1P,G14.5
     >  ,/,10X,' Z/T3   ',5X,   G14.5
     >  ,/,10X,' Z/Z3   ',5X,   G14.5
     >  ,/,10X,' Z/P3   ',5X,   G14.5,/)
        WRITE(NRES,126) T3(4,1),T3(4,2),T3(4,3),T3(4,4)
 126    FORMAT(
     >     10X,' P/Y3   ',5X,1P,G14.5
     >  ,/,10X,' P/T3   ',5X,   G14.5
     >  ,/,10X,' P/Z3   ',5X,   G14.5
     >  ,/,10X,' P/P3   ',5X,   G14.5)
 
        WRITE(NRES,101) IEX(1),(F(J,1),J=1,7)
  101   FORMAT(' TRAJ 1 IEX,D,Y,T,Z,P,S,time :',I3,1P,5G12.4,2G17.5)
      ENDIF

      RETURN

      ENTRY MATIM6(KWRI)
      KWRMAT = KWRI
      RETURN

 95   CALL ENDJOB('ERROR : no free unit # for '//FNAME,-99)
 96   KWRMAT = .FALSE.
      CALL ENDJOB('ERROR upon open  old  file '//FNAME,-99)
      RETURN
      END
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
C  Upton, NY, 11973
C  -------
      SUBROUTINE INITRA(I)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ***COORDONNEES DE LA TRAJECTOIRE I ,INITIALISATION
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      COMMON/PTICUL/ AM,Q,G,TO
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      COMMON/TRAJ/ Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT
      IT = I
      Y=F(2,I)
      T=F(3,I)*0.001D0
      Z=F(4,I)
      P=F(5,I)*0.001D0
      DP=F(1,I)
      QT = AMQ(2,I)
      QBR = Q*BORO*DP
      BRI = QT/QBR
      KEX=IEX(I)
      SAR= F(6,I)
      AMT = AMQ(1,I)
C----- AMQ(2,I) = Q/QE
      TAR = F(7,I)   *1.D5
      RETURN
      END
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
C  Brookhaven National Laboratory                    és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE FOCAL1(IRF,MX1,MX2,XI,YI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     -----------------------------------------------
C     CALCULE LA POSITION XI DU POINT DE FOCALISATION
C     APPELE PAR SPGM MATRIX.
C     -----------------------------------------------
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
 
      DIMENSION TTI(MXT)
 
      TTI(1)=TAN(F(3,IRF)*.001D0)
      STY=TTI(1)*F(2,IRF)
      ST2=TTI(1)*TTI(1)
      ST=TTI(1)
      SY=F(2,IRF)
      MXI=1
      DO 1 I=MX1,MX2
         IF(I .EQ. IRF) GOTO 1
         MXI=MXI+1
         TTI(I)=TAN(F(3,I)*1.D-3)
         STY=STY+TTI(I)*F(2,I)
         ST2=ST2+TTI(I)*TTI(I)
         ST=ST+TTI(I)
         SY=SY+F(2,I)
 1    CONTINUE
      XI=-(STY-(ST*SY)/MXI)/(ST2-(ST*ST)/MXI)
      YI =F(2,IRF)+XI*TTI(1)
 
      RETURN
      END
      SUBROUTINE CHAREF(EVNT,XC,YC,A)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL EVNT
C     -----------------------------------------------
C     CHANGEMENT DE REFERENCE PARTICULE PAR PARTICULE
C     -----------------------------------------------
      COMMON/AIM/ AE,AT,AS,RM,XI,XF,EN,EB1,EB2,EG1,EG2
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      COMMON/CONST2/ ZERO, UN
      INCLUDE "MAXTRA.H"
      COMMON/DESIN/ FDES(7,MXT),IFDES,KINFO,IRSAR,IRTET,IRPHI,NDES
     >,AMS,AMP,AM3,TDVM,TETPHI(2,MXT)
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      COMMON/GASC/ AI, DEN, KGA
      COMMON/MARK/ KART,KALC,KERK,KUASEX
      COMMON/PTICUL/ AM,Q,G,TO
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)
      COMMON/TRAJ/ Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT
      PARAMETER(I0=0)
      YO=Y
      Y=((Y-YC)*COS(T)+XC*SIN(T))/COS(T-A)
      T=T-A 
      XL=XC-Y*SIN(A)
      YL=YC-YO+Y*COS(A)
      DL=SQRT(XL*XL+YL*YL)
      DL=SIGN(DL,XL)
      DS = DL/COS(P)
      SAR= SAR+DS
      Z=Z+DL*TAN(P)
      QBRO = QBR*CL9
      DTAR = DS / (QBRO/SQRT(QBRO*QBRO+AMT*AMT)*CL9)
      TAR = TAR + DTAR
      RETURN 
      END
      SUBROUTINE MAJTRA(I)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C  Update coordinates of trajectory # I
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      COMMON/PTICUL/ AM,Q,G,TO
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      COMMON/TRAJ/ Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT

      F(2,I)=Y
      F(3,I)=T*1000.D0
      F(4,I)=Z
      F(5,I)=P*1000.D0
      DP=QBR/(Q*BORO)
      F(1,I)=DP
      IEX(I)=KEX
      F(6,I)= SAR
      F(7,I)= TAR    *1.D-5
      AMQ(1,I) = AMT
      AMQ(2,I) = QT
      RETURN
      END
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
C  Brookhaven National Laboratory                    és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE SYMPL(R)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(6,*) , T(6,6,*)

      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG

      DIMENSION S(14)

C----- First order. Cours INSTN G. Leleux, p5-14.
C      conditions equivalent to T(R).S.R = S
      S(1) = R(1,1)*R(2,2)-R(1,2)*R(2,1)+R(3,1)*R(4,2)-R(3,2)*R(4,1)-1D0
      S(2) = R(1,3)*R(2,4)-R(1,4)*R(2,3)+R(3,3)*R(4,4)-R(3,4)*R(4,3)-1D0
      S(3)= R(2,4)*R(1,1)-R(1,4)*R(2,1) + R(4,4)*R(3,1)-R(3,4)*R(4,1)
      S(4)= R(2,4)*R(1,2)-R(1,4)*R(2,2) + R(4,4)*R(3,2)-R(3,4)*R(4,2)
      S(5)=-R(2,3)*R(1,1)+R(1,3)*R(2,1) - R(4,3)*R(3,1)+R(3,3)*R(4,1)
      S(6)=-R(2,3)*R(1,2)+R(1,3)*R(2,2) - R(4,3)*R(3,2)+R(3,3)*R(4,2)
      
      IF(NRES.GT.0) WRITE(NRES,100) (S(J),J=1,6)
 100  FORMAT(/,5X,
     >' First order symplectic conditions (expected values = 0) :',
     >/,5X,1P,6G14.4)
      RETURN

      ENTRY SYMPL2(R,T)

C----- Second order. These J.L. Laclare, p18.
      S(1) =    R(1,1)*(T(2,1,2)+T(2,2,1)) + 2.D0*R(2,2)*T(1,1,1) 
     > - 2.D0*R(1,2)*T(2,1,1) -    R(2,1)*(T(1,1,2)+T(1,2,1))

      S(2) = 2.D0*R(1,1)*T(2,2,2) +    R(2,2)*(T(1,1,2)+T(1,2,1))
     > -    R(1,2)*(T(2,1,2)+T(2,2,1)) - 2.D0*R(2,1)*T(1,2,2)

      S(3) =    R(1,1)*(T(2,2,6)+T(2,6,2)) +  R(2,2)*(T(1,1,6)+T(1,6,1))
     > -    R(1,2)*(T(2,1,6)+T(2,6,1)) -    R(2,1)*(T(1,2,6)+T(1,6,2))

      S(4) = 2.D0*R(1,1)*T(2,3,3) +    R(4,3)*(T(3,1,3)+T(3,3,1)) 
     > -    R(3,3)*(T(4,1,3)+T(4,3,1)) - 2.D0*R(2,1)*T(1,3,3)

      S(5) =    R(1,1)*(T(2,3,4)+T(2,4,3)) +  R(4,3)*(T(3,1,4)+T(3,4,1))
     > -    R(3,3)*(T(4,1,4)+T(4,4,1)) -    R(2,1)*(T(1,3,4)+T(1,4,3))

      S(6) =    R(1,1)*(T(2,3,4)+T(2,4,3)) +  R(4,4)*(T(3,1,3)+T(3,3,1))
     > -    R(3,4)*(T(4,1,3)+T(4,3,1)) -    R(2,1)*(T(1,3,4)+T(1,4,3))

      S(7) = 2.D0*R(1,1)*T(2,4,4) +    R(4,4)*(T(3,1,4)+T(3,4,1))
     > -    R(3,4)*(T(4,1,4)+T(4,4,1)) - 2.D0*R(2,1)*T(1,4,4)

      S(8) =    R(1,2)*(T(2,3,4)+T(2,4,3)) +  R(4,4)*(T(3,2,3)+T(3,3,2))
     > -    R(3,4)*(T(4,2,3)+T(4,3,2)) -    R(2,2)*(T(1,3,4)+T(1,4,3))

      S(9) = 2.D0*R(1,2)*T(2,4,4) +    R(4,4)*(T(3,2,4)+T(3,4,2)) 
     > -    R(3,4)*(T(4,2,4)+T(4,4,2)) - 2.D0*R(2,2)*T(1,4,4)

C FM 07/97      S(10)= 2.D0*R(1,2)*T(2,3,3) +    R(2,1)*(T(3,2,3)+T(3,3,2)) 
      S(10)= 2.D0*R(1,2)*T(2,3,3) +    R(4,3)*(T(3,2,3)+T(3,3,2)) 
     > -    R(3,3)*(T(4,2,3)+T(4,3,2)) - 2.D0*R(2,2)*T(1,3,3)

C FM 07/97       S(11)= 2.D0*R(1,2)*(T(2,3,4)+T(2,4,3))+ R(2,1)*(T(3,2,4)+T(3,4,2)) 
      S(11)=    R(1,2)*(T(2,3,4)+T(2,4,3))+ R(4,3)*(T(3,2,4)+T(3,4,2)) 
     > -    R(3,3)*(T(4,2,4)+T(4,4,2)) -    R(2,2)*(T(1,3,4)+T(1,4,3))

      S(12)=    R(3,3)*(T(4,1,4)+T(4,4,1)) +  R(4,4)*(T(3,1,3)+T(3,3,1)) 
     > -    R(3,4)*(T(4,1,3)+T(4,3,1)) -    R(4,3)*(T(3,1,4)+T(3,4,1))

      S(13)=    R(3,3)*(T(4,2,4)+T(4,4,2)) +  R(4,4)*(T(3,2,3)+T(3,3,2)) 
     > -    R(3,4)*(T(4,2,3)+T(4,3,2)) -    R(4,3)*(T(3,2,4)+T(3,4,2))

      S(14)=    R(3,3)*(T(4,4,6)+T(4,6,4)) +  R(4,4)*(T(3,3,6)+T(3,6,3)) 
     > -    R(3,4)*(T(4,3,6)+T(4,6,3)) -    R(4,3)*(T(3,4,6)+T(3,6,4))

      IF(NRES.GT.0) WRITE(NRES,101) (S(J),J=1,14)
 101  FORMAT(/,5X,
     >' Second order symplectic conditions (expected values = 0) :',
     >/,1P,2(/,5X,7G14.4) )

      RETURN
      END
C  ZGOUBI, a program for computing the trajectories of charged particles
C  in electric and magnetic fields
C  Copyright (C) 1988-2007  FranÃ§ois Meot
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
C  FranÃ§ois Meot <fmeot@bnl.gov>
C  BNL
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE GO2END(LUN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ------------------------------------------------------
C     Extract substrings #1 up to #MSS, out of string STR. 
C     Strings are assumed spaced by (at least) one blank. 
C     They are saved in  array STRA, and their total number 
C     (possibly < mss) is NST.
C     ------------------------------------------------------
      CHARACTER TXT80*80

 1    CONTINUE
        READ(LUN,FMT='(A)',ERR=99,END=99) TXT80
        GOTO 1

 99   RETURN
      END
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
C  Brookhaven National Laboratory                    és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE ENDJOB(TXT,II)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*(*) TXT
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      NRES = ABS(NRES)
      IF(II.EQ.-99) THEN
        WRITE(   6,FMT='(/,A)') TXT
        WRITE(NRES,FMT='(/,A)') TXT
      ELSE
C        WRITE(   6,FMT='(/,A,1X,I8)') TXT,II
C        WRITE(NRES,FMT='(/,A,1X,I8)') TXT,II
        WRITE(   6,*) TXT,II
        WRITE(NRES,*) TXT,II
      ENDIF
      WRITE(NRES,FMT='(/,''End of job !'',//,''  '')')
      STOP
      END
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
C  Upton, NY, 11973
C  -------
      SUBROUTINE BEAMAT(R,
     >                    F0,PHY,PHZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(6,*),F0(6,*)
      COMMON/BEAM/ FI(6,6)
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M

      DIMENSION FII(6,6)

      F0(1,1) = R(1,1)*R(1,1)*FI(1,1)-2.D0*R(1,1)*R(1,2)*FI(2,1) +
     >  R(1,2)*R(1,2)*FI(2,2)
      F0(1,2) = -(  R(1,1)*R(2,1)*FI(1,1) -
     >  ( 1.D0 + 2.D0*R(1,2)*R(2,1) )*FI(2,1) +
     >  R(1,2)*R(2,2)*FI(2,2)  )
      F0(2,1) = F0(1,2) 
      F0(2,2) = R(2,1)*R(2,1)*FI(1,1) - 2.D0*R(2,2)*R(2,1)*FI(2,1) +
     >  R(2,2)*R(2,2)*FI(2,2)
      F0(3,3) = R(3,3)*R(3,3)*FI(3,3)-2.D0*R(3,3)*R(3,4)*FI(4,3) +
     >  R(3,4)*R(3,4)*FI(4,4)
      F0(3,4) = -(  R(3,3)*R(4,3)*FI(3,3) -
     >  ( 1.D0 + 2.D0*R(3,4)*R(4,3) )*FI(4,3) +
     >  R(3,4)*R(4,4)*FI(4,4)  )
      F0(4,3) = F0(3,4) 
      F0(4,4) = R(4,3)*R(4,3)*FI(3,3) - 2.D0*R(4,4)*R(4,3)*FI(4,3) +
     >  R(4,4)*R(4,4)*FI(4,4)
      F0(1,6) = R(1,1)*FI(1,6) +R(1,2)*FI(2,6) +R(1,3)*FI(3,6) +
     >  R(1,4)*FI(4,6) +R(1,5)*FI(5,6) +R(1,6)*FI(6,6) 
      F0(2,6) = R(2,1)*FI(1,6) +R(2,2)*FI(2,6) +R(2,3)*FI(3,6) + 
     >  R(2,4)*FI(4,6) +R(2,5)*FI(5,6) +R(2,6)*FI(6,6) 
      F0(3,6) = R(3,1)*FI(1,6) +R(3,2)*FI(2,6) +R(3,3)*FI(3,6) + 
     >  R(3,4)*FI(4,6) +R(3,5)*FI(5,6) +R(3,6)*FI(6,6) 
      F0(4,6) = R(4,1)*FI(1,6) +R(4,2)*FI(2,6) +R(4,3)*FI(3,6) + 
     >  R(4,4)*FI(4,6) +R(4,5)*FI(5,6) +R(4,6)*FI(6,6) 

C Betatron phase advance
c        PHY = atan2(R(1,2) , ( R(1,1)*F0(1,1) - R(1,2)*F0(1,2)))
c        PHZ = atan2(R(3,4) , ( R(3,3)*F0(3,3) - R(3,4)*F0(3,4)))
      PHY = atan2(R(1,2) , ( R(1,1)*FI(1,1) - R(1,2)*FI(1,2)))
       IF(PHY.LT.0.D0) PHY = 2.D0*PI + PHY
      PHZ = atan2(R(3,4) , ( R(3,3)*FI(3,3) - R(3,4)*FI(3,4)))
       IF(PHZ.LT.0.D0) PHZ = 2.D0*PI + PHZ

      RETURN

      ENTRY BEAMA1(FII)
      DO 2 J=1,6
        DO 2 I=1,6 
          FI(I,J) = FII(I,J)
 2    CONTINUE
      RETURN
      END
      SUBROUTINE OPTICC(LNOPTI,NOEL,KOPTIP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(6,6),F0(6,6)

      IORD = 1
      IFOC = 0
      KWR = 0

      CALL MATRIC(IORD,IFOC,KWR)
      CALL MATRI1(
     >            R)
      CALL BEAMAT(R, 
     >              F0,PHY,PHZ)
      CALL BEAIMP(F0,PHY,PHZ)

c Store in betaFromPlt.out
      IF(KOPTIP.EQ.1) CALL OPTIMP(LNOPTI,NOEL,F0,PHY,PHZ,R)

      RETURN
      END
      SUBROUTINE OPTIMP(LUN,NOEL,F0,PHY,PHZ,R) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(6,6),F0(6,6)
      PARAMETER (PI = 4.D0*ATAN(1.D0))

      CALL SCUMR(
     >            XL,SCUM,TCUM)

 104  FORMAT(1P,13(E13.5,1X),1X,I5,4(E13.5,1X))
      WRITE(LUN,104) F0(1,2), F0(1,1), F0(3,4), F0(3,3) 
     >, F0(5,6), F0(5,5)
     >, F0(1,6), F0(2,6), F0(3,6), F0(4,6)
     >, PHY/(2.D0*PI), PHZ/(2.D0*PI),SCUM,NOEL
     >, R(1,6),R(2,6),R(3,6),R(4,6)
      RETURN
      END
      SUBROUTINE SCUMUL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      COMMON/PTICUL/ AM,Q,G,TO
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      INCLUDE "MAXCOO.H"
      COMMON/UNITS/ UNIT(MXJ)

      SAVE XL, SCUM, TCUM

      DATA XL, SCUM, TCUM / 0.D0, 0.D0, 0.D0 / 

      ENTRY SCUMW(XLI)
      XL=XLI
      SCUM = SCUM + XL
C----- Compute cumulative time. Default is for proton. 
      AAM = AM
      QQ = Q
      IF(AAM .EQ. 0.D0) AAM = AMPROT
      IF(QQ .EQ. 0.D0) QQ = QE
C      PREF = (BORO*DPREF) *CL*1.D-9*QQ/QE
C      PREF = (BORO*DPREF) *CL*1.D-9*QQ
      PREF = (BORO*DPREF) *CL9*QQ
      BTA = PREF / SQRT(PREF*PREF+AAM*AAM)
C----- XL is in centimeters, TCUM is seconds
      DT = XL / (CL*BTA) *UNIT(5)
      TCUM = TCUM + DT
      RETURN

      ENTRY SCUMR(
     >            XLO,SCUMO,TCUMO)
      XLO=XL
      SCUMO = SCUM
      TCUMO = TCUM
      RETURN

      ENTRY SCUMS(SETS)
      SCUM = SETS
      RETURN

      ENTRY SCUMT(SETT)
      TCUM = SETT
      RETURN
      END
      SUBROUTINE MATRIC(JORD,JFOC,KWR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ------------------------------------
C     Compute transfer matrix coefficients
C     ------------------------------------
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
C      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
C      COMMON/DON/ A(09876,99),IQ(09876),IP(09876),NB,NOEL
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
 
C------         R_ref    +dp/p     -dp/p
      DIMENSION R(6,6), RPD(6,6), RMD(6,6) 
      DIMENSION T(6,6,6)
      DIMENSION T3(5,6) , T4(5,6)
      SAVE R,T, T3,       T4

C------        Beam_ref    +dp/p     -dp/p
      DIMENSION F0(6,6), F0PD(6,6), F0MD(6,6) 

      LOGICAL KWRMAT

      LOGICAL PRDIC

      DIMENSION RO(6,6)
      character(140) BUFFER
      integer DEBSTR, finstr
      DATA KWRMAT / .FALSE. /

      IF(NRES.LE.0) RETURN

      IF(.NOT. (KOBJ.EQ.5 .OR. KOBJ.EQ.6)) THEN
        WRITE(NRES,FMT='('' Matrix  cannot  be  computed :  need "OBJET" 
     >  with  KOBJ=5 or 6'')')
        RETURN
      ENDIF

C      IORD = A(NOEL,1)
      IORD = JORD
      IF(IORD .EQ. 0) THEN
        WRITE(NRES,FMT='(/,9X,'' Matrix  not  computed : IORD = 0'',/)')
        CALL IMPTRA(1,IMAX,NRES)
        RETURN
      ENDIF
      IF(KOBJ .EQ. 5) THEN
        IORD=1
      ELSEIF(KOBJ .EQ. 6) THEN
        IORD=2
      ENDIF 

C      IFOC = A(NOEL,2) 
      IFOC = JFOC
      PRDIC = IFOC .GT. 10
 
C      KWRMAT = NINT(A(NOEL,3)) .EQ. 1
      KWRMAT = KWR .EQ. 1
      IF(KWRMAT) CALL MATIM6(KWRMAT)

      IF    (IORD .EQ. 1) THEN
        CALL OBJ51(
     >             NBREF)
        IREF = 0
 1      CONTINUE
          IREF = IREF + 1

          IT1 = 1 + 11 * (IREF-1)
          IT2 = IT1+3
          IT3 = IT1+4

C FM, Nov. 2008
C          CALL REFER(1,IORD,IFOC,IT1,IT2,IT3)
          IFC = IFOC
          IF(PRDIC) IFC = 0
          CALL REFER(1,IORD,IFC,IT1,IT2,IT3)
          CALL MAT1(R,T,IT1)
          CALL MKSA(IORD,R,T,T3,T4)
          CALL MATIMP(R)
          IF(PRDIC) CALL TUNES(R,F0,IFOC-10,IERY,IERZ,.TRUE.,
     >                                                YNU,ZNU,CMUY,CMUZ)
C FM, Nov. 2008
C          CALL REFER(2,IORD,IFOC,IT1,IT2,IT3)
          CALL REFER(2,IORD,IFC,IT1,IT2,IT3)

          IF(IREF.LT.NBREF) GOTO 1
          
      ELSEIF(IORD .EQ. 2) THEN

C FM, Nov. 2008
C        CALL REFER(1,IORD,IFOC,1,6,7)
        IFC = IFOC
        IF(PRDIC) IFC = 0
        CALL REFER(1,IORD,IFC,1,6,7)
        CALL MAT2(R,T,T3,T4)
        CALL MKSA(IORD,R,T,T3,T4)
        CALL MATIMP(R)
        IF(PRDIC) CALL TUNES(R,F0,IFOC-10,IERY,IERZ,.TRUE.,
     >                                                YNU,ZNU,CMUY,CMUZ)
        CALL MATIM2(R,T,T3)
        IF(PRDIC) THEN 
          CALL MAT2P(RPD,DP)
          CALL MKSA(IORD,RPD,T,T3,T4)
          CALL MATIMP(RPD)
          CALL TUNES(RPD,F0PD,IFOC-10,IERY,IERZ,.TRUE.,
     >                                              YNUP,ZNUP,CMUY,CMUZ)
          CALL MAT2M(RMD,DP)
          CALL MKSA(IORD,RMD,T,T3,T4)
          CALL MATIMP(RMD)
          CALL TUNES(RMD,F0MD,IFOC-10,IERY,IERZ,.TRUE.,
     >                                              YNUM,ZNUM,CMUY,CMUZ)
C Momentum detuning
          NUML = 1
C          DNUYDP = (YNUP-YNUM)/2.D0/A(NUML,25)
C          DNUZDP = (ZNUP-ZNUM)/2.D0/A(NUML,25)
          DNUYDP = (YNUP-YNUM)/2.D0/DP
          DNUZDP = (ZNUP-ZNUM)/2.D0/DP
          IF(NRES .GT. 0) WRITE(NRES,FMT='(/,34X,'' Chromaticities : '',
     >      //,30X,''dNu_y / dp/p = '',G15.8,/, 
     >         30X,''dNu_z / dp/p = '',G15.8)') DNUYDP, DNUZDP
C             write(nres,*) dp, a(numl,25)
        ENDIF
C        CALL REFER(2,IORD,IFOC,1,6,7)
        CALL REFER(2,IORD,IFC,1,6,7)
      ENDIF

c---------------------------------------------------------------------------------
c                 Exportation of the matrix coefficients (Fred)
c---------------------------------------------------------------------------------
       return
C      if(.not. prdic) return

c----------------------------------------------
c Exportation of the matrix coefficients (Fred)
c----------------------------------------------
c      CALL SYSTEM('cp transfertM.dat transfertM_save.dat')
c      OPEN(11,FILE='transfertM_save.dat',STATUS='UNKNOWN',IOSTAT=IOS1)
      OPEN(12,FILE='transfertM.dat',STATUS='UNKNOWN',IOSTAT=IOS2)
c----------------------------------------------

 3    CONTINUE
           
c     Reading of 1turnM.dat
      READ(11,FMT='(a)',END=99,ERR=98) BUFFER
c     Writing of the new transfertM.dat
      WRITE(12,FMT='(a)') BUFFER(DEBSTR(BUFFER):FINSTR(BUFFER))
      GOTO 3
 98   CONTINUE
     
      WRITE(*,FMT='(/,/,''ERROR while reading transfertm.dat'',/,/)')
 99   CONTINUE

      WRITE(12,FMT='(/,/)')
      WRITE(12,FMT='(''TRANSFERT MATRIX:'')')
      DO I=1,4
         WRITE(12,FMT='(4(F15.8,X))') (R(I,J),J=1,4)
      ENDDO

      CLOSE(11,IOSTAT=IOS1)
      CLOSE(12,IOSTAT=IOS2)
c----------------------------------------------

      call system('/home/meot/zgoubi/struct/tools/ETparam/ETparam')

      RETURN

      ENTRY MATRI1(
     >             RO)
      DO IB = 1, 6
        DO IA = 1, 6
          RO(IA,IB)  = R(IA,IB)
        ENDDO
      ENDDO
      RETURN
      END
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
C  Upton, NY, 11973
C  -------
      SUBROUTINE IMPTRA(IMAX1,IMAX2,NRES)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      COMMON/DESIN/ FDES(7,MXT),IFDES,KINFO,IRSAR,IRTET,IRPHI,NDES
     >,AMS,AMP,AM3,TDVM,TETPHI(2,MXT)
C     1,AMS ,AMP,ENSTAR,BSTAR,TDVM,TETPHI(2,MXT)
      INCLUDE "MAXCOO.H"
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      CHARACTER LET
      COMMON/FAISCT/ LET(MXT)
      LOGICAL ZSYM
      COMMON/OPTION/ KORD,KFLD,MG,LC,ML,ZSYM
      COMMON/PTICUL/ AM,Q,G,TO
      COMMON/REBELO/ NPASS,IPASS,KWRT,NNDES,STDVM
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      COMMON/UNITS/ UNIT(MXJ) 

      DIMENSION SIG(4,4)
      CHARACTER TXT*10, TXT2*2 

      WRITE(NRES,100) IMAX2-IMAX1+1
 100  FORMAT('0',45X,'TRACE DU FAISCEAU',//,45X,I4,' TRAJECTOIRES',//
     >,35X,'OBJET',50X,'FAISCEAU',//,2(10X,'D',7X,'Y(CM)',5X,'T(MR)'
     >,5X,'Z(CM)',5X,'P(MR)',4X,'S(CM)'),/)
C     1,35X,'OBJET',50X,'FAISCEAU',//,2(10X,'D',6X,'Y(CM)',4X,'T(MR)'
C     2,4X,'Z(CM)',4X,'P(MR)',3X,'S(CM)'),/)
 
      DO 1 I=IMAX1,IMAX2
        WRITE(NRES,101) LET(I),IEX(I),(FO(J,I),J=1,6)
     >  ,(F(J,I),J=1,6),I
C 101    FORMAT(' ',A1,1X,I2,1X,F8.4,5F10.3,6X,F8.4,4F9.3,1X,F10.3,1X,I6)
 101    FORMAT(A1,1X,I3,1X,F8.4,5F10.3,5X,F8.4,4F9.3,1X,F12.3,1X,I5)
        IF(AM .NE. 0D0) THEN
          IF(IFDES.EQ.1) THEN
            WRITE(NRES,FMT='(15X,''Time of flight (mus) :'',
     >      1P,G16.8,'' mass (MeV/c2) :'',G14.6, ''    decay at (m) :'',
     >           G14.6)') F(7,I),AMQ(1,I),FDES(6,I)*UNIT(5)
          ELSE
            WRITE(NRES,FMT='(15X,''Time of flight (mus) :'',
     >      1P,G16.8,'' mass (MeV/c2) :'', G14.6)') F(7,I),AMQ(1,I)
          ENDIF
        ENDIF
 1    CONTINUE

Compute rms ellipse
      WRITE(NRES,FMT='(//,''  Beam  characteristics '', 1X
     >,'' (EMIT,ALP,BET,XM,XPM,NLIV,NINL,RATIN) : '',/)')
      TXT = 'B-Dim '
      PI4 = 4.D0 *      4.D0 * ATAN(1.D0)
      DO 10 JJ = 1, 3
        CALL LPSFIT(JJ, 
     >                           EMIT,ALP,BET,XM,XPM)
Compute number of particles alive and numberinside ellipse
        CALL CNTINL(JJ,PI4*EMIT,ALP,BET,XM,XPM,
     >                                        NLIV,NINL)
        RATIN = DBLE(NINL)/DBLE(IMAX)
        WRITE(TXT2,FMT='(I1)') JJ
        WRITE(NRES,110)
     >       PI4*EMIT,ALP,BET,XM,XPM,NLIV,NINL,RATIN,TXT//TXT2,IPASS
 110    FORMAT(1P,5(1X,G12.4),2I8,1X,G12.4,A,I6)
 10   CONTINUE

Compute 4-D sigma matrix
      WRITE(NRES,FMT='(//,''  Beam  characteristics '', 1X
     >,'' SIGMA(4,4) : '',/)')
      CALL LPSFI4( 
     >             sqx,sqz,SIG)

      WRITE(NRES,fmt='(10X,1P,A,2E14.6)') ' Ex, Ez =', sqx,sqz
      WRITE(NRES,fmt='(10X,1P,A,2E14.6)') ' AlpX, BetX =', 
     >                  sig(1,2)/sqx, sig(1,1)/sqx
      WRITE(NRES,fmt='(10X,1P,A,2E14.6)') ' AlpZ, BetZ =', 
     >                  sig(3,4)/sqZ, sig(3,3)/sqx
      WRITE(NRES,120) ((SIG(I,J),J=1,4),I=1,4)
 120  FORMAT(/,1P,4E14.6)

      RETURN
      END
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
C  Upton, NY, 11973
C  -------
      SUBROUTINE TUNES(R,F0,NMAIL,IERY,IERZ,OKPR,
     >                                           YNU,ZNU,CMUY,CMUZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(6,*), F0(6,*)
      LOGICAL OKPR

      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M

      CHARACTER*15 TXTYNU, TXTZNU

      IERY=0
      IERZ=0
      CALL RAZ(F0, 6*6)

      IF(OKPR) THEN
        IF(NMAIL.LE.0) THEN
          IF(NRES.GT.0) 
     >    WRITE(NRES,*) '  NUMBER OF PERIODS = IFOC-10 = ERRONEOUS !'
          RETURN
        ELSE
          IF(NRES.GT.0) WRITE(NRES,106) NMAIL
 106      FORMAT(//,15X,' TWISS  parameters,  periodicity  of',
     >           I4,'  is  assumed :')
        ENDIF
      ENDIF
 
C  HORIZONTAL
      CMUY = .5D0 * (R(1,1)+R(2,2))
      IF(CMUY*CMUY .GE. 1.D0) THEN
        IERY=-1
        F0(1,1)=0.D0
        F0(1,2)=0.D0
        F0(2,1)=0.D0
        F0(2,2)=0.D0
        F0(1,6)=0.D0
        F0(2,6)=0.D0
        YNU = 0.D0
      ELSE

        COSMU=CMUY
CBETA style
        SINMU=SIGN( SQRT(1.0D0-COSMU*COSMU) , R(1,2) )
        YNU=ACOS(COSMU)*DBLE(NMAIL)/(2.D0 * PI) 
        IF (R(1,2) .LT. 0.0D0) YNU=DBLE(NMAIL)-YNU
        YNU=YNU-AINT(YNU)
C         write(*,*) ' beta style, sinmu, ynu = ',sinmu,ynu,nmail
CMAD style. Optimized precisison close to 1/2 ? 
        SINMU=SIGN(SQRT(-R(1,2)*R(2,1)-.25D0*(R(1,1)-R(2,2))**2),R(1,2))
        YNU = SIGN(ATAN2(SINMU,COSMU) /(2.D0 * PI) ,R(1,2))
        IF (R(1,2) .LT. 0.0D0) YNU=DBLE(NMAIL)+YNU
C         write(*,*) ' mad style, sinmu, ynu = ',sinmu,ynu

        BX0=R(1,2) / SINMU
        AX0=(R(1,1) - R(2,2)) / (2.D0*SINMU)
        F0(1,1)=BX0
        F0(1,2)=-AX0
        F0(2,1)=-AX0
        F0(2,2)=(1.D0+AX0*AX0)/BX0
        S2NU = 2.D0 - R(1,1) - R(2,2)
        F0(1,6) = ( (1.D0 - R(2,2))*R(1,6) + R(1,2)*R(2,6) ) / S2NU
        F0(2,6) = ( R(2,1)*R(1,6) + (1.D0 - R(1,1))*R(2,6) ) / S2NU

      ENDIF

C  VERTICAL
      CMUZ = .5D0 * (R(3,3)+R(4,4))
      IF(CMUZ*CMUZ .GE. 1.D0) THEN
        IERZ=-1
        F0(3,3)=0.D0
        F0(3,4)=0.D0
        F0(4,3)=0.D0
        F0(4,4)=0.D0
        F0(3,6)=0.D0
        F0(4,6)=0.D0
        ZNU = 0.D0
      ELSE
        COSMU = CMUZ
CBETA style
            ZNU=ACOS(COSMU)*DBLE(NMAIL)/(2.D0 * PI) 
            IF (R(3,4) .LT. 0.0D0) ZNU=DBLE(NMAIL)-ZNU
            ZNU=ZNU-AINT(ZNU)
            SINMU=SIGN( SQRT(1.0D0-COSMU*COSMU) , R(3,4) )
C         write(*,*) ' beta style, sinmu, znu = ',sinmu,znu,nmail
CMAD style. Optimized precisison close to 1/2 ? 
        SINMU=SIGN(SQRT(-R(3,4)*R(4,3)-.25D0*(R(3,3)-R(4,4))**2),R(3,4))
        ZNU = SIGN(ATAN2(SINMU,COSMU) /(2.D0 * PI) ,R(3,4))
        IF (R(3,4) .LT. 0.0D0) ZNU=DBLE(NMAIL)+ZNU
C         write(*,*) ' mad style, sinmu, znu = ',sinmu,znu
C         write(*,*) '                                tunes...' 

        BZ0=R(3,4)/SINMU
        AZ0=(R(3,3)-R(4,4))/2.D0/SINMU
        F0(3,3)=BZ0
        F0(3,4)=-AZ0
        F0(4,3)=-AZ0
        F0(4,4)=(1.D0+AZ0*AZ0)/BZ0
        S2NU = 2.D0 - R(3,3) - R(4,4)
        F0(3,6) = ( (1.D0 - R(4,4))*R(3,6) + R(3,4)*R(4,6) ) / S2NU
        F0(4,6) = ( R(4,3)*R(3,6) + (1.D0 - R(3,3))*R(4,6) ) / S2NU

      ENDIF

      IF(OKPR) THEN

        IF(NRES.GT.0) THEN
          WRITE(NRES,103)
 103      FORMAT(/,6X,
     >    ' Beam  matrix  (beta/-alpha/-alpha/gamma)',
     >    ' and  periodic  dispersion  (MKSA units)',/)
          WRITE(NRES,104) (( F0(IA,IB) , IB=1,6) , IA=1,6)
 104      FORMAT(6X,6F13.6)
          WRITE(NRES,FMT='(/,35X,''Betatron  tunes'',/)') 
          WRITE(TXTYNU,FMT='(A)') 'undefined'
          WRITE(TXTZNU,FMT='(A)') 'undefined'

         IF    (ABS(CMUY).LT.1.D0 .AND. ABS(CMUZ).LT.1.D0) THEN
           WRITE(TXTYNU,FMT='(G15.8)') YNU
           WRITE(TXTZNU,FMT='(G15.8)') ZNU
         ELSEIF(ABS(CMUY).LT.1.D0 .OR. ABS(CMUZ).LT.1.D0) THEN
           IF(CMUY*CMUY .LT. 1.D0) WRITE(TXTYNU,FMT='(G15.8)') YNU
           IF(CMUZ*CMUZ .LT. 1.D0) WRITE(TXTZNU,FMT='(G15.8)') ZNU
         ENDIF

         WRITE(NRES,FMT='(15X,2(5X,A,A))') 
     >              'NU_Y = ', TXTYNU, 'NU_Z = ', TXTZNU

        ENDIF
      ENDIF
      RETURN
      END
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
C  Upton, NY, 11973
C  -------
      SUBROUTINE MAT2(R,T,T3,T4)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(6,*) , T(6,6,*)
      DIMENSION  T3(5,*), T4(5,*)
C     ------------------------------------
C     Option  IORD = 2 :
C       Matrix ordre 1, coeff ordre 2, & 
C        coeffs > 2; Taylor series ordre 5
C     ------------------------------------
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
C      COMMON/DON/ A(09876,99),IQ(09876),IP(09876),NB,NOEL
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
 
      DIMENSION RPD(6,6), RMD(6,6), RPDO(6,6), RMDO(6,6) 
      SAVE RPD, RMD, DP

      CALL RAZ(R,6*6)
      R(5,5) = 1D0
      R(6,6) = 1D0
      CALL RAZ(T,5*6*6)
      S1 = 2.D0 * F(6,1)
 
C      ++++ Ordre  1  ,  &  ordres 2  to  4  uncoupled
C           +++ utilise  traj. 1-21
 
      DP = ( FO(1,18)-FO(1,19) )/( FO(1,18)+FO(1,19) )
      DP2 = DP*DP
      DO 10 J=2,5
        J1=J-1
        QP = F(J,18)
        QM = F(J,19)
        QPP= F(J,20)
        QMM= F(J,21)
        R(J1,6) =  ( 8.D0*(QP-QM) -(QPP-QMM) )/12.D0/DP
        T(J1,6,6)= (16.D0*(QP+QM) -(QPP+QMM) )/24.D0/DP2
        T3(J1,6)=-( 2.D0*(QP-QM) -(QPP-QMM) )/12.D0/DP2/DP
        T4(J1,6)=-( 4.D0*(QP+QM) -(QPP+QMM) )/24.D0/DP2/DP2
        DO 10 I=1,4
          II = 4*I - 2
          UO = FO(I+1,II)
          UO2 = UO*UO
          QP = F(J,II)
          QM = F(J,II+1)
          QPP= F(J,II+2)
          QMM= F(J,II+3)
          R(J1,I) =  ( 8.D0*(QP-QM) -(QPP-QMM) ) /12.D0 / UO
          T(J1,I,I)= (16.D0*(QP+QM) -(QPP+QMM) ) /24.D0 / UO2
          T3(J1,I)=-( 2.D0*(QP-QM) -(QPP-QMM) ) /12.D0 / UO2/UO
          T4(J1,I)=-( 4.D0*(QP+QM) -(QPP+QMM) ) /24.D0 /  UO2/UO2
          IF(J .EQ. 5) THEN
            R(5,I) = ( 8.D0*( F(6,II) - F(6,II+1) )
     >        - ( F(6,II+2) - F(6,II+3) )          ) /12.D0/UO
            T(5,I,I)=( 16.D0*( F(6,II) +F(6,II+1) - S1 )
     >        - ( F(6,II+2) + F(6,II+3) - S1 )     ) /24.D0/UO2
          ENDIF
 10   CONTINUE
      QP = F(6,18)
      QM = F(6,19)
      QPP= F(6,20)
      QMM= F(6,21)
      R(5,6)  = ( 8.D0*( QP-QM   ) - ( QPP-QMM   )  )/12.D0/DP
      T(5,6,6)= (16.D0*( QP+QM-S1) - ( QPP+QMM-S1)  )/24.D0/DP2
      T3(5,6)=-( 2.D0*( QP-QM   ) - ( QPP-QMM   )  )/12.D0/DP2/DP
      T4(5,6)=-( 4.D0*( QP+QM-S1) - ( QPP+QMM-S1)  )/24.D0/DP2/DP2
 
      CALL RAZ(RPD,4*4)
      DO 20 J=2,5
        J1=J-1
          UO = FO(2,34) - FO(2,37)
          RPD(J1,1) = (F(J,34) - F(J,37)) / UO
          UO = FO(3,46) - FO(3,49)
          RPD(J1,2) = (F(J,46) - F(J,49)) / UO
          UO = FO(4,54) - FO(4,57)
          RPD(J1,3) = (F(J,54) - F(J,57)) / UO
          UO = FO(5,58) - FO(5,61)
          RPD(J1,4) = (F(J,58) - F(J,61)) / UO
 20   CONTINUE
 
      CALL RAZ(RMD,4*4)
      DO 21 J=2,5
        J1=J-1
          UO = FO(2,35) - FO(2,36)
          RMD(J1,1) = (F(J,35) - F(J,36)) / UO
          UO = FO(3,47) - FO(3,48)
          RMD(J1,2) = (F(J,47) - F(J,48)) / UO
          UO = FO(4,55) - FO(4,56)
          RMD(J1,3) = (F(J,55) - F(J,56)) / UO
          UO = FO(5,59) - FO(5,60)
          RMD(J1,4) = (F(J,59) - F(J,60)) / UO
 21    CONTINUE
 
C      ++++ Ordre 2  coupled
C           +++ utilise  traj. >=  #22
 
      III=4
C     ++++ COUPLAGE ./YT
      II=22
      J1=1
      J2=2
      UO = FO(J1+1,II)
      VO = FO(J2+1,II)
C     ... J/J1.J2 : Y/YT  ... L/YT
      DO 12 J=1,5
        U=UO
        V=VO
        QPP = F(J+1,II)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPP=QPP-F(6,1)
        U=UO
        V=-VO
        QPM = F(J+1,II+1)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPM=QPM-F(6,1)
        U=-UO
        V=-VO
        QMM = F(J+1,II+2)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMM=QMM-F(6,1)
        U=-UO
        V=VO
        QMP = F(J+1,II+3)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMP=QMP-F(6,1)
        T(J,J1,J2)=(QPP+QMM-QPM-QMP)/(4.D0*UO*VO)   *0.5D0
 12   CONTINUE
 
C     ++++ COUPLAGE ./YZ
      II=II+III
      J1=1
      J2=3
      UO = FO(J1+1,II)
      VO = FO(J2+1,II)
C     ... J/J1.J2 : Y/YZ  ... L/YZ
      DO 13 J=1,5
        U=UO
        V=VO
        QPP = F(J+1,II)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPP=QPP-F(6,1)
        U=UO
        V=-VO
        QPM = F(J+1,II+1)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPM=QPM-F(6,1)
        U=-UO
        V=-VO
        QMM = F(J+1,II+2)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMM=QMM-F(6,1)
        U=-UO
        V=VO
        QMP = F(J+1,II+3)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMP=QMP-F(6,1)
        T(J,J1,J2)=(QPP+QMM-QPM-QMP)/(4.D0*UO*VO)    *0.5D0
 13   CONTINUE
 
C     ++++ COUPLAGE ./YP
      II=II+III
      J1=1
      J2=4
      UO = FO(J1+1,II)
      VO = FO(J2+1,II)
C     ... J/J1.J2 : Y/YP  ... L/YP
      DO 14 J=1,5
        U=UO
        V=VO
        QPP = F(J+1,II)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPP=QPP-F(6,1)
        U=UO
        V=-VO
        QPM = F(J+1,II+1)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPM=QPM-F(6,1)
        U=-UO
        V=-VO
        QMM = F(J+1,II+2)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMM=QMM-F(6,1)
        U=-UO
        V=VO
        QMP = F(J+1,II+3)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMP=QMP-F(6,1)
        T(J,J1,J2)=(QPP+QMM-QPM-QMP)/(4.D0*UO*VO)    *0.5D0
 14   CONTINUE
 
C     ++++ COUPLAGE ./YD
      II=II+III
      J1=1
      J2=6
      UO = FO(J1+1,II)
      VO = FO(   1,II)-FO(1,1)
C     ... J/J1.J2 : Y/YD  ... L/YD
      DO 16 J=1,5
        U=UO
        V=VO
        QPP = F(J+1,II)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPP=QPP-F(6,1)
        U=UO
        V=-VO
        QPM = F(J+1,II+1)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPM=QPM-F(6,1)
        U=-UO
        V=-VO
        QMM = F(J+1,II+2)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMM=QMM-F(6,1)
        U=-UO
        V=VO
        QMP = F(J+1,II+3)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMP=QMP-F(6,1)
        T(J,J1,J2)=(QPP+QMM-QPM-QMP)/(4.D0*UO*VO)    *0.5D0
 16   CONTINUE
 
C     ++++ COUPLAGE ./TZ
      II=II+III
      J1=2
      J2=3
      UO = FO(J1+1,II)
      VO = FO(J2+1,II)
C     ... J/J1.J2 : Y/TZ  ... L/TZ
      DO 23 J=1,5
        U=UO
        V=VO
        QPP = F(J+1,II)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPP=QPP-F(6,1)
        U=UO
        V=-VO
        QPM = F(J+1,II+1)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPM=QPM-F(6,1)
        U=-UO
        V=-VO
        QMM = F(J+1,II+2)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMM=QMM-F(6,1)
        U=-UO
        V=VO
        QMP = F(J+1,II+3)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMP=QMP-F(6,1)
        T(J,J1,J2)=(QPP+QMM-QPM-QMP)/(4.D0*UO*VO)    *0.5D0
 23   CONTINUE
 
C     ++++ COUPLAGE ./TP
      II=II+III
      J1=2
      J2=4
      UO = FO(J1+1,II)
      VO = FO(J2+1,II)
C     ... J/J1.J2 : Y/TP  ... L/TP
      DO 24 J=1,5
        U=UO
        V=VO
        QPP = F(J+1,II)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPP=QPP-F(6,1)
        U=UO
        V=-VO
        QPM = F(J+1,II+1)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPM=QPM-F(6,1)
        U=-UO
        V=-VO
        QMM = F(J+1,II+2)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMM=QMM-F(6,1)
        U=-UO
        V=VO
        QMP = F(J+1,II+3)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMP=QMP-F(6,1)
        T(J,J1,J2)=(QPP+QMM-QPM-QMP)/(4.D0*UO*VO)    *0.5D0
 24   CONTINUE
 
C     ++++ COUPLAGE ./TD
      II=II+III
      J1=2
      J2=6
      UO = FO(J1+1,II)
      VO = FO(   1,II)-FO(1,1)
C     ... J/J1.J2 : Y/TD  ... L/TD
      DO 26 J=1,5
        U=UO
        V=VO
        QPP = F(J+1,II)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPP=QPP-F(6,1)
        U=UO
        V=-VO
        QPM = F(J+1,II+1)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPM=QPM-F(6,1)
        U=-UO
        V=-VO
        QMM = F(J+1,II+2)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMM=QMM-F(6,1)
        U=-UO
        V=VO
        QMP = F(J+1,II+3)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMP=QMP-F(6,1)
        T(J,J1,J2)=(QPP+QMM-QPM-QMP)/(4.D0*UO*VO)    *0.5D0
 26   CONTINUE
 
C     ++++ COUPLAGE ./ZP
      II=II+III
      J1=3
      J2=4
      UO = FO(J1+1,II)
      VO = FO(J2+1,II)
C     ... J/J1.J2 : Y/ZP ...  L/ZP
      DO 34 J=1,5
        U=UO
        V=VO
        QPP = F(J+1,II)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPP=QPP-F(6,1)
        U=UO
        V=-VO
        QPM = F(J+1,II+1)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPM=QPM-F(6,1)
        U=-UO
        V=-VO
        QMM = F(J+1,II+2)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMM=QMM-F(6,1)
        U=-UO
        V=VO
        QMP = F(J+1,II+3)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMP=QMP-F(6,1)
        T(J,J1,J2)=(QPP+QMM-QPM-QMP)/(4.D0*UO*VO)    *0.5D0
 34   CONTINUE
 
C     ++++ COUPLAGE ./ZD
      II=II+III
      J1=3
      J2=6
      UO = FO(J1+1,II)
      VO = FO(   1,II)-F(1,1)
C     ... J/J1.J2 : Y/ZD  ... L/ZD
      DO 36 J=1,5
        U=UO
        V=VO
        QPP = F(J+1,II)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPP=QPP-F(6,1)
        U=UO
        V=-VO
        QPM = F(J+1,II+1)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPM=QPM-F(6,1)
        U=-UO
        V=-VO
        QMM = F(J+1,II+2)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMM=QMM-F(6,1)
        U=-UO
        V=VO
        QMP = F(J+1,II+3)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMP=QMP-F(6,1)
        T(J,J1,J2)=(QPP+QMM-QPM-QMP)/(4.D0*UO*VO)    *0.5D0
 36   CONTINUE
 
C     ++++ COUPLAGE ./PD
      II=II+III
      J1=4
      J2=6
      UO = FO(J1+1,II)
      VO = FO(   1,II)-F(1,1)
C     ... J/J1.J2 : Y/PD  ... L/PD
      DO 46 J=1,5
        U=UO
        V=VO
        QPP = F(J+1,II)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPP=QPP-F(6,1)
        U=UO
        V=-VO
        QPM = F(J+1,II+1)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QPM=QPM-F(6,1)
        U=-UO
        V=-VO
        QMM = F(J+1,II+2)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMM=QMM-F(6,1)
        U=-UO
        V=VO
        QMP = F(J+1,II+3)-
     >  ( ( R(J,J1)+(T(J,J1,J1)+(T3(J,J1)+T4(J,J1)*U)*U)*U)*U
     >   +( R(J,J2)+(T(J,J2,J2)+(T3(J,J2)+T4(J,J2)*V)*V)*V)*V  )
        IF(J .EQ. 5) QMP=QMP-F(6,1)
        T(J,J1,J2)=(QPP+QMM-QPM-QMP)/(4.D0*UO*VO)    *0.5D0
 46   CONTINUE
 
      RETURN

      ENTRY MAT2P(RPDO,DPO)
      DO 66 J=1,6
        DO 66 I=1,6
 66       RPDO(I,J) = RPD(I,J) 
      DPO = DP
      RETURN

      ENTRY MAT2M(RMDO,DPO)
      DO 67 J=1,6
        DO 67 I=1,6
 67        RMDO(I,J) = RMD(I,J) 
      DPO = DP
      RETURN

      END
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
C  Brookhaven National Laboratory                    és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      FUNCTION DEBSTR(STRING)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DEBSTR
      CHARACTER * (*) STRING

C     --------------------------------------
C     RENVOIE DANS DEBSTR LE RANG DU
C     1-ER CHARACTER NON BLANC DE STRING,
C     OU BIEN 0 SI STRING EST VIDE ou BLANC.
C     --------------------------------------

      DEBSTR=0
      LENGTH=LEN(STRING)

      IF(LENGTH.EQ.0) RETURN

1     CONTINUE
        DEBSTR=DEBSTR+1
        IF (STRING(DEBSTR:DEBSTR) .EQ. ' ') THEN
          IF(DEBSTR .GE. LENGTH) THEN
            DEBSTR = 0
            RETURN
          ELSE
            GOTO 1
          ENDIF
        ENDIF

      RETURN
      END
      SUBROUTINE BEAIMP(F0,PHY,PHZ) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/BEAM/ FI(6,6)
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      DIMENSION F0(6,*)

      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      IF(NRES .LT. 0) RETURN
      WRITE(NRES,103) 'INITIAL'
 103  FORMAT(//,18X,'BEAM  MATRIX (beta/alpha/alpha/gamma, D,D''), 
     >       ',A,/)
      WRITE(NRES,104) (( FI(IA,IB) , IB=1,6) , IA=1,6)
 104  FORMAT(6X,1P,6G16.6)
      WRITE(NRES,103) 'FINAL' 
      WRITE(NRES,104) (( F0(IA,IB) , IB=1,6) , IA=1,6)
      WRITE(NRES,FMT='(/,18X,''Betatron phase advances (fractional),  ''
     >,''phi_y/2pi,'',''  phi_z/2pi :''
     >,//,18X,1P,2(5X,E14.6))') PHY/(2.D0*PI), PHZ/(2.D0*PI)
      RETURN
      END
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
C  Brookhaven National Laboratory                    és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE LPSFIT(JJ, 
     >                               SQ,A,B,XM,XPM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CONST/ CL9,CEL,PI,RAD,DEG,QE ,AMPROT, CM2M
      COMMON/CONST2/ ZERO, UN
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
C----- CONVERSION DES COORD. (CM,MRD) -> (M,RD)
      COMMON/UNITS/ UNIT(MXJ)

        J1 = 2*JJ
        J2 = J1+1
        UNIT1=  UNIT(J1-1)
        UNIT2=  UNIT(J2-1)
        IF(J1.EQ.6) THEN
C--------- Time-momentum
          J1=7
          J2=1
          UNIT1=  1.D0
          UNIT2=  1.D0
        ENDIF
        XM=0.D0
        XPM=0.D0
        NPTS = 0
        DO 21 I=1,IMAX
          IF(IEX(I) .LT. -1) GOTO 21
          NPTS = NPTS + 1
          X   = F(J1,I)*UNIT1
          XP  = F(J2,I)*UNIT2
          IF(J2.EQ.1) THEN
C--------- To get kineticE
              P = BORO*CL9 *XP * AMQ(2,I)
              XP= SQRT(P*P + AMQ(1,I)*AMQ(1,I))- AMQ(1,I)
          ENDIF
          XM  = XM + X
          XPM = XPM + XP
 21     CONTINUE
        XM = XM/NPTS
        XPM = XPM/NPTS
        YM=XM
        YPM=XPM

        J1 = 2*JJ
        J2 = J1+1
        UNIT1=  UNIT(J1-1)
        UNIT2=  UNIT(J2-1)
        IF(J1.EQ.6) THEN
C--------- Time-momentum
          J1=7
          J2=1
          UNIT1=  1.D0
          UNIT2=  1.D0
        ENDIF
        X2=0.D0
        XP2=0.D0
        XXP=0.D0
        DO 26 I=1,IMAX
          IF(IEX(I) .LT. -1) GOTO 26
          X   = F(J1,I)*UNIT1
          XP  = F(J2,I)*UNIT2
          IF(J2.EQ.1) THEN
C--------- To get kineticE
              P = BORO*CL9 *XP * AMQ(2,I)
              XP= SQRT(P*P + AMQ(1,I)*AMQ(1,I))- AMQ(1,I)
          ENDIF
          X2  = X2 + (X-YM)**2
          XP2 = XP2 + (XP-YPM)**2
          XXP = XXP + (X-YM)*(XP-YPM)
 26     CONTINUE
        X2  = X2/NPTS
        XP2 = XP2/NPTS
        XXP = XXP/NPTS

C G. Leleux : surface de l'ellipse S=4.pi.sqrt(DELTA)
C Soit d11=X2/sqrt(DELTA), d12=XXP/sqrt(DELTA), d22=XP2/sqrt(DELTA), alors 
C d22.x^2-2.d12.x.x'+d11.x'^2=S/pi=4sqrt(DELTA), ce qui permet d'ecrire 
C gamma=d22=XP2/sqrt(DELTA),-alpha=d12=XXP/sqrt(DELTA),beta=d11=X2/sqrt(DELTA).
C En outre, par definition des dij, 
C     2.sigma_x=sqrt(d11.S/pi),  2.sigma_x'=sqrt(d22.S/pi). 
C En outre, frontiere : 
C          <x^2>_frontiere=2.(sigma_x)^2,    <x'^2>_frontiere=2.(sigma_x')^2

C------- Courant invariant at 1 sigma is U=4.sqrt(DELTA)=Eps/pi (consistant with zgoubi !!) :
C Eps=ellipse surface
        SQ = SQRT(X2*XP2-XXP*XXP) 
        IF(SQ .GT. 0.D0) THEN
          B=  X2/SQ
          A=  -XXP/SQ
        ENDIF

      RETURN
      END
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
C  Brookhaven National Laboratory                    és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE CNTINL(JJ,EPSPI,ALP,BET,XM,XPM,
     >                                                 LIVE,NINL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
C----- CONVERT COORD. (CM,MRD) -> (M,RD)
      COMMON/UNITS/ UNIT(MXJ)

        J1 = 2*JJ
        J2 = J1+1
        UNIT1=  UNIT(J1-1)
        UNIT2=  UNIT(J2-1)
        IF(J1.EQ.6) THEN
C--------- Time-momentum
          J1=7
          J2=1
          UNIT1=  1.D0
          UNIT2=  1.D0
        ENDIF
        GAM = (1.D0+ALP*ALP) / BET

      LIVE = 0
      NINL = 0
      DO 1 I=1,IMAX
          IF(IEX(I) .LT. -1) GOTO 1
          LIVE = LIVE + 1
  
          IF    (JJ.EQ.1) THEN
C----------- Horizontal  
            Y2 = F(2,I)*UNIT(1) - XM
            T2 = F(3,I)*UNIT(2) - XPM
            YT = Y2*T2
            Y2 = Y2*Y2
            T2 = T2*T2 
          ELSEIF(JJ.EQ.2) THEN
C----------- Vertical  
            Y2 = F(4,I)*UNIT(3) - XM
            T2 = F(5,I)*UNIT(4) - XPM
            YT = Y2*T2
            Y2 = Y2*Y2
            T2 = T2*T2
          ELSEIF(JJ.EQ.3) THEN
C----------- Time-kineticE
CCCCC----------- Time-momentum
            Y2 = F(7,I) - XM
CCCCC            T2 = F(1,I)*UNIT(6) - XPM
            P = BORO*CL9 *F(1,I) * AMQ(2,I)
            T2 = SQRT(P*P + AMQ(1,I)*AMQ(1,I))- AMQ(1,I) - XPM
            YT = Y2*T2
            Y2 = Y2*Y2
            T2 = T2*T2
          ENDIF

          IF( GAM*Y2+2.D0*ALP*YT+BET*T2 .LE. EPSPI ) NINL = NINL + 1

 1      CONTINUE
      RETURN
      END 
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
C  Upton, NY, 11973
C  -------
      FUNCTION STRCON(STR,STR2,
     >                         IS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL STRCON
      CHARACTER STR*(*), STR2*(*)
C     ---------------------------------------------------------------
C     .TRUE. if the string STR contains the string STR2 at least once
C     IS = position of first occurence of STR2 in STR 
C     (i.e.,STR(IS:IS+LEN(STR2)-1)=STR2)
C     ---------------------------------------------------------------
      INTEGER DEBSTR,FINSTR
      LNG2 = LEN(STR2(DEBSTR(STR2):FINSTR(STR2)))
      IF(LEN(STR).LT.LNG2) GOTO 1
      DO I = DEBSTR(STR), FINSTR(STR)-LNG2+1
        IF( STR(I:I+LNG2-1) .EQ. STR2 ) THEN
          IS = I 
          STRCON = .TRUE.
          RETURN
        ENDIF
      ENDDO
 1    CONTINUE
      STRCON = .FALSE.
      RETURN
      END
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
C  Brookhaven National Laboratory                    és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE LPSFI4(
     >                  SQX,SQZ,S)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION S(4,4)

      COMMON/CONST/ CL9,CEL,PI,RAD,DEG,QE ,AMPROT, CM2M
      COMMON/CONST2/ ZERO, UN
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
C----- CONVERSION DES COORD. (CM,MRD) -> (M,RD)
      COMMON/UNITS/ UNIT(MXJ)

        UNIT1=  UNIT(1)
        UNIT2=  UNIT(2)
        UNIT3=  UNIT(3)
        UNIT4=  UNIT(4)

        XM=0.D0
        XPM=0.D0
        ZM=0.D0
        ZPM=0.D0
        NPTS = 0
        DO 21 I=1,IMAX
          IF(IEX(I) .LT. -1) GOTO 21
          NPTS = NPTS + 1
          X   = F(2,I)*UNIT1
          XP  = F(3,I)*UNIT2
          Z   = F(4,I)*UNIT3
          ZP  = F(5,I)*UNIT4
          XM  = XM  + X
          XPM = XPM + XP
          ZM  = ZM  + Z
          ZPM = ZPM + ZP
 21     CONTINUE
        XM =  XM/NPTS
        XPM = XPM/NPTS
        ZM =  ZM/NPTS
        ZPM = ZPM/NPTS
        UM= XM
        UPM=XPM
        VM= ZM
        VPM=ZPM

        X2  = 0.D0
        XP2 = 0.D0
        XXP = 0.D0
        Z2  = 0.D0
        ZP2 = 0.D0
        ZZP = 0.D0
        XZ  = 0.D0
        XZP = 0.D0
        XPZ = 0.D0
        XPZP= 0.D0
        DO 26 I=1,IMAX
          IF(IEX(I) .LT. -1) GOTO 26
          X   = F(2,I)*UNIT1
          XP  = F(3,I)*UNIT2
          Z   = F(4,I)*UNIT3
          ZP  = F(5,I)*UNIT4
          X2  = X2  + (X- UM)**2
          XP2 = XP2 + (XP-UPM)**2
          XXP = XXP + (X- UM)*(XP-UPM)
          Z2  = Z2  + (Z- VM)**2
          ZP2 = ZP2 + (ZP-VPM)**2
          ZZP = ZZP + (Z- VM)*(ZP-VPM)
          XZ  = XZ  + (X- UM)*(Z-VM)
          XZP = XZP + (X- UM)*(ZP-VPM)
          XPZ = XPZ + (XP- UPM)*(Z-VM)
          XPZP = XPZP + (XP- UPM)*(ZP-VPM)
 26     CONTINUE
        S(1,1) = X2/NPTS
        S(1,2) = XXP/NPTS
        S(1,3) = XZ /NPTS
        S(1,4) = XZP/NPTS
        S(2,1) = S(1,2)
        S(2,2) = XP2/NPTS
        S(2,3) = XPZ/NPTS
        S(2,4) = XPZP/NPTS
        S(3,1) = S(1,3)
        S(3,2) = S(2,3)
        S(3,3) = Z2/NPTS
        S(3,4) = ZZP/NPTS
        S(4,1) = S(1,4)
        S(4,2) = S(2,4)
        S(4,3) = S(3,4)
        S(4,4) = ZP2/NPTS

C G. Leleux : surface de l'ellipse S=4.pi.sqrt(DELTA)
C Soit d11=X2/sqrt(DELTA), d12=XXP/sqrt(DELTA), d22=XP2/sqrt(DELTA), alors 
C d22.x^2-2.d12.x.x'+d11.x'^2=S/pi=4sqrt(DELTA), ce qui permet d'ecrire 
C gamma=d22=XP2/sqrt(DELTA),-alpha=d12=XXP/sqrt(DELTA),beta=d11=X2/sqrt(DELTA).
C En outre, par definition des dij, 
C     2.sigma_x=sqrt(d11.S/pi),  2.sigma_x'=sqrt(d22.S/pi). 
C En outre, frontiere : 
C          <x^2>_frontiere=2.(sigma_x)^2,    <x'^2>_frontiere=2.(sigma_x')^2

C------- Courant invariant at 1 sigma is U=4.sqrt(DELTA)=Eps/pi (consistant with zgoubi !!) :
C Eps=ellipse surface
        SQX = SQRT(S(1,1)*S(2,2)-S(1,2)*S(1,2)) 
        SQZ = SQRT(S(3,3)*S(4,4)-S(3,4)*S(3,4)) 

      RETURN
      END
