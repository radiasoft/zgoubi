c      SUBROUTINE SPEANA(YM,BORNE,NC0,
c     >                               YNU,SPEC,PMAX)
      implicit double precision (a-h,o-z)
      DIMENSION YM(3), BORNE(2)
      parameter(mxpts=600000)
      dimension spin(3,mxpts)
      parameter(mxc=20000)
      dimension SPEC(mxc), q(mxc)
      parameter(mxspec=mxpts/100)
      dimension Qs(mxspec)
      PARAMETER ( PI=3.1415926536 , DEUXPI=2.0*PI )
      logical idluni

      PARAMETER (mxj=7,MXS=4)
      DIMENSION FO(MXJ),F(MXJ),SO(MXS),SF(MXS)

      PARAMETER (LBLSIZ=8)
      CHARACTER*(LBLSIZ) LBL1, LBL2

      PARAMETER (KSIZ=10)
      CHARACTER KLEY*(KSIZ) 
      CHARACTER TX1*1
      CHARACTER*1 LET

      DATA KSX, KSY, KSZ / 1, 2, 3 /

      logical okend 
      character*700 txt700
      logical exs

      dimension omga(3), b(3), c(3), u1(3), u2(3)
      dimension smi(3), sma(3)

      data nc0 / 2500 /
      data okend / .false. /
      data borne / 0.48d0, .52d0 / 

      data KPA, Ksmpl / 1000, 300 /
      data itraj /  1  / 

      write(*,*) ' '
      write(*,*) '----------------------------------------------------'
      write(*,*) 'NOW RUNNING PGM spinTuneFromFai... '
      write(*,*) '----------------------------------------------------'
      write(*,*) ' '

C Range of turns to be considered may be specified using spinTuneFromFai.In
      INQUIRE(FILE='spinTuneFromFai.In',exist=EXS)

      IF (IDLUNI(NLU)) THEN
        open(unit=nlu,file='spinTuneFromFai.In')
      ELSE
        stop 'Pgm spinTuneFromFai :   No idle unit number ! '
      ENDIF

      if(.NOT.exs) then
        write(*,*)'WARNING : File spinTuneFromFai.In does not exist'
        write(*,*)'Pgm spinTuneFromFai creates one from default values'

        write(nlu,fmt='(2(i6,1x),t60,t60,a)')  kpa, ksmpl
     >  ,' ! kpa, ksmpl : Fourier transf. ksmpl turns starting from kpa'
        write(nlu,fmt='(i6,1x,t60,a)')  itraj
     >  ,' ! itraj : number of the particle to be Fourier''ed'
        write(nlu,fmt='(2(f8.4,1x),t60,a)')  (borne(i),i=1,2)
     >  ,' ! Q_1, Q_2 : x/y/l  spectrum range'
        write(nlu,fmt='(i6,1x,t60,a)')  nc0
     >  ,' ! nbin : x/y/l # of bins in spectrum range'
      endif

      rewind(nlu)

      read(nlu,*,err=11,end=11) kpa, ksmpl
      read(nlu,*,err=11,end=11) itraj
      read(nlu,*,err=11,end=11) borne(1), borne(2) 
      read(nlu,*,err=11,end=11) nc0

      close(nlu)

      KPB = KPA + ksmpl -1

      if(idluni(lunR)) then
        open(unit=lunR,file='b_zgoubi.fai',FORM='UNFORMATTED')
            CALL HEADER(lunR,4,.TRUE.,*98)
c        open(unit=lunR,file='zgoubi.fai')
c            CALL HEADER(lunR,4,.FALSE.,*98)
      else
        stop ' PGM SPINTUNEFROMFAI : no idle unit'
      endif

      if(idluni(lunW)) 
     >   open(unit=lunW,file='spinTuneFromFai.Out')
      if(idluni(lunW2)) 
     >   open(unit=lunW2,file='spinTuneFromFai_spectrum.Out')

      ii = 1
 7    continue

      do jc=1, 3
        smi(jc) = 1d10
        sma(jc) = -1d10
      enddo

      i = 1
 1    continue
            READ(lunR,ERR=3,END=3) 
     >      KEX,(FO(J),J=1,7),
     >      (F(J),J=1,7), 
     >      (SO(J),J=1,4),spin(1,i),spin(2,i),spin(3,i),dum,
     >      ENEKI, ENERG, 
     >      IT, IREP, SORT, AMQ1,AMQ2,AMQ3,AMQ4,AMQ5, RET, DPR, PS,
     >      BORO, IPASS, NOEL ,KLEY,LBL1,LBL2,LET
c            READ(LUNr,110,ERR=99,END=3)
c     >      KEX,(FO(J),J=1,7),
c     >      (F(J),J=1,7), 
c     >      (SO(J),J=1,4),(SF(J),J=1,4),
c     >      ENEKI, ENERG, 
c     >      IT, IREP, SORT, AMQ1,AMQ2,AMQ3,AMQ4,AMQ5, RET, DPR,  PS,
c     >      BORO, IPASS,NOEL, 
c     >      TX1,KLEY,TX1,TX1,LBL1,TX1,TX1,LBL2,TX1,TX1,LET,TX1
C        write(*,*) (SF(J),J=1,4),it,i 
       IF(IT.NE.itraj) goto 1
       IF(ipass .lt. kpa) goto 1
c        write(*,*) ' spectrum,  i, it, itraj :',i,it,itraj
c        write(*,*) ' spin(3,i) :',spin(3,i),i,ksmpl
        if(i.ge.Ksmpl) goto 2
        if(i.ge.mxpts) goto 2
    
        do jc = 1, 3
          if(smi(jc) .gt. spin(jc,i)) smi(jc)=spin(jc,i)
          if(sma(jc) .lt. spin(jc,i)) sma(jc)=spin(jc,i)
        enddo
c        write(*,*)  smi,  sma

        i = i + 1
      goto 1

 3    i = i-1
      okend = .true.
 2    continue
      npt = i

C omga is ~ the local precession direction
      omga(1) = (sma(1) + smi(1)) / 2.d0
      omga(2) = (sma(2) + smi(2)) / 2.d0
      omga(3) = (sma(3) + smi(3)) / 2.d0
C u1 and u2 define the plan normal to omga
      nt = 1, npt
      b(1) = spin(1,nt)
      b(2) = spin(2,nt)
      b(3) = spin(3,nt)
      call pvect(omga,b,
     >                  c)
      u1(1) = c(1)
      u1(2) = c(2)
      u1(3) = c(3)
      call pvect(u1,omga,
     >                   c)
      u2(1) = c(1)
      u2(2) = c(2)
      u2(3) = c(3)

      write(*,*) 'u1.u2 = ',u1(1)*u2(1)+u1(2)*u2(2)+u1(3)*u2(3)
      write(*,*) 'u1.u2 = ',u1(1)*u2(1)+u1(2)*u2(2)+u1(3)*u2(3)
      write(*,*) 'u1.u2 = ',u1(1)*u2(1)+u1(2)*u2(2)+u1(3)*u2(3)
      if(npt.le.2) goto 22
c      write(*,fmt='(a,3(1x,i6,a))') 'Fourier-transform of pass # '
c     >,kpa,' to pass # ', kpb,'  (',ksmpl,' points)'
      
        ANUI = BORNE(1)
        ANUF = BORNE(2)
        DELNU=(ANUF - ANUI) / NC0
        PAS=DEUXPI * DELNU
        VAL=DEUXPI *(ANUI - 0.5d0 * DELNU)
        PMAX=0.D0
        PMIN=1.D12
        DO NC=1,NC0
          VAL=VAL+PAS
          SR=0.D0
          SI=0.D0
          DO NT=1,npt
            SR=SR + u1(ksx,nt)*cos(nt*val)
            SI=SI + u1(ksx,nt)*sin(nt*val)
          ENDDO
          PP=SR*SR+SI*SI
          IF(PP.GT. PMAX) THEN
            PMAX=PP
            KMAX=NC
          ELSEIF(PP.LT. PMIN) THEN
            PMIN=PP
          ENDIF
          SPEC(NC)=PP
          Q(NC)=val/(2.d0*pi)
        ENDDO

        IF (PMAX .GT. PMIN) THEN
          IF (KMAX .LT. mxc) THEN
            DEC=0.5D0 * (SPEC(KMAX-1)-SPEC(KMAX+1))
     >      /(SPEC(KMAX-1) - 2.D0 *SPEC(KMAX)+SPEC(KMAX+1))
          ELSE
            DEC=0.5D0 
          ENDIF
          YNU= ANUI + (DBLE(KMAX)+DEC-0.5D0) * DELNU
        ELSE
           YNU = 0.D0
        ENDIF

      write(*,*) '---------------------------------'
      write(*,*) '---------------------------------'
      write(*,*) ' Number of channels : ',nc0
      write(*,fmt='(3(a,i6))') 
     >  ' Number of turns analyzed : ',npt,', from ',kpa,' to ',kpb
      write(*,*) ' Qs / 1-Qs = ',ynu,' / ',1.d0-ynu
      write(*,*) '---------------------------------'
      write(*,*) '---------------------------------'
      
      Qs(1) = ynu  
      write(lunW,fmt='(a,1x,1p,4(1x,e14.5),2(1x,i),5x,a)') 
     >'  ',Qs(1),omga(1), omga(2), omga(3), kpa, kpb
     >,' ! Qs at max amplitude ; omga_x,y,z ; turn# range'
      write(lunW2,fmt='(a)') '# Spectrum (Qs / amp. / bin#) : '
      write(lunW2,fmt='(1p,2(e16.8,1x),i6)') (q(j),spec(j),j,j=1,nc0)


 22   continue
      write(*,*) ' Job completed !'
      close(lunW)
      close(lunW2)
      close(lunR)
      stop
 
 11   continue
      stop 'Error during read from spinTuneFromFai.In'
 98   continue
      stop 'Error during read header in zgoubi.fai'
 99   continue
      stop 'Error during read in zgoubi.fai'

 110     FORMAT(1X,

C1               KEX,      XXXO,(FO(J,IT),J=2,MXJ)
     >  1P,      1X,I2,    7(1X,E16.8)

C2              XXX,Y,T*1.D3,
     >         ,3(1X,E24.16)

C3       Z,P*1.D3,SAR,TAR 
     >   ,4(1X,E24.16)

C4       SXo, SYo, SZo, So, SX, SY, SZ, S
     >   ,8(1X,E15.7)

C5        ENEKI,ENERG
     >   ,2(1X,E16.8)

C6       IT,IREP(IT),   SORT(IT), (AMQ(J,I),J=1,5), RET(IT), DPR(IT), PS
     >   ,2(1X,I6),   9(1X,E16.8)

C7         BORO,   IPASS, NOEL,    
     >   ,1X,E16.8,  2(1X,I6)

C8            'KLEY',    ('LABEL(NOEL,I)',I=1,2),           'LET(IT)'
     >    ,1X,A1,A10,A1,   2(1X,A1,A8,A1),                   1X,3A1)
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
      SUBROUTINE HEADER(NL,N,BINARY,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL BINARY
      CHARACTER*80 TITL
      CHARACTER*270 TXT
      WRITE(6,FMT='(/,A)') ' File header : '
      IF(.NOT.BINARY) THEN
        READ(NL,FMT='(A)',ERR=99,END=99) TXT
        WRITE(6,FMT='(A)') TXT
        READ(NL,FMT='(A)',ERR=99,END=99) TITL
        WRITE(6,FMT='(A)') TITL
      ELSE
        READ(NL,ERR=99,END=89) TXT
        WRITE(6,FMT='(A)') TXT
        READ(NL,ERR=99,END=89) TITL
        WRITE(6,FMT='(A)') TITL
      ENDIF
      IF(.NOT.BINARY) THEN
        DO 1 I=3, N
           READ(NL,FMT='(A)',ERR=99,END=89) TXT
           WRITE(6,FMT='(A)') TXT
 1      CONTINUE
      ELSE
        DO 2 I=3, N
           READ(NL,          ERR=99,END=89) TXT
           WRITE(6,FMT='(A)') TXT
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
      subroutine pvect(a,b,
     >                     c)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension a(3),b(3),c(3)
      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)
      return
      end
      
