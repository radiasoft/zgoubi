C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL OKOPN, CHANGE
      COMMON/LUN/ NDAT,NRES,NPLT,NFAI,NMAP,NSPN

      CHARACTER*20 NOMFIC
      logical exs
      LOGICAL IDLUNI, ok

      write(*,*) ' '
      write(*,*) '------------------------- '
      write(*,*) 'Pgm averageS.  '
      write(*,*) '------------------------- '
      write(*,*) ' '

      call block 
      call INIGR(
     >           NLOG, LM, NOMFIC)
      nl = nfai
      okopn = .false.
      change = .true.

      INQUIRE(FILE='averageS.in',EXIST=EXS)
      if(exs) then
         ok = IDLUNI(lunin)         
         OPEN(UNIT=lunin,FILE='averageS.in')
         write(*,*) 'Data read from averageS.in :'
         read(lunin,*) nx,ny
         read(lunin,*) nomfic
         close(lunin)
      else
        nx = 59
c        ny = 23
        ny = 25
      endif

      write(*,*) ' '
      write(*,*) '------------------------- '
      write(*,*) 'nx, ny : ',nx,ny
      write(*,*) 'file : ',nomfic
      write(*,*) '------------------------- '
      write(*,*) ' '

      nl = nfai
      OPEN(UNIT=nfai,FILE=nomfic)

      i12 = 1
      call avrg(NLOG,i12,NL,LM,OKOPN,CHANGE,nx,ny)

      if (ny.eq.25) then
        i12 =2 
        call avrg(NLOG,i12,NL,LM,OKOPN,CHANGE,nx,ny)
      endif

      close(nfai)
      close(lunout)
      stop
      end

      SUBROUTINE avrg(NLOG,I12,NL,LM,OKOPN,CHANGE,nx,ny)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL OKOPN, CHANGE
      PARAMETER (NCANAL=2500)
      COMMON/SPEDF/BORNE(6),SPEC(NCANAL,3),PMAX(3),NC0(3)
      COMMON/TRACKM/ NPTS,NPTR

      DIMENSION YM(3), YPM(3), U(3), A(3), B(3), YNU(3)
      DIMENSION YMX(6), YPMX(6)
 
      LOGICAL OKECH
      CHARACTER REP, NOMFIC*20
      LOGICAL BINARY, BINARF
      CHARACTER HVL(3)*12

      SAVE NT
      DATA NT / -1 /

      LOGICAL OPN
      DATA OPN / .FALSE. /
      LOGICAL IDLUNI

      parameter (mxpass=999999)
      dimension xvar(mxpass), sum(mxpass,10), nbtraj(mxpass)

      DATA HVL / 'Horizontal', 'Vertical', 'Longitudinal' /
    
      call iraz(nbtraj,mxpass)
      if(i12.eq.1) call raz(sum,mxpass*10)

      NPTR=NPTS

      OKECH = .FALSE.

      IF(NT.EQ.-1) THEN
        CALL READC6B(1,NPTS)
      ELSE
        CALL READC6B(NT,NT)        
      ENDIF

c           write(*,*) ' npts ',npts,nt
c               read(*,*)


 6    CONTINUE
C          OPN = .FALSE.
        IF(.NOT. OPN) THEN
          IF(NT.EQ.-1) THEN
            IF (IDLUNI(IUN)) THEN
              OPEN(UNIT=IUN,FILE='averageS.out',ERR=699)
              OPN = .TRUE.
              write(iun,fmt='(a)')
     >   '# yzxb(nx), av sx sy sz, av |S|, av ctta stta, av tta tta^2, '
     >        //'sigTta,  # of traj at that pass'
            ELSE
              GOTO 698
            ENDIF
          ENDIF
        ENDIF

            CALL STORCO(I12,NL,LM,1  ,BINARY,nx,ny,
     >                              xvar, sum,nbtraj,npass)

        do ipass = 1, npass
          if(ny.eq.25) then 
            if(i12.eq.1) then 
              sum(ipass,1)= sum(ipass,1)/dble(nbtraj(ipass)) 
              sum(ipass,2)= sum(ipass,2)/dble(nbtraj(ipass)) 
              sum(ipass,3)= sum(ipass,3)/dble(nbtraj(ipass)) 
              sqrts = sqrt(
     >        sum(ipass,1)**2+sum(ipass,2)**2+sum(ipass,3)**2)

              sum(ipass,1)= sum(ipass,1)/sqrts
              sum(ipass,2)= sum(ipass,2)/sqrts
              sum(ipass,3)= sum(ipass,3)/sqrts

            else              

              sum(ipass,4)= sum(ipass,4)/dble(nbtraj(ipass))  ! avrg cos tta
              sum(ipass,5)= sum(ipass,5)/dble(nbtraj(ipass))  ! avrg sin tta
              sum(ipass,6)= sum(ipass,6)/dble(nbtraj(ipass))  ! avrg tta
              sum(ipass,7)= sum(ipass,7)/dble(nbtraj(ipass))  ! avrg tta^2
              xnrm=sqrt(sum(ipass,1)**2+sum(ipass,2)**2+sum(ipass,3)**2)
              write(iun,
     >          fmt='(1p,e17.8,1x,9e17.8,1x,i6,1x,i6,1x,e17.8)')
c    >          fmt='(1p,e17.8,1x,10e17.8,1x,i6,1x,i6,1x,e17.8,a)')
     >        abs(xvar(ipass)), 
     >        sum(ipass,1), sum(ipass,2), sum(ipass,3),               ! avrg S_x,_y,_z
     >        xnrm,  ! |av S|
     >        sum(ipass,4), sum(ipass,5),                       ! cos_tta, sin_tta
     >        sum(ipass,6), sum(ipass,7),                       ! av tta, av tta^2
     >        sqrt(sum(ipass,7) -   (sum(ipass,6))**2 )  ,      ! sig_tta
     >        nbtraj(ipass), ipass,
     >        sqrt(sum(ipass,7) -   (sum(ipass,6))**2 ) 
     >         * 180./(4.* atan(1.d0)) 
c     >         * 11.8393845 / (8.80238350E-02/(4.d0*atan(1.d0))*180.)
c     >        ,' yzxb(nx), av sx sy sz, av ctta stta, av tta tta^2, '
c     >        //'sigTta, sin2+cos2, # of traj at that pass'
            endif
          else
            write(iun,fmt='(1p,e17.8,1x,e17.8,1x,i6,1x,i6,a)')
     >      abs(xvar(ipass)), sum(ipass,1)/dble(nbtraj(ipass)), 
     >      nbtraj(ipass), ipass,
     >         ' yzxb(nx), avrge, # of traj at that pass'
          endif
        enddo

         if(ny.eq.25 ) then 
           if( i12.eq.2) then 
             CLOSE(IUN)
           endif
         else
           CLOSE(IUN)
         endif

      RETURN

 698  WRITE(6,*) ' *** Problem : No idle unit for averageS.out '
      GOTO 99
 699  WRITE(6,*) ' *** Problem at OPEN averageS.out '
      GOTO 99

 99   RETURN

      END

      SUBROUTINE STORCO(i12,NL,LM,KPS,BINARY,nx,ny,
     >                            xvar,sum,nbtraj,npass)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ---------------------------------------------------
C     Read coordinates from zgoubi output file, and store  
C     ---------------------------------------------------
      parameter (mxpass=999999)
      dimension xvar(mxpass), sum(mxpass,10), nbtraj(mxpass)

      LOGICAL BINARY
      COMMON/TRACKM/ NPTS,NPTR

      CHARACTER LET 

      INCLUDE 'MXVAR.H'
      DIMENSION YZXB(MXVAR),NDX(5)
      dimension a(3),b(3),vv(3)

      CALL REWIN2(NL,*96)
      WRITE(6,*) '  READING  AND  STORING  COORDINATES...'

      NOC=0
      NRBLT = -1 
C----- BOUCLE SUR READ FICHIER NL 
 44   CONTINUE

        CALL READCO(NL,
     >                 KART,LET,YZXB,NDX,*10,*44)

C----- NDX: 1->KEX, 2->IT, 3->IREP, 4->IMAX

        NOC=NOC+1

        IF(NINT(YZXB(39)) .GE. NRBLT+1) NRBLT = NINT(YZXB(39)) -1

        ipass = nint(yzxb(39))
        xvar(ipass) = yzxb(nx)
c         write(*,*) yzxb(nx), nx, yzxb(21), yzxb(22), yzxb(23), ny,i12
        if(ipass .gt. mxpass) stop ' Increase mxpass !! '

        if(ny.eq.25) then   !  average spin vector
          if(i12.eq.1) then 
            sum(ipass,1) = sum(ipass,1) + yzxb(21)
            sum(ipass,2) = sum(ipass,2) + yzxb(22)
            sum(ipass,3) = sum(ipass,3) + yzxb(23)
c             write(*,*) ' Sbr storco. i12 = 1; ipass, |s_i|^2 : ',   ! this quantities IS OK
c     >         ipass,
c     >      yzxb(21)**2 + 
c     >      yzxb(22)**2 + 
c     >      yzxb(23)**2
c                  read(*,*) 
          else
C cos angle between avrge \vec S and current \vec S
            a(1) = sum(ipass,1)
            a(2) = sum(ipass,2)
            a(3) = sum(ipass,3)
            xa = sqrt(pscal(a,a))

            b(1) =  yzxb(21)
            b(2) =  yzxb(22)
            b(3) =  yzxb(23)
            xb = sqrt(pscal(b,b))

            ctta = pscal(a,b)/xa/xb
            call vvect(a,b,vv)
            stta = sqrt(pscal(vv,vv))/xa/xb

c        write(*,*) ' |a|, |b| :',xa,xb,ctta**2+stta**2

            sum(ipass,4) = sum(ipass,4) + ctta
            sum(ipass,5) = sum(ipass,5) + stta
            tta = atan2(stta,ctta)
            tta = acos(ctta)
            sum(ipass,6) = sum(ipass,6) + tta
            sum(ipass,7) = sum(ipass,7) + tta*tta

c      write(*,*)' tta, stta, ctta, atg, acos, sum_tta, sum_tta2 : ',
c     >tta,stta, ctta, atan(stta/ctta),acos(ctta),(ctta**2 + stta**2)
c                  read(*,*)
          endif 
        else
          sum(ipass,1) = sum(ipass,1) + yzxb(ny)
        endif
        nbtraj(ipass) = nbtraj(ipass) + 1

C        IF(NOC.EQ. NPTR) GOTO 10

      GOTO 44             
C     ----------------------------------

 99   CONTINUE
      WRITE(6,*) ' *** Coordinates  storage  stopped: error during',
     > ' read of event # ',NOC+1
      GOTO 11

 10   CONTINUE
      WRITE(6,*) ' READ  OK; END  OF  FILE  ENCOUNTERED'

 11   CONTINUE
      NPASS = NRBLT + 1
      NPTR=NOC
      CALL READC5(KT1,KT2)
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

      SUBROUTINE READCO(NL,
     >                     KART,LET,YZXB,NDX,*,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ----------------------------------------------------
C     Look for and read coordinates, etc. of particle # NT
C     ----------------------------------------------------
      CHARACTER*1 LET
      INCLUDE 'MXVAR.H'
      DIMENSION YZXB(MXVAR),NDX(5)

      INCLUDE 'MXLD.H'
      COMMON/LABCO/ ORIG(MXL,6) 
      COMMON/LUN/ NDAT,NRES,NPLT,NFAI,NMAP,NSPN
      COMMON/TRACKM/NPTS,NPTR
      INCLUDE 'MAXCOO.H'
      COMMON/UNITS/ UNIT(MXJ-1) 
      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB

      PARAMETER (MXS=4)
      DIMENSION FO(MXJ),F(MXJ),SI(MXS),SF(MXS)
      
      PARAMETER (LBLSIZ=8)
      CHARACTER*(LBLSIZ) LBL1, LBL2

      PARAMETER (KSIZ=10)
      CHARACTER KLEY*(KSIZ) 
      CHARACTER TX1*1

      LOGICAL BINARY,BINAR,OKKP,OKKT,OKKL

      CHARACTER*1 KLET, KLETO, KLETI

      SAVE KP1, KP2, KP3, BINARY
      SAVE KL1, KL2
      SAVE KT1, KT2
      SAVE KKEX, KLET

      SAVE MOD, RFR, RFR2 
      SAVE NOEL1, NOC

      DATA MOD / 0 /
      DATA RFR, RFR2 / 0.D0, 0.D0 /

      DATA KP1, KP2, KP3 / 1, 999999, 1 /
      parameter (MXT=999999)
      DATA KT1, KT2 / 1, 999999 /
      DATA KL1, KL2 / 1, 999999 /
      DATA KKEX, KLET / 1, '*' / 

      DATA NOEL1, NOC / -1, 0 /

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
     >      BORO, IPASS, NOEL ,KLEY,LBL1,LBL2,LET

            jpass = ipass
            IF(.NOT. OKKT(KT1,KT2,IT,KEX,LET,
     >                                IEND)) GOTO 222

            IF(.NOT. OKKP(KP1,KP2,KP3,jPASS,
     >                                IEND))  THEN
C              IF(IEND.EQ.0) THEN
                                             GOTO 222
C              ELSEIF(IEND.EQ.1) THEN
C                GOTO 91
C              ENDIF
            ENDIF
            IF(.NOT. OKKL(KL1,KL2,NOEL,
     >                                IEND)) GOTO 222
            IF(IEND.EQ.1) GOTO 91

          ELSE
 21         READ(NL,110,ERR=99,END=10)
     >      KEX,(FO(J),J=1,7),
     >      (F(J),J=1,7), 
     >      (SI(J),J=1,4),(SF(J),J=1,4),
     >      ENEKI, ENERG, 
     >      IT, IREP, SORT, AMQ1,AMQ2,AMQ3,AMQ4,AMQ5, RET, DPR,  PS,
     >      BORO, IPASS,NOEL, 
     >      TX1,KLEY,TX1,TX1,LBL1,TX1,TX1,LBL2,TX1,TX1,LET,TX1

            INCLUDE "FRMFAI.H"

             jpass = ipass
            IF(.NOT. OKKT(KT1,KT2,IT,KEX,LET,
     >                                IEND)) GOTO 21
            IF(.NOT. OKKP(KP1,KP2,KP3,jPASS,
     >                                IEND))  THEN

C              IF(IEND.EQ.0) THEN
                                             GOTO 21
C              ELSEIF(IEND.EQ.1) THEN
C                GOTO 91
C              ENDIF
            ENDIF
            IF(.NOT. OKKL(KL1,KL2,NOEL,
     >                                IEND)) GOTO 21
            IF(IEND.EQ.1) GOTO 91

          ENDIF
C                write(*,*) ' averageS ',ENEKI, ENERG,ipass,jpass
               ipass = jpass
        ELSEIF(NL .EQ. NPLT) THEN
C--------- read in zgoubi.plt type storage file

          IMAX = 0
          IF(BINARY) THEN
 232         CONTINUE
            READ(NL,ERR=99,END=10) 
     >      KEX,(FO(J),J=1,7),
     >      (F(J),J=1,7), BTI, DS, 
     >      KART, IT, IREP, SORT, XX, BX, BY, BZ, RET, DPR, PS,
     >      (SI(J),J=1,4),(SF(J),J=1,4),
     >      EX, EY, EZ, BORO, IPASS, NOEL, KLEY,LBL1,LBL2,LET
C            IF(LM .NE. -1) THEN
C              IF(LM .NE. NOEL) GOTO 232
C            ENDIF

            IF(.NOT. OKKT(KT1,KT2,IT,KEX,LET,
     >                                IEND)) GOTO 232
            IF(.NOT. OKKP(KP1,KP2,KP3,IPASS,
     >                                IEND))  THEN
C              IF(IEND.EQ.0) THEN
                                             GOTO 232
C              ELSEIF(IEND.EQ.1) THEN
C                GOTO 91
C              ENDIF
            ENDIF
            IF(.NOT. OKKL(KL1,KL2,NOEL,
     >                                IEND)) GOTO 232
            IF(IEND.EQ.1) GOTO 91

          ELSE
 31         READ(NL,100,ERR=99,END=10)
     >      KEX,(FO(J),J=1,MXJ),
     >      (F(J),J=1,MXJ), BTI, DS,
     >      KART, IT, IREP, SORT, XX, BX, BY, BZ, RET, DPR, PS,
     >      (SI(J),J=1,4),(SF(J),J=1,4),
     >      EX, EY, EZ, BORO, IPASS, NOEL,
     >      TX1,KLEY,TX1,TX1,LBL1,TX1,TX1,LBL2,TX1,TX1,LET,TX1
            INCLUDE "FRMPLT.H"
CCCCCCCCCCC           if(it.eq.1) yref = f(2)

C            IF(LM .NE. -1) THEN
C              IF(LM .NE. NOEL) GOTO 31
C            ENDIF

            IF(.NOT. OKKT(KT1,KT2,IT,KEX,LET,
     >                                IEND)) GOTO 31
            IF(.NOT. OKKP(KP1,KP2,KP3,IPASS,
     >                                IEND))  THEN
C              IF(IEND.EQ.0) THEN
                                             GOTO 31
C              ELSEIF(IEND.EQ.1) THEN
C                GOTO 91
C              ENDIF
            ENDIF
            IF(.NOT. OKKL(KL1,KL2,NOEL,
     >                                IEND)) GOTO 31
            IF(IEND.EQ.1) GOTO 91

          ENDIF
        ENDIF        !NL = NFAI, NPLT

C------- dp/p
        J = 1
        JU = 6
        YZXB(J)   =   F(J)   * UNIT(JU)     
C        YZXB(J)   =   1.D0 + F(J)  
        YZXB(J+10) =  FO(J)   * UNIT(JU)        ! dp/p_initial

        DO J=2,MXJ
          JU = J-1
C------- J=2,7 : Y, T, Z, P, S, Time
          YZXB(J)   =  F(J)   * UNIT(JU)     
CCCCCCC          if(j.eq.2) YZXB(J) = (f(j)-yref) * UNIT(JU)
C------- J=2,7 : Y_0, ..., Time_0
          YZXB(J+10) = FO(J)  * UNIT(JU) 
         ENDDO

C------- KART=1 : Cartesian coordinates, X is current x-coordinate (normally 
C        ranging in [XI,XF] as defined in quasex.f)
C        KART=2 : Cylindrical coordinates, X is current angle (normally 
C        ranging in [XI,XF] as defined in aimant.f)
C        write(*,*) '  * * * * * * * * ',xx
        YZXB(8) = XX

        IF(KART .EQ. 1) THEN
          YZXB(8) = YZXB(8) * UNIT(5)
        ELSE
C Rustine  RACCAM pour plot avec ffag-spi
          if(noel.ne.noel1) then
            if(kley .eq. 'FFAG-SPI') noc = noc+1
            noel1 = noel
          endif
          nbCell = 10
          pnCell = 4.d0 * atan(1.d0) / DBLE(nbCell)
c            write(*,*) '  noc, yzxb(8) ', noc, yzxb(8)
          YZXB(8) = YZXB(8) + pnCell * (2.d0*noc -1.d0) 
C          YZXB(2) = YZXB(2) +   DY * UNIT(5) 
        ENDIF

C         step size :
        YZXB(9) = DS       * UNIT(5)
C         r = sqrt(y^2+z^2) :
        YZXB(10) = SQRT(YZXB(2)*YZXB(2) + YZXB(4)*YZXB(4))
        YZXB(18) = RET
C------- (p_ps)/ps
        YZXB(19) = DPR            
C-------- momentum
C        YZXB(19) = BORO * (1.D0+F(1))*0.299792458D0   
        YZXB(20) = ENEKI
        YZXB(21) = SF(1)
        YZXB(22) = SF(2)    
        YZXB(23) = SF(3)    
        YZXB(24) = SQRT(SF(1)**2 + SF(2)**2 + SF(3)**2)
C        YZXB(24) = SF(4)
C         convert B from kG to T
        YZXB(30) = BX      * .1D0
        YZXB(31) = BY      * .1D0
        YZXB(32) = BZ      * .1D0
        YZXB(33) = SQRT(BX*BX + BY*BY +  BZ*BZ) * .1D0

        YZXB(34) = EX
        YZXB(35) = EY     !!!/YZXB(2)
        YZXB(36) = EZ 
        YZXB(37) = SQRT(EX*EX + EY*EY +  EZ*EZ)

C AMAG is magnyfying factor for plotting of element synoptic and trajectories
C        CALL INSY1(
C     >             AMAG)        
        PI2 = 2.D0*ATAN(1.D0)
        YINL = F(2)* UNIT(1)  
        ZINL = F(4)* UNIT(3)  
C FM, Dec. 05       XINL = XX* UNIT(5) - ORIG(NOEL,5)
        IF    (KART .EQ. 1) THEN
          XINL = XX* UNIT(5) + ORIG(NOEL,5)
        ELSEIF(KART .EQ. 2) THEN
C--------- Plot from tracking in polar frame
          XINL = XX
          X = XINL
          Y = YINL
          IF( (KX .EQ. 48 .AND. KY .EQ. 42) ) THEN
C----------- Plot in Lab X-Y
            IF    (MOD.NE.0) THEN
C             mod=22 for 150MeV FFAG field map
C             mod=20 for RACCAM spiral field map
C             write(*,*) mod,rfr,'  ploter'
              IF(MOD.NE.20) THEN
                Y = Y + RFR
                TEMP = X
                X = Y * SIN(TEMP) 
                Y = Y * COS(TEMP) - RFR2
              ENDIF
            ELSEIF(MOD.EQ.0) THEN
C Example : SPES3 using DIPOLE-M
C           FFAG-SPI
              Y = Y + RFR
              TEMP = X 
C-------------------------
C Rustine  RACCAM pour plot avec ffag-spi
              if(noel.ne.noel1) then
                if(kley .eq. 'FFAG-SPI') noc = noc+1
                noel1 = noel
              endif
              nbCell = 10
              pnCell = 4.d0 * atan(1.d0) / DBLE(nbCell)
C               write(*,*) '  noc, yzxb(8) ', noc, yzxb(8)
              temp = temp + pnCell * (2.d0*noc -1.d0) 
C-------------------------
              X = Y * SIN(TEMP) 
              Y = Y * COS(TEMP)  - RFR2
            ENDIF
            XINL = X - ORIG(NOEL,5)
            YINL = Y
          ENDIF
        ENDIF
        YZXB(44) = ZINL 
        PHI = ORIG(NOEL,6)  
        CT = COS(ORIG(NOEL,4)+PHI) 
        ST = SIN(ORIG(NOEL,4)+PHI)
        YZXB(48) = ( XINL*CT - YINL*ST) + ORIG(NOEL,1)
        YZXB(42) = ( XINL*ST + YINL*CT) + ORIG(NOEL,2)

        CONTINUE

      ENDIF ! NL = NSPN, NFAI, NPLT

C      Location about where particle was lost
      YZXB(38) = SORT * 1.D-2
      YZXB(39) = IPASS 
      YZXB(57) = NOEL
      YZXB(58) = IT
      YZXB(59) = ENERG / AMQ1 *AMQ3   ! G.gamma

C        write(*,*) ipass, ENERG , AMQ1 ,AMQ3   
            

C- For RACCAM design --------------------------
c       nCell = 8
c        pi = 4.d0 *atan(1.d0)
c      if(noel.ne.noel1) noc = noc+1
c      YZXB(62) = (Y+RFR) * SIN(XX + 2.d0 * pi / nCell * DBLE(noc-1))
c      YZXB(68) = (Y+RFR) * COS(XX + 2.d0 * pi / nCell * DBLE(noc-1))
C-----------------------------------------------

      NDX(1)=KEX
      NDX(2)=IT
      NDX(3)=IREP
      NDX(4)=IMAX
      NDX(5)=NOEL

      RETURN

 91   CONTINUE
      write(*,*) ipass,' AT 91,   readco'
      write(*,*) ipass,' AT 91,   readco'
      write(*,*) ipass,' AT 91,   readco'
      write(*,*) ipass,' AT 91,   readco'
      write(*,*) ipass,' AT 91,   readco'
      write(*,*) ipass,' AT 91,   readco'
      RETURN 1

C------------------ Pass # KP1 to KP2, ipass-modulo KP3
      ENTRY READC1(
     >              KP1O,KP2O,KP3O)
C------- Read pass #  KP1 to KP2 step KP3
      KP1O=KP1
      KP2O=KP2
      KP3O=KP3
      RETURN
C--
      ENTRY READC2(LN)
C------- Write pass #,  KP1 to KP2, ipass-modulo KP3
 12   WRITE(6,FMT='(''  Option status is now : KP1='',I6
     >  ,'',   KP2='',I6,'', ipass-modulo ='',I6)') KP1, KP2, KP3
        WRITE(6,FMT='(''    Expected data : '',
     >  /,10X,'' KP1>0,  KP2>=KP1 : will plot in range [KP1,KP2]''
     >  ,'',  ipass-modulo  KP3>0''
C     >  /,10X,''(ii) KP1=-1, KP2 > 0 : will plot all ipass-modulo KP2 ''
     >  )')
        WRITE(6,FMT='(/
     >  ,'' Enter desired values KP1>0, KP2>KP1, modulo-KP3>0''
     >  )')
        READ(LN,*,ERR=12) KP1NEW, KP2NEW, KP3NEW
        IF(KP1NEW.LE.0 .OR. KP2NEW.LE.0 .OR. KP3NEW.LE.0) GOTO 12
      GOTO 11
C--
      ENTRY READC2B(KP1W,KP2W,KP3W)
C------- Pass KP1 to pass KP2, KP-modulo
        KP1NEW=KP1W
        KP2NEW=KP2W
        KP3NEW=KP3W
 11     CONTINUE
        IF(KP1NEW.NE.0) KP1=KP1NEW
        IF(KP2NEW.NE.0) KP2=KP2NEW
        IF(KP3NEW.NE.0) KP3=KP3NEW
C          IF(KP1.NE.-1) THEN
C            CALL PLOT31        ! OKS set to .FALSE.
            WRITE(6,*)
            WRITE(6,FMT='(
     >      '' Warning : although possibly requested (Menu-7/3/11),'')')
            WRITE(6,FMT='(''      reset of S coordinate upon '')')
            WRITE(6,FMT='(''      multiturn plot is NOW inhibited'')')
C            WRITE(6,FMT='(''      (only compatible with KP1=-1) '')')
C          ENDIF
      RETURN
C----------------------------

C------------------ Element #, KL1 to KL2
      ENTRY READC3(
     >             KL1O,KL2O)
C------- Read lmnt #,  KL1 to KL2
      KL1O=KL1
      KL2O=KL2
      RETURN
C--
      ENTRY READC4(LN)
        KLA = KL1
        KLB = KL2
C------- Write lmnt #,  KL1 to KL2
        WRITE(6,FMT='('' Observation now at elements  KL1='',I6,
     >   ''to   KL2='',I6)') KL1, KL2
        WRITE(6,FMT='(''     Available options for elment are : '',
     >   /,10X,'' 0 < KL1 < KL2  : will plot within range [KL1,KL2]'', 
     >   /,10X,'' KL1=-1, KL2 > 0 : will plot every KL2 other pass '')')
        WRITE(6,FMT='(/,
     >        '' Enter desired values KL1, KL2  : '')')
        READ(LN,fmt='(2I9)',ERR=33) KL1NEW, KL2NEW
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
      ENTRY READC5(
     >             KT1O,KT2O)
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

      ENTRY READC9(
     >             KKEXO,KLETO)
        KKEXO=KKEX
        KLETO=KLET
      RETURN

      ENTRY READCA(KKEXI,KLETI)
        KKEX=KKEXI
        KLET=KLETI
      RETURN      

      ENTRY READCC(MODI,RFRI,RFR2I)
        MOD = MODI
        RFR = RFRI
        RFR2 = RFR2I
      RETURN

 10   continue
        noc = 0
        noel1 = 0
        RETURN 1
 99   continue
        noc = 0
        noel1 = 0
        RETURN 2      

      END

      SUBROUTINE BLOCK
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

C----- NUMERO DES UNITES LOGIQUES D'ENTREES-SORTIE
      COMMON/LUN/ NDAT,NRES,NPLT,NFAI,NMAP,NSPN
 
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

      DATA NDAT,NRES,NPLT,NFAI,NMAP,NSPN / 3, 4, 1, 2, 8, 9 /
 
      DATA CL, PI, DPI, RAD, DEG, QE, AH /
     >2.99792458D8 , 3.141592653589D0, 6.283185307178D0,
     > .01745329252D0 , 57.29577951D0, 1.60217733D-19, 6.626075D-34 /
 
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

      COMMON/TRACKM/ NPTS,NPTR

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
      DATA NC0/ 2000, 2000, 2000 /

      DATA NPTS / 999999 /

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

C----- zpop log unit
      NLOG = 30
C----- Input data file name
C      Normally, default is zgoubi.fai, .plt, .spn, .map...
C      NOMFIC = 'b_zgoubi.fai'
      NOMFIC = 'zgoubi.fai'

      RETURN
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
      CALL HEADER(NL,4,BINAR,*99)
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
      SUBROUTINE IRAZ(ITAB,N)
      DIMENSION ITAB(1)
      DO 1 I=1,N
 1      ITAB(I) = 0
      RETURN
      END
      SUBROUTINE HEADER(NL,N,BINARY,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL BINARY
      CHARACTER*270 TXT80
      WRITE(6,FMT='(/,A,/)') ' File header : '
      IF(.NOT.BINARY) THEN
        READ(NL,FMT='(A)',ERR=99,END=99) TXT80
        WRITE(6,FMT='(A)') TXT80
        READ(NL,FMT='(A)',ERR=99,END=99) TXT80
        WRITE(6,FMT='(A)') TXT80
      ELSE
        READ(NL,ERR=99,END=89) TXT80
        WRITE(6,FMT='(A)') TXT80
        READ(NL,ERR=99,END=89) TXT80
        WRITE(6,FMT='(A)') TXT80
      ENDIF
      IF(.NOT.BINARY) THEN
        DO 1 I=3, N
           READ(NL,FMT='(A)',ERR=99,END=89) TXT80
           WRITE(6,FMT='(A)') TXT80
 1      CONTINUE
      ELSE
        DO 2 I=3, N
           READ(NL,          ERR=99,END=89) TXT80
           WRITE(6,FMT='(A)') TXT80
 2      CONTINUE
      ENDIF
C      WRITE(6,*) ' Header has been read,  ok'
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
      SUBROUTINE VVECT(A,B,
     >                     C)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(*),B(*),C(*)
      C(1) = A(2)*B(3) - A(3)*B(2)
      C(2) = A(3)*B(1) - A(1)*B(3)
      C(3) = A(1)*B(2) - A(2)*B(1)
      RETURN
      END
      FUNCTION PSCAL(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(*),Y(*)
 
      PSCAL=X(1)*Y(1)+X(2)*Y(2)+X(3)*Y(3)
 
      RETURN
      END
