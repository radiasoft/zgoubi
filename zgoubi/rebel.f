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
      SUBROUTINE REBEL(READAT,KLE,LABEL,
     >                                  REBFLG,NOELRB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER(*) KLE(*)
      INCLUDE 'MXLD.H'
      CHARACTER(*) LABEL(MXL,2)
      LOGICAL READAT

      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      COMMON/CHAMBR/ LIMIT,IFORM,YLIM2,ZLIM2,SORT(MXT),FMAG,BMAX
     > ,YCH,ZCH
      PARAMETER (MXPUD=9,MXPU=1000)
      COMMON/CO/ FPU(MXPUD,MXPU),KCO,NPU,NFPU,IPU
      COMMON/CONST2/ ZERO, UN
      COMMON/DESIN/ FDES(7,MXT),IFDES,KINFO,IRSAR,IRTET,IRPHI,NDES
     >,AMS,AMP,AM3,TDVM,TETPHI(2,MXT)
C     >,AMS ,AMP,ENSTAR,BSTAR,TDVM ,TETPHI(2,MXT)
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER(80) TA
      PARAMETER (MXTA=45)
      COMMON/DONT/ TA(MXL,MXTA)
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IMAXD,IMAXT
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)
      COMMON/SYNRA/ KSYN
      COMMON/UNITS/ UNIT(MXJ)

      DIMENSION SSP(4,MXT)
      SAVE SSP
      CHARACTER(9) HMS
      CHARACTER(108) TXTBUF

      SAVE KREB3, KREB31, KREB4

      SAVE KWRI6 
      LOGICAL REBFLG

      LOGICAL OKPCKP
      SAVE OKPCKP

      PARAMETER (MXPRM=MXD)
      DIMENSION PARAM(MXPRM)
      SAVE PARAM

      save noela, noelb

      PARAMETER (LBLSIZ=10)
      CHARACTER(LBLSIZ) LBL1, LBL2
      PARAMETER (KSIZ=10)
      CHARACTER(KSIZ) KLEI
      LOGICAL FOUND
      LOGICAL EMPTY
      integer debstr, finstr

      DATA KREB3, KREB31, kreb4 / 0, 0, 0 /
      DATA OKPCKP / .FALSE. /
      DATA PARAM / MXPRM*0.D0 /
      data noela, noelb / 1, mxl  /

      NRBLT = NINT(A(NOEL,1))
C----- Switch for print into zgoubi.res :
      KWRT = INT(A(NOEL,2)) 
C----- Switch for print to standard output :
      KWRI6=NINT((A(NOEL,2)-KWRT)*10)
      KREB3 = NINT(A(NOEL,3))
      KREB4 = NINT(A(NOEL,4))
C----- For multiturn injection
C----- If A(NOEL,3)=99.xx, then KREB31=xx. For instance, KREB3=99.15 -> KREB31=15 for 16-turn injection
      KREB31 = NINT(100*(A(NOEL,3)-KREB3))

      IF(KREB4 .EQ. 1) THEN
        IF(IPASS .EQ.1 ) THEN
C            write(*,*) ' rebel empty ',EMPTY(TA(NOEL,1)),TA(NOEL,1)
          IF(.NOT. EMPTY(TA(NOEL,1))) THEN
            klei = TA(NOEL,1)
            nlbl = 0
            IF(.NOT. EMPTY(TA(NOEL,2))) THEN
              LBL1 = TA(NOEL,2)
              nlbl = 1
              IF(.NOT. EMPTY(TA(NOEL,3))) THEN
                LBL2 = TA(NOEL,3)
                nlbl = 2
              ENDIF
            ENDIF
            FOUND = .FALSE.
            NEL = 1
            DO WHILE (.NOT. FOUND .AND. NEL .LE. Noel)  
              IF(KLE(IQ(NEL))(debstr(KLE(IQ(NEL))):finstr(KLE(IQ(NEL))))
     >             .EQ. KLEI(debstr(KLEI):finstr(KLEI))) THEN
                IF    (NLBL .EQ. 0) THEN                
                  FOUND = .TRUE.
                ENDIF
              ENDIF
              nel = nel + 1
            ENDDO
            IF(.NOT. FOUND) THEN
              CALL ENDJOB('Sbr rebel. No such keyword in sequence',-99)
            ENDIF 
          ENDIF
        ENDIF
        KLM = A(NOEL,10)
        KPRM = A(NOEL,11)
        DO I = 1, NRBLT
          J = 20 + 20*((I-1)/20)
          PARAM(I) = A(NOEL,J +I -1)
          IF(I .GT. MXPRM) STOP ' SBR REBEL : Too many parameters under 
     >    ''REBELOTE''. Reduce it or increase MXPRM.'
        ENDDO

        A(KLM,KPRM) = PARAM(IPASS)

        IF(IPASS .EQ. NRBLT) THEN
C------- Last but one pass through structure
          IF(KWRT .NE. 1) THEN
C--------- reactive WRITE
            IF(NRES.LT.0) NRES=-NRES
          ENDIF

        ELSEIF(IPASS .EQ. NRBLT+1) THEN
C------- Last pass has been completed
C Now last occurence of REBELOTE => carry on beyond REBELOTE
          LUN=ABS(NRES)
          IF(LUN.GT.0) THEN
            WRITE(LUN,101) IPASS
 
C            IF(KSPN .EQ. 1) THEN
C              IF(KREB3 .EQ. 99) THEN
C                WRITE(LUN,126)
C                WRITE(LUN,125) ( I,( SSP(J,I)/IPASS,J=1,4 ) ,I=1,IMAX)
C              ENDIF
C            ENDIF
 
            CALL CNTOUR(
     >                NOUT)
            IF(NOUT.GT. 0) WRITE(LUN,107) NOUT,IMX
            CALL CNTNRR(
     >                NRJ)
            IF(NRJ .GT. 0) WRITE(LUN,108) NRJ,IMX
          ENDIF
          READAT = .TRUE.
          stop
        ENDIF
      ENDIF

C Will stop at element # NOELB when doing last turn
      REBFLG = NOELB .LT. NOEL

      IF(KWRI6 .NE. 0) THEN
        CALL TIME2(HMS)
        II = 10**(KWRI6-1)
        IF(IPASS .EQ. 1 .OR. II*(IPASS/II).EQ.IPASS) THEN
          CALL CNTSTO(
     >                NSTOP)
          WRITE(TXTBUF,FMT='(A20,I8,A1,I8,A34,I7,A1,I7,2A9)') 
     >    ' Pass #/Requested : ',IPASS,'/', NRBLT+1,
     >    '.  Particles remaining/launched = ',IMAX-NSTOP,'/',IMAX,
     >    '.  Time :', HMS
          CALL ARRIER(TXTBUF)
        ENDIF
        IF(IPASS.EQ.1) WRITE(6,FMT='(/)')
        IF(IPASS.EQ.NRBLT) WRITE(6,FMT='(/)')
      ENDIF
 
C-----  structure length ------------
      CALL SCUMS(ZERO)
C------------------------------------------------

C-----  Count # of particles last gone through optical structure ------------
C       Initialisation is done in OBJECT routines
C      IF(IPASS .EQ. 1) THEN
C        CALL CNTMXW(IMAX)
C      ELSE
      IF(IPASS .GT. 1) THEN
C------- If not multiturn tracking : 
        IF(KREB3.NE.99) CALL CNTMXT(IMAX)
      ENDIF
C-----  endif particle count ------------

C----- Pick-up signal ---------------------------------------- 
      IF(KCO .EQ. 1) THEN
        CALL AVORPR(NFPU,1)
        IF(NRES .GT. 0) CALL AVORPR(NRES,2)
      ENDIF

C------- In-flight decay ----------------------------------
      IF(IFDES .EQ. 1) THEN
        IF(IPASS .EQ. 1) THEN
C--------- 1st pass at REBELOTE
          STDVM = 0.D0
          NNDES = 0
        ENDIF
        STDVM = STDVM + TDVM*IMAX
        IF(KREB3.NE.99) NNDES = NNDES + NDES
      ENDIF
C------- endif In-flight decay ----------------------------------

C--------- spin tracking ----------------------------------
      IF(KSPN .EQ. 1) THEN
        IF(KREB3 .EQ. 99) THEN
C--------- multiturn
          IF(IPASS .EQ. 1) THEN
C----------- 1-er pass at REBELOTE
            DO J=1,4
              DO I=1,IMAX
                SSP(J,I) = 0D0
              ENDDO
            ENDDO
          ENDIF
          DO J=1,4
            DO I=1,IMAX
              SSP(J,I)  = SSP(J,I) +SF(J,I)
            ENDDO
          ENDDO
        ENDIF
      ENDIF
C--------- endif spin tracking ----------------------------------
         
C--------- SR loss ----------------------------------
      IF(KSYN .EQ. 1) THEN
        IF(KREB3.NE.99) THEN
c          DUM=SCAL0W(1.D0)
c          CALL SYNPA0
        ENDIF
      ENDIF
C--------- endif SR loss ----------------------------------

      IF( IPASS .LT. NRBLT ) THEN

        LUN=ABS(NRES) 
        IF(LUN.GT.0 .AND. KWRT.EQ.2) THEN
          WRITE(LUN,100) IPASS
 100      FORMAT(/,30X,'  -----  REBELOTE  -----',//
     >    ,5X,'End of pass # ',I8,' through the optical structure ',/)
          CALL CNTMXR(
     >                IMX)
          WRITE(LUN,103) IMX
          IF(KREB3.NE.99) THEN
            KNDES = NNDES
          ELSE
            KNDES = NDES
          ENDIF
          IF(IFDES .EQ. 1) WRITE(LUN,105) STDVM*UNIT(5)/IMX/IPASS, KNDES
C          IF(KSPN .EQ. 1) THEN
C            IF(KREB3 .EQ. 99) THEN
C              WRITE(LUN,126)
C              WRITE(LUN,125) (I,(SSP(J,I)/IPASS,J=1,4 ) ,I=1,IMAX)
C            ENDIF
C          ENDIF
          CALL CNTOUR(
     >                NOUT)
          IF(NOUT.GT. 0) WRITE(LUN,107) NOUT
          CALL CNTNRR(
     >                NRJ)
          IF(NRJ .GT. 0) WRITE(LUN,108) NRJ
        ENDIF
 
        IF(IPASS .EQ. 1) THEN
          IF(NRES .GT. 0) WRITE(NRES,FMT='(/,5X,
     >    ''Multiple pass, '', /, 
     >    10X,''from element # '',I5,'' : '',
     >    A,''/label1='',A,''/label2='',A,
     >    '' to REBELOTE '',''/label1='',A,''/label2='',A,/, 
     >    10X,''ending at pass # '',I7,'' at element # '',I5,'' : '',
     >    A,''/label1='',A,''/label2='',A,/)') 
     >    NOELA,KLE(IQ(NOELA)),LABEL(NOELA,1),LABEL(NOELA,2),
     >    LABEL(NOELA,1),LABEL(NOELA,2), NRBLT+1,
     >    NOELB,KLE(IQ(NOELB)),LABEL(NOELB,1),LABEL(NOELB,2)
C          WRITE(NRES,FMT='(/,5X,
C     >    ''Total nuber of passes will be : '',I7,/)') NRBLT+1

          IF(NRES .GT. 0)  THEN
            IF(KREB4 .EQ. 1) THEN
              WRITE(NRES,FMT='(/,5X,2(A,1x,I4),A)') 'Parameter #',KPRM
     >        ,' in element #',KLM,' will be modified at each pass. '
              WRITE(NRES,FMT='(/,5X,A)') 'List of parameter values :'
              WRITE(NRES,FMT='(15X,i3,1p,e17.8)')(i,param(i),i=1, nrblt)
            ENDIF
          ENDIF

          IF(NRBLT.GT.1) THEN
            IF(KWRT .NE. 1) THEN
C------------- inihibit WRITE if KWRT.NE.1 and more than 1 pass
              IF(KREB4 .NE. 1)  THEN 
                IF(NRES .GT. 0) NRES =-NRES
              ENDIF
            ENDIF
CC----------- If not FIT :
C            IF(KREB3.NE.22) READAT = .FALSE.
            READAT = .FALSE.
          ENDIF
          IF(REBFLG) NOELRB = NOEL
        ENDIF
 
C        JJJ = 0
C        DO III = 1, IMAX
C           IF(IEX(III).LT.0) JJJ = JJJ+1
C        ENDDO
C           WRITE(LUN,*) '    SUM OVER IEX : ',JJJ

        IPASS=IPASS+1
        NOEL=NOELA-1

        IF(OKPCKP) CALL PCKUP3(NOELA)

C--------- SR loss ----------------------------------
        IF(KSYN .EQ. 1) THEN
          WRITE(LUN,FMT='(/,2X,
     >    '' * Theoretical S.R. parameters in BEND and MULTIPOL ''
     >    ,''*dipole*  field :'')')
          CALL SYNPA3(LUN,
     >                  SMELPP,EE)
          WRITE(LUN,FMT='(5X,
     >    '' Particle E / Radiated energy per turn : '',1P,G16.8,/,
     >    '',   time : '',G16.8,'' mu_sec'')') 
     >    EE/SMELPP,EE/SMELPP*F(7,1)
        ENDIF
C--------- endif SR loss ----------------------------

      ELSEIF(IPASS .EQ. NRBLT) THEN
C------- Last but one pass through structure
        IF(KWRT .NE. 1) THEN
C--------- reactive WRITE
          IF(NRES.LT.0) NRES=-NRES
        ENDIF

        LUN=ABS(NRES)
        IF(LUN.GT.0) THEN
          WRITE(LUN,100) IPASS
          CALL CNTMXR(
     >                IMX)
          WRITE(LUN,103) IMX
          IF(KREB3.NE.99) THEN
            KNDES = NNDES
          ELSE
            KNDES = NDES
          ENDIF
          IF(IFDES .EQ. 1) WRITE(LUN,105) STDVM*UNIT(5)/IMX/IPASS, KNDES
C          IF(KSPN .EQ. 1) THEN
C            IF(KREB3 .EQ. 99) THEN
C              WRITE(LUN,126)
C              WRITE(LUN,125) (I,( SSP(J,I)/IPASS,J=1,4) ,I=1,IMAX)
C            ENDIF
C          ENDIF
          CALL CNTOUR(
     >                NOUT)
          IF(NOUT.GT. 0) WRITE(LUN,107) NOUT
          CALL CNTNRR(
     >                NRJ)
          IF(NRJ .GT. 0) WRITE(LUN,108) NRJ
 
          WRITE(LUN,104)
 104      FORMAT(/,128('*'),//,128('*'),//,128('*'))
          WRITE(LUN,102) NRBLT+1
 102      FORMAT(//,5X,' Next  pass  is  #',I6
     >    ,' and  last  pass  through  the  optical  structure',/)

        ENDIF
 
        IPASS=IPASS+1
        NOEL=NOELA-1
        IF(OKPCKP) CALL PCKUP3(NOELA)
 
      ELSEIF(IPASS .EQ. NRBLT+1) THEN
C------- Last pass has been completed
C Now last occurence of REBELOTE => carry on beyond REBELOTE
        LUN=ABS(NRES)
        IF(LUN.GT.0) THEN
          WRITE(LUN,101) IPASS
 101      FORMAT(/,25X,'****  End  of  ''REBELOTE''  procedure  ****',//
     >     ,5X,' There  has  been ',I10,
     >              '  passes  through  the  optical  structure ',/)
 
          CALL CNTMXR(
     >                IMX)
          WRITE(LUN,103) IMX
 103      FORMAT(20X,' Total of ',I10,' particles have been launched')
 
          IF(KREB3.NE.99) THEN
            KNDES = NNDES
          ELSE
            KNDES = NDES
          ENDIF
          IF(IFDES .EQ. 1) WRITE(LUN,105) STDVM*UNIT(5)/IMX/NRBLT, KNDES
 105      FORMAT(20X,
     >    ' Average  life  distance  from  Monte Carlo :',F10.3,' m',/,
     >    20X,' Number  of  decays  in  flight  :',3I10)
C 105      FORMAT(20X,' LIBRE   TEMPS   DE   VOL   MOYEN      :',F10.3
C     >    ,' CM',/,20X,' NOMBRE  DE  DESINTEGRATIONS  EN  VOL  :',I10)
 
c          IF(KSPN .EQ. 1) THEN
c            IF(KREB3 .EQ. 99) THEN
c              WRITE(LUN,126)
C 126          FORMAT(/,20X,' Average values of spin components :'
C     >        ,//,24X,'<SX>',T37,'<SY>',T49,'<SZ>',T61,'<S>')
c              WRITE(LUN,125) ( I,( SSP(J,I)/IPASS,J=1,4 ) ,I=1,IMAX)
C 125          FORMAT(15X,I5,2X,1P,4(1X,G12.4))
c            ENDIF
c          ENDIF
 
          CALL CNTOUR(
     >                NOUT)
          IF(NOUT.GT. 0) WRITE(LUN,107) NOUT,IMX
 107      FORMAT(/,5X,' *** # of particles out of acceptance  :',
     >          I10,'/',I10)
          CALL CNTNRR(
     >                NRJ)
          IF(NRJ .GT. 0) WRITE(LUN,108) NRJ,IMX
 108      FORMAT(/,5X,' Number of particles stopped :',I10,'/',I10)
        ENDIF

        READAT = .TRUE.

C REBELOTE should be usable within FIT -> under developement. 
        IPASS = 1
        IF(OKPCKP) CALL PCKUP3(NOEL)

      ENDIF

      IF(KREB4 .EQ. 1) READAT = .FALSE.

      RETURN

      ENTRY REBELR(
     >             KREB3O,KREB31O,KREB4O)
      KREB3O = KREB3
      KREB31O = KREB31
      KREB4O = KREB4
      RETURN

      ENTRY REBEL2(KPCKUP)
      OKPCKP = KPCKUP .EQ. 1
      RETURN

      ENTRY REBEL5(
     >             PARAMO)
      PARAMO = PARAM(IPASS-1)
      if(nres.gt.0)  write(nres,*) ipass, param
      RETURN

      ENTRY REBEL6(NLAI, NLBI)
      NOELA = NLAI
      NOELB = NLBI
C      if(nres.gt.0)
C     >  write(nres,*) ' rebel6 noela, noelb ',noela,noelb
      RETURN

      ENTRY REBEL7(
     >             NLBO)
      NLBO = NOELB 
      RETURN

      END
