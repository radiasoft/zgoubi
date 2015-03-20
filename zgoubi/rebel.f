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

      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      INCLUDE "C.CHAMBR.H"     ! COMMON/CHAMBR/ LIMIT,IFORM,YLIM2,ZLIM2,SORT(MXT),FMAG,YCH,ZCH
 
      PARAMETER (MXPUD=9,MXPU=5000)
       INCLUDE "C.CO.H"     ! COMMON/CO/ FPU(MXPUD,MXPU),KCO,NPU,NFPU,IPU
      INCLUDE "C.CONST2.H"     ! COMMON/CONST2/ ZERO, UN
      INCLUDE "C.DESIN.H"     ! COMMON/DESIN/ FDES(7,MXT),IFDES,KINFO,IRSAR,IRTET,IRPHI,NDES
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
      INCLUDE "C.INTEG.H"     ! COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      INCLUDE "C.OBJET.H"     ! COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      INCLUDE "C.REBELO.H"   ! COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      INCLUDE "C.SPIN.H"     ! COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)
      INCLUDE "C.SYNRA.H"     ! COMMON/SYNRA/ KSYN
      INCLUDE "C.UNITS.H"     ! COMMON/UNITS/ UNIT(MXJ)

      DIMENSION SSP(4,MXT)
      SAVE SSP
      CHARACTER(9) HMS
      CHARACTER(108) TXTBUF

      SAVE KREB3, KREB31, KREB4

      SAVE KWRI6 
      LOGICAL REBFLG

      LOGICAL OKPCKP
      SAVE OKPCKP

      SAVE NOELA, NOELB

      PARAMETER (MXPRM=10, MXLST=4000)
      PARAMETER (LBLSIZ=10)
      CHARACTER(LBLSIZ) LBL1(MXPRM), LBL2(MXPRM)
      SAVE LBL1, LBL2
      PARAMETER (KSIZ=10)
      CHARACTER(KSIZ) KLEI(MXPRM)
      LOGICAL FOUND
      LOGICAL EMPTY
      INTEGER DEBSTR, FINSTR
      DIMENSION PARAM(MXPRM,MXLST), PARAMI(MXPRM,MXLST)
      CHARACTER(KSIZ) TPRM(MXPRM,3), TPRMI(MXPRM,3)
      SAVE PARAM, TPRM
      DIMENSION KLM(MXPRM), KPRM(MXPRM)
      SAVE KLM, KPRM
      DIMENSION NLBL(MXPRM)
      SAVE NLBL

      PARAMETER (MXPL= MXPRM*MXLST, MXPRM3= MXPRM*3)
      LOGICAL OKLSTP
      dimension aold(MXPRM)

C      logical fiting

      DATA KREB3, KREB31, KREB4 / 0, 0, 0 /
      DATA OKPCKP / .FALSE. /
      DATA PARAM / MXPL*0.D0 /
      DATA TPRM / MXPRM3*' ' /
      DATA NOELA, NOELB / 1, MXL /
      DATA OKLSTP / .FALSE. /
      DATA NPRM / 0 /

      NRBLT = NINT(A(NOEL,1))
C----- Switch for print into zgoubi.res :
      KWRT = INT(A(NOEL,2)) 
C----- Switch for print to standard output :
      KWRI6 = NINT((A(NOEL,2)-KWRT)*10)
      KREB3 = NINT(A(NOEL,3))
C----- For multiturn injection : 
C      If A(NOEL,3)=99.xx, then KREB31=xx. For instance, KREB3=99.15 -> KREB31=15 for 16-turn injection
      KREB31 = NINT(100*(A(NOEL,3)-KREB3))
C----- KREB4=1 allows changing parameter values prior to rebelote
      KREB4 = NINT(A(NOEL,4))

      IF(KREB4 .EQ. 1) THEN
        NPRM =  A(NOEL,10)
        IF(NPRM .GT. MXPRM) 
     >    CALL ENDJOB('SBR REBEL. Too many parameters, has to be .le. '
     >    ,MXPRM)
        IF(NRBLT .GT. MXLST) 
     >  CALL ENDJOB('SBR REBEL. Parameter list too large. Has to be < '
     >  //' NRBLT=',NRBLT)

        DO IPRM = 1, NPRM
          KLM(IPRM)   = NINT(A(NOEL,20+10*(IPRM-1)))
          KPRM(IPRM)  = NINT(A(NOEL,21+10*(IPRM-1)))
          AOLD(IPRM) = A(KLM(IPRM),KPRM(IPRM))
          A(KLM(IPRM),KPRM(IPRM)) = PARAM(IPRM,IPASS)
        ENDDO

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
C------- If not periodic sequence type of tracking (i.e., final coord .ne. initial coord) : 
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
        IF(LUN .GT. 0) THEN

          CALL OPTIO1(
     >                KWROF)
          IF(KWROF .NE. 1) THEN
            WRITE(LUN,100) IPASS
 100        FORMAT(/,30X,'  -----  REBELOTE  -----',//
     >      ,5X,'End of pass # ',I8,' through the optical structure ',/)
            CALL CNTMXR(
     >                  IMX)
            WRITE(LUN,103) IMX
            IF(KREB3.NE.99) THEN
              KNDES = NNDES
            ELSE
              KNDES = NDES
            ENDIF
            IF(IFDES .EQ. 1) 
     >        WRITE(LUN,105) STDVM*UNIT(5)/IMX/IPASS, KNDES
            CALL CNTOUR(
     >                NOUT)
            IF(NOUT.GT. 0) WRITE(LUN,107) NOUT
            CALL CNTNRR(
     >                  NRJ)
            IF(NRJ .GT. 0) WRITE(LUN,108) NRJ
          ENDIF
        ENDIF

        IF(IPASS .EQ. 1) THEN

          LUN=ABS(NRES)
          IF(LUN .GT. 0) THEN
            WRITE(LUN,FMT='(/,5X,
     >      ''Multiple pass, '', /, 
     >      10X,''from element # '',I5,'' : '',
     >      A,''/label1='',A,''/label2='',A,
     >      '' to REBELOTE '',''/label1='',A,''/label2='',A,/, 
     >      10X,''ending at pass # '',I7,'' at element # '',I5,'' : '',
     >      A,''/label1='',A,''/label2='',A,/)') 
     >      NOELA,KLE(IQ(NOELA)),LABEL(NOELA,1),LABEL(NOELA,2),
     >      LABEL(NOELA,1),LABEL(NOELA,2), NRBLT+1,
     >      NOELB,KLE(IQ(NOELB)),LABEL(NOELB,1),LABEL(NOELB,2)
          ENDIF

          IF(KREB4 .EQ. 1) THEN
C----- Will first change parameter values in zgoubi.dat, prior to rebelote.
            DO IPRM = 1, NPRM
              IF(.NOT. EMPTY(TPRM(IPRM,1))) THEN
                KLEI(IPRM) = TPRM(IPRM,1)
                NLBL(IPRM) = 0
                IF(.NOT. EMPTY(TPRM(IPRM,2))) THEN
                  LBL1(IPRM) = TPRM(IPRM,2)
                  NLBL(IPRM) = 1
                  IF(.NOT. EMPTY(TPRM(IPRM,3))) THEN
                    LBL2(IPRM) = TPRM(IPRM,3)
                    NLBL(IPRM) = 2
                  ENDIF
                ENDIF
              ENDIF
              FOUND = .FALSE.
              NEL = 1
              DO WHILE (.NOT. FOUND .AND. NEL .LE. NOEL)  
                IF(
     >           KLE(IQ(NEL))(DEBSTR(KLE(IQ(NEL))):FINSTR(KLE(IQ(NEL))))
     >           .EQ. KLEI(IPRM)(DEBSTR(KLEI(IPRM)):FINSTR(KLEI(IPRM))))
     >          THEN
                  IF    (NLBL(IPRM) .EQ. 0) THEN                
                    FOUND = .TRUE.
                  ENDIF
                ENDIF
                NEL = NEL + 1
              ENDDO
              IF(.NOT. FOUND) THEN
                CALL ENDJOB(
     >          'Sbr rebel. No such keyword in optical sequence',-99)
              ENDIF 
            ENDDO

            LUN=ABS(NRES)
            IF(LUN .GT. 0) THEN
              DO IPRM = 1, NPRM
                WRITE(LUN,FMT='(/,5X,2(A,1x,I4),A)') 
     >          'Parameter #',KPRM(iprm),' in element #',
     >          KLM(IPRM),' will be modified at each pass. '
                OKLSTP = NRBLT .LE.20    ! Just to avoid big listing 
                IF(OKLSTP) THEN
                  WRITE(LUN,FMT='(5X,A)') 'list of requested '
     >            //'parameter values :'
                  WRITE(LUN,FMT='(15X,I5,1P,E17.8)')
     >            (I,PARAM(IPRM,I),I=1, NRBLT)
                ENDIF
              ENDDO
            ENDIF

          ENDIF

C          IF(NRBLT.GT.1) READAT = .FALSE.

          IF(REBFLG) NOELRB = NOEL

        ENDIF  !IF(IPASS .EQ. 1) 
 
        LUN = NRES
        IF(LUN .GT. 0) THEN
          DO IPRM = 1, NPRM
            WRITE(6,fmt='(/,'' SBR rebel. At pass # '',I4,''/'',
     >      I4,''.  In element # '',I4,
     >      '',  parameter #'',I3,''  changed  to  '',
     >      1P,E16.8,''   (was  '',E16.8,'')'',/)') IPASS,NRBLT+1,
     >      KLM(IPRM),KPRM(IPRM),PARAM(IPRM,IPASS),AOLD(IPRM)
            WRITE(LUN,FMT='(/,'' Pgm rebel. At pass # '',I4,''/'',
     >      I4,''.  In element # '',I4,
     >      '',  parameter #'',I3,''  changed to  '',
     >      1P,E16.8,''   (was  '',E16.8,'')'')') IPASS,NRBLT+1,
     >      KLM(IPRM),KPRM(IPRM),PARAM(IPRM,IPASS),AOLD(IPRM)
          ENDDO
        ENDIF

        IPASS=IPASS+1
        NOEL=NOELA-1

        IF(OKPCKP) CALL PCKUP3(NOELA)

C--------- SR loss ----------------------------------
        IF(KSYN .EQ. 1) THEN
          LUN = NRES
          IF(LUN.GT.0) THEN 
            WRITE(LUN,FMT='(/,2X,
     >      '' * Theoretical S.R. parameters in BEND and MULTIPOL ''
     >      ,''*dipole*  field :'')')
            CALL SYNPA3(LUN,
     >                      SMELPP,EE)
            WRITE(LUN,FMT='(5X,
     >      '' Particle E / Radiated energy per turn : '',1P,G16.8,/,
     >      '',   time : '',G16.8,'' mu_sec'')') 
     >      EE/SMELPP,EE/SMELPP*F(7,1)
          ENDIF
        ENDIF
C--------- endif SR loss ----------------------------

        READAT = .FALSE.

      ELSEIF(IPASS .EQ. NRBLT) THEN
C------- Last but one pass through structure HAS JUST BEEN completed
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
          CALL CNTOUR(
     >                NOUT)
          IF(NOUT.GT. 0) WRITE(LUN,107) NOUT
          CALL CNTNRR(
     >                NRJ)
          IF(NRJ .GT. 0) WRITE(LUN,108) NRJ
 
          WRITE(LUN,102) NRBLT+1
 102      FORMAT(//,5X,' Next  pass  is  #',I6
     >    ,' and  last  pass  through  the  optical  structure',/)
          DO IPRM = 1, NPRM
            WRITE(6,fmt='(/,'' SBR rebel. At pass # '',I4,''/'',
     >      I4,''.  In element # '',I4,
     >      '',  parameter #'',I3,''  changed  to  '',
     >      1P,E16.8,''   (was  '',E16.8,'')'',/)') IPASS,NRBLT+1,
     >      KLM(IPRM),KPRM(IPRM),PARAM(IPRM,IPASS),AOLD(IPRM)
            WRITE(LUN,FMT='(/,'' Pgm rebel. At pass # '',I4,''/'',
     >      I4,''.  In element # '',I4,
     >      '',  parameter #'',I3,''  changed to  '',
     >      1P,E16.8,''   (was  '',E16.8,'')'')') IPASS,NRBLT+1,
     >      KLM(IPRM),KPRM(IPRM),PARAM(IPRM,IPASS),AOLD(IPRM)
          ENDDO

c          WRITE(LUN,104)
c 104      FORMAT(/,128('*'),//,128('*'),//,128('*'))

        ENDIF
 
        READAT = .FALSE.

        IPASS=IPASS+1
        NOEL=NOELA-1
        IF(OKPCKP) CALL PCKUP3(NOELA)
 
      ELSEIF(IPASS .EQ. NRBLT+1) THEN
C------- Last pass HAS JUST BEEN completed
C This is the last occurence of REBELOTE. Wiil carry on beyond REBELOTE

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

        NOEL = NOELB
 
C REBELOTE should be usable within FIT -> under developement. 

        IPASS = ipass+1
        IF(OKPCKP) CALL PCKUP3(NOEL)

      ENDIF

      NXTPAS = IPASS
C      IF(NXTPAS .EQ. 2) THEN
      IF(NXTPAS .LE. NRBLT) THEN
        IF(KWRT .NE. 1) THEN
          IF(NRES .GT. 0) THEN
            WRITE(NRES,FMT='(/,5X,''WRITE statements''
     >      ,'' to zgoubi.res are inhibeted from now on.'')')
            WRITE(NRES,201)
 201        FORMAT(/,132('*'))
            NRES =-NRES
          ENDIF
        ENDIF
      ELSEIF(NXTPAS .EQ. NRBLT+1) THEN
C------- Next pass is last pass through structure
        IF(KWRT .NE. 1) THEN
C--------- reactive WRITE
          KWRT = 1
          IF(NRES .GT. 0) WRITE(NRES,FMT='(/,5X,''WRITE statements''
     >    ,'' to zgoubi.res are re-established from now on.'')')
          IF(NRES .LT. 0) NRES=-NRES
        ENDIF
      ENDIF

      RETURN

      ENTRY REBELR(
     >             KREB3O,KREB31O,KREB4O)
      KREB3O = KREB3
      KREB31O = KREB31
      KREB4O = KREB4
      RETURN

      ENTRY REBEL1(
     >             KWRTO)
      KWRTO = KWRT
      RETURN

      ENTRY REBEL2(KPCKUP)
      OKPCKP = KPCKUP .EQ. 1
      RETURN

      ENTRY REBEL4(PARAMI,TPRMI)
      PARAM = PARAMI
      TPRM = TPRMI
      RETURN

      ENTRY REBEL6(NLAI, NLBI)
      NOELA = NLAI
      NOELB = NLBI
      RETURN

      ENTRY REBEL7(
     >             NLBO)
      NLBO = NOELB 
      RETURN

      END
