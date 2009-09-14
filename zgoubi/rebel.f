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
C  François Méot <meot@lpsc.in2p3.fr>
C  Service Accélerateurs
C  LPSC Grenoble
C  53 Avenue des Martyrs
C  38026 Grenoble Cedex
C  France
      SUBROUTINE REBEL(READAT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL READAT

      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      COMMON/CHAMBR/ LIMIT,IFORM,YLIM2,ZLIM2,SORT(MXT),FMAG,BMAX
     > ,YCH,ZCH
      PARAMETER (MXPUD=9,MXPU=1000)
      COMMON/CO/ FPU(MXPUD,MXPU),KCO,NPU,NFPU,IPU
      COMMON/DESIN/ FDES(7,MXT),IFDES,KINFO,IRSAR,IRTET,IRPHI,NDES
     >,AMS,AMP,AM3,TDVM,TETPHI(2,MXT)
C     >,AMS ,AMP,ENSTAR,BSTAR,TDVM ,TETPHI(2,MXT)
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
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
      CHARACTER*9 HMS
      CHARACTER*104 TXTBUF

      SAVE KREB3, KREB31

      SAVE KWRI6 

      DATA KREB3, KREB31 / 0, 0 /

      NRBLT = NINT(A(NOEL,1))
C----- Switch for print into zgoubi.res :
      KWRT = INT(A(NOEL,2)) 
C----- Switch for print to standard output :
      KWRI6=NINT((A(NOEL,2)-KWRT)*10)

C----- For multiturn injection
      KREB3 = NINT(A(NOEL,3))
C----- If A(NOEL,3)=99.xx, then KREB31=xx. For instance, KREB3=99.15 -> KREB31=15 for 16-turn injection
      KREB31 = NINT(100*(A(NOEL,3)-KREB3))

      IF(KWRI6 .NE. 0) THEN
        CALL TIME2(HMS)
        II = 10**(KWRI6-1)
        IF(II*(IPASS/II).EQ.IPASS) THEN
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
        CALL FLUSH(6)
      ENDIF
 
C-----  structure length ------------
      CALL SCUMS(0.D0)
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

C----- Average orbit---------------------------------------- 
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
            DO 21 J=1,4
              DO 21 I=1,IMAX
                SSP(J,I) = 0D0
 21         CONTINUE
          ENDIF
          DO 22 J=1,4
            DO 22 I=1,IMAX
              SSP(J,I)  = SSP(J,I) +SF(J,I)
 22       CONTINUE
        ENDIF
      ENDIF
C--------- endif spin tracking ----------------------------------
         
C--------- SR loss ----------------------------------
      IF(KSYN .EQ. 1) THEN
        IF(KREB3.NE.99) THEN
          DUM=SCAL0W(1.D0)
          CALL SYNPA0
        ENDIF
      ENDIF
C--------- endif SR loss ----------------------------------

      IF( IPASS .LT. NRBLT ) THEN
        LUN=ABS(NRES) 
        IF(LUN.GT.0) THEN
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
          IF(KSPN .EQ. 1) THEN
            IF(KREB3 .EQ. 99) THEN
              WRITE(LUN,126)
              WRITE(LUN,125) (I,(SSP(J,I)/IPASS,J=1,4 ) ,I=1,IMAX)
            ENDIF
          ENDIF
          CALL CNTOUR(
     >                NOUT)
          IF(NOUT.GT. 0) WRITE(LUN,107) NOUT
          CALL CNTNRR(
     >                NRJ)
          IF(NRJ .GT. 0) WRITE(LUN,108) NRJ
        ENDIF
 
        IF(IPASS .EQ. 1) THEN
          IF(NRBLT.GT.1) THEN
            IF(KWRT .EQ. 0) THEN
C------------- inihibit WRITE if KWRT=0 and more than 1 pass
              IF(NRES .GT. 0) NRES =-NRES
            ENDIF
            READAT = .FALSE.
          ENDIF
        ENDIF
 
        JJJ = 0
        DO III = 1, IMAX
           IF(IEX(III).LT.0) JJJ = JJJ+1
        ENDDO
           WRITE(LUN,*) '    SUM OVER IEX : ',JJJ

        IPASS=IPASS+1
        NOEL=0 

        RETURN
 
      ELSEIF(IPASS .EQ. NRBLT) THEN
C------- Last pass through structure will occur
        IF(KWRT .EQ. 0) THEN
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
          IF(KSPN .EQ. 1) THEN
            IF(KREB3 .EQ. 99) THEN
              WRITE(LUN,126)
              WRITE(LUN,125) (I,( SSP(J,I)/IPASS,J=1,4) ,I=1,IMAX)
            ENDIF
          ENDIF
          CALL CNTOUR(
     >                NOUT)
          IF(NOUT.GT. 0) WRITE(LUN,107) NOUT
          CALL CNTNRR(
     >                NRJ)
          IF(NRJ .GT. 0) WRITE(LUN,108) NRJ
 
          WRITE(LUN,104)
 104      FORMAT(/,128(1H*),//,128(1H*),//,128(1H*))
          WRITE(LUN,102) NRBLT+1
 102      FORMAT(//,5X,' Next  pass  is  #',I6
     >    ,' and  last  pass  through  the  optical  structure',/)
        ENDIF
 
        IPASS=IPASS+1
        NOEL=0 
        RETURN
 
      ELSEIF(IPASS .GT. NRBLT) THEN
C------- Last pass through REBELOTE 
        LUN=ABS(NRES)
        IF(LUN.GT.0) THEN
          WRITE(LUN,101) IPASS
 101      FORMAT(/,25X,'****  End  of  ''REBELOTE''  procedure  ****',//
     >     ,5X,' There  has  been ',I10,
     >              '  pass  through  the  optical  structure ',/)
 
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
 
          IF(KSPN .EQ. 1) THEN
            IF(KREB3 .EQ. 99) THEN
              WRITE(LUN,126)
 126          FORMAT(/,20X,' Average values of spin components :'
     >        ,//,24X,'<SX>',T37,'<SY>',T49,'<SZ>',T61,'<S>')
              WRITE(LUN,125) ( I,( SSP(J,I)/IPASS,J=1,4 ) ,I=1,IMAX)
 125          FORMAT(15X,I3,2X,1P,4G12.4)
            ENDIF
          ENDIF
 
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

      ENDIF

      RETURN

      ENTRY REBELR(
     >             KREB3O,KREB31O)
      KREB3O = KREB3
      KREB31O = KREB31
      RETURN

      END
