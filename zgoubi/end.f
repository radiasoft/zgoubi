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
      SUBROUTINE END(
     >               READAT,NOEL,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL READAT
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      PARAMETER (MXPUD=9,MXPU=1000)
      COMMON/CO/ FPU(MXPUD,MXPU),KCO,NPU,NFPU,IPU
      INCLUDE 'MXSTEP.H'
      INCLUDE "MAXTRA.H"
      INCLUDE 'CSR.H'
      COMMON/DESIN/ FDES(7,MXT),IFDES,KINFO,IRSAR,IRTET,IRPHI,NDES
     >,AMS,AMP,AM3,TDVM,TETPHI(2,MXT)
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      CHARACTER LET
      COMMON/FAISCT/ LET(MXT)
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      COMMON/UNITS/ UNIT(MXJ)
C      CHARACTER * 9   HMS

C----- Coherent synchrotron radiation
      IF(KCSR.EQ.1) THEN
        IF(IPASS .EQ. 1) THEN
          IF(NRES.GT.0) WRITE(NRES,100) IPASS, NRBLT+1
 100      FORMAT(/,30X,'  -----  CSR procedure  -----',//,5X,
     >     'End of pass # ',I1,'/',I1,
     >     ' through the optical structure',/)
          READAT = .FALSE.
          IPASS=IPASS+1
          NOEL=0 
          CALL SCUMS(0.D0)
          RETURN 1
        ELSEIF(IPASS .EQ. 2) THEN
C--------- CSR calculations done

          IF(NRES.GT.0) WRITE(NRES,101) IPASS          
 101      FORMAT(/,25X,' ****  End  of  ''CSR''  procedure  ****',//
     >     ,5X,' There  has  been ',I10,
     >          '  pass  through  the  optical  structure ',/)
            write(nres,fmt='(i6,1p,g12.4)') (it,dwc(it),it=1,imax)
            write(  88,fmt='(i6,1p,g12.4)') (it,dwc(it),it=1,imax)
        ENDIF
      ENDIF

C----- Average orbit---------------------------------------- 
      IF(KCO .EQ. 1) THEN
        CALL AVORPR(NFPU,1)
        IF(NRES .GT. 0) CALL AVORPR(NRES,2)
      ENDIF

        IF(IPASS.EQ.1) THEN
          CALL CNTMXR(
     >                IMX)
          WRITE(NRES,103) IMX
 103      FORMAT(/,20X,I10,' particles have been launched')
          NP = 0
          NS = 0
          NSURV = 0
          N1 = 0
          DO 10, I=1,IMX
            IF(IEX(I).GE.-1) THEN
              NSURV = NSURV + 1
              IF(LET(I).EQ.'S') THEN
                NS = NS+1
              ELSE
                NP = NP+1
              ENDIF
              IF(IEX(I).EQ.-1)  N1 = N1+1
            ENDIF
 10       CONTINUE
          WRITE(NRES,FMT='(20X,'' Made  it  to  the  end : '',I6)')NSURV
          IF(N1.NE.0) WRITE(NRES,FMT='(30X,
     >             '' (including '',I6,'' with KEX=-1'')') N1
          IF(IFDES .EQ. 1) THEN
            WRITE(NRES,FMT='(/,20X,'' ditributed  as  follows  :  '',
     >      I6,''  primaries,     '',I6,''  secondaries'')') NP,NS
            WRITE(NRES,105) TDVM*UNIT(5), NDES
 105        FORMAT(20X,' Average  life  distance  :',F10.3
     >      ,' m,    yielded ',I7,'  decays  in  flight')
          ENDIF
          CALL CNTOUR(
     >                NOUT)
          IF(NOUT.GT. 0) WRITE(NRES,107) NOUT,IMX
 107      FORMAT(/,5X,'     # of particles stopped by collimation :',
     >          I10,'/',I10)

          CALL CNTNRR(
     >                NRJ)
          IF(NRJ .GT. 0) WRITE(NRES,108) NRJ,IMX
 108      FORMAT(/,5X,'     # of particles stopped during ',
     >      'integration in field :', I10,'/',I10) 
        ENDIF
      RETURN
      END
