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
        subroutine bb
C S. White & F. Meot, Jan. 2012
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      INCLUDE "C.REBELO.H"   ! COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      INCLUDE "C.UNITS.H"     ! COMMON/UNITS/ UNIT(MXJ)

        PARAMETER (EPSILON0 = 8.854187817D-12)
        PARAMETER (R0 = 1.535E-18)
        PARAMETER (GS = 1.79285D0)

        DOUBLE PRECISION  CHARGE,MASS,ENERGEV,INTENSITY
        DOUBLE PRECISION  ALFX,BETX,EPSNX,EPSX
        DOUBLE PRECISION  ALFY,BETY,EPSNY,EPSY
        DOUBLE PRECISION  CIRC,ALFMOM,DPP,SIGZ
        DOUBLE PRECISION  TUNEX,TUNEY,TUNEZ
        DOUBLE PRECISION  AMPX,AMPY,AMPZ
C        DIMENSION  SPINCOORD(3)
        DIMENSION  SIGMA(6)
        DOUBLE PRECISION  COEF

        PARAMETER (MXPT=10000)
        DIMENSION  PTSL(9,MXPT)
        DOUBLE PRECISION LINEARMAP(6,6) 

        DOUBLE PRECISION  GA, BET

        SAVE NPT,COEF,SIGMA,PTSL,GA
        SAVE LINEARMAP

        LOGICAL OK, IDLUNI

C        OPEN(1,FILE="ZGOUBBI.IN",STATUS="OLD")
C        READ(1,*)NTURNS
C        READ(1,*)Npt  
c        if(Npt.gt.mxpt) stop ' Npt > max allowed'
c        read(1,*)charge,mass,energev,intensity
c        read(1,*)alfx,betx,epsnx
c        read(1,*)alfy,bety,epsny
c        read(1,*)sigz,dpp
c        read(1,*)circ,alfmom
c        read(1,*)tunex,tuney,tunez
c        read(1,*)spincoord(1:3)
c        read(1,*)ampx,ampy,ampz
c        close(1)

      IF(NINT(A(NOEL,1)).EQ.0) THEN
          WRITE(NRES,*) ' BEAM-BEAM is OFF'
         
      ELSE
        IF(IPASS.EQ.1) THEN

       WRITE(NRES,*) '   intensity= ',A(NOEL,2)
       WRITE(NRES,*)  ' alfx= ',A(NOEL,10)
       WRITE(NRES,*)  ' betx= ',A(NOEL,11)
       WRITE(NRES,*)  ' epsnx= ',A(NOEL,12)
       WRITE(NRES,*)  ' alfy= ',A(NOEL,20)
       WRITE(NRES,*)  ' bety= ',A(NOEL,21)
       WRITE(NRES,*)  ' epsny= ',A(NOEL,22)
       WRITE(NRES,*)  ' sigz= ',A(NOEL,30)
       WRITE(NRES,*)  ' dpp= ',A(NOEL,31)
       WRITE(NRES,*)  ' circ= ',A(NOEL,40)
       WRITE(NRES,*)  ' alfmom= ',A(NOEL,41)
       WRITE(NRES,*)  ' tunex= ',A(NOEL,50)
       WRITE(NRES,*)  ' tuney= ',A(NOEL,51)
       WRITE(NRES,*)  ' tunez= ',A(NOEL,52)
       WRITE(NRES,*)  ' ampx= ',A(NOEL,60)
       WRITE(NRES,*)  ' ampy= ',A(NOEL,61)
       WRITE(NRES,*)  ' ampz= ',A(NOEL,62)

          MASS=AMQ(1,1)
          CHARGE=AMQ(2,1)
          P0 = BORO*CL9*CHARGE
          ENERGEV=.001D0 * SQRT(P0*P0 + MASS*MASS)
          INTENSITY=A(NOEL,2)
          ALFX=A(NOEL,10)
          BETX=A(NOEL,11)
          EPSNX=A(NOEL,12)
          ALFY=A(NOEL,20)
          BETY=A(NOEL,21)
          EPSNY=A(NOEL,22)
          SIGZ=A(NOEL,30)
          DPP=A(NOEL,31)
          CIRC=A(NOEL,40)
          ALFMOM=A(NOEL,41)
          TUNEX=A(NOEL,50)
          TUNEY=A(NOEL,51)
          TUNEZ=A(NOEL,52)
****        SPINCOORD(1)
          AMPX=A(NOEL,60)
          AMPY=A(NOEL,61)
          AMPZ=A(NOEL,62)

          GA = SQRT(ENERGEV**2+MASS**2)/MASS
          BET = SQRT(1.0-1.0/GA/GA)
          EPSX = EPSNX/GA/BET
          EPSY = EPSNY/GA/BET

          GX = (1.0+ALFX*ALFX)/BETX
          GY = (1.0+ALFY*ALFY)/BETY

          SIGMA(1) = BETX*EPSX
          SIGMA(2) = -ALFX*EPSX
          SIGMA(3) = GX*EPSX
          SIGMA(4) = BETY*EPSY
          SIGMA(5) = -ALFY*EPSY
          SIGMA(6) = GY*EPSY

          COEF = CHARGE*QE*INTENSITY*(1.0+BET**2) 
     >          /(GA*BET*(BET+BET)*2*PI*EPSILON0*MASS)

          CALL BBLMAP(LINEARMAP,TUNEX,TUNEY,TUNEZ,
     >              ALFX,BETX,ALFY,BETY,ALFMOM,GA,BET,CL,CIRC)

        ENDIF
      ENDIF

      RETURN

      ENTRY BBKICK

      CALL BBKCK(NPT,COEF/2.0D0,SIGMA,PTSL,GA,GS)
      CALL BBKLM(PTSL,NPT,LINEARMAP)
        
      IF(IPASS.EQ.1) THEN

        OK = IDLUNI(
     >              LUNB)
        OPEN(LUNB,FILE="ptcl.out",STATUS='UNKNOWN',POSITION='APPEND',
     >                                        FORM='FORMATTED')
          WRITE(LUNB,*)NPT
      ENDIF

        DO I = 1, NPT
C           WRITE(LUNA,FMT='(2I6,6(1X,E15.6))') NPT,I,(PTSL(J,I),J=1,6)
          F(1,I) = PTSL(6,I)
          F(2,I) = PTSL(1,I)*1D2
          F(3,I) = PTSL(2,I)*1D3
          F(4,I) = PTSL(3,I)*1D2
          F(5,I) = PTSL(4,I)*1D3
CCCCCCC          F(6,I) = PTSL(5,I)*1D6
        ENDDO
C        CLOSE(LUNA)


        WRITE(LUNB,*) IPASS
        DO I =1, NPT
          SNORM = SQRT(PTSL(7,I)**2+PTSL(8,I)**2+PTSL(9,I)**2)
          WRITE(LUNB,222) PTSL(1:9,I),SNORM,I,IPASS
        ENDDO

C        CLOSE(LUNB)
222     FORMAT(10(1X,E17.9),2(1X,I6))
        RETURN
        END
