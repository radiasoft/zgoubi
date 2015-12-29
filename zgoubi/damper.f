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
      SUBROUTINE DAMPER
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      INCLUDE "C.CHAMBR.H"     ! COMMON/CHAMBR/ LIMIT,IFORM,YLIM2,ZLIM2,SORT(MXT),FMAG,YCH,ZCH
 
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
C      PARAMETER (LNTA=132) ; CHARACTER(LNTA) TA
C      PARAMETER (MXTA=45)
      INCLUDE "C.DONT.H"     ! COMMON/DONT/ TA(MXL,MXTA)
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      INCLUDE "C.MARK.H"     ! COMMON/MARK/ KART,KALC,KERK,KUASEX
C      LOGICAL ZSYM
      INCLUDE "C.TYPFLD.H"     ! COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      INCLUDE "C.PTICUL.H"     ! COMMON/PTICUL/ AM,Q,G,TOO
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      INCLUDE 'MXFS.H'
      INCLUDE 'MXSCL.H'
      INCLUDE "C.SCAL.H"     ! COMMON/SCAL/ SCL(MXF,MXS,MXSCL),TIM(MXF,MXS),NTIM(MXF),KSCL
C      COMMON/SCAL/SCL(MXF,MXS),TIM(MXF,MXS),NTIM(MXF),JPA(MXF,MXP),KSCL
      INCLUDE "C.SPIN.H"     ! COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)
 
      CHARACTER(10) TYP(2)
      PARAMETER (LBLSIZ=10)
      CHARACTER(LBLSIZ) NAMPU(2)
      DIMENSION BTAD(2), BTAPU(2,2)

      integer debstr

      DATA TYP / 'HORIZONTAL',' VERTICAL ' /
 
C     ... FACTEUR D'ECHELLE DES ChampS. UTILISE PAR 'SCALING'
      SCAL = SCAL0()
      IF(KSCL .EQ. 1) SCAL = SCAL0()*SCALER(IPASS,NOEL,
     >                                                 ZERO)
C     ...........................................................
C        KHV = 1/2 for x/y plane
C     ...........................................................
 
      
      KHV = NINT(A(NOEL,10))
      GAIN = A(NOEL,11)
      BTAD(1) = A(NOEL,12)
      BTAD(2) = A(NOEL,13)
      NPU = NINT(A(NOEL,20))
      DO IPU = 1, NPU
        NAMPU(I) = 
     >   TA(NOEL,IPU)(DEBSTR(TA(NOEL,IPU)):DEBSTR(TA(NOEL,IPU))+9)
      ENDDO

      BTAPU(1,1) = A(NOEL,30)
      BTAPU(1,2) = A(NOEL,31)
      BTAPU(2,1) = A(NOEL,32)
      BTAPU(2,2) = A(NOEL,33)

      call picku1(
     >            KPU)

      if(kpu.eq.0 .and. khv .ne. 0) 
     >  call endjob(' SBR damper : need switch PICKUP on ',-99)
       

      IF(KHV .NE. 0) THEN
        IF(NRES.GT.0) THEN
          WRITE(NRES,110) TYP(KHV), Gain, BTAD(1), BTAD(2)
 110      FORMAT(/,25X,' ******  DAMPER - zero length element ******'
     >          ,/,40X,A
     >          ,/,30X,' Gain                      = ',G12.5,'    '
     >          ,/,30X,' Local beta function, Y    = ',G12.5,'  m '
     >          ,/,30X,' Local beta function, Z    = ',G12.5,'  m ')
          do ipu = 1, npu 
            WRITE(NRES,111) ipu, NAMPU(ipu)
 111        FORMAT(/,25X,' Pick-up # ',I1,' : ',A)
            WRITE(NRES,112) btaPU(IPU,1:2)
 112        FORMAT(30X,' Beta function at PU       = ',G12.5,'  m ')
          enddo
        ENDIF 
 
      ELSEIF(KHV .EQ. 0.D0) THEN
 
        IF(NRES.GT.0) WRITE(NRES,109) 
 109    FORMAT(//,25X,' +++++++++++ DAMPER OFF ++++++++++++',//)
 
      ENDIF

c          write(*,*) ' damper '
c           read(*,*)


 
      IF(KHV .NE. 0) THEN
 
        DO 10 IT = 1,IMAX
           Y = F(2,IT)    ! x
           T = F(3,IT)    ! xp
           Z = F(4,IT)    ! y
           P = F(5,IT)    ! yp
            

               F(2,IT) =  y  
               F(3,IT) =  t   
               F(4,IT) =  z  
               F(5,IT) =  p   
            
 10     CONTINUE
 
      ELSEIF(IJ .EQ. 0) THEN
C------- DAMPER OFF
 
      ENDIF
 
      RETURN
      END
