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
      SUBROUTINE IMPCSR(DS,Y,T,Z,P,X,SAR,TAR,IT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      COMMON/CHAMBR/ LIMIT,IFORM,YLIM2,ZLIM2,SORT(MXT),FMAG,BMAX
     > ,YC,ZC
      COMMON/CHAVE/ B(5,3),V(5,3),E(5,3)
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXSTEP.H'
      INCLUDE 'CSR.H'
C      COMMON/CSR/ KTRA,KCSR,YZXB(MXSTEP,41,36),DWC(MXT)
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      CHARACTER LET
      COMMON/FAISCT/ LET(MXT)
      COMMON/MARK/ KART,KALC,KERK,KUASEX
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      LOGICAL ZSYM
      COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      COMMON/PTICUL/ AM,Q,G,TOO
      COMMON/REBELO/ NPASS,IPASS,KWRT,NNDES,STDVM
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      INCLUDE 'MXFS.H'
      PARAMETER (LBLSIZ=10)
      PARAMETER (KSIZ=10)
      CHARACTER FAM*(KSIZ),LBF*(LBLSIZ)
      COMMON/SCALT/ FAM(MXF),LBF(MXF,MLF)
      COMMON/SYNCH/ RET(MXT), DPR(MXT),PS
      COMMON/UNITS/ UNIT(MXJ)

C----- READ STEP #
      CALL INTEG4(NSTEP)
      K = LSTK()

C----- Case  IL=2*10**n with n>1 (IL=20, 200, 2000...)
      IF(1+K*((NSTEP-1)/K) .NE. NSTEP) RETURN

      YZXB(NSTEP,IT,9) = DS   * UNIT(5)      ! Size of next step (m)
      YZXB(NSTEP,IT,6) = SAR  * UNIT(5)        ! Path length (m)
      YZXB(NSTEP,IT,30) = B(1,1)     * .1D0           ! Bx, Tesla
      YZXB(NSTEP,IT,31) = B(1,2)     * .1D0           ! By
      YZXB(NSTEP,IT,32) = B(1,3)     * .1D0           ! Bz
      YZXB(NSTEP,IT,3) = T * UNIT(2)
      YZXB(NSTEP,IT,5) = P * UNIT(4)
      YZXB(NSTEP,IT,7) = X   * UNIT(5)         ! m
      YZXB(NSTEP,IT,2) = Y   * UNIT(1)         ! m
      YZXB(NSTEP,IT,4) = Z   * UNIT(3)         ! m
      YZXB(NSTEP,IT,36) = BORO * 1.D-3       ! T.m
      YZXB(NSTEP,IT,1) = -1.D0 + QBR/(Q*BORO)     ! -1+p/p0
      YZXB(NSTEP,IT,8) = TAR  * 1.D-11       ! s
C            PP = BORO*CL9 *F(1,IT) * AMQ(2,IT)
C            ENERG = SQRT(PP*PP + AMQ(1,IT)*AMQ(1,IT))
C            ENEKI = ENERG - AMQ(1,IT)
      RETURN
      END
