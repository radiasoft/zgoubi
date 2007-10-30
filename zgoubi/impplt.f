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
      SUBROUTINE IMPPLT(LN,DS,Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-------- Store into zgoubi.plt. 
C         Post-processing of stored data possible with zpop. 
      COMMON/AIM/ AE,AT,AS,RM,XI,XF,EN,EB1,EB2,EG1,EG2
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      COMMON/CHAMBR/ LIMIT,IFORM,YLIM2,ZLIM2,SORT(MXT),FMAG,BMAX
     > ,YC,ZC
      COMMON/CHAVE/ B(5,3),V(5,3),E(5,3)
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),IMAX,IEX(MXT),IREP(MXT)
      CHARACTER LET
      COMMON/FAISCT/ LET(MXT)
      COMMON/MARK/ KART,KALC,KERK,KUASEX
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      LOGICAL ZSYM
      COMMON/OPTION/ KFLD,MG,LC,ML,ZSYM
      COMMON/PTICUL/ AM,Q,G,TOO
      COMMON/REBELO/ NPASS,IPASS,KWRT,NNDES,STDVM
      COMMON/RIGID/ BORO,DPREF,DP,BR
      INCLUDE 'MXFS.H'
      CHARACTER FAM*8,LBF*8,KLEY*10,LABEL*8
      COMMON/SCALT/ FAM(MXF),LBF(MXF,2),KLEY,LABEL(MXL,2)
      COMMON/SYNCH/ RET(MXT), DPR(MXT),PS
 
      LOGICAL BINARY,BINARI
      SAVE BINARY

      CALL INTEG4(NSTEP)
      K = LSTK()

C----- Case  IL=2*10**n with n>1 (IL=20, 200, 2000...)
      IF(1+K*((NSTEP-1)/K) .NE. NSTEP) RETURN

      ENTRY IMPPLA(LN,DS,Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT)

C------- Bx,y,z (kG)
      BX=B(1,1)  *BR 
      BY=B(1,2)  *BR
      BZ=B(1,3)  *BR

      EX=E(1,1)  *BR
      EY=E(1,2)  *BR
      EZ=E(1,3)  *BR

C-----------------------------------------------
C  Assign XXX

C       P0 = BR*CL*1.D-9 *Q/QE
C       w=sqrt(p0*p0+am*am) -am
C       XXX=W

C      XXX=BR/BORO
      XXX= -1.D0 + BR/BORO
      XXXO=-1.D0 + FO(1,IT)
C----- Tests installation miroirs elec
C      CALL ENRGY(XXX)

C-----------------------------------------------
      DY = 0.D0
CCCC test spiral injection      IF(KART.EQ.2) DY = RM
   
      GOTO 10

      ENTRY IMPPLB(LN,DS,Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT)

      BX=0.D0
      BY=0.D0
      BZ=0.D0
      EX=0.D0
      EY=0.D0
      EZ=0.D0
      XXX =-1.D0 + F(1,IT)
      XXXO=-1.D0 + FO(1,IT)
      DY = 0.D0

 10   CONTINUE

      IF(BINARY) THEN
        WRITE(LN)
     1   LET(IT),KEX,      XXXO,(FO(J,IT),J=2,MXJ),
     2   XXX,Y-DY,T*1.D3,
     3   Z,P*1.D3,SAR,     TAR,     DY, 
C     3   Z,P*1.D3,SAR,     TAR,     DS, 
     4   KART, IT,IREP(IT),SORT(IT),X, BX,BY,BZ, RET(IT), DPR(IT),
     5   EX,EY,EZ,  BORO,  IPASS, KLEY,(LABEL(NOEL,I),I=1,2),NOEL
      ELSE
        WRITE(LN,100)
     1   LET(IT),KEX,      XXXO,(FO(J,IT),J=2,MXJ),
C       Initial coordinates: D,  Y,T,Z,P,X,Time

     2   XXX,Y-DY,T*1.D3,
     3   Z,P*1.D3,SAR,     TAR,     DY, 
C     3   Z,P*1.D3,SAR,     TAR,     DS, 
C     current coordinates  time    Step 
C                          mu_s    size 

     4   KART,  IT,IREP(IT),SORT(IT),X, BX,BY,BZ, RET(IT), DPR(IT),
C        Cart.              Path out     - kG -     (S)   dp/p_Synchro
C         or                 CHAMBR                Synchrotron
C        Polar                                       motion

     5   EX,EY,EZ,   BORO,  IPASS, KLEY,(LABEL(NOEL,I),I=1,2),NOEL
C        - eV/cm -  reference         
C                 rigidity(kG.cm)

        INCLUDE "FRMPLT.H"
      ENDIF 
C      CALL FLUSH2(LN,.FALSE.)

      RETURN

      ENTRY IMPPL2(BINARI)
      BINARY=BINARI
      RETURN
      END
