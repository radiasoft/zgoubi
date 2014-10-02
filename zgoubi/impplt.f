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
      SUBROUTINE IMPPLT(LN,Y,T,Z,P,X,SAR,TAR,DS,AMT,QT,KEX,IT)
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
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      CHARACTER(1) LET
      COMMON/FAISCT/ LET(MXT)
      PARAMETER (LBLSIZ=10)
      CHARACTER(LBLSIZ) LABEL
      COMMON /LABEL/ LABEL(MXL,2)
      COMMON/MARK/ KART,KALC,KERK,KUASEX
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      LOGICAL ZSYM
      COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      COMMON/PTICUL/ AM,Q,G,TOO
      COMMON/REBELO/ NPASS,IPASS,KWRT,NNDES,STDVM
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      PARAMETER (KSIZ=10)
      INCLUDE 'MXFS.H'
      CHARACTER FAM*(KSIZ),LBF*(LBLSIZ)
      COMMON/SCALT/ FAM(MXF),LBF(MXF,MLF)
      COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)
      COMMON/SYNCH/ RET(MXT), DPR(MXT),PS
 
      CHARACTER(KSIZ) KLEY

      CHARACTER TX1*1
      PARAMETER (TX1='''')

      LOGICAL BINARY,BINARI
      SAVE BINARY
      LOGICAL IDLUNI, OPN
      SAVE LST2, LUN

      CALL INTEG4(NSTEP)
      K = LSTK()

C----- Case  IL=2*10**n with n>1 (IL=20, 200, 2000...)
      IF(1+K*((NSTEP-1)/K) .NE. NSTEP) RETURN

      ENTRY IMPPLA(LN,Y,T,Z,P,X,SAR,TAR,DS,AMT,QT,KEX,IT)

C------- Bx,y,z (kG)
      BX=B(1,1)  /BRI
      BY=B(1,2)  /BRI
      BZ=B(1,3)  /BRI

      EX=E(1,1)  /BRI
      EY=E(1,2)  /BRI
      EZ=E(1,3)  /BRI

C-----------------------------------------------
C  Assign XXX

C       P0 = BR*CL*1.D-9 *Q/QE
C       w=sqrt(p0*p0+am*am) -am
C       XXX=W

C      XXX=BR/BORO
      XXX= -1.D0 + QBR/(QT*BORO)
      XXXO=-1.D0 + FO(1,IT)
C----- Tests installation miroirs elec
C      CALL ENRGY(XXX)

C-----------------------------------------------
      DY = 0.D0
CCCC test spiral injection      IF(KART.EQ.2) DY = RM
   
      GOTO 10

      ENTRY IMPPLB(LN,Y,T,Z,P,X,SAR,TAR,DS,AMT,QT,KEX,IT)

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

          PPI = BORO*CL9*QT
          EI = SQRT(PPI*PPI+AMT*AMT)
          BTI = PPI / EI

      CALL ZGKLEY(
     >            KLEY)

      IF(BINARY) THEN
        WRITE(LN)
     >   KEX,      XXXO, (FO(J,IT),J=2,MXJ),
     >   XXX, Y-DY, T*1.D3, Z, P*1.D3, SAR, TAR, BTI, DS, 
     >   KART, IT, IREP(IT), SORT(IT), X, BX,BY,BZ, RET(IT),DPR(IT), PS,
     >   (SI(J,IT),J=1,4),(SF(J,IT),J=1,4),
     >   EX,EY,EZ, BORO, IPASS,NOEL,KLEY,(LABEL(NOEL,I),I=1,2),LET(IT)
      ELSE
        WRITE(LN,100)
     1   KEX,      XXXO,(FO(J,IT),J=2,MXJ),
C       Initial coordinates: D,Y,T,Z,P,X,Time

     2   XXX,Y-DY,T*1.D3,Z,P*1.D3,SAR,   TAR,     BTI,   DS,
C       current coordinates              time     v/c    Step 
C                                        mu_s            size 

     4   KART,  IT,IREP(IT), SORT(IT),X, BX,BY,BZ, RET(IT), DPR(IT), PS,
C        Cart.               Path out     - kG -     (S)   dp/p_Synchro
C         or                  CHAMBR                Synchrotron
C        Polar                                       motion

     5   (SI(J,IT),J=1,4),(SF(J,IT),J=1,4),
C        spin

     6   EX,EY,EZ,   BORO,  IPASS, NOEL,  
C        - eV/cm -  refrnce       lmntt#   keywrd   2 labels at keyword
C                  rigdty(kG.cm)

     7  TX1,KLEY,TX1,TX1,LABEL(NOEL,1),TX1,TX1,LABEL(NOEL,2),TX1,
     8                                            TX1,LET(IT),TX1
C        keywrd   2 labels at keyword
C         

        INCLUDE "FRMPLT.H"
      ENDIF 
C      CALL FLUSH2(LN,.FALSE.)

      if(lst2.gt.0) write(lun,fmt='(8e12.4,I6,4e12.4)')
     2   XXX,Y-DY,T*1.D3,
     3   Z,P*1.D3,SAR,     TAR,     DY, 
     4   IT,X, BX,BY,BZ


      RETURN

      ENTRY IMPPL2(BINARI)
      BINARY=BINARI
      RETURN

      ENTRY IMPPL3(LST2I)
      LST2=LST2I
      IF(LST2.GT.0) THEN
        INQUIRE(FILE='zgoubi.impplt',OPENED=OPN)
        IF(.NOT.OPN) THEN
          IF(IDLUNI(
     >              LUN)) THEN
              OPEN(UNIT=LUN,FILE='zgoubi.impplt',ERR=99)
          ENDIF
        ENDIF
 99     WRITE(6,*) '**** SBR impplt : '
        WRITE(6,*) '         error upon open zgoubi.impplt'
      ENDIF
      RETURN
      END
