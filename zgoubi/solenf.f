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
C  Upton, NY, 11973,  USA
C  -------
      SUBROUTINE SOLENF(XX,Y,Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.AIM.H"     ! COMMON/AIM/ BO,RO,FG,GF,XI,XF,EN,EB1,EB2,EG1,EG2
      PARAMETER(MCOEF=6)
      INCLUDE "C.CHAFUI.H"     ! COMMON/CHAFUI/ XE,XS,CE(MCOEF),CS(MCOEF),QCE(MCOEF),QCS(MCOEF)
      INCLUDE "C.CHAVE.H"     ! COMMON/CHAVE/ B(5,3),V(5,3),E(5,3)
      INCLUDE "C.DDBXYZ.H"     ! COMMON/DDBXYZ/ DB(3,3),DDB(3,3,3)
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL,PI,RAD,DEG,QEL,AMPROT,CM2M
      INCLUDE "C.INTEG.H"     ! COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,HDPRF,DP,QBR,BRI

      DOUBLE PRECISION K2L,KP2L,K2LR,KL
      DOUBLE PRECISION K2R,KP2R,K2RR,KR

      DIMENSION BC(2),DBC(2,2),DDBC(2,2,2)
      PARAMETER (MDX=6)
      DIMENSION BX(MDX)

      SAVE MODL, BO2, RO2

      DATA MODL / 1 /

      GOTO(1,2) MODL
      CALL ENDJOB('SBR SOLENF. No such model # ',MODL)

C Axial field model
 1    CONTINUE

      XL = XS - XE
      X = XX - XE
      CALL SOLBAX(XL,BO2*BRI,RO2,X,
     >                              BX)
      R2  =Y*Y + Z*Z
      R   =SQRT(R2)

c         write(*,*) ' solenf XL,BO2, BRI,RO2,X : ',XL,BO2,BRI,RO2,X
c         write(*,*) ' solenf  Y,Z,R,r2,Bx ',Y,Z,R,r2,Bx

      CALL BAXBXR(BX,R,R2,
     >                    BC,DBC,DDBC)

c         write(*,*) ' solenf ',BC,DBC,DDBC,Y,Z,R
c               read(*,*)

      CALL BXRXYZ(BC,DBC,DDBC,Y,Z,R,2,
     >                                B,DB,DDB)

      RETURN

C Elliptic integral model
 2    CONTINUE

      X=XX-(XS-XE)/2.D0-XE

      XL  =-((XS-XE)/2.D0 +X)
      XR  =(XS-XE)/2.D0 -X

      Y2=Y*Y
      Z2=Z*Z
      R2  =Z2 + Y2
      R   =SQRT(R2)
      DR=2.D0*R
      ROR =RO+R
      RL  =SQRT(ROR*ROR + XL*XL)
      RR  =SQRT(ROR*ROR + XR*XR)
      R2R =RR*RR
      R2L =RL*RL
      C2  =4.D0*RO*R/ROR/ROR
      CP2 =1.D0 - C2
      K2L =4.D0*RO*R/(RL*RL)
      K2R =4.D0*RO*R/(RR*RR)
      KP2L=1.D0-K2L
      KP2R=1.D0-K2R
C            write(88,*) ' solenf bo',bo
      BBC  =BO/(4.D0*PI)*BRI

CCCCCCC Modified because the approximation (GOTO 1000) causes dirac
C pics (in an apparently stochastic maner...) in the field values in
C the uniform region within solenoid. To be debugged
C FM March/02          IF (R .LE. (RO*1.D-4)) GOTO 1000
C       IF (R .LE. (RO*1.D-8)) GOTO 1000
       IF (R .EQ. 0.D0) GOTO 1000

      CALL ELLIP(K2L,KP2L,C2,CP2,KL,EL,PL)
      CALL ELLIP(K2R,KP2R,C2,CP2,KR,ER,PR)

      KL=KL*BBC
      EL=EL*BBC
      PL=PL*BBC
      KR=KR*BBC
      ER=ER*BBC
      PR=PR*BBC

      AL =2.D0*(KL-EL)-K2L*KL
      AR =2.D0*(KR-ER)-K2R*KR
      GL  =RL/R
      GR  =RR/R

      BRL=GL*AL
      BRR=GR*AR
      BR =BRR-BRL
      BSR=BR/R


      BXR=(4.D0*RO*XR/(ROR*RR))*(KR+(RO-R)*(PR-KR)/(2.D0*RO))
      BXL=(4.D0*RO*XL/(ROR*RL))*(KL+(RO-R)*(PL-KL)/(2.D0*RO))

      BC(1)=(BXR-BXL)
      BC(2)=BR

      K2LR =SQRT(K2L)
      K2RR =SQRT(K2R)
      XL   =-XL
      RRR  =R*RR
      RRL  =R*RL
      SRO  =SQRT(RO)
      SR   =SQRT(R)
      R3R  =R2R*RR
      R3L  =R2L*RL


      DGRX =-XR/RRR
      DGLX = XL/RRL
      DGRR =-RR/R2+ROR/RRR
      DGLR =-RL/R2+ROR/RRL
      D2GRX= (1.D0-XR*XR/R2R)/RRR
      D2GLX= (1.D0-XL*XL/R2L)/RRL
      D2GRR= (2.D0*R2R/R2-2.D0*ROR/R+1.D0-ROR*ROR/R2R)/RRR
      D2GLR= (2.D0*R2L/R2-2.D0*ROR/R+1.D0-ROR*ROR/R2L)/RRL
      DGRXR= (XR/R+ROR*XR/R2R)/RRR
      DGLXR=-(XL/R+ROR*XL/R2L)/RRL

      DKRX = 2.D0*SRO*SR*XR/R3R
      DKLX =-2.D0*SRO*SR*XL/R3L
      DKRR = ( 1.D0-DR*ROR/R2R)*SRO/SR/RR
      DKLR = ( 1.D0-DR*ROR/R2L)*SRO/SR/RL
      D2KRX= (-1.D0+3.D0*XR*XR/R2R)*2.D0*SRO*SR/R3R
      D2KLX= (-1.D0+3.D0*XL*XL/R2L)*2.D0*SRO*SR/R3L
      D2KRR= SRO/R3R/SR*(-ROR*ROR/DR-XR*XR/DR-2.D0*ROR-DR
     >                   +6.D0*ROR*ROR*R/R2R)
      D2KLR= SRO/R3L/SR*(-ROR*ROR/DR-XL*XL/DR-2.D0*ROR-DR
     >                   +6.D0*ROR*ROR*R/R2L)
      DKRXR= (1.D0-6.D0*R*ROR/R2R)*SRO*XR/SR/R3R
      DKLXR=-(1.D0-6.D0*R*ROR/R2L)*SRO*XL/SR/R3L

      DAR  =K2RR*(ER/KP2R-KR)
      DAL  =K2LR*(EL/KP2L-KL)
      D2AR =((1.D0+2.D0*K2R/KP2R)*ER-KR)/KP2R
      D2AL =((1.D0+2.D0*K2L/KP2L)*EL-KL)/KP2L

      DBRRX =DGRX*AR+GR*DKRX*DAR
      DBRLX =DGLX*AL+GL*DKLX*DAL
      DBRRR =DGRR*AR+GR*DKRR*DAR
      DBRLR =DGLR*AL+GL*DKLR*DAL
      D2BRRX=D2GRX*AR+DAR*(GR*D2KRX+2*DGRX*DKRX)+
     >       GR*DKRX*DKRX*D2AR
      D2BRLX=D2GLX*AL+DAL*(GL*D2KLX+2*DGLX*DKLX)+
     >       GL*DKLX*DKLX*D2AL
      D2BRRR=D2GRR*AR+DAR*(GR*D2KRR+2*DGRR*DKRR)+
     >       GR*DKRR*DKRR*D2AR
      D2BRLR=D2GLR*AL+DAL*(GL*D2KLR+2*DGLR*DKLR)+
     >       GL*DKLR*DKLR*D2AL
      DBRRXR=DGRXR*AR+DAR*(GR*DKRXR+DGRR*DKRX+DGRX*DKRR)+
     >       GR*DKRR*DKRX*D2AR
      DBRLXR=DGLXR*AL+DAL*(GL*DKLXR+DGLR*DKLX+DGLX*DKLR)+
     >       GL*DKLR*DKLX*D2AL

      DBC(2,2)=DBRRR-DBRLR
      DBC(1,1)=-BC(2)/R-DBC(2,2)
      DBC(2,1)=DBRRX-DBRLX

      DDBC(1,1,1)=-((DBRRX-DBRLX)/R +DBRRXR-DBRLXR)
      DDBC(2,1,1)= D2BRRX-D2BRLX
      DDBC(2,2,1)= DBRRXR-DBRLXR
      DDBC(2,2,2)= D2BRRR-D2BRLR

      GOTO 2000

CCCCCCCCCCCCCCCCCCC Seems to not work well...CCCCCCCCCCC
1000  XL=-XL
      U =(RO*RO +XR*XR)
      VV =(RO*RO +XL*XL)
      SQU=SQRT(U)
      SQV=SQRT(VV)
      U3=SQU*U
      V3=SQV*VV
      U5=U3*U
      V5=V3*VV

      BC(1)=(XR/SQU+XL/SQV)
      BC(2)=0.D0

      DBC(1,1) =-1.D0/SQU+1/SQV+XR*XR/U3-XL*XL/V3
CCC DBC(2,2)=0?????
      DBC(2,2)=0.D0
      DBC(2,1)=0.D0

      DDBC(1,1,1)=3.D0*(-XR/U3-XL/V3+XR*XR*XR/U5+XL*XL*XL/V5)
      DDBC(2,1,1)= 0.D0
CCC DDBC(2,2,1)=0?????
      DDBC(2,2,1)=0.D0
      DDBC(2,2,2)= 0.D0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 2000 CONTINUE

      CALL BXRXYZ(BC,DBC,DDBC,Y,Z,R,2,
     >                                B,DB,DDB)
      RETURN

      ENTRY SOLEN2(MODLI,BO2I,RO2I)
      MODL = MODLI
      BO2 = BO2I
      RO2 = RO2I
      RETURN

      END
