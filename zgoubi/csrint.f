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
      SUBROUTINE CSRINT(DS,IMAX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXSTEP.H'
      INCLUDE 'MAXTRA.H'
      INCLUDE 'CSR.H'
C      COMMON/CSR/ KTRA,KCSR,YZXB(MXSTEP,41,36),DWC(MXT)
      INCLUDE "MAXCOO.H"
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      COMMON/RIGID/ BORO,DPREF,DP,BR
      COMMON/TRAJ/ Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT
      COMMON/UNITS/ UNIT(MXJ)

      SAVE FE, NPART, BUNL, XMI

C----- Receiver position & speed, E-field at receiver
      DIMENSION BET(3), EVB(3), ON(3), BE(3)

C IT=1 is the front particle in the bunch. Cannot be victim of its own CSR
      IF(IT.EQ.IMAX) RETURN

      STP = DS * UNIT(5)

C----- Velocity of observer
      P0 = BR*CL9*QT
      BTA = P0 / SQRT( P0*P0 + AMT*AMT )  
      CP = COS(P)
      BET(1) = BTA * COS(T) * CP
      BET(2) = BTA * SIN(T) * CP
      BET(3) = BTA * SIN(P)

C------- Loop on all  positions (tagged KSTP) upstream of IT, including current step
        CALL INTEG4(NSTP)

C------- Neglect (avoid!) action of bunch on itself
        IF(NSTP.EQ.1) RETURN
        KSTP = NSTP-1

 1      CONTINUE

C------- Loop on all emittors within bunch at step KSTP, except current particle itself
          ITE = IT+1
 3        CONTINUE
           
C--------- Compute propagation time on cord from ITE to IT. 
C           Direction ON and distance RR from ITE at step KSTP (tagged KSTP) to current particle position, 
C              BE is velocity of ITE
            CALL SREF1(KSTP,ITE,
     >                          RR,ON,BE,BTE)
          
          TCORD = RR / CL 
C--------- Time of departure of signal from ITE 
C                   YZXB(KSTP,ITE,8)=time of flight of ITE at KSTP
C                   FO(MXJ,MXT)/(BE*CL) accounts for position of ITE in bunch (front 
C                                        particles i.e. with lower ITE have FO>0)
          TITE =  YZXB(KSTP,ITE,8) + FO(6,ITE) * UNIT(5)  / (BTE*CL)
C--------- Time of arrival of signal at IT 
          TIT = TITE + TCORD
C--------- Time at observer. FO(6,IT)/(BTA*CL) accounts for position of IT in bunch 
          TARO = TAR*1.D-11 + FO(6,IT) * UNIT(5) /(BTA*CL)
            write(89,fmt='(4I4,1P,4G14.6,A)') 
     >       ITE,kstp,it,nstp, TITE,TCORD,TIT,TARO, 
     >       '   ITE,kstp,IT,nstp,TITE,TCORD,TIT,TARO'
            DTARO = STP / (BTA*CL)
            write(89,fmt='(4I4,1P,6G12.4,A)') 
     >       ITE,kstp, it,nstp, YZXB(KSTP,ITE,6), 
     >                  YZXB(KSTP,ITE,8),tcord,TIT,sar/1.D2,TARO,
     >       ' ITE,kstp,it,nstp,sare,TARE,tcord, TIT, saro, TARO'
          IF(TIT .GE. TARO) THEN     !!!!!!!!!!!!!!!!!! .AND. TIT.LT.TARO+DTARO) THEN
            CALL SREF(KSTP,FE,ITE,
     >                        EVB,*97)

            DDWC = QT * VSCAL(EVB,BET,3) * DTARO
            IF(KTRA .NE. 0) THEN
              X0 = DBLE(IMAX)/2.D0
              SIG = DBLE(IMAX)/DBLE(KTRA)
              DDWC =  DDWC  *  EXP(-(ITE-X0)**2/(2.D0*SIG**2))   
            ENDIF

            DWC(IT)=DWC(IT) + DDWC

               write(87,fmt='(4I4,1P,5G12.4,I4)') 
     >                  it,nstp, ITE,kstp,ddwc,tITE,tcord,tit,taro,ktra

              KSTP = KSTP-1
              ITE=ITE+1
              IF(ITE.GT.IMAX) GOTO 10 
              IF(KSTP.EQ.0) GOTO 10
                GOTO 3
          ELSE
            KSTP = KSTP-1
              IF(KSTP.EQ.0) GOTO 10 
              GOTO 3

          ENDIF
C             ITE=ITE-1
C               IF(ITE.EQ.0) THEN
C                 KSTP = KSTP-1
C                 IF(KSTP.EQ.0) GOTO 10 
C                 GOTO 1
C               ENDIF
C             GOTO 3

 10   CONTINUE
      RETURN

      ENTRY CSRIN1(QM,NPI)
      FE = QM * CL * 1.D-13     ! = Q/4.Pi.Epsilon0.c;  gives  E in MeV/m
      NPART=NPI
      RETURN

      ENTRY CSRIN2(XMII,BUNLI)
      XMI=XMII
      BUNL = BUNLI
      RETURN

 97   CONTINUE
      WRITE(6,*) ' Problem in CSRINT :  in sref '
      stop
      RETURN
 98   CONTINUE
      WRITE(6,*) 
     >  'Radiator is still ahead of receiver -> goto next receiver step'
      RETURN
      END
