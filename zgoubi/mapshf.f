C  ZGOUBI, a program for computing the trajectories of charged particles
C  in electric and magnetic fields
C  Copyright (C) 1988-2007  François Mot
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
C  François Mot <fmeot@bnl.gov>
C  BNL
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE MAPSHF(HC,XH,YH,DY,IXMA,JYMA,
     >                                        HCB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C------------------------------------------------------
C     Compute new field map, HCB, in D-translated frame
C------------------------------------------------------
      INCLUDE 'PARIZ.H'
      DIMENSION HC(ID,MXX,MXY,IZ), HCB(ID,MXX,MXY,IZ)
      DIMENSION XH(MXX),YH(MXY)

      PARAMETER (IRD=2, KART=1)
      PARAMETER (DX = 0.D0)

      DIMENSION BMESH(5,5)

          call raz(hcb,ID*MXX*MXY*IZ)


      DA=XH(2)-XH(1)
      DO IA1=1,IXMA
        A1 = XH(IA1)+DX
        IAC = IA1
        JAC=(A1-XH(1))/DA+1.5D0

        IF(IRD .EQ. 2) THEN
C---------  2-D/3*3 points or 3-D/3*3*3 points grid
          IAC=MAX0(IAC,2)
          IAC=MIN0(IAC,IXMA-1)
        ELSE
C---------  2-D  5*5 points
          IAC=MAX0(IAC,3)
          IAC=MIN0(IAC,IXMA-2)
        ENDIF
        A=A1-XH(IAC)
                   
        DR=YH(2)-YH(1)           

C           write(88,*) dy,yh(1),yh(jyma),dr
        DO IR1 = 1,JYMA
          R1 = YH(IR1)+DY
          IRC=(R1-YH(1))/DR+1.5D0
          JRC=IR1

          IF    (IRC.LT.1 .OR. IRC.GT.JYMA) THEN
C            ... POINT HORS CARTE DE Champ
C             write(*,*) ' POINT HORS CARTE DE Champ'
             KERK=1
          ELSE
C            ... Particle is inside field map
             KERK=0
          ENDIF

          IF(IRD .EQ. 2) THEN
C          +++ 2-D  3*3 points or 3-D  3*3*3 points grid
            IRC=MAX0(IRC,2)
            IRC=MIN0(IRC,JYMA-1)
          ELSE
C          +++  2-D grid,  5*5 points
            IRC=MAX0(IRC,3)
            IRC=MIN0(IRC,JYMA-2)
          ENDIF

          R=R1-YH(IRC)
            
C         .... 2-D mid-plane field maps
          IF (IRD.EQ.4) GOTO 21
          IF (IRD.GT.4) GOTO 22
 
          CONTINUE
C          WRITE(*,*) ' .... ORDRE 2 ,  3*3 points grid'
 
          F1= HC(ID,IAC-1,IRC-1,1)
          F2= HC(ID,IAC  ,IRC-1,1)
          F3= HC(ID,IAC+1,IRC-1,1)
          F4= HC(ID,IAC-1,IRC  ,1)
          F5= HC(ID,IAC  ,IRC  ,1)
          F6= HC(ID,IAC+1,IRC  ,1)
          F7= HC(ID,IAC-1,IRC+1,1)
          F8= HC(ID,IAC  ,IRC+1,1)
          F9= HC(ID,IAC+1,IRC+1,1)

          C1=F1+F2+F3+F4+F5+F6+F7+F8+F9
          C4=F1+F3+F4+F6+F7+F9
          C5=F1+F2+F3+F7+F8+F9

          BZ=(20.D0*C1-12.D0*C4-12.D0*C5)/36.D0   

C         write(88,*) ' mapshf  r,r1,irc,jrc,bz : ',irc, jrc

          GOTO 30
 
 21       CONTINUE
C           .... ORDRE 4 , GRILLE A 5*5 POINTS
 
C         *** CALCUL DES 15 COEFFS DU POLYNOME DE DEGRE 4 :
C          BZ=A0+A10*X+A11*Y+A20*X2+A21*XY+A22*Y2+...+A42*X2Y2+A43*XY3+A44*Y4
          S1   =0.D0
          SI   =0.D0
          SJ   =0.D0
          SI2  =0.D0
          SIJ  =0.D0
          SJ2  =0.D0
          SI3  =0.D0
          SI2J =0.D0
          SIJ2 =0.D0
          SJ3  =0.D0
          SI4  =0.D0
          SI3J =0.D0
          SI2J2=0.D0
          SIJ3 =0.D0
          SJ4  =0.D0
          DO 1 J=1,5
             JR=3-J
             RJ=DBLE(JR)
             RJ2=RJ*RJ
             RJ3=RJ2*RJ
             IRCJR=IRC+JR
             DO 2 I=1,5
                IA=I-3
                BIAJR= HC(ID,IAC+IA,IRCJR,1)
                BMESH(I,J) = BIAJR
                AI=DBLE(IA)
                AI2=AI*AI*BIAJR
                AI3=AI2*AI
                AI4=AI3*AI
                AI =AI*BIAJR
                S1   =S1   +         BIAJR
                SI   =SI   +AI
                SJ   =SJ   +RJ      *BIAJR
                SI2  =SI2  +AI2
                SIJ  =SIJ  +AI*RJ
                SJ2  =SJ2  +RJ2     *BIAJR
                SI3  =SI3  +AI3
                SI2J =SI2J +AI2*RJ
                SIJ2 =SIJ2 +AI *RJ2
                SJ3  =SJ3  +RJ3     *BIAJR
                SI4  =SI4  +AI4
                SI3J =SI3J +AI3*RJ
                SI2J2=SI2J2+AI2*RJ2
                SIJ3 =SIJ3 +AI *RJ3
                SJ4  =SJ4  +RJ3*RJ  *BIAJR
2            CONTINUE
1         CONTINUE
 
C         *** CALCUL DE BZ ET SES DERIVEES AU POINT DE MAILLAGE (IAC,IRC)
          DA1  =1.D0/DA
          DA12 =DA1 *DA1
          DA13 =DA12*DA1
          DR1  =1.D0/DR
          DR12 =DR1 *DR1
          DR13 =DR12*DR1
          BZXXX = (-3.4D0 *SI + SI3  )*F1S12
          BZXYY = (-2.D0  *SI + SIJ2 )*F1S70
          BZX   = (.02D0 *SI - 3.4D0*BZXXX*F1S6 - BZXYY) *DA1        
          BZXXX = BZXXX                                  *DA13       
          BZXYY = BZXYY                                  *DA1 *DR12  
          BZXXY = (-2.D0  *SJ + SI2J )*F1S70
          BZYYY = (-3.4D0 *SJ + SJ3  )*F1S12
          BZY   = (.02D0 *SJ - BZXXY - 3.4D0*BZYYY*F1S6)      *DR1   
          BZXXY = BZXXY                                  *DA12*DR1   
          BZYYY = BZYYY                                       *DR13  
          BZX3Y = (-3.4D0 *SIJ+ SI3J )*F1S24
          BZXY3 = (-3.4D0 *SIJ+ SIJ3 )*F1S24
          BZXY  = (.01D0*SIJ- 3.4D0*(BZX3Y + BZXY3)*F1S6)*DA1 *DR1   
          BZX3Y = BZX3Y                                  *DA13*DR1   
          BZXY3 = BZXY3                                  *DA1 *DR13  
          BZX2Y2= (    SI2J2-2.D0*(SI2+SJ2)+ 4.D0 *S1)*F1S49
          BZY4  = (7.D0 *SJ4  -31.D0 *SJ2    +14.4D0 *S1)*F1S12
          BZYY  = ((SJ2-2.D0*S1-310.D0*F1S24*BZY4)*F1S35-BZX2Y2)
          BZX4  = (7.D0 *SI4  -31.D0 *SI2    +14.4D0 *S1)*F1S12
          BZXX  = (SI2-2.D0*S1-155.D0*BZX4*F1S12)*F1S35-BZX2Y2
          BZ =(.01D0*SI2J2 - 1.7D0*(BZXX+BZYY)-13.D0*(BZX4+BZY4)*F1S24 - 
     >    2.89D0*BZX2Y2) 
          BZXX  = BZXX   *DA12       
          BZYY  = BZYY   *DR12       
          BZX4  = BZX4   *DA13*DA1   
          BZX2Y2= BZX2Y2 *DA12*DR12  
          BZY4  = BZY4   *DR13*DR1   
 
C  CALCUL DE BZ ET SES DERIVEES AU POINT COURANT A1,R1
          A2=A *A
          A3=A2*A
          R2=R *R
          R3=R2*R
          BZ=BZ+BZX*A+BZY*R+0.5D0*BZXX*A2+BZXY*A*R+0.5D0*BZYY*R2
     >    +F1S6*BZXXX*A3+.5D0*BZXXY*A2*R 
     >     +.5D0*BZXYY*A *R2+F1S6*BZYYY*R3
     >    +(0.25D0*BZX4*A+BZX3Y*R)*F1S6*A3+.25D0*BZX2Y2*A2*R2
     >    +(BZXY3*A+0.25D0*BZY4*R)*F1S6*R3
 
          GOTO 99
 
 22       CONTINUE
C--------- ORDRE 2 , GRILLE A 5*5 POINTS
 
C     --- CALCUL DES 6 COEFFS DU POLYNOME DE DEGRE 2 :
C       BZ=A0+A10*X+A11*Y+A20*X2+A21*XY+A22*YY
          A0 =0D0
          A10=0D0
          A11=0D0
          A20=0D0
          A21=0D0
          A22=0D0
          DO 3 J=1,5
             JR=3-J
             RJ=DBLE(JR)
             RJ2=RJ*RJ
             IRCJR=IRC+JR
             DO 4 I=1,5
                IA=I-3
                AI=DBLE(IA)
                AI2=AI*AI
                BIAJR= HC(ID,IAC+IA,IRCJR,1)
                A0 =A0 +(27.D0-5.D0*(AI2+RJ2))*BIAJR
                A10=A10+         AI       *BIAJR
                A11=A11+             RJ   *BIAJR
                A20=A20+       ( AI2-2.D0 ) *BIAJR
                A21=A21+         AI *RJ   *BIAJR
                A22=A22+       (-2.D0 +RJ2) *BIAJR
 4           CONTINUE
 3        CONTINUE

          FAC=1.D0/70D0
          A0 =A0 *FAC*.4D0

C  CALCUL DE BZ
          BZ=A0                 
 
 30       CONTINUE
C  +++ CALCUL DE BZ ET SES DERIVEES AU POINT COURANT A1,R1
              bzav = bz
c              write(88,*) ' bz avant ',bz
          BZ=BZ+BZX*A+BZY*R+0.5D0*BZXX*A*A+0.5D0*BZYY*R*R+BZXY*A*R
              write(88,*) ' irc,jrc,bzav-bz   ',irc,jrc,a,r,bz-bzav
 
 99       CONTINUE

          HCB(ID,jAC,jRC,1) = BZ

        ENDDO ! R1 
      ENDDO ! A1 

               do jjj=1,JYMA
                 write(89,*) ' y ',yh(jjj)
                 do iii=1,IXMA
                   write(89,*) ' xh  ',xh(iii)
                     write(89,*) ' bz av/ap         '
     >          ,irc,jrc,HCB(id,iii,jjj,1)-HC(id,iii,jjj,1)
                 enddo
               enddo
c              stop
      RETURN
      END

