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
C  Brookhaven National Laboratory                    és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE DEVTRA
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     -----------------------------------------
C     CALCULE u ET SES DERIVEES dnu/dsn, n=1,5,
C     i.e., calcul des coeffs des series de Taylor
C     (eqs. 2.2.4-5) 
C     -----------------------------------------
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CHAVE/ B(5,3),V(5,3),E(5,3)
      COMMON/D3BXYZ/ D3BX(27), D3BY(27), D3BZ(27)
      COMMON/D4BXYZ/ D4BX(81),D4BY(81),D4BZ(81)
      COMMON/DDEXYZ/ DE(9),DDE(27)
      COMMON/D3EXYZ/ D3EX(27), D3EY(27), D3EZ(27)
      COMMON/D4EXYZ/ D4EX(81), D4EY(81), D4EZ(81)
      COMMON/DDBXYZ/ DB(9),DDB(27)
      COMMON/MARK/ KART,KALC,KERK,KUASEX
      LOGICAL ZSYM
      COMMON/OPTION/ KFLD,MG,LC,ML,ZSYM
      COMMON/ORDRES/ KORD,IRD,IDS,IDB,IDE,IDZ
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      COMMON/TRAJ/ R1,T1,Z1,P1,A1,SAR,TAR,KEX,IT,AMT,QT
      COMMON/VITES/ U(6,3),DQBR(6),DDT(6)
 
      DIMENSION EU(5)
 
      PARAMETER(IMAX=3, IJMAX=9, IJLMAX=27, IM2 = IMAX+IJMAX)
      PARAMETER(IM3 = IM2+IJLMAX)

C      SAVE CL, CL9, CQ, AM2, QEV
      SAVE CL, CL9
 
C-------------------------------------------------------------
C     CALCUL D'ADRESSE : ADRESSE(I,J) = I + (J-1)*IMAX
C     EXEMPLES:
C      DB(2,1) = DB(2)
C      DB(3,3) = DB(9)
C     CALCUL D'ADRESSE : ADRESSE(I,J,K) = I+(J-1)*IMAX+(K-1)*IMAX*JMAX
C     EXEMPLES:
C      DDB(3,1,1) = DDB(3)
C      DDB(3,3,1) = DDB(3*3**0+2*3**1+0*3**2)=DDB(9)
C-------------------------------------------------------------

CALCUL u=DX/DS
      CP=COS(P1)
      U(1,1)=CP*COS(T1)
      U(1,2)=CP*SIN(T1)
      U(1,3)=SIN(P1)
 
C      umod1 = sqrt(U(1,1)*U(1,1) +U(1,2)*U(1,2) +U(1,3)*U(1,3))

      GOTO(1,2,3) KFLD
 
C----- MAGNETIC FIELD
 1    CONTINUE
 
CALCUL u'=uxB
      CALL PVECT(1,1,1)
      DO 10 K=1,3
 10     U(2,K)=V(1,K)
 
C      umod2 = sqrt(U(2,1)*U(2,1) +U(2,2)*U(2,2) +U(2,3)*U(2,3))

      IF(IDS.LE.2) GOTO 99
 
C      B'=dB/dS
      KIM = -IMAX
      DO 12 K=1,3
        TP=0.D0
        KIM = KIM + IMAX
        IF(IDB.GE.1) THEN
          DO 11 I=1,3
            IK = I + KIM
            TP=TP+U(1,I)*DB(IK)
 11       CONTINUE
        ENDIF
        B(2,K)=TP
 12   CONTINUE
 
CALCUL u''=u'xB + uxB'
      CALL PVECT(1,1,2)
      CALL PVECT(2,2,1)
      DO 13 K=1,3
 13     U(3,K)=V(1,K)+V(2,K)
 
C      umod3 = sqrt(U(3,1)*U(3,1) +U(3,2)*U(3,2) +U(3,3)*U(3,3))

      IF(IDS.LE.3) GOTO 99
 
CALCUL B''=d2B/ds2
      KIJM = -IJMAX-IMAX
      KIM = -IMAX
      DO 15 K=1,3
        TP=0.D0
        IF(IDB.GE.1) THEN
          KIJM = KIJM + IJMAX
          KIM = KIM + IMAX
          DO 14 I=1,3
            IK = I + KIM
            TP=TP+U(2,I)*DB(IK)
C
C           symetrie
C
            IF(IDB.GE.2) THEN
              IKIJ = I + KIJM
              DO 141 J=I,3
                IF(I .EQ. J) THEN
                  XMUL=1.0D0
                ELSE
                  XMUL=2.0D0
                ENDIF
                IJK = IKIJ + J*IMAX
                TP=TP+XMUL*U(1,I)*U(1,J)*DDB(IJK)
 141          CONTINUE
            ENDIF
 14       CONTINUE
        ENDIF
        B(3,K)=TP
 15   CONTINUE
 
CALCUL u'''=u''xB + 2*u'B' + uxB''
      CALL PVECT(1,1,3)
      CALL PVECT(2,2,2)
      CALL PVECT(3,3,1)
      DO 16 K=1,3
 16     U(4,K)=V(1,K)+2.D0*V(2,K)+V(3,K)
 
C      umod4 = sqrt(U(4,1)*U(4,1) +U(4,2)*U(4,2) +U(4,3)*U(4,3))

      IF(IDS.LE.4) GOTO 99
 
C      B'''=d3B/ds3
      KIJM = -IJMAX-IMAX
      KIM = -IMAX
      DO 18 K=1,3
        TP=0.D0
        IF(IDB.GE.1) THEN
          KIJM = KIJM + IJMAX
          KIM = KIM + IMAX
          DO 17 I=1,3
            IK = I + KIM
            TP=TP+U(3,I)*DB(IK)
            IF(IDB.GE.2) THEN
              IMIM2 = I - IM2
              IKIJ = I + KIJM
              DO 172 J=1,3
                IJK = IKIJ + J*IMAX
                TP=TP+3.D0*U(2,I)*U(1,J)*DDB(IJK)
C
C               symetrie
C
                IF(J.GE.I) THEN
                  IF(IDB.GE.3) THEN
                    IJL2 = IMIM2 + J*IMAX
                    DO 171 L = J,3
                      IF(J .EQ. I) THEN
                        IF(L .EQ. J) THEN
                          XMUL=1.0D0
                        ELSE
                          XMUL=3.0D0
                        ENDIF
                      ELSE
                        IF(L .EQ. J) THEN
                          XMUL=3.0D0
                        ELSE
                          XMUL=6.0D0
                        ENDIF
                      ENDIF
                      IJL = IJL2 + L*IJMAX
                      IF    ( K .EQ. 1 ) THEN
                        TP = TP+XMUL*U(1,I)*U(1,J)*U(1,L)*D3BX(IJL)
                      ELSEIF( K .EQ. 2 ) THEN
                        TP = TP+XMUL*U(1,I)*U(1,J)*U(1,L)*D3BY(IJL)
                      ELSEIF( K .EQ. 3 ) THEN
                        TP = TP+XMUL*U(1,I)*U(1,J)*U(1,L)*D3BZ(IJL)
                      ENDIF
 171                CONTINUE
                  ENDIF
                ENDIF
 172          CONTINUE
            ENDIF
 17     CONTINUE
        ENDIF
        B(4,K)=TP
 18   CONTINUE
 
CALCUL u''''=u'''xB + 3*u''xB' + 3*u'xB''+ B'''
      CALL PVECT(1,1,4)
      CALL PVECT(2,2,3)
      CALL PVECT(3,3,2)
      CALL PVECT(4,4,1)
      DO 19 K=1,3
 19     U(5,K)=V(1,K)+3.D0*V(2,K)+3.D0*V(3,K)+V(4,K)
 
C      umod5 = sqrt(U(5,1)*U(5,1) +U(5,2)*U(5,2) +U(5,3)*U(5,3))

      IF (IDS.LE.5) GOTO 99
 
C     B''''=d4B/ds4
      KIJM = -IJMAX-IMAX
      KIM = -IMAX
      DO 44 K=1,3
        TP=0.D0
        IF(IDB.GE.1) THEN
          KIJM = KIJM + IJMAX
          KIM = KIM + IMAX
          DO 43 I=1,3
            IK = I + KIM
            TP=TP+U(4,I)*DB(IK)
            IF(IDB.GE.2) THEN
              IMIM3 = I - IM3
              IMIM2 = I - IM2
              IKIJ = I + KIJM
              DO 432 J=1,3
                IJK = IKIJ + J*IMAX
                TP=TP+(4.D0*U(3,I)*U(1,J)+3.D0*U(2,I)*U(2,J))*DDB(IJK)
                IF(IDB.GE.3) THEN
                  JIM = J*IMAX
                  IJIM = IMIM3 + JIM
                  IJL  = IMIM2 + JIM
                  DO 433 L = 1,3
                    IJL = IJL + IJMAX
                    IF    (K .EQ. 1) THEN
                      TP = TP+6.D0*U(2,I)*U(1,J)*U(1,L)*D3BX(IJL)
                    ELSEIF(K .EQ. 2) THEN
                      TP = TP+6.D0*U(2,I)*U(1,J)*U(1,L)*D3BY(IJL)
                    ELSEIF(K .EQ. 3) THEN
                      TP = TP+6.D0*U(2,I)*U(1,J)*U(1,L)*D3BZ(IJL)
                    ENDIF
                    IF(J.GE.I) THEN
                      IF(L.GE.J) THEN
                        IF(IDB.GE.4) THEN
                          IJLM2 = IJIM + L*IJMAX
                          DO 431 M=L,3
                            IF(J .EQ. I) THEN
                              IF(L .EQ. J) THEN
                                IF(M .EQ. L) THEN
                                  XMUL=1.0D0
                                ELSE
                                  XMUL=4.0D0
                                ENDIF
                              ELSE
                                IF(M .EQ. L) THEN
                                  XMUL=6.0D0
                                ELSE
                                  XMUL=12.0D0
                                ENDIF
                              ENDIF
                            ELSE
                              IF(L .EQ. J) THEN
                                IF(M .EQ. L) THEN
                                  XMUL=4.0D0
                                ELSE
                                  XMUL=12.0D0
                                ENDIF
                              ELSE
                                XMUL=12.0D0
                              ENDIF
                            ENDIF
                            IJLM = IJLM2 + M*IJLMAX
                            IF    (K .EQ. 1) THEN
                              TP = TP+U(1,I)*U(1,J)*U(1,L)*U(1,M)*
     S                             D4BX(IJLM)*XMUL
                            ELSEIF(K .EQ. 2) THEN
                              TP = TP+U(1,I)*U(1,J)*U(1,L)*U(1,M)*
     S                             D4BY(IJLM)*XMUL
                            ELSEIF(K .EQ. 3) THEN
                              TP = TP+U(1,I)*U(1,J)*U(1,L)*U(1,M)*
     S                             D4BZ(IJLM)*XMUL
                            ENDIF
 431                      CONTINUE
                        ENDIF
                      ENDIF
                    ENDIF
 433              CONTINUE
                ENDIF
 432          CONTINUE
            ENDIF
 43       CONTINUE
        ENDIF
        B(5,K)=TP
 44   CONTINUE
CALCUL   u'''''=u''''xB + 4*u'''xB' + 6*u''xB''+ 4*u'xB''' + B''''
      CALL PVECT(1,1,5)
      CALL PVECT(2,2,4)
      CALL PVECT(3,3,3)
      CALL PVECT(4,4,2)
      CALL PVECT(5,5,1)
      DO 45 K=1,3
 45     U(6,K)=V(1,K)+4.D0*(V(2,K)+V(4,K))+6.D0*V(3,K)+V(5,K)
 
C      umod6 = sqrt(U(6,1)*U(6,1) +U(6,2)*U(6,2) +U(6,3)*U(6,3))

      GOTO 99
 
C----- ELECTRIC FIELD
C----- E=e/Bro
 2    CONTINUE
 
      AM2 = AMT*AMT
      BRCQ = QBR*CL9
C----- Wt=sqrt(p2+m2), p=beta.Wt/c
      CSV = SQRT(1.D0 + AM2/(BRCQ*BRCQ))
Calcul (e.u)/Bro, Bro'
      EU(1)= E(1,1)*U(1,1)+E(1,2)*U(1,2)+E(1,3)*U(1,3)
      VL=CL9/CSV
      DDT(1) = 1.D0/VL
      DBSB=EU(1)/VL
      DQBR(1)=DBSB*QBR
Calcul u'=E/v + uxB - uBro'/Bro
      DO 21 K=1,3
 21     U(2,K)=E(1,K)/VL  - U(1,K)*DBSB
 
      IF(IDS.LE.2) GOTO 99

Calcul e'/Bro = (e/Bro)'|Bro=cste = E'|Bro=cste
      KIM = -IMAX
      DO 23 K=1,3
        TP=0.D0
        IF(IDE.GE.1) THEN
          KIM = KIM + IMAX
          DO 22 I=1,3
            IK = I + KIM
            TP=TP+U(1,I)*DE(IK)
 22       CONTINUE
        ENDIF
        E(2,K)=TP
 23   CONTINUE
Calcul (e.u)'/Bro, Bro''
      EU(2)= E(2,1)*U(1,1)+E(2,2)*U(1,2)+E(2,3)*U(1,3)
     >+ E(1,1)*U(2,1)+E(1,2)*U(2,2)+E(1,3)*U(2,3)
      DCSV=QT*EU(1)/BRCQ-CSV*DBSB
      D1SV=DCSV/CL9
      DDT(2) =D1SV
      D2BSB=D1SV*EU(1)+EU(2)/VL
      DQBR(2)=D2BSB*QBR
Calcul u''=(1/v)'E + (e'/Bro)/v - 2Bro'/Brou' - Bro''/Brou
      DO 26 K=1,3
        U(3,K)=E(1,K)*D1SV + E(2,K)/VL
     >  -2.D0*U(2,K)*DBSB - U(1,K)*D2BSB
 26   CONTINUE
 
      IF(IDS.LE.3) GOTO 99
 
Calcul e''/Bro = E''|Bro=cste            
      KIJM = -IJMAX-IMAX
      KIM = -IMAX
      DO 25 K=1,3
        TP=0.D0
        IF(IDE.GE.1) THEN
          KIJM = KIJM + IJMAX
          KIM = KIM + IMAX
          DO 24 I=1,3
            IK = I + KIM
            TP=TP+U(2,I)*DE(IK)
            IF(IDE.GE.2) THEN
              IKIJ = I + KIJM
C
C               symmetry
C
              DO 241 J=I,3
                IF(I .EQ. J) THEN
                  XMUL=1.0D0
                ELSE
                  XMUL=2.0D0
                ENDIF
                IJK = IKIJ + J*IMAX
                TP=TP+XMUL*U(1,I)*U(1,J)*DDE(IJK)
 241          CONTINUE
            ENDIF
 24       CONTINUE
        ENDIF
        E(3,K)=TP
 25   CONTINUE
Calcul (e.u)''/Bro, Bro'''
      EU(3)= E(3,1)*U(1,1)+E(3,2)*U(1,2)+E(3,3)*U(1,3)
     >+ 2.D0*  ( E(2,1)*U(2,1)+E(2,2)*U(2,2)+E(2,3)*U(2,3))
     >+        E(1,1)*U(3,1)+E(1,2)*U(3,2)+E(1,3)*U(3,3)
      D2CSV= QT*EU(2)/BRCQ - 2.D0*DBSB*DCSV - D2BSB*CSV
      D21SV=D2CSV/CL9
      DDT(3) = D21SV
      D3BSB=D21SV*EU(1)+2.D0*D1SV*EU(2)+EU(3)/VL
      DQBR(3)=D3BSB*QBR
Calcul u'''=(1/v)''E + 2(1/v)'(e'/Bro) + (e''/Bro)/v
C           - 3Bro'/Brou'' - 3Bro''/Brou' - Bro'''/Brou
      DO 27 K=1,3
        U(4,K)=E(1,K)*D21SV + 2.D0*E(2,K)*D1SV + E(3,K)/VL
     >  -3.D0*( DBSB*U(3,K) + D2BSB*U(2,K) ) - D3BSB*U(1,K)
 27   CONTINUE
 
      IF(IDS.LE.4) GOTO 99
 
Calcul e'''/Bro = E'''|Bro=cste
      KIJM = -IJMAX-IMAX
      KIM = -IMAX
      DO 28 K=1,3
        TP=0.D0
        IF(IDE.GE.1) THEN
          KIJM = KIJM + IJMAX
          KIM = KIM + IMAX
          DO 29 I=1,3
            IK = I + KIM
            TP=TP+U(3,I)*DE(IK)
            IF(IDE.GE.2) THEN
              IMIM2 = I - IM2
              IKIJ = I + KIJM
              DO 292 J=1,3
                IJK = IKIJ + J*IMAX
                TP=TP+3.D0*U(2,I)*U(1,J)*DDE(IJK)
C
C               symmetry
C
                IF(J.GE.I) THEN
                  IF(IDE.GE.3) THEN
                    IJL2 = IMIM2 + J*IMAX
                    DO 291 L = J,3
                      IF(J .EQ. I) THEN
                        IF(L .EQ. J) THEN
                          XMUL=1.0D0
                        ELSE
                          XMUL=3.0D0
                        ENDIF
                      ELSE
                        IF(L .EQ. J) THEN
                          XMUL=3.0D0
                        ELSE
                          XMUL=6.0D0
                        ENDIF
                      ENDIF
                      IJL = IJL2 + L*IJMAX
                      IF    ( K .EQ. 1 ) THEN
                        TP = TP+XMUL*U(1,I)*U(1,J)*U(1,L)*D3EX(IJL)
                      ELSEIF( K .EQ. 2 ) THEN
                        TP = TP+XMUL*U(1,I)*U(1,J)*U(1,L)*D3EY(IJL)
                      ELSEIF( K .EQ. 3 ) THEN
                        TP = TP+XMUL*U(1,I)*U(1,J)*U(1,L)*D3EZ(IJL)
                      ENDIF
 291                CONTINUE
                  ENDIF
                ENDIF
 292          CONTINUE
            ENDIF
 29       CONTINUE
        ENDIF
        E(4,K)=TP
 28   CONTINUE
Calcul (e.u)'''/Bro, Bro''''
      EU(4)= E(4,1)*U(1,1)+E(4,2)*U(1,2)+E(4,3)*U(1,3)
     >+ 3.D0* (  E(3,1)*U(2,1)+E(3,2)*U(2,2)+E(3,3)*U(2,3))
     >+ 3.D0* (  E(2,1)*U(3,1)+E(2,2)*U(3,2)+E(2,3)*U(3,3))
     >+        E(1,1)*U(4,1)+E(1,2)*U(4,2)+E(1,3)*U(4,3)
      D3CSV= QT*EU(3)/BRCQ-3.D0*(DBSB*D2CSV + D2BSB*DCSV) - D3BSB*CSV
      D31SV=D3CSV/CL9
      DDT(4) =D31SV
      D4BSB=D31SV*EU(1)+3.D0*(D21SV*EU(2)+D1SV*EU(3))+EU(4)/VL
      DQBR(4)=D4BSB*QBR
Calcul u''''=(1/v)'''E+3(1/v)''(e'/Bro)+3(1/v)'(e''/Bro)+(1/v)E'''/Bro
C           - 4Bro'/Brou''' - 6Bro''/Brou'' - 4Bro'''/Brou' - Bro''''/Brou
      DO 240 K=1,3
        U(5,K)=E(1,K)*D31SV + 3.D0*E(2,K)*D21SV + 3.D0*E(3,K)*D1SV
     >  + E(4,K)/VL 
     >  -4.D0*(DBSB*U(4,K)+D3BSB*U(2,K))-6.D0*D2BSB*U(3,K)-D4BSB*U(1,K)
 240  CONTINUE
 
 
      IF(IDS.LE.5) GOTO 99
 
Calcul e''''/Bro = E''''|Bro=cste
      KIJM = -IJMAX-IMAX
      KIM = -IMAX
      DO 828 K=1,3
        TP=0.D0
        IF(IDE.GE.1) THEN
          KIJM = KIJM + IJMAX
          KIM = KIM + IMAX
          DO 829 I=1,3
            IK = I + KIM
            TP=TP+U(4,I)*DE(IK)
            IF(IDE.GE.2) THEN
              IMIM3 = I - IM3
              IMIM2 = I - IM2
              IKIJ = I + KIJM
              DO 892 J=1,3
                IJK = IKIJ + J*IMAX
                TP=TP+(4.D0*U(3,I)*U(1,J)+3.D0*U(2,I)*U(2,J))*DDE(IJK)
                IF(IDE.GE.3) THEN
                  JIM = J*IMAX
                  IJIM = IMIM3 + JIM
                  IJL  = IMIM2 + JIM
                  DO 833 L = 1,3
                    IJL = IJL + IJMAX
                    IF    (K .EQ. 1) THEN
                      TP = TP+6.D0*U(2,I)*U(1,J)*U(1,L)*D3EX(IJL)
                    ELSEIF(K .EQ. 2) THEN
                      TP = TP+6.D0*U(2,I)*U(1,J)*U(1,L)*D3EY(IJL)
                    ELSEIF(K .EQ. 3) THEN
                      TP = TP+6.D0*U(2,I)*U(1,J)*U(1,L)*D3EZ(IJL)
                    ENDIF
                    IF(J.GE.I) THEN
                      IF(L.GE.J) THEN
                        IF(IDE.GE.4) THEN
                          IJLM2 = IJIM + L*IJMAX
                          DO 831 M=L,3
                            IF(J .EQ. I) THEN
                              IF(L .EQ. J) THEN
                                IF(M .EQ. L) THEN
                                  XMUL=1.0D0
                                ELSE
                                  XMUL=4.0D0
                                ENDIF
                              ELSE
                                IF(M .EQ. L) THEN
                                  XMUL=6.0D0
                                ELSE
                                  XMUL=12.0D0
                                ENDIF
                              ENDIF
                            ELSE
                              IF(L .EQ. J) THEN
                                IF(M .EQ. L) THEN
                                  XMUL=4.0D0
                                ELSE
                                  XMUL=12.0D0
                                ENDIF
                              ELSE
                                XMUL=12.0D0
                              ENDIF
                            ENDIF
                            IJLM = IJLM2 + M*IJLMAX
                            IF    (K .EQ. 1) THEN
                              TP = TP+U(1,I)*U(1,J)*U(1,L)*U(1,M)*
     S                             D4EX(IJLM)*XMUL
                            ELSEIF(K .EQ. 2) THEN
                              TP = TP+U(1,I)*U(1,J)*U(1,L)*U(1,M)*
     S                             D4EY(IJLM)*XMUL
                            ELSEIF(K .EQ. 3) THEN
                              TP = TP+U(1,I)*U(1,J)*U(1,L)*U(1,M)*
     S                             D4EZ(IJLM)*XMUL
                            ENDIF
 831                      CONTINUE
                        ENDIF
                      ENDIF
                    ENDIF
 833              CONTINUE
                ENDIF
 892          CONTINUE
            ENDIF
 829      CONTINUE
        ENDIF
        E(5,K)=TP
 828  CONTINUE
Calcul (e.u)''''/Bro, Bro'''''
      EU(5)=     E(5,1)*U(1,1)+E(5,2)*U(1,2)+E(5,3)*U(1,3)
     >+ 4.D0* (  E(4,1)*U(2,1)+E(4,2)*U(2,2)+E(4,3)*U(2,3))
     >+ 6.D0* (  E(3,1)*U(3,1)+E(3,2)*U(3,2)+E(3,3)*U(3,3))
     >+ 4.D0* (  E(2,1)*U(4,1)+E(2,2)*U(4,2)+E(2,3)*U(4,3))
     >+          E(1,1)*U(5,1)+E(1,2)*U(5,2)+E(1,3)*U(5,3)
      D4CSV= QT*EU(4)/BRCQ-(4.D0*DBSB*D3CSV + 6.D0*D2BSB*D2CSV +4.D0*
     >   D3BSB*DCSV + D4BSB *CSV)
      D41SV=D4CSV/CL9
      DDT(5) = D41SV
      D5BSB=D41SV*EU(1)+4.D0*D31SV*EU(2)+6.D0*D21SV*EU(3)+
     >       4.D0*D1SV*EU(4) + EU(5)/VL
      DQBR(5)=D5BSB*QBR
Calcul u'''''
      DO 840 K=1,3
        U(6,K)=E(1,K)*D41SV + 4.D0*E(2,K)*D31SV + 6.D0*E(3,K)*D21SV +
     >  4.D0*E(4,K)*D1SV + E(5,K)/VL 
     >  -(5.D0*DBSB*U(5,K)+10.D0*D2BSB*U(4,K)+10.D0*D3BSB*U(3,K)+
     >  5.D0*D4BSB*U(2,K)+D5BSB*U(1,K))
 840  CONTINUE

      GOTO 99

C----- MAGNETIC+ELECTRIC FIELD
C----- B=b/Bro, E=e/Bro
 3    CONTINUE

      AM2 = AMT*AMT
      BRCQ = QBR*CL9
C----- Wt=sqrt(p2+m2), p=beta.Wt/c
      CSV = SQRT(1.D0 + AM2/(BRCQ*BRCQ))
 
Calcul (e.u)/Bro, Bro'
      EU(1)= E(1,1)*U(1,1)+E(1,2)*U(1,2)+E(1,3)*U(1,3)
      VL=CL9/CSV
      DDT(1) = 1.D0/VL
      DBSB=EU(1)/VL
      DQBR(1)=DBSB*QBR
CALCUL u'=uxB
      CALL PVECT(1,1,1)
Calcul u'=E/v + uxB - uBro'/Bro
      DO 31 K=1,3
 31     U(2,K)=E(1,K)/VL + V(1,K) - U(1,K)*DBSB
 
      IF(IDS.LE.2) GOTO 99
 
Calcul e'/Bro = (e/Bro)'|Bro=cste = E'|Bro=cste
      KIM = -IMAX
      DO 33 K=1,3
        TP=0.D0
        IF(IDE.GE.1) THEN
          KIM = KIM + IMAX
          DO 32 I=1,3
            IK = I + KIM
            TP=TP+U(1,I)*DE(IK)
 32       CONTINUE
        ENDIF
        E(2,K)=TP
 33   CONTINUE
Calcul (e.u)'/Bro, Bro''
      EU(2)= E(2,1)*U(1,1)+E(2,2)*U(1,2)+E(2,3)*U(1,3)
     >+ E(1,1)*U(2,1)+E(1,2)*U(2,2)+E(1,3)*U(2,3)
      DCSV=QT*EU(1)/BRCQ-CSV*DBSB
      D1SV=DCSV/CL9
      DDT(2) = D1SV
      D2BSB=D1SV*EU(1)+EU(2)/VL
      DQBR(2)=D2BSB*QBR
C      b'/Bro = B'|Bro=cste
      KIM = -IMAX
      DO 312 K=1,3
        TP=0.D0
        IF(IDB.GE.1) THEN
          KIM = KIM + IMAX
          DO 311 I=1,3
            IK = I + KIM
            TP=TP+U(1,I)*DB(IK)
 311      CONTINUE
        ENDIF
        B(2,K)=TP
 312  CONTINUE
CALCUL (uxb)'/Bro
      CALL PVECT(1,1,2)
      CALL PVECT(2,2,1)
Calcul u''=(1/v)'E + (e'/Bro)/v + (uxb)'/Bro - 2Bro'/Brou' - Bro''/Brou
      DO 36 K=1,3
        U(3,K)=E(1,K)*D1SV + E(2,K)/VL + V(1,K)+V(2,K)
     >  -2.D0*U(2,K)*DBSB - U(1,K)*D2BSB
 36   CONTINUE
 
      IF(IDS.LE.3) GOTO 99
 
Calcul e''/Bro = E''|Bro=cste
      KIJM = -IJMAX-IMAX
      KIM = -IMAX
      DO 35 K=1,3
        TP=0.D0
        IF(IDE.GE.1) THEN
          KIJM = KIJM + IJMAX
          KIM = KIM + IMAX
          DO 34 I=1,3
            IK = I + KIM
            TP=TP+U(2,I)*DE(IK)
            IF(IDE.GE.2) THEN
              IKIJ = I + KIJM
C
C             symetrie
C
              DO 341 J=I,3
                IF(I .EQ. J) THEN
                  XMUL=1.0D0
                ELSE
                  XMUL=2.0D0
                ENDIF
                IJK = IKIJ + J*IMAX
                TP=TP+XMUL*U(1,I)*U(1,J)*DDE(IJK)
 341          CONTINUE
            ENDIF
 34       CONTINUE
        ENDIF
        E(3,K)=TP
 35   CONTINUE
Calcul (e.u)''/Bro, Bro'''
      EU(3)= E(3,1)*U(1,1)+E(3,2)*U(1,2)+E(3,3)*U(1,3)
     >+ 2.D0*  ( E(2,1)*U(2,1)+E(2,2)*U(2,2)+E(2,3)*U(2,3))
     >+        E(1,1)*U(3,1)+E(1,2)*U(3,2)+E(1,3)*U(3,3)
      D2CSV= QT*EU(2)/BRCQ - 2.D0*DBSB*DCSV - D2BSB*CSV
      D21SV=D2CSV/CL9
      DDT(3) = D21SV
      D3BSB=D21SV*EU(1)+2.D0*D1SV*EU(2)+EU(3)/VL
      DQBR(3)=D3BSB*QBR
CALCUL b''/Bro = B''|Bro=cste
      KIJM = -IJMAX-IMAX
      KIM = -IMAX
      DO 315 K=1,3
        TP=0.D0
        IF(IDB.GE.1) THEN
          KIJM = KIJM + IJMAX
          KIM = KIM + IMAX
          DO 314 I=1,3
            IK = I + KIM
            TP=TP+U(2,I)*DB(IK)
C
C           symetrie
C
            IF(IDB.GE.2) THEN
              IKIJ = I + KIJM
              DO 3141 J=I,3
                IF(I .EQ. J) THEN
                  XMUL=1.0D0
                ELSE
                  XMUL=2.0D0
                ENDIF
                IJK = IKIJ + J*IMAX
                TP=TP+XMUL*U(1,I)*U(1,J)*DDB(IJK)
 3141         CONTINUE
            ENDIF
 314      CONTINUE
        ENDIF
        B(3,K)=TP
 315  CONTINUE
Calcul (uxb)''/Bro
      CALL PVECT(1,1,3)
      CALL PVECT(2,2,2)
      CALL PVECT(3,3,1)
Calcul u'''=(1/v)''E + 2(1/v)'(e'/Bro) + (e''/Bro)/v + (uxb)''/Bro
C           - 3Bro'/Brou'' - 3Bro''/Brou' - Bro'''/Brou
      DO 37 K=1,3
        U(4,K)=E(1,K)*D21SV + 2.D0*E(2,K)*D1SV + E(3,K)/VL
     >  +V(1,K)+2.D0*V(2,K)+V(3,K)
     >  -3.D0*(DBSB*U(3,K) + D2BSB*U(2,K)) - D3BSB*U(1,K)
 37   CONTINUE
 
      IF(IDS.LE.4) GOTO 99
 
Calcul e'''/Bro = E'''|Bro=cste
      KIJM = -IJMAX-IMAX
      KIM = -IMAX
      DO 38 K=1,3
        TP=0.D0
        IF(IDE.GE.1) THEN
          KIJM = KIJM + IJMAX
          KIM = KIM + IMAX
          DO 39 I=1,3
            IK = I + KIM
            TP=TP+U(3,I)*DE(IK)
            IF(IDE.GE.2) THEN
              IMIM2 = I - IM2
              IKIJ = I + KIJM
              DO 392 J=1,3
                IJK = IKIJ + J*IMAX
                TP=TP+3.D0*U(2,I)*U(1,J)*DDE(IJK)
C
C               symetrie
C
                IF(J.GE.I) THEN
                  IF(IDE.GE.3) THEN
                    IJL2 = IMIM2 + J*IMAX
                    DO 391 L = J,3
                      IF(J .EQ. I) THEN
                        IF(L .EQ. J) THEN
                          XMUL=1.0D0
                        ELSE
                          XMUL=3.0D0
                        ENDIF
                      ELSE
                        IF(L .EQ. J) THEN
                          XMUL=3.0D0
                        ELSE
                          XMUL=6.0D0
                        ENDIF
                      ENDIF
                      IJL = IJL2 + L*IJMAX
                      IF    ( K .EQ. 1 ) THEN
                        TP = TP+XMUL*U(1,I)*U(1,J)*U(1,L)*D3EX(IJL)
                      ELSEIF( K .EQ. 2 ) THEN
                        TP = TP+XMUL*U(1,I)*U(1,J)*U(1,L)*D3EY(IJL)
                      ELSEIF( K .EQ. 3 ) THEN
                        TP = TP+XMUL*U(1,I)*U(1,J)*U(1,L)*D3EZ(IJL)
                      ENDIF
 391                CONTINUE
                  ENDIF
                ENDIF
 392          CONTINUE
            ENDIF
 39       CONTINUE
        ENDIF
        E(4,K)=TP
 38   CONTINUE
Calcul (e.u)'''/Bro, Bro''''
      EU(4)= E(4,1)*U(1,1)+E(4,2)*U(1,2)+E(4,3)*U(1,3)
     >+ 3.D0* (  E(3,1)*U(2,1)+E(3,2)*U(2,2)+E(3,3)*U(2,3))
     >+ 3.D0* (  E(2,1)*U(3,1)+E(2,2)*U(3,2)+E(2,3)*U(3,3))
     >+        E(1,1)*U(4,1)+E(1,2)*U(4,2)+E(1,3)*U(4,3)
      D3CSV= QT*EU(3)/BRCQ-3.D0*(DBSB*D2CSV + D2BSB*DCSV) - D3BSB*CSV
      D31SV=D3CSV/CL9
      DDT(4) = D31SV
      D4BSB=D31SV*EU(1)+3.D0*(D21SV*EU(2)+D1SV*EU(3))+EU(4)/VL
      DQBR(4)=D4BSB*QBR
CALCUL b'''/Bro = B'''|Bro=cste
      KIJM = -IJMAX-IMAX
      KIM = -IMAX
      DO 318 K=1,3
        TP=0.D0
        IF(IDB.GE.1) THEN
          KIJM = KIJM + IJMAX
          KIM = KIM + IMAX
          DO 317 I=1,3
            IK = I + KIM
            TP=TP+U(3,I)*DB(IK)
            IF(IDB.GE.2) THEN
              IMIM2 = I - IM2
              IKIJ = I + KIJM
              DO 3172 J=1,3
                IJK = IKIJ + J*IMAX
                TP=TP+3.D0*U(2,I)*U(1,J)*DDB(IJK)
C
C               symetrie
C
                IF(J.GE.I) THEN
                  IF(IDB.GE.3) THEN
                    IJL2 = IMIM2 + J*IMAX
                    DO 3171 L = J,3
                      IF(J .EQ. I) THEN
                        IF(L .EQ. J) THEN
                          XMUL=1.0D0
                        ELSE
                          XMUL=3.0D0
                        ENDIF
                      ELSE
                        IF(L .EQ. J) THEN
                          XMUL=3.0D0
                        ELSE
                          XMUL=6.0D0
                        ENDIF
                      ENDIF
                      IJL = IJL2 + L*IJMAX
                      IF    ( K .EQ. 1 ) THEN
                        TP = TP+XMUL*U(1,I)*U(1,J)*U(1,L)*D3BX(IJL)
                      ELSEIF( K .EQ. 2 ) THEN
                        TP = TP+XMUL*U(1,I)*U(1,J)*U(1,L)*D3BY(IJL)
                      ELSEIF( K .EQ. 3 ) THEN
                        TP = TP+XMUL*U(1,I)*U(1,J)*U(1,L)*D3BZ(IJL)
                      ENDIF
 3171               CONTINUE
                  ENDIF
                ENDIF
 3172         CONTINUE
            ENDIF
 317      CONTINUE
        ENDIF
        B(4,K)=TP
 318  CONTINUE
Calcul (uxb)'''/Bro
      CALL PVECT(1,1,4)
      CALL PVECT(2,2,3)
      CALL PVECT(3,3,2)
      CALL PVECT(4,4,1)
Calcul u''''=(1/v)'''E+3(1/v)''(e'/Bro)+3(1/v)'(e''/Bro)+(1/v)E'''/Bro
C           + (uxb)'''/Bro
C           - 4Bro'/Brou''' - 6Bro''/Brou'' - 4Bro'''/Brou' - Bro''''/Brou
      DO 40 K=1,3
        U(5,K)=E(1,K)*D31SV + 3.D0*E(2,K)*D21SV + 3.D0*E(3,K)*D1SV
     >  + E(4,K)/VL + V(1,K)+3.D0*V(2,K)+3.D0*V(3,K)+V(4,K)
     >  -4.D0*(DBSB*U(4,K)+D3BSB*U(2,K))-6.D0*D2BSB*U(3,K)-D4BSB*U(1,K)
 40   CONTINUE
 
 99   CONTINUE
      IF(LST .GE. 1) CALL IMPDEV
      RETURN

      ENTRY DEVTRW(CEL)
      CL = CEL
      CL9=CL*1.D-9
      RETURN
      END
