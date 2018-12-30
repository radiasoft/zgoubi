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
      SUBROUTINE SPNTRK(DS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.CHAVE_2.H"     ! COMMON/CHAVE/ B(5,3),V(5,3),E(5,3)
      INCLUDE "C.CONST_2.H"     ! COMMON/CONST/ CL9,CL,PI,RAD,DEG,QEL,AMPROT,CM2M
C      LOGICAL ZSYM
      INCLUDE "C.TYPFLD.H"     ! COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      INCLUDE "C.PTICUL.H"     ! COMMON/PTICUL/ AMASS,Q,G,TO
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,HDPRF,DP,QBR,BRI
      INCLUDE "MAXTRA.H"
      INCLUDE "C.SPIN.H"     ! COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)
      INCLUDE "C.SPTRK.H"     ! COMMON/SPTRK/ EU(5),CSV,DCSV,D2CSV,D3CSV,D4CSV,DBSB,D2BSB,D3BSB,D4BSB,D5BSB
C     >,DBSB,D2BSB,D3BSB,D4BSB,D5BSB
      INCLUDE "C.TRAJ.H"     ! COMMON/TRAJ/ Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT
      INCLUDE "C.VITES.H"     ! COMMON/VITES/ U(6,3),DQBR(6),DDT(6)

Correction, FM, 02/98
CCC   DIMENSION BU(4,3), BP(4,3), O(4,3), S(5,3)
C      DIMENSION BU(4,4), BP(4,3), O(4,3), S(5,3)
      DIMENSION BU(5,5), BP(5,3), O(5,3), S(6,3)
      DIMENSION GA1(5), GAC(5), OE(5,3), DSE(6,3), VEU(5,3)

      DIMENSION SPMI(4,MXT), SPMA(4,MXT)
      DIMENSION SMI(4,MXT), SMA(4,MXT)

      SAVE SMI, SMA
      PARAMETER (ICMXT=4*MXT)
      DATA SMI, SMA / ICMXT*1D10, ICMXT* -1D10 /

C////// debug...
c      call ZGNOEL(
c     >             NOEL)
c       if(noel.eq.20) then
c           write(*,*) 'spntrk IN : ',noel
c           write(*,*) ' QBR, DP, DS : ',QBR,DP,DS
c           write(*,*) ' S_X,Y,Z, |S| : ',it,(sf(i,it),i=1,4)
c         endif
C///////////////////////

      PM = QBR*CL9/AMT
      GG = G * SQRT( 1.D0+PM*PM )
      GP = 1.D0+ GG
      GM = G-GG

C////// debug...
c       if(noel.eq.20) then
c           write(*,*) ' PM, GG, GP, GM : ',PM, GG, GP, GM ,kfld
c         endif
C////// debug...

      GOTO(1,2,1) KFLD

C----- Mag pure (KFLD=1) or Mag+Elc (KFLD=3)
 1    CONTINUE

CALCUL B.U, B.U', B.U'',B.U''', B'.U, B'.U',B'U'', B''.U, B''.U', B'''.U, ...
      DO IDB=1,5
        DO IDU=1,6-IDB
          TEMP = 0.D0
          DO I=1,3
            TEMP = TEMP + B(IDB,I)*U(IDU,I)
          ENDDO
          BU(IDB,IDU) = TEMP
        ENDDO
      ENDDO

CALCUL d^nB||/ds^n  = d^n((B.U)U)/ds^n
      DO I=1,3
        BP(1,I)=BU(1,1)* U(1,I)
        DBU  =    BU(2,1) +        BU(1,2)
        BP(2,I)  =  DBU *U(1,I) +      BU(1,1)*U(2,I)
        D2BU =    BU(3,1) + 2.D0*  BU(2,2) +       BU(1,3)
        BP(3,I)  =  D2BU*U(1,I) +     2.D0*DBU*U(2,I) +
     >  BU(1,1)*U(3,I)
        D3BU =    BU(4,1) + 3.D0*( BU(3,2) +       BU(2,3) )
     >   + BU(1,4)
        BP(4,I) =  D3BU*U(1,I) +   3.D0*(D2BU*U(2,I)  +
     >      DBU*U(3,I) ) +  BU(1,1)*U(4,I)

c       aaaa =  D3BU*U(1,I) +   3.D0*(D2BU*U(2,I)  +
c     >      DBU*U(3,I) ) +  BU(1,1)*U(4,I)
c         BP(4,I)  = aaaa
C////// debug...
c       if(noel.eq.20) then
c           write(*,*) ' i, BP(4,i) : ',i, bp(4,i),aaaa
c           write(*,*) D3BU, D2BU, DBU, BU(1,1)
c           write(*,*) U(1,I), U(2,I),  U(3,I),  U(4,I)
c           write(*,*) D3BU*U(1,I) +   3.D0*(D2BU*U(2,I)  +
c     >      DBU*U(3,I) ) +  BU(1,1)*U(4,I)
c                read(*,*)
c         endif
C////// debug...

        D4BU =    BU(5,1) + 4.D0*  BU(4,2) + 6.D0* BU(3,3) +
     >     4.D0 * BU(2,4) + BU(1,5)


        BP(5,I)  =  D4BU*U(1,I) +    4.D0*D3BU*U(2,I)  +
     >   6.D0*D2BU*U(3,I)   + 4.D0*DBU*U(4,I)  + BU(1,1)*U(5,I)


      ENDDO



CALCUL Omega
      DO ID=1,5
        O(ID,1) = GP*B(ID,1) + GM*BP(ID,1)
        O(ID,2) = GP*B(ID,2) + GM*BP(ID,2)
        O(ID,3) = GP*B(ID,3) + GM*BP(ID,3)

c       if(noel.eq.20) then
c           write(*,*) ' O : ', id, (o(id,ii),ii=1,3)
c           write(*,*) ' B : ', id, (b(id,ii),ii=1,3)
c           write(*,*) ' BP : ', id, (bp(id,ii),ii=1,3)  ! bp(4,1) is NaN
c         endif

      ENDDO


C Pure B field
      IF(KFLD.EQ.1) GOTO 10


C----- Elc pure (KFLD=2) or Mag+Elc (KFLD=3)
 2    CONTINUE
      FAC = 1.d0/AMT
      GM1 = - GM/G         ! gamma-1
      CSVM = CSV-1.D0
      CSVP = CSV+1.D0
      QM = QT/AMT

      GA1(1)= CSVM*CSVP*GM1
      GA1(2)= DCSV*CSVP*GM1 + CSVM*DCSV*GM1 + CSVM*CSVP*QM*EU(1)
      GA1(3)= D2CSV*CSVP*GM1 + CSVM*D2CSV*GM1+  CSVM*CSVP*QM*EU(2) +
     >2.D0*(DCSV*DCSV*GM1 + DCSV*CSVP*QM*EU(1) + CSVM*DCSV*QM*EU(1))
      GA1(4)=
     >D3CSV*CSVP*GM1 +
     >CSVM*D3CSV*GM1 +
     >CSVM*CSVP*QM*EU(3) +
     >3.D0*(
     >     D2CSV*DCSV*GM1 + DCSV*D2CSV*GM1 +
     >D2CSV*CSVP*QM*EU(1) + 2.D0*DCSV*DCSV*QM*EU(1)+ DCSV*CSVP*QM*EU(2)+
     >CSVM*D2CSV*QM*EU(1) + CSVM*DCSV*QM*EU(2)  )

C Compute the   beta*gamma/c * d^n(G + 1/(1+gamma))/ds^n   series, starting w n=0
CCC Erreur ?      GAC(1) = FAC * QBR * GA1(1)
      GGA1 = (G + GA1(1))
      GAC(1) = FAC * QBR * GGA1
      GAC(2) = FAC * (DQBR(1) * GGA1 + QBR * GA1(2))
      GAC(3) = FAC * (
     >    DQBR(2) * GGA1 + 2.D0 * DQBR(1) * GA1(2) +
     >    QBR * GA1(3))
      GAC(4) = FAC * (
     >    DQBR(3) * GGA1 + 3.D0 * DQBR(2) * GA1(2) +
     >    3.D0 * DQBR(1) * GA1(3) +   QBR * GA1(4))

c      write(*,fmt='(1a,1x,20e12.4)') ' spntrk fac : ',gg,gm,gp,g
c          write(*,fmt='(1a,1x,20e12.4)') ' spntrk dqbr : ',qbr,dqbr
c         write(*,fmt='(1a,1x,20e12.4)') ' spntrk ga1 : ',1.d0/csv,ga1(1)
c          write(*,fmt='(1a,1x,20e12.4)') ' spntrk gac : ',gac

      CALL PEU(E,U,
     >             VEU)

c          write(*,fmt='(1a,1x,20e12.4)') ' spntrk   ex,y,z : ',
c     >               e(1,1),e(1,2),e(1,3)
c          write(*,fmt='(1a,1x,20e12.4)') ' spntrk   bx,y,z : ',
c     >               b(1,1),b(1,2),b(1,3)
c          write(*,fmt='(1a,1x,20e12.4)') ' spntrk   u : ',u
c          write(*,fmt='(1a,1x,20e12.4)') ' spntrk   veu : ',veu
c                 read(*,*)

      DO I = 1, 3
        OE(1,I) =GAC(1)*VEU(1,I)
        OE(2,I) =GAC(2)*VEU(1,I)+     GAC(1)*VEU(2,I)
        OE(3,I) =GAC(3)*VEU(1,I)+2.D0*GAC(2)*VEU(2,I) +GAC(1)*VEU(3,I)
        OE(4,I) =GAC(4)*VEU(1,I)  +
     >                         3.D0*( GAC(3)*VEU(2,I) +GAC(2)*VEU(3,I))
     >                                                +GAC(1)*VEU(4,I)
        OE(5,I) =GAC(5)*VEU(1,I)+4.D0*GAC(4)*VEU(2,I)+
     >  6.d0*GAC(3)*VEU(3,I) +4.d0*GAC(2)*VEU(4,I) +   GAC(1)*VEU(5,I)
      ENDDO

c          write(*,fmt='(1a,1x,20e12.4)') ' spntrk oe : ',oe
c                 read(*,*)

C Initial value of \vec S (d^0S(I)/ds^0 == S(I))
      DO I = 1, 3
        DSE(1,I) = SF(I,IT)
      ENDDO

c          write(*,fmt='(1a,1x,20e12.4)') ' spntrk dse : ',dse
c                 read(*,*)
CALCUL S'=SxO
      DSE(2,1) = DSE(1,2)*OE(1,3) - DSE(1,3)*OE(1,2)
      DSE(2,2) = DSE(1,3)*OE(1,1) - DSE(1,1)*OE(1,3)
      DSE(2,3) = DSE(1,1)*OE(1,2) - DSE(1,2)*OE(1,1)
C      DO I = 1, 3
C        DSE(2,I) = DSE(2,I) - DBSB * DSE(1,I)
C      ENDDO
CALCUL S''=S'xO + SxO'
      DSE(3,1) = DSE(2,2)*OE(1,3) - DSE(2,3)*OE(1,2) +
     >           DSE(1,2)*OE(2,3) - DSE(1,3)*OE(2,2)
      DSE(3,2) = DSE(2,3)*OE(1,1) - DSE(2,1)*OE(1,3) +
     >           DSE(1,3)*OE(2,1) - DSE(1,1)*OE(2,3)
      DSE(3,3) = DSE(2,1)*OE(1,2) - DSE(2,2)*OE(1,1) +
     >           DSE(1,1)*OE(2,2) - DSE(1,2)*OE(2,1)
      DO I = 1, 3
C        DSE(3,I) = DSE(3,I) - D2BSB * DSE(1,I) - 2.D0*DBSB * DSE(2,I)
        DSE(3,I) = DSE(3,I) - DBSB * DSE(2,I)
      ENDDO


          goto 222



CALCUL S'''=S''xO + 2*S'xO' + SxO''
      DSE(4,1) = DSE(3,2)*OE(1,3) - DSE(3,3)*OE(1,2) +
     >                                  2.D0*(DSE(2,2)*OE(2,3)
     >  - DSE(2,3)*OE(2,2)) +       DSE(1,2)*OE(3,3) - DSE(1,3)*OE(3,2)
      DSE(4,2) = DSE(3,3)*OE(1,1) - DSE(3,1)*OE(1,3) +
     >                                  2.D0*(DSE(2,3)*OE(2,1)
     >  - DSE(2,1)*OE(2,3)) +       DSE(1,3)*OE(3,1) - DSE(1,1)*OE(3,3)
      DSE(4,3) = DSE(3,1)*OE(1,2) - DSE(3,2)*OE(1,1) +
     >                                  2.D0*(DSE(2,1)*OE(2,2)
     >  - DSE(2,2)*OE(2,1)) +       DSE(1,1)*OE(3,2) - DSE(1,2)*OE(3,1)
      DO I = 1, 3
C        DSE(4,I) = DSE(4,I) - D3BSB * DSE(1,I)
C     >             -  3.D0*( D2BSB * DSE(2,I)+ DBSB * DSE(3,I) )
        DSE(4,I) = DSE(4,I) - (D2BSB * DSE(2,I)+ 2.D0*DBSB * DSE(3,I) )
      ENDDO
CALCUL S''''=S'''xO + 3*S''xO' + 3*S'xO'' + SxO'''
      DSE(5,1) = DSE(4,2)*OE(1,3) - DSE(4,3)*OE(1,2) +
     >                                  3.D0*(DSE(3,2)*OE(2,3)
     >- DSE(3,3)*OE(2,2)) + 3.D0*(DSE(2,2)*OE(3,3) - DSE(2,3)*OE(3,2)) +
     >DSE(1,2)*OE(4,3) - DSE(1,3)*OE(4,2)
      DSE(5,2) = DSE(4,3)*OE(1,1) - DSE(4,1)*OE(1,3) +
     >                                  3.D0*(DSE(3,3)*OE(2,1)
     >- DSE(3,1)*OE(2,3)) + 3.D0*(DSE(2,3)*OE(3,1) - DSE(2,1)*OE(3,3)) +
     >DSE(1,3)*OE(4,1) - DSE(1,1)*OE(4,3)
      DSE(5,3) = DSE(4,1)*OE(1,2) - DSE(4,2)*OE(1,1) +
     >                                  3.D0*(DSE(3,1)*OE(2,2)
     >- DSE(3,2)*OE(2,1)) + 3.D0*(DSE(2,1)*OE(3,2) - DSE(2,2)*OE(3,1)) +
     >DSE(1,1)*OE(4,2) - DSE(1,2)*OE(4,1)
      DO I = 1, 3
C        DSE(5,I) = DSE(5,I) - D4BSB * DSE(1,I) - 4.D0*D3BSB * DSE(2,I)
C     >               - 6.D0*D2BSB * DSE(3,I)   - 4.D0*DBSB * DSE(4,I)
        DSE(5,I) = DSE(5,I) - D4BSB * DSE(1,I) - (D3BSB * DSE(2,I)
     >               + 3.D0*D2BSB * DSE(3,I) + 2.D0*DBSB * DSE(4,I))
      ENDDO
CALCUL S'''''=S''''xO + 4*S'''xO' + 6*S''xO'' + 4*S'xO''' + SxO''''
      DSE(6,1) = DSE(5,2)*OE(1,3) - DSE(5,3)*OE(1,2) +
     >                                  4.D0*(DSE(4,2)*OE(2,3)
     >- DSE(4,3)*OE(2,2)) + 6.D0*(DSE(3,2)*OE(3,3) - DSE(3,3)*OE(3,2)) +
     >4.D0*(DSE(2,2)*OE(4,3) -DSE(2,3)*OE(4,2)) +
     >                                 DSE(1,2)*OE(5,3)-DSE(1,3)*OE(5,2)
      DSE(6,2) = DSE(5,3)*OE(1,1) - DSE(5,1)*OE(1,3) +
     >                                  4.D0*(DSE(4,3)*OE(2,1)
     >- DSE(4,1)*OE(2,3)) + 6.D0*(DSE(3,3)*OE(3,1) - DSE(3,1)*OE(3,3)) +
     >4.D0*(DSE(2,3)*OE(4,1) -DSE(2,1)*OE(4,3)) +
     >                                 DSE(1,3)*OE(5,1)-DSE(1,1)*OE(5,3)
      DSE(6,3) = DSE(5,1)*OE(1,2) - DSE(5,2)*OE(1,1) +
     >                                  4.D0*(DSE(4,1)*OE(2,2)
     >- DSE(4,2)*OE(2,1)) + 6.D0*(DSE(3,1)*OE(3,2) - DSE(3,2)*OE(3,1)) +
     >4.D0*(DSE(2,1)*OE(4,2) -DSE(2,2)*OE(4,1)) +
     >                                 DSE(1,1)*OE(5,2)-DSE(1,2)*OE(5,1)
      DO I = 1, 3
C        DSE(6,I) = DSE(6,I) - D5BSB * DSE(1,I) - 5.D0*D4BSB * DSE(2,I)
C     >  -10.D0*(D3BSB *DSE(3,I) + D2BSB* DSE(4,I))- 5.D0*DBSB * DSE(5,I)
        DSE(6,I) = DSE(6,I) - D5BSB * DSE(1,I) - (D4BSB * DSE(2,I)
     >  +4.D0*D3BSB *DSE(3,I) + 2.D0*D2BSB* DSE(4,I)
     >          + 2.D0*DBSB * DSE(5,I))
      ENDDO

 222  continue

C Pure E field
      IF(KFLD.EQ.2) THEN

        DO I = 1, 3
          DO ID = 1, 6
            S(ID,I) = DSE(ID,I)
          ENDDO
        ENDDO

        GOTO 20

      ENDIF


C----- MAGNETIC+ELECTRIC FIELD
C----- B=b/Bro, E=e/Bro

Compute s' = s x omega and derivatives
 10   CONTINUE

      S(1,1) = SF(1,IT)
      S(1,2) = SF(2,IT)
      S(1,3) = SF(3,IT)
CALCUL S'=SxO
      S(2,1) = S(1,2)*O(1,3) - S(1,3)*O(1,2)
      S(2,2) = S(1,3)*O(1,1) - S(1,1)*O(1,3)
      S(2,3) = S(1,1)*O(1,2) - S(1,2)*O(1,1)
      IF(KFLD.EQ.3) THEN
        DO I = 1, 3
          S(2,I) = S(2,I) + DSE(2,I)
        ENDDO
      ENDIF

CALCUL S''=S'xO + SxO'
      S(3,1) = S(2,2)*O(1,3) - S(2,3)*O(1,2) +
     >       S(1,2)*O(2,3) - S(1,3)*O(2,2)
      S(3,2) = S(2,3)*O(1,1) - S(2,1)*O(1,3) +
     >       S(1,3)*O(2,1) - S(1,1)*O(2,3)
      S(3,3) = S(2,1)*O(1,2) - S(2,2)*O(1,1) +
     >       S(1,1)*O(2,2) - S(1,2)*O(2,1)
      IF(KFLD.EQ.3) THEN
        DO I = 1, 3
          S(3,I) = S(3,I) + DSE(3,I)
        ENDDO
      ENDIF
CALCUL S'''=S''xO + 2*S'xO' + SxO''
      S(4,1) = S(3,2)*O(1,3) - S(3,3)*O(1,2) + 2.D0*(S(2,2)*O(2,3)
     >  - S(2,3)*O(2,2)) +       S(1,2)*O(3,3) - S(1,3)*O(3,2)
      S(4,2) = S(3,3)*O(1,1) - S(3,1)*O(1,3) + 2.D0*(S(2,3)*O(2,1)
     >  - S(2,1)*O(2,3)) +       S(1,3)*O(3,1) - S(1,1)*O(3,3)
      S(4,3) = S(3,1)*O(1,2) - S(3,2)*O(1,1) + 2.D0*(S(2,1)*O(2,2)
     >  - S(2,2)*O(2,1)) +       S(1,1)*O(3,2) - S(1,2)*O(3,1)
      IF(KFLD.EQ.3) THEN
        DO I = 1, 3
          S(4,I) = S(4,I) + DSE(4,I)
        ENDDO
      ENDIF
CALCUL S''''=S'''xO + 3*S''xO' + 3*S'xO'' + SxO'''
      S(5,1) = S(4,2)*O(1,3) - S(4,3)*O(1,2) + 3.D0*(S(3,2)*O(2,3)
     >- S(3,3)*O(2,2)) + 3.D0*(S(2,2)*O(3,3) - S(2,3)*O(3,2)) +
     >S(1,2)*O(4,3) - S(1,3)*O(4,2)
      S(5,2) = S(4,3)*O(1,1) - S(4,1)*O(1,3) + 3.D0*(S(3,3)*O(2,1)
     >- S(3,1)*O(2,3)) + 3.D0*(S(2,3)*O(3,1) - S(2,1)*O(3,3)) +
     >S(1,3)*O(4,1) - S(1,1)*O(4,3)
      S(5,3) = S(4,1)*O(1,2) - S(4,2)*O(1,1) + 3.D0*(S(3,1)*O(2,2)
     >- S(3,2)*O(2,1)) + 3.D0*(S(2,1)*O(3,2) - S(2,2)*O(3,1)) +
     >S(1,1)*O(4,2) - S(1,2)*O(4,1)
      IF(KFLD.EQ.3) THEN
        DO I = 1, 3
          S(5,I) = S(5,I) + DSE(5,I)
        ENDDO
      ENDIF
CALCUL S'''''=S''''xO + 4*S'''xO' + 6*S''xO'' + 4*S'xO''' + SxO''''
      S(6,1) = S(5,2)*O(1,3) - S(5,3)*O(1,2) + 4.D0*(S(4,2)*O(2,3)
     >- S(4,3)*O(2,2)) + 6.D0*(S(3,2)*O(3,3) - S(3,3)*O(3,2)) +
     >4.D0*(S(2,2)*O(4,3) -S(2,3)*O(4,2)) + S(1,2)*O(5,3)-S(1,3)*O(5,2)
      S(6,2) = S(5,3)*O(1,1) - S(5,1)*O(1,3) + 4.D0*(S(4,3)*O(2,1)
     >- S(4,1)*O(2,3)) + 6.D0*(S(3,3)*O(3,1) - S(3,1)*O(3,3)) +
     >4.D0*(S(2,3)*O(4,1) -S(2,1)*O(4,3)) + S(1,3)*O(5,1)-S(1,1)*O(5,3)
      S(6,3) = S(5,1)*O(1,2) - S(5,2)*O(1,1) + 4.D0*(S(4,1)*O(2,2)
     >- S(4,2)*O(2,1)) + 6.D0*(S(3,1)*O(3,2) - S(3,2)*O(3,1)) +
     >4.D0*(S(2,1)*O(4,2) -S(2,2)*O(4,1)) + S(1,1)*O(5,2)-S(1,2)*O(5,1)
      IF(KFLD.EQ.3) THEN
        DO I = 1, 3
          S(6,I) = S(6,I) + DSE(6,I)
        ENDDO
      ENDIF



 20   CONTINUE

C Now push the spin one step further.
CALCUL S(s+ds) = S(s) +S'(s)ds +S''(s)ds^2/2 +S'''(s)ds^3/6 +S''''(s)ds4/24...
      SF1IT=S(1,1)+(S(2,1)+(S(3,1)/2.D0+(S(4,1)/6.D0+(S(5,1)/24.D0
     > +S(6,1)/120.D0*DS)*DS)*DS)*DS)*DS
      SF2IT=S(1,2)+(S(2,2)+(S(3,2)/2.D0+(S(4,2)/6.D0+(S(5,2)/24.D0
     > +S(6,2)/120.D0*DS)*DS)*DS)*DS)*DS
      SF3IT=S(1,3)+(S(2,3)+(S(3,3)/2.D0+(S(4,3)/6.D0+(S(5,3)/24.D0
     > +S(6,3)/120.D0*DS)*DS)*DS)*DS)*DS

      AN=SQRT(SF1IT*SF1IT+SF2IT*SF2IT+SF3IT*SF3IT)
      SF4IT = AN

C----- NORMALISATION (may be useless... check step size instead...)
      SF1IT=SF1IT/AN
      SF2IT=SF2IT/AN
      SF3IT=SF3IT/AN
      SF4IT = SQRT(SF1IT*SF1IT+SF2IT*SF2IT+SF3IT*SF3IT)
C--------------------------------------------------------------------

      SF(1,IT) = SF1IT
      SF(2,IT) = SF2IT
      SF(3,IT) = SF3IT
      SF(4,IT) = SF4IT

      DO ICOO = 1, 4
        IF(SMI(ICOO,IT).GT.SF(ICOO,IT)) SMI(ICOO,IT) = SF(ICOO,IT)
        IF(SMA(ICOO,IT).LT.SF(ICOO,IT)) SMA(ICOO,IT) = SF(ICOO,IT)
      ENDDO

C////// debug...
c      call ZGNOEL(
c     >             NOEL)
c       if(noel.eq.20) then
c           write(*,*) 'spntrk OUT '
c           write(*,*) ' S_X,Y,Z, |S| : ',it,(sf(i,it),i=1,4)
c         endif
C///////////////////////

      RETURN

      ENTRY SPNTR2(IMAX)
      DO IIT = 1, IMAX
        DO ICOO = 1, 4
          SMI(ICOO,IIT) = +1.D10
          SMA(ICOO,IIT) = -1.D10
        ENDDO
      ENDDO
      RETURN

      ENTRY SPNTR3(IMAX,
     >                  SPMI,SPMA)
      DO IIT = 1, IMAX
        DO ICOO = 1, 4
          SPMI(ICOO,IIT) = SMI(ICOO,IIT)
          SPMA(ICOO,IIT) = SMA(ICOO,IIT)
        ENDDO
      ENDDO
      RETURN

      ENTRY SPNTR4(IMAX,SPMI,SPMA)
      DO IIT = 1, IMAX
        DO ICOO = 1, 4
          SMI(ICOO,IIT) = SPMI(ICOO,IIT)
          SMA(ICOO,IIT) = SPMA(ICOO,IIT)
        ENDDO
      ENDDO
      RETURN
      END
