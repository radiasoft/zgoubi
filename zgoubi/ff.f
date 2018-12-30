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
      FUNCTION FF()
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.CONST_2.H"     ! COMMON/CONST/ CL9,CL,PI,RAD,DEG,QEL,AMPROT,CM2M
      PARAMETER (MXV=60)
      INCLUDE "C.CONTR.H"     ! COMMON/CONTR/ VAT(MXV),XI(MXV)
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      INCLUDE "C.FINIT.H"     !COMMON/FINIT/ FINI
      INCLUDE "C.OBJET.H"     ! COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT,KZOB
      INCLUDE "C.SPIN.H"     ! COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)
      INCLUDE "C.VARY.H"  ! COMMON/VARY/ NV,IR(MXV),NC,I1(MXV),I2(MXV),V(MXV),IS(MXV),W(MXV),
                          !     >IC(MXV),IC2(MXV),I3(MXV),XCOU(MXV),CPAR(MXV,27)

      DIMENSION U(6,6), T(6,6,6)
      CHARACTER(1) BLANC
      CHARACTER(10) AST
      PARAMETER (BLANC=' ')
C      PARAMETER (AST  ='*')
      PARAMETER (AST  ='improving!')
      PARAMETER (NMAIL  = 1)

      DIMENSION F0P(6,6), F0PD(6,6), F0MD(6,6)
      PARAMETER (MXRF = 7)
      DIMENSION YNUI(MXRF), ZNUI(MXRF)

C      DIMENSION IC2(MXV)
C      SAVE IC2

      DIMENSION RPD(6,6), RMD(6,6)

      LOGICAL READAT
      PARAMETER (I0=0,ZRO=0.D0)

      DIMENSION SMAT(3,3)
      DIMENSION TR(3)

      M=1
      CALL REMPLI(M)
      Z=0D0
      KK=0


      DO 3 I=1,NC
         ICONT=IC(I)
         ICONT2=IC2(I)
         K=I1(I)
         L=I2(I)

         IF(I3(I) .NE. KK) THEN
            KK=I3(I)
C----------- Compute present value of constraint
C            FITING = .FALSE.
C            CALL FITSTA(6,FITING)
            READAT = .FALSE.

            CALL ZGOUBI(1,KK,READAT,
     >                              NBEL)

         ENDIF

         IF     (ICONT .EQ. 0) THEN
C----------- Constrain BEAM matrix [sigma_ij]. Normally used with OBJECT, Kobj=5or 6
C 11=beta_y, 12=21=-alpha_y, 22=gamma_y ; 3-4 for ._z ; 5-6 for dp-s
           IORD=1
           IF(KOBJ .EQ. 6) IORD=2

           IF    (ICONT2 .EQ. 0) THEN
             CALL COEFFS(1,IORD,U,T,1,
     >                                F0P,Cstrn)
             VAL= F0P(K,L)

           ELSEIF(ICONT2 .GE. 1 .AND. ICONT2 .LE. 9) THEN
C------------ Periodic case :  Twiss coefficients or tunes
C  R16=periodic dispersion D_y, R26= its derivative D'_y ; 36=D_z, 46=D'_z
C ICONT2 uses the reference trajectory as defined in OBJET/KOBJ=5.2-6

             CALL COEFFS(0,IORD,U,T,ICONT2,
     >                                     F0P,Cstrn)
             CALL TUNES(U,F0P,NMAIL,IERY,IERZ,.FALSE.,
     >                                                YNU,ZNU,CMUY,CMUZ)
             IF(K .LE. 6 .AND. L .LE. 6 ) THEN
               VAL= F0P(K,L)
             ELSEIF( K .EQ. 7 ) THEN
C--------------Constraint  Y  tune
               VAL= YNU
             ELSEIF( K .EQ. 8 ) THEN
               VAL= ZNU
             ELSEIF( K .EQ. 9 ) THEN
C--------------Constraint  cos(muY) = trace/2
               VAL= CMUY
             ELSEIF( K .EQ. 10 ) THEN
               VAL= CMUZ
             ENDIF
           ENDIF

         ELSEIF(ICONT .EQ. 1) THEN

           IORD=1

           IF    (ICONT2 .EQ. 0) THEN
C-----------Contraints are first order transport coeffs

             IF(KOBJ .EQ. 6) IORD=2
             CALL COEFFS(0,IORD,U,T,1,
     >                                F0P,Cstrn)
             IF(K .LE. 6 .AND. L .LE. 6 ) THEN
               VAL= U(K,L)

             ELSEIF( K .EQ. 7 ) THEN
               IF(.NOT. L .EQ. 8 ) THEN
C---------------Contraint determinant H
                 VAL= U(1,1)*U(2,2)-U(1,2)*U(2,1)
               ELSE
C---------------Contraint determinant HV
                 VAL= U(1,3)*U(2,4)-U(1,4)*U(2,3)
               ENDIF
             ELSEIF( K .EQ. 8 ) THEN
               IF(.NOT. L .EQ. 7 ) THEN
C---------------Contraint determinant V
                 VAL= U(3,3)*U(4,4)-U(3,4)*U(4,3)
               ELSE
C---------------Contraint determinant VH
                 VAL= U(3,1)*U(4,2)-U(3,2)*U(4,1)
               ENDIF
             ENDIF

           ELSEIF(ICONT2 .EQ. 1) THEN
C-----------Chromaticity

             CALL OBJ51(
     >                  NBREF)
             IREF = 0
 1           CONTINUE
               IREF = IREF + 1
               IF(IREF.GT.MXRF) STOP ' SBR FF '
               IT1 = 1 + 11 * (IREF-1)
               IT2 = IT1+3
               IT3 = IT1+4
               IFC = 0
               CALL REFER(1,IORD,IFC,IT1,IT2,IT3)
               CALL MAT1(IT1,F,IMAX,
     >                              U,T)
               CALL REFER(2,IORD,IFC,IT1,IT2,IT3)
               CALL TUNES(U,F0P,NMAIL,IERY,IERZ,.FALSE.,
     >                                                YNU,ZNU,CMUY,CMUZ)
               YNUI(IREF) = YNU
               ZNUI(IREF) = ZNU

             IF(IREF.LT.NBREF) GOTO 1

             IF    (NBREF .EQ. 1) THEN
               STOP ' SBR FF : not enough MATRIX blocks  '
             ELSEIF(NBREF .EQ. 2) THEN
               DP1 =  F(1,1)
               DP2 =  F(1,12)
               DNUYDP = (YNUI(2) - YNUI(1)) / (DP2-DP1)
               DNUZDP = (ZNUI(2) - ZNUI(1)) / (DP2-DP1)

             ELSEIF(NBREF .EQ. 3) THEN
               STOP ' SBR FF : too many MATRIX blocks  '
C              q' = (q+ - q-) / 2dp
C              q'' = (q+ + q- -2q0) / dp^2
             ELSE
               STOP ' SBR FF : too many MATRIX blocks  '
             ENDIF
             IF( K .EQ. 7 ) THEN
C-------------Contraint dNu_Y/dpp
               VAL= DNUYDP
             ELSEIF( K .EQ. 8 ) THEN
C-------------Contraint dNu_Z/dpp
               VAL= DNUZDP
             ENDIF

           ENDIF

         ELSEIF(ICONT .EQ. 2) THEN
C-----------Contraints are second order transport coeffs

           IORD=2

           IF    (ICONT2 .EQ. 0) THEN

             CALL COEFFS(0,IORD,U,T,1,
     >                                F0P,Cstrn)
             L1=L/10
             L2=L-10*L1
             VAL= T(K,L1,L2)

           ELSEIF(ICONT2 .EQ. 1) THEN

             CALL COEFFS(0,IORD,U,T,1,
     >                                F0P,Cstrn)
             CALL MAT2P(RPD,DP)
             CALL TUNES(RPD,F0PD,NMAIL,IERY,IERZ,.TRUE.,
     >                                             YNUP,ZNUP,CMUY,CMUZ)
             CALL MAT2M(RMD,DP)
             CALL TUNES(RMD,F0MD,NMAIL,IERY,IERZ,.TRUE.,
     >                                             YNUM,ZNUM,CMUY,CMUZ)
             DNUYDP = (YNUP-YNUM)/2.D0/DP
             DNUZDP = (ZNUP-ZNUM)/2.D0/DP

             IF( K .EQ. 7 ) THEN
C-------------Contraint dNu_Y/dpp
               VAL= DNUYDP
             ELSEIF( K .EQ. 8 ) THEN
C-------------Contraint dNu_Z/dpp
               VAL= DNUZDP
             ENDIF

           ENDIF

        ELSE IF(ICONT .EQ. 3) THEN
C----------- Constraints on particle coordinates or bundle

           IF    (ICONT2.EQ.0) THEN
             IF(K .GT. 0) THEN
C-------------- Constraint is value of coordinate L of particle K
               VAL=F(L,K)
             ELSEIF(K.EQ.-1) THEN
C------------ Constraint on beam : average value of coordinate L
               if    (NINT(CPAR(I,1)) .eq. 0) then
                 VAL=FITLAV(L,1,imax,ZRO)
               elseif(NINT(CPAR(I,1)) .eq. 2) then
                 VAL=FITLAV(L,NINT(CPAR(I,2)),NINT(CPAR(I,3)),ZRO)
               endif
             ELSEIF(K.EQ.-2) THEN
C------------- Constraint on beam : max value of coordinate L
               VAL=FITLMA(L)
             ELSEIF(K.EQ.-3) THEN
C------------ Constraint : minimze distance between paticles for coord. L.
C             List of particles concerned is entered via parameter list
               CALL DIST2(L,
     >                      VAL)
             ELSEIF(K.EQ.-4) THEN
C------------ Constraint : bring distance between bunch centroids at various pickups closest to VAL.
C             PUs concerned are those included in the range NOELA - NOELB.
C             What that does : get PU signals for a given coordinate, then
C             dist3 computes the sum of the absolute values of the differences between these averages.
               CALL DIST3(L,nint(CPAR(I,2)),NINT(CPAR(I,3)),
     >                                                      VAL)
             ENDIF
           ELSEIF(ICONT2.EQ.1) THEN
C------------ Constraint on closed orbit :
C             e.g., particle #K has equal values for coordinate L,
C                 at ends of cell
C                   (hence expected constraint value in zgoubi.dat is 0).
             VAL=ABS(F(L,K) - FO(L,K))

           ELSEIF(ICONT2.EQ.2) THEN
C------------ Constraint on closed orbit :
C             or opposite angle values (L=3,5)
C                 at ends of cell
C                   (hence expected constraint value in zgoubi.dat is 0)
             VAL=ABS(F(L,K) + FO(L,K))

c           ELSEIF(ICONT2.EQ.3) THEN
cC------------ Constraint on min/max value (MIMA=1/2) of coordinate L reached inside optical element KK
cC             MIMA = 1
c             MIMA = NINT(CPAR(I,2))
c             CALL FITMM1(K,L,KK,MIMA,
c     >                             VAL1)
c             if(mima .eq. 1) mima=2
c             if(mima .eq. 2) mima=1
c             CALL FITMM1(K,L,KK,MIMA,
c     >                             VAL2)
c             val = val1 + val2

           ELSEIF(ICONT2.EQ.3) THEN
C------------ Same coordinate L, two different particles K, K2
             K2 = NINT(CPAR(I,2))
             VAL=ABS( F(L,K) + F(L,K2) )

           ELSEIF(ICONT2.EQ.4) THEN
C------------ Same coordinate L, two different particles K, K2
             K2 = NINT(CPAR(I,2))
             VAL=ABS( F(L,K) - F(L,K2) )

           ELSEIF(ICONT2.EQ.5) THEN
C------------ Same coordinate L, two different particles K, K2
             K2 = NINT(CPAR(I,2))
             VAL=F(L,K) / F(L,K2) -1.D0

           ENDIF

         ELSE IF(ICONT .EQ. 4) THEN
C----------- Constraint ellipse parameters

            IF    (K.LE.2) THEN
              K2=1
            ELSEIF(K.LE.4) THEN
              K2=2
            ELSEIF(K.LE.6) THEN
              K2=3
            ENDIF
            CALL LPSFIT(K2,
     >                     EMIT,ALP,BET,XM,XPM)
            IF(K.EQ.L) THEN
              IF(K.EQ.2*K2) THEN
                VAL = (1.D0+ALP*ALP)/BET
              ELSE
                VAL = BET
              ENDIF
            ELSE
              VAL = ALP
            ENDIF

         ELSE IF(ICONT .EQ. 5) THEN

            IF( K .EQ. -1 ) THEN
C------------- Ratio of iex>0 particles,  (IMAX-Nstopped)/IMAX
              NLIV = 0
              DO 5 IT=1,IMAX
 5              IF(IEX(IT).GE.-1) NLIV = NLIV+1
              NJUMP = NLIV/100
              IF(NJUMP.EQ.0) NJUMP=1
              VAL = ((DBLE(NJUMP)*(NLIV/NJUMP))/DBLE(IMAX))

            ELSE IF( K .GE. 1) THEN

              IF( K .LE. 3) THEN
C--------------- Ratio N_InLips/IMAX of particles within ph-space ellipse ;
C                ellipse surface is entered as parameter in constraint data list
                CALL ACCEN(K,
     >                       RATIN)
                VAL=RATIN
C                VAL = 1.D0/VAL
CCCCCCC     ?? remplace VAL = 1.D0/VAL par VAL = 1.D0-VAL dans FF ????

              ELSEIF( K .LE. 6) THEN
C--------------- Maximum ratio N_InLips/IMAX  of particles within a
C                ph-space ellipse which is found at best in region of rms ellipse ;
C                ellipse surface is entered as parameter in constraint data list
                CALL ACCEP(K-3,
     >                         EMIT,ALP,BET,XM,XPM,NLIV,MXINL)
C     >                         EMIT,ALP,BET,XM,XPM,NLIV,NINL)

C--------------- For a follow-up of the optimized ellipse, using zpop (Menu8/14/7/6) ------
c                WRITE(88,FMT='(1P,I7,5G12.4,A13)')
c     >              MXINL,EMIT,ALP,BET,XM,XPM,'  Function FF'
c                BACKSPACE(88)
c                READ(88,*) DUM
C-----------------------------------------------------------------------------------------

                VAL=DBLE(MXINL)/DBLE(IMAX)
C                VAL=DBLE(NINL)/DBLE(IMAX)

              ELSE
                STOP ' *** Error, FCT FF -> No such constraint 5.K'
              ENDIF

            ELSE
              STOP ' *** Error, FCT FF -> No such constraint 5.K'
            ENDIF

         ELSE IF(ICONT .EQ. 6) THEN
C    constraint rms  emittance

            CALL LPSFIT(K,
     >                    EMIT,ALP,BET,XM,XPM)
            VAL=EMIT

        ELSE IF(ICONT .EQ. 7) THEN
C----------- Constraints on coordinates and fields *inside* optical elements

           IF(ICONT2.EQ.1) THEN
C------------ Constraint on min or max value (MIMA=1 or 2) of coordinate L reached inside optical element KK
             MIMA = NINT(CPAR(I,2))
             IF(MIMA .NE. 1 .AND.
     >          MIMA .NE. 2)
     >          CALL KSTOP(' FF, MIMA should = 1 or 2. ',-99)

             CALL FITMM1(K,L,KK,MIMA,icont2,
     >                                    VAL)

           ELSEIF(ICONT2.EQ.2) THEN
C------------ Constraint on |min-max| value of coordinate L reached inside optical element KK
             MIMA = 1
             CALL FITMM1(K,L,KK,MIMA,icont2,
     >                                    VAL1)
             MIMA=2
             CALL FITMM1(K,L,KK,MIMA,icont2,
     >                                    VAL2)
             VAL = VAL2 - VAL1

           ELSEIF(ICONT2.EQ.3) THEN
C------------ Constraint on  min+max  value of coordinate L reached inside optical element KK
             MIMA = 1
             CALL FITMM1(K,L,KK,MIMA,icont2,
     >                                    VAL1)
             MIMA=2
             CALL FITMM1(K,L,KK,MIMA,icont2,
     >                                    VAL2)
             VAL = VAL1 + VAL2

           ELSEIF(ICONT2.EQ.6) THEN
C------------ Constraint on  min or max value  (MIMA=1 or 2) of field L-component, across optical element KK
             MIMA = NINT(CPAR(I,2))
             IF(MIMA .NE. 1 .AND.
     >          MIMA .NE. 2)
     >          CALL KSTOP(' FF, MIMA should = 1 or 2. ',-99)

             CALL FITMM1(K,L,KK,MIMA,icont2,
     >                                    VAL)

           ELSEIF(ICONT2.EQ.7) THEN
C------------ Constraint on  |min-max| value of field L-component, across optical element KK
             MIMA = 1
             CALL FITMM1(K,L,KK,MIMA,icont2,
     >                                    VAL1)
             MIMA = 2
             CALL FITMM1(K,L,KK,MIMA,icont2,
     >                                    VAL2)
             VAL = VAL2 - VAL1

           ELSEIF(ICONT2.EQ.8) THEN
C------------ Constraint on  min+max  value of field L-component, across optical element KK
             MIMA = 1
             CALL FITMM1(K,L,KK,MIMA,ICONT2,
     >                                    VAL1)
             MIMA = 2
             CALL FITMM1(K,L,KK,MIMA,ICONT2,
     >                                    VAL2)
             VAL = VAL1 + VAL2

           ELSEIF(ICONT2.EQ.9) THEN
C------------ Constraint on integral of field L-component for particle K, along optical element KK
             MIMA = 0
             CALL FITMM1(K,L,KK,MIMA,ICONT2,
     >                                      VAL)

           ELSEIF(ICONT2.EQ.10) THEN
C------------ Constraint on value of
C             min. coordinate L of particle K reached inside optical element KK, and
C             max. coordinate L2 of particle K2 inside optical element KK.
C Ex. of use: in CBEAT linear FFAG cell, allows centering the orbits in the QF quad such that
C extreme excursions (particle K,L (resp. K2,L2) has extreme negative (resp. positiv) excur.) have equal abs. value / opposite sign.

             MIMA=1
             CALL FITMM1(K,L,KK,MIMA,ICONT2,
     >                                    VAL1)
             K2 = NINT(CPAR(I,2))
             L2 = NINT(CPAR(I,3))
             MIMA=2
             CALL FITMM1(K2,L2,KK,MIMA,ICONT2,
     >                                      VAL2)
             VAL = VAL2 + VAL1

           ELSE
             CALL ENDJOB(' SBR ff.f : no such FIT option 7.',ICONT2)

           ENDIF

         ELSE IF(ICONT .EQ. 10) THEN
C-----------Contraints on spin

           IF    (ICONT2.EQ.0) THEN

             VAL=SF(L,K)

           ELSEIF(ICONT2.EQ.1) THEN
C------------ Constraint on spin components,
C             e.g., equal spin values  (L=2,4 or other), at ends of cell
C                   (hence expected value for the constraint is 0)
             VAL=ABS(SF(L,K) - SI(L,K))

           ELSEIF(ICONT2.EQ.2) THEN
C------------ Constraint on spin rotation angle of momentum group #K.
C             Requires OBJET/KOBJ=2, w/ groups of 3 particles,
C             all particles in a group have same momenta and respective spins in direction X, Y, Z

             CALL SPNMAT(K,
     >                     SMAT,TRM, SROT,TR(1),TR(2),TR(3),QS)
             VAL = SROT

           ELSEIF(ICONT2.EQ.3) THEN
C------------ Constraint on spin rotation axis of momentum group #K.
C             Requires OBJET/KOBJ=2, w/ groups of 3 particles,
C             all particles in a group have same momenta and respective spins in direction X, Y, Z

             CALL SPNMAT(K,
     >                     SMAT,TRM, SROT,TR(1),TR(2),TR(3),QS)
             VAL = TR(L)

           ELSEIF(ICONT2.EQ.4) THEN
C------------ Constraint on spin rotation axis of momentum group #K.
C             Requires OBJET/KOBJ=2, w/ groups of 3 particles,
C             all particles in a group have same momenta and respective spins in direction X, Y, Z

             VAL = 0.D0
             jj = 0
             DO  II = 1, IMAX
                 jj = jj + 1
               VAL = VAL + SF(3,II)
             ENDDO
             val = val /dble(jj)

           ELSE
             CALL ENDJOB(' SBR ff.f : no such FIT option 10.',ICONT2)
           ENDIF

         ENDIF

         Z=Z+((VAL-V(I))/W(I))**2
         VAT(I)=VAL

 3    CONTINUE

      IF(FINI.LE.Z) THEN
C        WRITE(6,100) CHAR(13),IV,FINI,Z,BLANC
C100      FORMAT(1H+,A1,I3,1P,2E20.5,1X,A10,$)
      ELSE
C         WRITE(6,100) CHAR(13),IV,FINI,Z,AST
         FINI=Z
      ENDIF
      FF=Z
      RETURN

      END
