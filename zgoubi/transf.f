C  ZGOUBI, a program for computing the trajectories of charged particles
C  in electric and magnetic fields
C  Copyright (C) 1988-2007  Fran�ois M�ot
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
C  Fran�ois M�ot <fmeot@bnl.gov>
C  Brookhaven National Laboratory 
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  -------
      SUBROUTINE TRANSF(QSHROE,VSHROE,QSHROS,VSHROS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ----------------------------------------------------------------
C     APPELE PAR AIMANT ET QUASEX .
C     APPELLE INTEGR POUR LE CALCUL DES IMAX TRAJECTOIRES DU FAISCEAU.
C     ----------------------------------------------------------------
      PARAMETER (MSR=8)
      CHARACTER(2) QSHROE(MSR),QSHROS(MSR)
      DIMENSION VSHROE(MSR),    VSHROS(MSR)

      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      PARAMETER(MCOEF=6)
      INCLUDE "C.CHAFUI.H"     ! COMMON/CHAFUI/ XE,XS,CE(MCOEF),CS(MCOEF),QCE(MCOEF),QCS(MCOEF)
      INCLUDE "MAXTRA.H"
      INCLUDE "C.CHAMBR.H"     ! COMMON/CHAMBR/ LIMIT,IFORM,YLIM2,ZLIM2,SORT(MXT),FMAG,YCH,ZCH
 
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE "C.CONST2.H"     ! COMMON/CONST2/ ZERO, UN
      INCLUDE 'MXSTEP.H'
      INCLUDE 'CSR.H'
      INCLUDE "C.DESIN.H"     ! COMMON/DESIN/ FDES(7,MXT),IFDES,KINFO,IRSAR,IRTET,IRPHI,NDES
C     >,AMS,AMP,AM3,TDVM,TETPHI(2,MXT)
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "C.DROITE_2.H"     ! COMMON/DROITE/ AM(9),BM(9),CM(9),IDRT
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      CHARACTER(1) LET
      INCLUDE "C.FAISCT.H"     ! COMMON/FAISCT/ LET(MXT)
      INCLUDE "C.GASC.H"     ! COMMON/GASC/ AI, DEN, KGA
      INCLUDE "C.MARK.H"     ! COMMON/MARK/ KART,KALC,KERK,KUASEX
      INCLUDE "C.OBJET.H"     ! COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT,KZOB
C      LOGICAL ZSYM
      INCLUDE "C.TYPFLD.H"     ! COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      INCLUDE "C.PTICUL_2.H"     ! COMMON/PTICUL/ AAM,Q,G,TO
      INCLUDE "C.SPIN.H"     ! COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)
      INCLUDE 'MXFS.H'
      PARAMETER (LBLSIZ=20)
      PARAMETER (KSIZ=10)
      CHARACTER(KSIZ) FAM ; CHARACTER(LBLSIZ) LBF
      INCLUDE "C.SCALT.H"     ! COMMON/SCALT/ FAM(MXF),LBF(MXF,MLF)
      INCLUDE "C.SYNRA.H"     ! COMMON/SYNRA/ KSYN
      INCLUDE "C.TRAJ.H"     ! COMMON/TRAJ/ Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT
      INCLUDE "C.TRNSF.H"     ! COMMON/TRNSF/ XFE,XFS

      LOGICAL CHGRFE,CHGRFS,CHPFE,CHPFS,EVNT,BACKW,DRT
     
      PARAMETER(MPOL=10)
      
      LOGICAL MIRROR, MIRRIN, BACKIN, MIRROU, BACKOU
      SAVE MIRROR, BACKW

      LOGICAL CONSTY, CONSTI
      SAVE CONSTY
      
      DATA MIRROR / .FALSE. /
      DATA CONSTY / .FALSE. /
      DATA BACKW / .FALSE. /

      IF(LST .GE. 1 ) CALL CTRLB(1)

      CHGRFE= NINT(VSHROE(MSR)).NE.0
      CHGRFS= NINT(VSHROS(MSR)).NE.0
 
      CHPFE = ( KART .EQ. 1 .AND. XFE .NE. ZERO )
      CHPFS = ( KART .EQ. 1 .AND. XFS .NE. ZERO )
 
C----- Trajectory deviations T > Pi/2 are allowed in the field map of SPES3
C      BACKW = ( KALC .EQ. 2 .AND. KUASEX .EQ. 3 )
C      also with CARTEMES, 
C     >     .OR. ( KALC .EQ. 2 .AND. KUASEX .EQ. 1 )

C----- Events, such as spin tracking, in-flight decay, etc...
      EVNT = KSPN .EQ. 1 .OR. IFDES .EQ. 1 .OR. KGA .EQ. 1 .OR.
     >   LIMIT .EQ. 1  .OR. KSYN.GE.1  .OR. KCSR.EQ.1
 
C----- Droite de coupure entree
      DRT = IDRT .EQ. -1 .OR. IDRT .GE. 2

C------ Some initialisations in SBR DEVTRA...
      IF(KFLD.GE.LC) CALL DEVTRW(CL)

      XO=X
      DO 1 IT=1,IMAX
 
C-------- IEX<-1 <=> Particle stopped
        IF(IEX(IT) .LT. -1) GOTO 1
 
        IF(IT .EQ. IREP(IT) .OR. .NOT. ZSYM
     >     .OR. CONSTY                      ) THEN

          CALL INITRA(IT)
          X=XO

          IF( CHGRFE ) CALL CHANRF(EVNT,QSHROE,VSHROE)
          IF( CHPFE  ) CALL CHAREF(EVNT,XFE,ZERO,ZERO)

          IF( DRT  ) CALL DRTENT
 
          CALL INTEGR(EVNT,BACKW,MIRROR,KFLD)
 
          IF( CHPFS  ) CALL CHAREF(EVNT,XFS,ZERO,ZERO)

          IF( CHGRFS ) THEN
C--------print phi and theta in mrad, and time in ns--------------------
            IF((NRES.GT.0) .AND. (IT.EQ.1)) WRITE(NRES,199)
            IF(NRES.GT.0) WRITE(NRES,100) KEX,(FO(J,IT),J=1,5),
     > X, Y, T*1000.D0, Z, P*1000.D0, SAR, TAR*1.0D-2, IT

 199  FORMAT(2X, '  KPOS  DP         Y(cm)   T(mrad)     ',
     >  'Z(cm)   P(mrad)   |', 2X, ' X(cm)         Y(cm)        ',
     >  'T(mrdd)      Z(cm)        P(mrad)     S(cm)         Time(ns)',
     >   4X, '  I')
 100  FORMAT(2X,'A',2X,I3,F8.4,4F10.3, '   |', F12.6, 6F13.6, 2X, I5)

              CALL CHANRF(EVNT,QSHROS,VSHROS)
          ENDIF

          CALL MAJTRA(IT)

        ELSE

          CALL DEJACA(IT)

        ENDIF
 
 1    CONTINUE
 
      IF(LST .GE. 1 ) CALL CTRLB(2)
 
      MIRROR = .FALSE.
      BACKW = .FALSE.
      RETURN

      ENTRY TRANSW(MIRRIN,BACKIN)
      MIRROR = MIRRIN
      BACKW = BACKIN
      RETURN

      ENTRY TRANSR(
     >             MIRROU,BACKOU)
      MIRROU=MIRROR
      BACKOU=BACKW
      RETURN

      ENTRY TRANS2(CONSTI)
      CONSTY = CONSTI
      RETURN
      END
