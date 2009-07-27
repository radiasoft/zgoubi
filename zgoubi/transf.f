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
      SUBROUTINE TRANSF
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ----------------------------------------------------------------
C     APPELE PAR AIMANT ET QUASEX .
C     APPELLE INTEGR POUR LE CALCUL DES IMAX TRAJECTOIRES DU FAISCEAU.
C     ----------------------------------------------------------------
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      PARAMETER(MCOEF=6)
      COMMON/CHAFUI/ XE,XS,CE(MCOEF),CS(MCOEF),QCE(MCOEF),QCS(MCOEF)
C      COMMON/CHAFUI/ XE,XS,CE(6),CS(6),QCE(6),QCS(6)
      INCLUDE "MAXTRA.H"
      COMMON/CHAMBR/ LIMIT,IFORM,YL2,ZL2,SORT(MXT),FMAG,BMAX
     > ,YCH,ZCH
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      COMMON/CONST2/ ZERO, UN
      INCLUDE 'MXSTEP.H'
      INCLUDE 'CSR.H'
C      COMMON/CSR/ KTRA,KCSR,YZXB(MXSTEP,41,36),DWC(MXT)
      COMMON/DESIN/ FDES(7,MXT),IFDES,KINFO,IRSAR,IRTET,IRPHI,NDES
     >,AMS,AMP,AM3,TDVM,TETPHI(2,MXT)
C     1,AMS  ,AMP,ENSTAR,BSTAR,TDVM,TETPHI(2,MXT)
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      COMMON/DROITE/ AM(9),BM(9),CM(9),IDRT
      INCLUDE "MAXCOO.H"
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),IMAX,IEX(MXT),IREP(MXT)
      CHARACTER LET
      COMMON/FAISCT/ LET(MXT)
      COMMON/GASC/ AI, DEN, KGA
      COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      COMMON/MARK/ KART,KALC,KERK,KUASEX
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      LOGICAL ZSYM
      COMMON/OPTION/ KFLD,MG,LC,ML,ZSYM
      COMMON/PTICUL/ AAM,Q,G,TO
      COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)
      INCLUDE 'MXFS.H'
      CHARACTER FAM*8,LBF*8,KLEY*10,LABEL*8
      COMMON/SCALT/ FAM(MXF),LBF(MXF,2),KLEY,LABEL(MXL,2)
      COMMON/SYNRA/ KSYN
      COMMON/TRAJ/ Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT
      COMMON/TRNSF/ XFE,XFS

      LOGICAL CHGRFE,CHGRFS,CHPFE,CHPFS,EVNT,BACKW,DRT
     
      PARAMETER(MPOL=10)
      PARAMETER(I0=0)
      
      LOGICAL MIRROR, MIRRIN, BACKIN, MIRROU, BACKOU

      SAVE MIRROR, BACKW
        
      DATA MIRROR / .FALSE. /

      IF(LST .GE. 1 ) CALL CTRLB(1)

C----- Change ref. at entrance or exit 
C-----   (Including Automatic positionning, i.e. case KPOS=3 : MULTIPOL with non 
C                                       zero B(1), BEND, etc.)
C Tests, Jan 06, FM
C      CHGRFE= ( KP .GE. 1 .AND. XCE*XCE + YCE*YCE + ALE*ALE .GT. ZERO )
C      CHGRFS= ( KP .GE. 1 .AND. XCS*XCS + YCS*YCS + ALS*ALS .GT. ZERO )
      CHGRFE= ( KP .GE. 2 .AND. XCE*XCE + YCE*YCE + ALE*ALE .GT. ZERO )
      CHGRFS= ( KP .GE. 2 .AND. XCS*XCS + YCS*YCS + ALS*ALS .GT. ZERO )
 
C----- Champ DE FUITE ENTREE OU SORTIE
      CHPFE = ( KART .EQ. 1 .AND. XE .NE. ZERO )
      CHPFS = ( KART .EQ. 1 .AND. XS .NE. ZERO )
 
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
 
        IF(IT .EQ. IREP(IT) .OR. .NOT. ZSYM) THEN

          CALL INITRA(IT)
          X=XO

          IF( CHGRFE ) CALL CHAREF(EVNT,XCE,YCE,ALE)
          IF( CHPFE  ) CALL CHAREF(EVNT,XFE,ZERO,ZERO)

          IF( DRT  ) CALL DRTENT
 
          CALL INTEGR(EVNT,BACKW,MIRROR)
 
          IF( CHPFS  ) CALL CHAREF(EVNT,XFS,ZERO,ZERO)
          IF( CHGRFS ) THEN
            IF(NRES.GT.0)
     >        WRITE(NRES,100)LET(IT),KEX,(FO(J,IT),J=1,5),X ,Y,T,Z,P,IT
 100        FORMAT(2X,A1,2X,I2,F8.4,4F10.3,8X,F9.3,4F10.3,8X,I4)
            CALL CHAREF(EVNT,XCS,YCS,ALS)
          ENDIF

          CALL MAJTRA(IT)

        ELSE

          CALL DEJACA(IT)

        ENDIF
 
 1    CONTINUE
 
      IF(LST .GE. 1 ) CALL CTRLB(2)
 
      MIRROR = .FALSE.
      BACKW = .FALSE.

C----- Unset coded step
      CALL CHXC1W(I0,I0)
      CALL DEPLAW(.FALSE.,I0)

      RETURN

      ENTRY TRANSW(MIRRIN,BACKIN)
      MIRROR = MIRRIN
      BACKW = BACKIN
      RETURN

      ENTRY TRANSR(MIRROU,BACKOU)
      MIRROU=MIRROR
      BACKOU=BACKW
      RETURN
      END
