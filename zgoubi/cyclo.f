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
C  USA
C  -------
      SUBROUTINE CYCLO(SCAL,
     >                      DSREF,IRD,IDB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C---------------------------------------------------------
C     CYCLOTRON sector with several dipoles.
C     Each dipole has entrance and exit field boundaries,
C     possibly a third lateral one, much like DIPOLES.
C     Up to NMAG dipoles.
C     The total angular extent of the field region is given
C     by AT, each dipole is positionned within AT by  ACENT
C     angle and by its RM radial positionniing.
C     For each dipole, the three faces are positionned wrt.
C     the dipole's ACENT value.
C---------------------------------------------------------
      INCLUDE "C.AIM_2.H"     ! COMMON/AIM/ AE,AT,AS,RM,XI,XF,EN,EB1,EB2,EG1,EG2
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL,PI,RAD,DEG,QEL,AMPROT,CM2M
      INCLUDE "C.CONST2.H"     ! COMMON/CONST2/ ZERO, UN
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "C.DROITE.H"     ! COMMON/DROITE/ CA(9),SA(9),CM(9),IDRT
      INCLUDE "C.REBELO.H"   ! COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      INCLUDE "C.SPIRALE.H"     ! COMMON/spiral_ent/UMEG,ASP0,ASP1,ASP2,ASP3
      INCLUDE "C.SPIRALX.H"     ! COMMON/spiral_ext/UMEGs,ASPS0,ASPS1,ASPS2,ASPS3
      INCLUDE "C.RADIALS.H"     ! COMMON/radial_sec/aen,ben,cen,aex,bex,cex

      DIMENSION FTAB(5,5)

      PARAMETER (NM=5)
      DIMENSION NBFACE(NM)
      DIMENSION CE(NM,8),CS(NM,8),CC(NM,6)
      DIMENSION ACN(NM),DRM(NM),RRM(NM),HNORM(NM),Rref(NM)
      DIMENSION H1(NM),H2(NM),H3(NM),H4(NM), FAC(NM)
      DIMENSION COEFN(NM)
      DOUBLE PRECISION LAMBDE, LAMBDS, LAMBD3
      DOUBLE PRECISION NORMS
      DIMENSION LAMBDE(NM),QAPPAE(NM),NCOEFE(NM),RNME(NM)
      DOUBLE PRECISION G10,G11,G20,G21
      DIMENSION G10(NM),G11(NM),G20(NM),G21(NM)
      DIMENSION LAMBDS(NM),QAPPAS(NM),NCOEFS(NM),NORMS(NM)
      DIMENSION LAMBD3(NM),QAPPAL(NM),NCOEF3(NM),SHIFT3(NM),RM3(NM)
      DIMENSION UMEGA(NM),ASPIE0(NM),ASPIE1(NM),ASPIE2(NM),ASPIE3(NM)
      DIMENSION UMEGAS(NM),ASPIS0(NM),ASPIS1(NM),ASPIS2(NM),ASPIS3(NM)
      DIMENSION UMEGA3(NM),ASPI3(NM),R1L(NM),U13(NM),U23(NM),R2L(NM)
      DIMENSION aen1(NM), ben1(NM), cen1(NM), aex1(NM),bex1(NM),cex1(NM)

      SAVE ACN,DRM,HNORM,Rref,CE,CS,CC,H1,H2,H3,H4
      SAVE COEFN, FAC
      SAVE LAMBDE,QAPPAE,NCOEFE,RNME
      SAVE LAMBDS,QAPPAS,NCOEFS,NORMS
      SAVE LAMBD3,QAPPAL,NCOEF3,SHIFT3
      SAVE UMEGA,ASPIE0,ASPIE1,ASPIE2,ASPIE3
      SAVE UMEGAS,ASPIS0,ASPIS1,ASPIS2,ASPIS3
      SAVE NN, RESOL
      SAVE NBMAG, NBFACE
C      SAVE RO1,RO2,B1,B2,AMIN,AMAX,TTA1,TTA2
C FM Jan 2015
      SAVE RO1,RO2,AMIN,AMAX,TTA1,TTA2


C FM Jan 2015
C      DIMENSION BZ0(5,5)

      CHARACTER(14) TYPCAL(2)
C      CHARACTER(12) TYPGAP(2)
C FM Jan 2015
      CHARACTER(15) TYPGAP(2)
      SAVE TYPCAL, TYPGAP

      PARAMETER (PLIM=40.D0)

      LOGICAL SHARPE, SHARPS

      DATA TYPCAL / ' analytic', ' interpolation'/
      DATA TYPGAP / ' constant', ' g_0(1-r**2)**k' /

c      parameter(lunW=12)

c      open(unit=lunW,file='Distance.H')


C  NBMAG=number of magnets.  AT=total extent angle of field
      NP = 2
      NBMAG = NINT(A(NOEL,NP))
      NP=NP+1
      AT    = A(NOEL,NP)
      NP=NP+1
C The rm value in the spiral equation:
      RM    = A(NOEL,NP)
      NP=NP+1
      Typ   = A(NOEL,NP)
      IF(NRES.GT.0) WRITE(NRES,200) NBMAG,AT,RM
  200 FORMAT(20X,'Cyclotron N-tuple,  number  of  dipoles  N : ',I2,//,
     1 11X,' Total angular extent of the magnet : ',F6.2,' Degres',/,
     2 11X,' Reference geometrical radius R0  : ',F10.2,' cm',/)

C  HNORM=Champ MAX DANS LE DIPOLE.
C  COEFN=N=INDICE DE Champ, B=N', G=N''.
C  Rref= rayon de reference utilise uniquement dans la loi de champ: B(r)~ [1-(r/Rref)**2]**(-0.5)
C  AT=ANGLE TOTAL DE LA zone DE Champ, ACENT='ANGLE AU CENTRE',
C    RM=RAYON MOYEN DE LA zone DE Champ.
C  NBFACE(KMAG)=(2)3 : dipole limited by (2)3 field boundaries

C-----  a list of NBMAG.LE.NM magnets
      KMAG = 0
 10   CONTINUE
      KMAG = KMAG+1

      NP=NP+1
      ACENT = A(NOEL,NP)
      ACN(KMAG) = ACENT * RAD
      NP=NP+1
      DRM(KMAG)    = A(NOEL,NP)
      NP=NP+1
      FAC(KMAG) = A(NOEL,NP)
      NP=NP+1
      HNORM(KMAG) = A(NOEL,NP)*SCAL
      NP=NP+1
      COEFN(KMAG) = A(NOEL,NP)
      NP=NP+1
      Rref(KMAG)  = A(NOEL,NP)
      NP=NP+1
      H1(KMAG) = A(NOEL,NP)*SCAL
      NP=NP+1
      H2(KMAG) = A(NOEL,NP)*SCAL
      NP=NP+1
      H3(KMAG) = A(NOEL,NP)*SCAL
      NP=NP+1
      H4(KMAG) = A(NOEL,NP)*SCAL

      IF(NRES.GT.0) WRITE(NRES,100) KMAG,ACN(KMAG)/RAD,DRM(KMAG),
     >       HNORM(KMAG),COEFN(KMAG)
 100  FORMAT(5X,'Dipole # ',I1,/,1P,
     > 11X,' Positionning  angle ACENT : ',G12.4,' degrees',/,
     > 11X,' Positionning  wrt.  R0  : ',G10.2,' cm',/,
     5 11X,' B0 =',G12.4,' kGauss,',7X,'K =',G13.5)

      NP=NP+1
      LAMBDE(KMAG) = A(NOEL,NP)
      NP=NP+1
      QAPPAE(KMAG) = A(NOEL,NP)
      NP=NP+1
      G10(KMAG) = A(NOEL,NP)
      NP=NP+1
      G11(KMAG) = A(NOEL,NP)

      SHARPE=LAMBDE(KMAG) .LE. 0.D0
      IF(SHARPE) THEN
        GPE = -LAMBDE(KMAG)
        IF(NRES.GT.0)
     >  WRITE(NRES,FMT='(''Entrance hard edge is to be implemented'')')
        CALL INTEG1(ZERO,ZERO,GPE)
      ENDIF
      NP=NP+1
      NCOEFE(KMAG) = NINT(A(NOEL,NP))
      DO 227 I=1,8
        NP=NP+1
 227    CE(KMAG,I) = A(NOEL,NP)
      NP=NP+1
      RNME(KMAG) = A(NOEL,NP)
!      write(*,*) RNME(KMAG)  ! malek
      NP=NP+1
      UMEGA(KMAG)  = A(NOEL,NP)
      NP=NP+1
      ASPIE0(KMAG) = A(NOEL,NP)
      NP=NP+1
      ASPIE1(KMAG) = A(NOEL,NP)
      NP=NP+1
      ASPIE2(KMAG) = A(NOEL,NP)
      NP=NP+1
      ASPIE3(KMAG) = A(NOEL,NP)

      NP=NP+1
      aen1(KMAG) = A(NOEL,NP)
      NP=NP+1
      ben1(KMAG) = A(NOEL,NP)
      NP=NP+1
      cen1(KMAG) = A(NOEL,NP)

      IF(NRES.GT.0) THEN
        KGAP = 2
        IF(QAPPAE(KMAG).EQ.0) KGAP=1
        WRITE(NRES,106)'Entrance',LAMBDE(KMAG),TYPGAP(KGAP),QAPPAE(KMAG)
106     FORMAT (/,5X,A8,'  EFB',/,10X,'Fringe  field  :  gap at R0 is',
     >                                 F7.2,' cm,    type is : ',A,F5.2)
        WRITE(NRES,127) NCOEFE(KMAG),(CE(KMAG,I),I=1,6),RNME(KMAG)
127     FORMAT (10X,' COEFFICIENTS :',I3,6F10.5
     2  ,/,10X,' Shift  of  EFB  is ',G12.4,' cm',/)
        WRITE(NRES,103) UMEGA(KMAG),ASPIE0(KMAG)
C 103    FORMAT(10X,7HOMEGA =,F7.2,5X,17HANGLE  DE  FACE =,F7.2,/ ,
 103    FORMAT(10X,'OMEGA =',F7.2,' deg.',
     >   5X,'Spiral  angle  =',F7.2,' deg.')
      ENDIF

C Exit Fringe Field
      NP=NP+1
      LAMBDS(KMAG) = A(NOEL,NP)
      NP=NP+1
      QAPPAS(KMAG)   = A(NOEL,NP)
      NP=NP+1
      G20(KMAG) = A(NOEL,NP)
      NP=NP+1
      G21(KMAG) = A(NOEL,NP)
      SHARPS=LAMBDS(KMAG) .LE. 0.D0
      IF(SHARPS) THEN
        GPS = -LAMBDS(KMAG)
        IF(NRES.GT.0)
     >  WRITE(NRES,FMT='(''Exit hard edge is to be implemented'')')
        CALL INTEG2(ZERO,ZERO,GPS)
      ENDIF
      NP=NP+1
      NCOEFS(KMAG) = NINT(A(NOEL,NP))
      DO 228 I=1,8
        NP=NP+1
 228    CS(KMAG,I) = A(NOEL,NP)
      NP=NP+1
      NORMS(KMAG) = A(NOEL,NP)

      NP=NP+1
      UMEGAS(KMAG) = A(NOEL,NP)
      NP=NP+1
      ASPIS0(KMAG) = A(NOEL,NP)
      NP=NP+1
      ASPIS1(KMAG) = A(NOEL,NP)
      NP=NP+1
      ASPIS2(KMAG) = A(NOEL,NP)
      NP=NP+1
      ASPIS3(KMAG) = A(NOEL,NP)

      NP=NP+1
      aex1(KMAG) = A(NOEL,NP)
      NP=NP+1
      bex1(KMAG) = A(NOEL,NP)
      NP=NP+1
      cex1(KMAG) = A(NOEL,NP)

      IF(NRES.GT.0) THEN
        KGAP = 2
        IF(QAPPAS(KMAG).EQ.0) KGAP=1
        WRITE(NRES,106) 'Exit',LAMBDS(KMAG),TYPGAP(KGAP),QAPPAS(KMAG)
        WRITE(NRES,127) NCOEFS(KMAG),(CS(KMAG,I),I=1,6),NORMS(KMAG)
        WRITE(NRES,103) UMEGAS(KMAG),ASPIS0(KMAG)
      ENDIF

      NP=NP+1
      LAMBD3(KMAG) = A(NOEL,NP)
      NP=NP+1
      QAPPAL(KMAG)   = A(NOEL,NP)
      NP=NP+1
      NCOEF3(KMAG) = NINT(A(NOEL,NP))
      DO 229 I=1,6
        NP=NP+1
 229    CC(KMAG,I) = A(NOEL,49+I)
      NP=NP+1
      SHIFT3(KMAG) = A(NOEL,NP)

       NP=NP+1
       UMEGA3(KMAG) = A(NOEL,NP)
       NP=NP+1
       ASPI3(KMAG) = A(NOEL,NP)
       NP=NP+1
       R1L(KMAG)    = A(NOEL,NP)
       NP=NP+1
       U13(KMAG)    = A(NOEL,NP)
       NP=NP+1
       U23(KMAG)    = A(NOEL,NP)
       NP=NP+1
       R2L(KMAG)    = A(NOEL,NP)

       NBFACE(KMAG)=2
       IF(QAPPAL(KMAG) .GE. 0) NBFACE(KMAG)=3

       IF(NRES.GT.0) THEN
         IF(QAPPAL(KMAG).LT.0) THEN
            WRITE(NRES,*)
            WRITE(NRES,*) '        Lateral face :  unused'
            WRITE(NRES,*)
         ELSE
           STOP 'Lateral EFB is not implemented. Use Entrance/Exit only'
           KGAP = 2
           IF(QAPPAL(KMAG).EQ.0) KGAP=1
           WRITE(NRES,106) 'Lateral ',LAMBD3(KMAG),TYPGAP(KGAP),
     >                                                    QAPPAL(KMAG)
           WRITE(NRES,127) NCOEF3(KMAG),(CC(KMAG,I),I=1,6),SHIFT3(KMAG)
           WRITE(NRES,113) RM3(KMAG)
113        FORMAT(20X,' EFB is centred on direction ACENT+OMEGA,  at ',
     >     G12.4,' cm'/)
           WRITE(NRES,103)UMEGA3(KMAG),ASPI3(KMAG),R1L(KMAG),U13(KMAG),
     >        U23(KMAG),R2L(KMAG)
         ENDIF
       ENDIF

      IF(KMAG.LT.NBMAG) GOTO 10

C-----------------------------

C Get type of field & deriv. calculation
      NP=NP+1
      KIRD = NINT(A(NOEL,NP))
C Get resol, or idb
      NP=NP+1
      RESOL=A(NOEL,NP)
      IF    (KIRD.NE.0) THEN
C    interpolation
        IF(SHARPE .OR. SHARPS) CALL ENDJOB
     >    ('ERROR :  sharp edge not compatible with num. deriv.',-99)
        IRD = KIRD
        KIRD = 1
        IF    (IRD.EQ.2) THEN
          NN=3
        ELSEIF(IRD.EQ.25) THEN
          NN=5
        ELSEIF(IRD.EQ.4) THEN
          NN=5
        ELSE
          STOP ' *** ERROR - SBR FFAGI, WRONG VALUE IRD'
        ENDIF
      ELSEIF(KIRD.EQ.0) THEN
C    analytic
        IDB=NINT(RESOL)
        IF(IDB.NE.4) IDB=2
      ENDIF
      CALL CHAMC6(KIRD)

      AE=0.D0
      AS=0.D0
      AT = AT * RAD

      AMIN = -AT/2.D0 !- ACN(KMAG)
      AMAX =  AT/2.D0 !- ACN(KMAG)

c          write(*,*) ' sbr ffgspi,   amin,amax :', amin,amax
c                 pause

      XI = AMIN
      XF = AMAX

C--- Formule à revoir...
      DSREF = RM * (UMEGA(KMAG)-UMEGAS(KMAG))
C--------------------

      IF(NRES.GT.0) THEN
        WRITE(NRES,FMT='(/,5X,'' Field & deriv. calculation :'',A)')
     >  TYPCAL(KIRD+1)
        IF    (KIRD.NE.0) THEN
          IF(IRD .EQ. 2) THEN
            WRITE(NRES,121) '3*3', RESOL
 121        FORMAT(20X,A3,
     >        '-point  interpolation, size of flying mesh :  ',
     >        'integration step /',1P,G12.4)
          ELSEIF(IRD .GE. 4) THEN
C----------- IRD= 4 OR 25
            WRITE(NRES,121) '5*5', RESOL
          ENDIF
        ELSEIF(KIRD.EQ.0) THEN
          WRITE(NRES,FMT='(20X,
     >    ''Derivatives computed to order '',I1)') IDB
        ENDIF
      ENDIF

      RETURN

C------------------------------------------------------------
C  Compute CYCLOTRON field  and derivatives from flying field-mesh
      ENTRY CYCLOF(TTA,RO,
     >                   DTTA,DRO,FTAB)

C           write(*,*) ' tta, ro  ',tta, ro,acn(kmag),tta-amin,amin
      CALL INTEG5(
     >            STEP)

      CALL RAZ(FTAB,5*5)
      DRO = STEP/RESOL
      DTTA = DRO/RM

      KMAG = 0
 20   CONTINUE
      KMAG = KMAG+1
      RRM(KMAG) = RM +DRM(KMAG)

C------------------------------------------------------------

C Entrance EFB

      UMEG = UMEGA(KMAG) * RAD
      ASP0 = ASPIE0(KMAG) * RAD
      ASP1 = ASPIE1(KMAG) * RAD
      ASP2 = ASPIE2(KMAG) * RAD
      ASP3 = ASPIE3(KMAG) * RAD
C Exit EFB
      UMEGS = UMEGAS(KMAG) * RAD
      ASPS0 = ASPIS0(KMAG) * RAD
      ASPS1 = ASPIS1(KMAG) * RAD
      ASPS2 = ASPIS2(KMAG) * RAD
      ASPS3 = ASPIS3(KMAG) * RAD


C Entrance and exit EFBs DO NOT necessarily have the same parameters

C------------------------------------------------------------

C  Entrance and Exit radial face equation parameters: ax+by+c=0

      aen=aen1(KMAG)
      ben=ben1(KMAG)
      cen=cen1(KMAG)

      aex=aex1(KMAG)
      bex=bex1(KMAG)
      cex=cex1(KMAG)



C------------------------------------------------------------

C----- CALCUL LE Champ en X, Y
C      COORDONNEES DU POINT COURANT

      XACC=0.0001D0     ! accuracy for newton zero method

      RO1=RM

CCCCCCCCCCCCCCCCCCCCCCCCCCCC  MALEK  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IF(ASP .EQ. 0.D0) GOTO 13

C      ELSE
C      B1=1.D0/TAN(ASP0)
C      ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      RO2=RM


      AMIN =  -AT/2.D0 !- ACN(KMAG)
      AMAX =   AT/2.D0 !- ACN(KMAG)

c      TTA1 = UMEG
c      TTA2 = -UMEG
      TTA1 = acn(kmag) + UMEG
      TTA2 = acn(kmag) + UMEGS

      DO  1  JRO = 1,NN
        ROJ = RO + DRO * DBLE(NN-JRO-INT(NN/2))
!        tanxi1=TAN(ASP0+ASP1*ROJ+ASP2*ROJ**2+ASP3*ROJ**3)
!        B1= 1.D0/tanxi1
!        tanxi2=TAN(ASPS0+ASPS1*ROJ+ASPS2*ROJ**2+ASPS3*ROJ**3)
!        B2= 1.D0/tanxi2

        DO  1  ITTA = 1,NN
          TTAI = TTA - DTTA * DBLE(NN-ITTA-INT(NN/2))
          X = ROJ * COS(TTAI)
          Y = ROJ * SIN(TTAI)

C CALCUL DE LA DISTANCE DE (X,Y) A LA FACE D'ENTREE
          FA=1.0
          D=-DSTEFB2(X,Y,RO1,AMIN,AMAX,XACC,TTA1,YN,FA,Typ)


C          write(*,*) B1
C            IF( Y .GT. YN .OR. D .LE. 1.D-6) D = -D
            IF( Y .GT. YN ) D = -D
C            D=( D + SHIFTE(KMAG) )


C CALCUL DE FE

          IF(LAMBDE(KMAG) .EQ. 0.D0) THEN
            IF(D.LE.0.D0) THEN
              FE=1.D0
            ELSE
              FE=0.D0
            ENDIF
          ELSE
            IF(QAPPAE(KMAG).EQ.0) THEN
              GAP = LAMBDE(KMAG)
            ELSE
         GAP = G10(KMAG)+G11(KMAG)*ROJ/100.0

!        GAP = QAPPAE(KMAG)/((HNORM(KMAG)+H1(KMAG)*ROJ+H2(KMAG)*(ROJ)**2+
!     >         H3(KMAG)*(ROJ)**3+H4(KMAG)*(ROJ)**4))  !* (1-(ROJ/Rref(KMAG))**2)**QAPPAE(KMAG)          !  loi du gap pour satisfaire l'isochronisme
            ENDIF
            DNTR=D
            D = D/GAP

            Dent=D




            P=CE(KMAG,1)+(CE(KMAG,2)+(CE(KMAG,3)+(CE(KMAG,4)+
     >              (CE(KMAG,5)+(CE(KMAG,6)+CE(KMAG,7)*D)*D)*D)*D)*D)*D

            PNTR=P
            IF    (P .GE.  PLIM) THEN
              FE = 0.D0
            ELSEIF(P .LE. -PLIM) THEN
              FE = 1.D0
            ELSE
              FE = (1.D0/(1.D0+EXP(P)))
            ENDIF
          ENDIF

C CALCUL DE LA DISTANCE DE (X,Y) A LA FACE DE SORTIE
         FA=2.0
         D=DSTEFB2(X,Y,RO2,AMIN,AMAX,XACC,TTA2,YN,FA,Typ)

         IF( Y .GT. YN ) D = -D




C CALCUL DE FS

          IF(LAMBDS(KMAG) .EQ. 0.D0) THEN
            IF(D.LE.0.D0) THEN
              FS=1.D0
            ELSE
              FS=0.D0
            ENDIF
          ELSE
            IF(QAPPAS(KMAG).EQ.0) THEN
              GAP = LAMBDS(KMAG)
            ELSE
         GAP = G20(KMAG)+G21(KMAG)*ROJ/100.0

!       GAP = QAPPAS(KMAG)/((HNORM(KMAG)+H1(KMAG)*ROJ+H2(KMAG)*(ROJ)**2+
 !    >         H3(KMAG)*(ROJ)**3+H4(KMAG)*(ROJ)**4))  !* (1-(ROJ/Rref(KMAG))**2)**0.5
            ENDIF


            D = D/GAP

            Dexit=D




C            GAP = LAMBDS(KMAG) * (1-(ROJ/Rref(KMAG))**2)**0.5
            P=CS(KMAG,1)+(CS(KMAG,2)+(CS(KMAG,3)+(CS(KMAG,4)+
     >              (CS(KMAG,5)+(CS(KMAG,6)+CS(KMAG,7)*D)*D)*D)*D)*D)*D

            IF    (P .GE.  PLIM) THEN
              FS = 0.D0
            ELSEIF(P .LE. -PLIM) THEN
              FS = 1.D0
            ELSE
              FS = (1.D0/(1.D0+EXP(P)))
            ENDIF
          ENDIF


c          A1=-0.00584785
c          A2=1.13695053


          F = FE + FS -1
C             Calcul du champ B au point (ITTA,JRO) de la grille volante
          FTAB(ITTA,JRO)=FTAB(ITTA,JRO) + FAC(KMAG)*
     >                 F*(HNORM(KMAG)+H1(KMAG)*ROJ+H2(KMAG)*(ROJ)**2+
     >         H3(KMAG)*(ROJ)**3+H4(KMAG)*(ROJ)**4)


cccccccccccccccccccc   write the distance to a file  malek   ccccccccccccccccccccccccc
c          if ((JRO .EQ. NN) .and. (ITTA .EQ. NN) .and.
c     >         (KMAG .EQ. 1))  THEN
c        write(lunW,*) TTAI,Dent,Dexit, FTAB(ITTA,JRO), FE,FS,gap,F
c          endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    1 CONTINUE
C--------- end loop on   NN x tta  &  NN x ro

      IF(KMAG.LT.NBMAG) GOTO 20
C-----


c      CLOSE(lunW)

C-----------------------------------------------------------
C  Compute FFAG field  and derivatives from analytical model

      RETURN

      END
