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
      SUBROUTINE MULTPO(KUASEX,LMNT,KFL,SCAL,
     >          DEV,RT,XL,BM,DLE,DLS,DE,DS,XE,XS,CE,CS,BORO,DPREF,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(MPOL=10)
      CHARACTER(*) LMNT(*)
      DIMENSION RT(*),BM(*),DLE(*),DLS(*),DE(MPOL,*),DS(MPOL,*)
      PARAMETER(MCOEF=6)
      DIMENSION CE(MCOEF), CS(MCOEF)
 
      COMMON/AIM/ BO,RO,FG,GF,XI,XF,EN,EB1,EB2,EG1,EG2
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      COMMON/CHAMBR/ LIMIT,IFORM,YLIM2,ZLIM2,SORT(MXT),FMAG,BMAX
     > ,YCH,ZCH
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QEL,AMPROT, CM2M
      COMMON/CONST2/ ZERO, UN
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER(80) TA
      PARAMETER (MXTA=45)
      COMMON/DONT/ TA(MXL,MXTA)
      COMMON/DROITE/ CA(9),SA(9),CM(9),IDRT
      COMMON/EFBS/ AFB(2), BFB(2), CFB(2), IFB
      COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      PARAMETER (LBLSIZ=10)
      CHARACTER(LBLSIZ) LABEL
      COMMON /LABEL/ LABEL(MXL,2)
      LOGICAL ZSYM
      COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      COMMON/PTICUL/ AM,Q,G,TO
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      INCLUDE 'MXFS.H'
      INCLUDE 'MXSCL.H'
      COMMON/SCAL/ SCL(MXF,MXS,MXSCL),TIM(MXF,MXS),NTIM(MXF),KSCL
C      COMMON/SCAL/SCL(MXF,MXS),TIM(MXF,MXS),NTIM(MXF),JPA(MXF,MXP),KSCL
      COMMON/SYNRA/ KSYN

C----------- MIXFF = true if combined sharp edge multpole + fringe field multpole
      LOGICAL SKEW, MIXFF
     
C----- FM, Fermilab, 1996, For special simulation of b10 in LHC low-beta quads
      LOGICAL CASPI 

      CHARACTER(3) DIM(2)
      CHARACTER(2) BE(2)
      CHARACTER(80) TXT(10)

      DIMENSION  AREG(2),BREG(2),CREG(2)

      DATA DIM / 'kG ', 'V/m'/
      DATA BE / 'B-', 'E-'/
 
      DATA CASPI / .TRUE. /
      
      dimension AKS(6)
      SAVE AKS, sXL

      PARAMETER (MXERR=MXTA)
      CHARACTER(LBLSIZ) LBL1(MXERR), LBL2(MXERR)
      CHARACTER(LBLSIZ) LBL1i, LBL2i
      save lbl1, lbl2

      DIMENSION kpol(mxerr,mpol)
      character(2) TYPERR(mxerr,mpol)
      character(1) TYPAR(mxerr,mpol),TYPDIS(mxerr,mpol)
      character(2) TYPERI
      character(1) TYPAI,TYPDII
      DIMENSION ERRCEN(mxerr,mpol),ERRSIG(mxerr,mpol),ERRCUT(mxerr,mpol)
      SAVE TYPERR,TYPAR,TYPDIS,ERRCEN,ERRSIG,ERRCUT
      LOGICAL empty
      LOGICAL erron
      SAVE ERRON
C      Errors
      dimension db(MXL,mpol),dpos(MXL,mpol,3),tilt(MXL,mpol,3)
      save db, dpos, tilt
      logical ok
      logical fiting

c      data db / mpol*0.d0 /
c      data dpos / mpol*0.d0, mpol*0.d0, mpol*0.d0 /
c      data dtilt / mpol*0.d0, mpol*0.d0, mpol*0.d0 /
      data erron / .false. /

      IER = 0
      SKEW=.FALSE.
      CALL RAZ(BM,MPOL)

      GAP = 0.D0
      IF(KUASEX .LE. MPOL) THEN
C------- Single-pole, from QUAD (KUASEX=2) up to 20-POLE (KUASEX=10)

        XL =A(NOEL,10)
        RO =A(NOEL,11)
        GAP = RO/KUASEX
        BM(KUASEX) =A(NOEL,12)*SCAL
        XE =A(NOEL,20)
        DLE(KUASEX) =A(NOEL,21)

        CALL RAZ(CE,MCOEF)
        NCE = NINT(A(NOEL,30))
        DO 22 I=1, NCE
          CE(I) =A(NOEL,30+I)
 22     CONTINUE

        XLS =A(NOEL,40)
        DLS(KUASEX) =A(NOEL,41)

C        IF(XE+XLS.GE.XL) 
C     >   CALL ENDJOB('SBR MULTIP : fringe field extent too long',-99)

        CALL RAZ(CS,MCOEF)
        NCS = NINT(A(NOEL,50))
        DO 24 I=1,NCS
          CS(I) =A(NOEL,50+I)
 24     CONTINUE

C------- Multipole rotation
        RT(KUASEX)=ZERO
 
        NM0 = KUASEX
        NM = KUASEX        

      ELSEIF(KUASEX .EQ. MPOL+1) THEN
C-------  Mag MULTIPOLE, Elec Multipole ELMULT, Elec & Mag Multipole EBMULT
 
        IF    (KFLD .EQ. MG .OR. KFL .EQ. LC) THEN
C--------- Mag MULTIPOLE or electr ELMULT or Electric part of EBMULT...
          IA = 2
        ELSEIF(KFL .EQ. MG) THEN
C--------- ... or magnetic part of EBMULT
          IA = 60
        ENDIF
 
        XL =A(NOEL,IA)
        IA = IA + 1
        RO =A(NOEL,IA)
        DO IM=1,MPOL
          IA = IA + 1
          BM(IM) =A(NOEL,IA)*SCAL
          IF(RO .EQ. 0.D0) THEN
            IF(IM .GE. 2) THEN 
              IF(BM(IM) .NE. 0.D0) THEN
                IER = IER+1
                TXT(IER) = 'SBR  MULTPO - RO must be non-zero '
                CALL ENDJOB(TXT(IER),-99)
              ENDIF
            ENDIF
          ENDIF
          IF(GAP.EQ.0.D0) THEN
            IF(BM(IM).NE.0.D0) GAP = RO/IM
          ENDIF
        ENDDO

C------- If SR-loss switched on by procedure SRLOSS
        IF(KSYN.GE.1) THEN
          IF(KFL .EQ. MG) THEN
            IF(BM(1).NE.0.D0) CALL SYNPAR(BM(1),XL)
          ENDIF
        ENDIF

        IA = IA + 1
        XE =A(NOEL,IA)
        DO 31 IM=1,MPOL
          IA = IA + 1
          DLE(IM) =A(NOEL,IA)
 31     CONTINUE
        IA = IA + 1
        CALL RAZ(CE,MCOEF)
        NCE = NINT(A(NOEL,IA))
C        DO 32 I=1,6
        DO 32 I=1, NCE
          IA = IA + 1
          CE(I) =A(NOEL,IA)
 32     CONTINUE

        IA = IA + MCOEF - NCE + 1
        XLS =A(NOEL,IA)
        DO 33 IM = 1,MPOL
          IA = IA + 1
          DLS(IM) =A(NOEL,IA)
 33     CONTINUE
        IA = IA + 1
        CALL RAZ(CS,MCOEF)
        NCS = NINT(A(NOEL,IA))
         DO 34 I=1,NCS
          IA = IA + 1
          CS(I) =A(NOEL,IA)
 34     CONTINUE
C------- FRINGE FIELD NOT INSTALLED FOR
C        DECA, DODECA, ... 18-POLE
        DLE(5)=ZERO
        DLS(5)=ZERO
        DLE(6)=ZERO
        DLS(6)=ZERO
        DLE(7)=ZERO
        DLS(7)=ZERO
        DLE(8)=ZERO
        DLS(8)=ZERO
        DLE(9)=ZERO
        DLS(9)=ZERO
  
C----- Multipole rotation
        IA = IA + MCOEF - NCS 
        DO 35 IM=1,MPOL
          IA = IA + 1
          RT(IM)=A(NOEL,IA) 
          SKEW=SKEW .OR. RT(IM) .NE. ZERO
 35     CONTINUE
 
        NM0 = 1
        NM = MPOL
 
C------- FM, LHC purpose, Fermilab, 1996
        IF(NCE .EQ. 999 .OR .NCS .EQ. 999) CALL MULTI1(CASPI)

      ENDIF
C------------ KUASEX
 
C----- Magnetic MULTIPOL with non-zero dipole
C--------- Automatic positionning in SBR TRANSF
C           Dev normally identifies with the deviation
C             that would occur in a sharp edge dipole magnet.
        IF(BM(1) .NE. 0.D0) THEN
          DEV= 2.D0* ASIN(.5D0*XL*BM(1)/BORO)
        ELSE
          DEV= 0.D0
        ENDIF
C------------------------  TESTS for COSY
C         DEV = 2.D0 * PI /24.D0  
CC         DEV = 2.D0 * ASIN(.5D0 * XL * 14.32633183D0 / BORO )
C         DEV = ALE
C-----------------------------------------

C----- Case erron (errors)
      if(erron) then
        do irr = 1, mxerr 
          ok = (empty(lbl1(irr)) .or. lbl1(irr).eq.label(noel,1)) 
     >    .and.(empty(lbl2(irr)) .or. lbl2(irr).eq.label(noel,2)) 
c                 write(*,*) ' multpo ok err ',ok,
c     >         empty(lbl1(irr)), lbl1(irr).eq.label(noel,1),  
c     >      empty(lbl2(irr)),  lbl2(irr).eq.label(noel,2)
c                 write(*,*) '               ',
c     >         lbl1(irr), label(noel,1),lbl2(irr),label(noel,2)
c                      read(*,*)
          if(ok) then
            if(ipass.eq.1) then 
              CALL FITSTA(5,FITING)
              if(.not.fiting) 
     >        call mulerr(noel,irr,BM, 
     >        KPOL,TYPERR,TYPAR,TYPDIS,ERRCEN,ERRSIG,ERRCUT,
     >                                     DB,dpos,tilt)
c           if(noel.le.20)write(*,*) 
c     >     ' sbr multpo. fiting, db : ',fiting,noel, db(noel,1)
c           if(noel.le.20)write(66,*) 
c     >     ' sbr multpo. fiting, db : ',fiting,noel, db(noel,1)
c                      read(*,*)
            endif
            IF(KUASEX .LE. MPOL) THEN
              BM(KUASEX) = BM(KUASEX) + db(noel,kuasex)
            ELSEIF(KUASEX .EQ. MPOL+1) THEN
              DO IM=1,MPOL
                if(kpol(irr,im).eq.1) BM(IM) = BM(IM) + db(noel,im)
              enddo
            endif      
          endif      
        enddo
      endif

C----- MULTIPOLE
      DO IM = NM0+1,NM
        DLE(IM)  = DLE(NM0)*DLE(IM)
        DLS(IM)  = DLS(NM0)*DLS(IM)
      ENDDO

      SXL = XL

      AKS(1) =  BM(1)/1.D0/(BORO*DPREF)*1.D2
      AKS(2) =  BM(2)/(RO)/(BORO*DPREF)*1.D4
      AKS(3) =  BM(3)/(RO**2)/(BORO*DPREF)*1.D6

      IF(NRES.GT.0) THEN
        WRITE(NRES,100) LMNT(KUASEX),XL,RO
 100    FORMAT(/,5X,' -----  ',A10,'  : ', 1P
     >  ,/,15X,' Length  of  element  = ',G16.8,'  cm'
     >  ,/,15X,' Bore  radius      RO = ',G13.5,'  cm')
        WRITE(NRES,103) (BE(KFL),LMNT(IM),BM(IM),DIM(KFL),IM=NM0,NM)
 103    FORMAT(15X,2A,'  =',1P,E15.7,1X,A)
        IF(SKEW) WRITE(NRES,101) (LMNT(IM),RT(IM),IM=NM0,NM)
 101    FORMAT(15X,A,'  Skew  angle =',1P,E15.7,' rd')
        IF(XL .NE. 0.D0) THEN
          IF( (XL-DLE(NM)-DLS(NM)) .LT. 0.D0) WRITE(NRES,102)
 102      FORMAT(/,10X,'Entrance  &  exit  fringe  fields  overlap, ',
     >    /,10X,'  =>  computed  gradient  is ',' G = GE + GS - 1 ')
        ELSE
          GOTO 98
        ENDIF
      ENDIF
 
      DL0=0.D0
      SUM=0.D0

      DO IM=NM0,NM
        DL0=DL0+DLE(IM)+DLS(IM)
        SUM=SUM+BM(IM)*BM(IM)
C------- E converted to MeV/cm
        IF(KFL .EQ. LC) THEN
          IF(BM(IM).NE.0.D0) BM(IM) = 2.D0*BM(IM)/RO*1.D-6
          RT(IM) = RT(IM) + .5D0*PI/DBLE(IM)
        ENDIF
        IF(BM(IM).NE.0.D0) BM(IM) = BM(IM)/RO**(IM-1)
      ENDDO

      IF(SUM .EQ. 0.D0) KFLD=KFLD-KFL
      IF(DL0 .EQ. 0.D0) THEN
C-------- Sharp edge at entrance and exit
        FINTE = XE
        XE=0.D0
        FINTS = XLS
        XLS=0.D0
        IF(NRES.GT.0) THEN
          WRITE(NRES,105) 'Entrance/exit field models are '
 105      FORMAT(/,15X,A,'sharp edge')
          WRITE(NRES,FMT='(15X,''FINTE, FINTS, gap : '',
     >    1P,3(1X,E12.4))') FINTE,FINTS,GAP
        ENDIF
        IF(KFL .EQ. MG) THEN
          IF(KUASEX .EQ. MPOL+1) THEN
C------------- MULTIPOL case. Set entrance & exit wedge correction in SBR INTEGR
            IF(BM(1) .NE. 0.D0) THEN
C FM, Oct. 2011. Avoid wedge correction with short lens
              IF(XL .GT. 2.D0) THEN            
C FM, 2006
                CALL INTEG1(ZERO,FINTE,GAP)
                CALL INTEG2(ZERO,FINTS,GAP)
              ENDIF
            ENDIF
          ENDIF
        ENDIF
        
      ELSE
C-------- Gradient G(s) at entrance or exit
C-----    Let's see entrance first
        IF(NRES.GT.0) WRITE(NRES,104)
 104    FORMAT(/,15X,' Entrance  face  ')
C 104    FORMAT(/,15X,' FACE  D''ENTREE  ')
        DL0 = 0.D0
        DO 5 IM=NM0,NM
 5        DL0 = DL0+DLE(IM)

        IF(DL0 .EQ. 0.D0) THEN
          FINTE = XE
          XE=0.D0
          IF(NRES.GT.0) THEN
            WRITE(NRES,105) 'Entrance field model is '
            WRITE(NRES,FMT='(15X,''FINTE, gap : '',
     >      1P,2(1X,E12.4))') FINTE,GAP
          ENDIF
          IF(KFL .EQ. MG) THEN
            IF(KUASEX .EQ. MPOL+1) THEN
              IF(BM(1) .NE. 0.D0) THEN
C----------- Set entrance wedge correction in  SBR INTEGR
C FM, 2006
C FM, Oct. 2011. Avoid wedge correction with short lens
                IF(XL .GT. 2.D0) CALL INTEG1(ZERO,FINTE,GAP)
              ENDIF
            ENDIF
          ENDIF
        ELSE
C---------- IFB = 0 if no mixff
C----------     set to -1 if mixff at entrance
C----------     set to  1 if mixff at exit
C----------     set to  2 if mixff at entrance & exit
          IF( IFB .EQ. 0 .OR. IFB .EQ. 1 ) THEN
            MIXFF = .FALSE.
            DO 51 IM = NM0,NM
              IF(.NOT. MIXFF) THEN
C--------------- MIXFF = true if combined sharp edge multpole + fringe field multpole
                IF(DLE(IM) .EQ. 0.D0 .AND. BM(IM) .NE. 0.D0)
     >            MIXFF= .TRUE.
              ENDIF
 51         CONTINUE
            IF(MIXFF) THEN
              IF(IFB .EQ. 0) THEN
                IFB = -1
              ELSE
                IFB = 2
              ENDIF
            ENDIF
          ENDIF

          IF(NRES.GT.0) THEN
            WRITE(NRES,130) XE
 130        FORMAT(20X,' with  fringe  field :'
C 130        FORMAT(20X,' AVEC  Champ  DE  FUITE  :'
     >      ,/,20X,' DX  = ',F7.3,'  CM ')
            WRITE(NRES,131) ( LMNT(IM),DLE(IM) ,IM=NM0,NM)
 131        FORMAT(20X,' LAMBDA-',A,' =',F7.3,'  CM')
C            WRITE(NRES,132) (CE(I),I=1,6)
            WRITE(NRES,132) NCE, (CE(I),I=1,NCE)
 132        FORMAT(20X,I1,' COEFFICIENTS :',6F9.5)
          ENDIF
          DO 45 IM=NM0,NM
            IF(DLE(IM) .NE. 0.D0) THEN
              DE(IM,1)= -BM(IM)/DLE(IM)
C Error - Corrctn FM Nov. 2009
C              DO 44 I=2, 10 !MCOEF
              DO 44 I=2, MCOEF
                DE(IM,I)=-DE(IM,I-1)/DLE(IM)
 44           CONTINUE
            ENDIF
 45       CONTINUE
        ENDIF
 
C--------- Let's see exit, next
        IF(NRES.GT.0) WRITE(NRES,107)
 107    FORMAT(/,15X,' Exit  face  ')
C 107    FORMAT(/,15X,' FACE  DE  SORTIE  ')
 
        DL0 = 0.D0
        DO 6 IM=NM0,NM
 6        DL0 = DL0+DLS(IM)
        IF(DL0 .EQ. 0.D0) THEN
          FINTS = XLS
          XLS=0.D0
          IF(NRES.GT.0) THEN
            WRITE(NRES,105) 'Exit field model is '
            WRITE(NRES,FMT='(15X,''FINTS, gap : '',
     >      1P,2(1X,E12.4))') FINTS,GAP
          ENDIF
          IF(KFL .EQ. MG) THEN
            IF(KUASEX .EQ. MPOL+1) THEN
C------------- Set exit wedge correction in SBR INTEGR
              IF(BM(1) .NE. 0.D0) THEN
C FM, Oct. 2011. Avoid wedge correction with short lens
C FM, 2006
                IF(XL .GT. 2.D0) CALL INTEG2(ZERO,FINTS,GAP)
              ENDIF
            ENDIF
          ENDIF
        ELSE
          IF( IFB .EQ. 0 .OR. IFB .EQ. -1 ) THEN
            MIXFF = .FALSE.
            DO 61 IM = NM0,NM
              IF(.NOT. MIXFF) THEN
                IF(DLS(IM) .EQ. 0.D0 .AND. BM(IM) .NE. 0.D0) 
     >           MIXFF= .TRUE.
              ENDIF
 61         CONTINUE

            IF(MIXFF) THEN
              IF(IFB .EQ. 0) THEN
                IFB = 1

              ELSE
                IFB = 2

              ENDIF
            ENDIF
          ENDIF
          IF(NRES.GT.0) THEN
            WRITE(NRES,130) XLS
            WRITE(NRES,131) ( LMNT(IM),DLS(IM) ,IM=NM0,NM)
            WRITE(NRES,132) NCS,(CS(I),I=1,NCS)
          ENDIF
          DO 46 IM=NM0,NM
            IF(DLS(IM) .NE. 0.D0) THEN
              DS(IM,1)=  BM(IM)/DLS(IM)
C Error - Corrctn FM Nov. 2009
C              DO 461 I=2, 10 !MCOEF
              DO 461 I=2,MCOEF
                DS(IM,I)= DS(IM,I-1)/DLS(IM)
 461          CONTINUE
            ENDIF
 46       CONTINUE
        ENDIF
 
      ENDIF
C---------- end of test DLE or DLS=0
 
C----- Some more actions about Magnetic Dipole components :
      IF( KUASEX .EQ. MPOL+1 .AND. KFL .EQ. MG ) THEN
C----- MULTIPOL
        IF(XE .EQ. 0.D0) THEN
C------- Entrance sharp edge field model
C          IF(NM .EQ. 1 .AND. BM(1) .NE. 0.D0) THEN
          IF(BM(1) .NE. 0.D0) THEN
            IF(NRES.GT.0) 
     >      WRITE(NRES,FMT='(/,''  ***  Warning : sharp edge model, '',
     >      ''vertical wedge focusing approximated with '',
     >      ''first order kick. FINT at entrance = '',1P,2G12.4)') 
     >      FINTE
          ENDIF
        ENDIF
        IF(XLS .EQ. 0.D0) THEN
C------- Exit sharp edge field model
C          IF(NM .EQ. 1 .AND. BM(1) .NE. 0.D0) THEN
          IF(BM(1) .NE. 0.D0) THEN
            IF(NRES.GT.0) 
     >      WRITE(NRES,FMT='(/,''  ***  Warning : sharp edge model, '',
     >      '' vertical wedge focusing approximated with '',
     >      ''first order kick. FINT at exit = '',1P,2G12.4)') 
     >      FINTS
          ENDIF
        ENDIF
      ENDIF

C----- Case erron (errors)
      IF(NRES.GT.0) then
        do irr = 1, mxerr 
          ok = (empty(lbl1(irr)) .or. lbl1(irr).eq.label(noel,1)) 
     >    .and.(empty(lbl2(irr)) .or. lbl2(irr).eq.label(noel,2)) 
          if(ok) then
            do i = 1, mpol
              if(kpol(irr,i) .eq. 1) then 
                WRITE(NRES,FMT='(/,15x,
     >          ''ERRORS ARE SET, accounted for in the fields '', 
     >          ''above.'')')
                WRITE(NRES,FMT=
     >          '(15x,''Case of MULTIPOL with labels : '',4a,i4,a,i4)')
     >          lbl1(irr), ' / ',lbl2(irr),' /  error set # ',irr,
     >          ', # element = ',noel  
                WRITE(NRES,FMT='(15x,
     >          ''Pole#, error type, A/R, G/U : '',i1,3(2x,a))')
     >          i, typerr(irr,i), typar(irr,i), typdis(irr,i)
                WRITE(NRES,FMT='(15x,a,1p,3(e14.6,2x))') 
     >          'err_center, err_sigma, err_cutOff : ',
     >          errcen(irr,i),errsig(irr,i),errcut(irr,i)
              endif
            enddo
          ENDIF
        enddo
      ENDIF

      XI = 0.D0
      XLIM = XL + XE + XLS
      XF = XLIM
      XS = XL + XE
      SUM=0.D0

C----- Passage obligatoire sur les EFB's si
C       melange Mpoles-crenau + Mpoles-champ de fuite
      IF( IFB .EQ. -1 ) THEN
        AFB(1) = 1.D0
        BFB(1) = 0.D0
        CFB(1) = -XE         
      ELSEIF( IFB .EQ. 1 ) THEN
        AFB(2) = 1.D0
        BFB(2) = 0.D0
        CFB(2) = -XS       
      ELSEIF( IFB .EQ. 2 ) THEN
        AFB(1) = 1.D0
        BFB(1) = 0.D0
        CFB(1) = -XE       
        AFB(2) = 1.D0
        BFB(2) = 0.D0
        CFB(2) = -XS       
      ENDIF
C      IF( XE .GE. XL .OR. XLS .GE. XL ) THEN
C        IER = IER+1
C        TXT(IER) =  
C     >   'Overlapping of fringe fields is too large. Check XE, XS < XL'
C      ENDIF

cC----- Magnetic MULTIPOL with non-zero dipole
cC--------- Automatic positionning in SBR TRANSF
cC           Dev normally identifies with the deviation
cC             that would occur in a sharp edge dipole magnet.
c        IF(BM(1) .NE. 0.D0) THEN
c          DEV= 2.D0* ASIN(.5D0*XL*BM(1)/BORO)
c        ELSE
c          DEV= 0.D0
c        ENDIF
cC------------------------  TESTS for COSY
cC         DEV = 2.D0 * PI /24.D0  
cCC         DEV = 2.D0 * ASIN(.5D0 * XL * 14.32633183D0 / BORO )
cC         DEV = ALE
cC-----------------------------------------

      CALL CHXC1R(
     >            KPAS)
      IF(KPAS.GE.1) THEN
        AREG(1)=1.D0
        BREG(1)=0.D0
        CREG(1)=-2.D0*XE
        AREG(2)=1.D0
        BREG(2)=0.D0
        CREG(2)=-2.D0*XS+XLIM
        CALL INTEG6(AREG,BREG,CREG)
      ENDIF

      IF(IER.NE.0) GOTO 99

 98   RETURN

 99   CONTINUE
      LUN=NRES
      IF(LUN.LE.0) LUN = 6
      WRITE(LUN,FMT='(//)')
      WRITE(LUN,FMT='(10X,A11,I2,3X,A)')
     >               ('*** ERROR #',I,TXT(I),I=1,IER)
C----- Execution stopped :
      RETURN 1
      
      ENTRY MULTKL(
     >             AL, AK1)
C     >     aL, aK1, aK2, aK3)
      aL = sXL
      aK1 = AKS(1)
C      aK2 = AKS(2)
C      aK3 = AKS(3)
      RETURN
      
      ENTRY MULTP2(irri,iPOLI,TYPERI,TYPAI,TYPDII,
     >ERRCEI,ERRSII,ERRCUI,lbl1i,lbl2i)
      ERRON = .TRUE.
      irr = irri
      IPOL = IPOLI
      KPOL(irr,IPOL) = 1
      TYPERR(irr,IPOL)=      TYPERI
      TYPAR(irr,IPOL)=      TYPAI
      TYPDIS(irr,IPOL)=      TYPDII
      ERRCEN(irr,IPOL)=      ERRCEI
      ERRSIG(irr,IPOL)=      ERRSII
      ERRCUT(irr,IPOL)=      ERRCUI
      LBL1(irr) = LBL1I
      LBL2(irr) = LBL2I
      RETURN

      ENTRY MULTP4
      ERRON = .FALSE.
      RETURN      
      END
