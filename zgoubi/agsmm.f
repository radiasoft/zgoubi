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
      SUBROUTINE AGSMM(LMNT,KFL,MPOL,NPOL,SCAL,
     >  DEV,RT,XL,BM,DLE,DLS,DE,DS,XE,XS,CE,CS,BORO,DPREF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER LMNT(*)*(*)
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
      CHARACTER*80 TA
      COMMON/DONT/ TA(MXL,40)
      COMMON/DROITE/ CA(9),SA(9),CM(9),IDRT
      COMMON/EFBS/ AFB(2), BFB(2), CFB(2), IFB
      COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      LOGICAL ZSYM
      COMMON/OPTION/ KFLD,MG,LC,ML,ZSYM
      COMMON/PTICUL/ AM,Q,G,TO
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      INCLUDE 'MXFS.H'
      COMMON/SCAL/SCL(MXF,MXS),TIM(MXF,MXS),NTIM(MXF),KSCL
      PARAMETER (LBLSIZ=8)
      PARAMETER (KSIZ=10)
      CHARACTER FAM*(KSIZ),LBF*(LBLSIZ),KLEY*(KSIZ),LABEL*(LBLSIZ)
      COMMON/SCALT/ FAM(MXF),LBF(MXF,2),KLEY,LABEL(MXL,2)
      COMMON/SYNRA/ KSYN

C----------- MIXFF = true if combined sharp edge multpole + fringe field multpole
      LOGICAL SKEW, MIXFF
     
      CHARACTER DIM(2)*3, BE(2)*2, TXT(10)*80
      DIMENSION  AREG(2),BREG(2),CREG(2)

      DIMENSION AK1(6), AK2(6)
      DIMENSION WN(2), WA(2)
      PARAMETER (I3=3)

      DIMENSION DB(3)

      CHARACTER TXT30*80

      DATA DIM / 'kG ', 'V/m'/
      DATA BE / 'B-', 'E-'/
 
      DATA AK1, AK2 / 6*1.D0, 6*1.D0 /

      SKEW=.FALSE.
      CALL RAZ(BM,NPOL)
         
      GAP = 0.D0
      RO = 10.d0

      MOD = NINT(A(NOEL,10))
      MOD2 = NINT(10.D0*A(NOEL,10)) - 10*MOD
      DLL =A(NOEL,11)
      GAP =A(NOEL,12)
      IF(GAP.EQ.0) GAP = RO
      NBLW = NINT(A(NOEL,20))
      IF(NBLW .GT. 2) STOP ' SBR agsmm, NBLW cannot exceed 2'
      DO I = 1, NBLW
        WN(I) = A(NOEL,20+(2*I-1))
        WA(I) = A(NOEL,20+ 2*I)
      ENDDO

      CALL AGSKS(BORO*DPREF*CL9/1.D3,
     >                         AK1,AK2)
      CALL AGSK12(NOEL,RO,AK1,AK2,MOD,
     >                                XL,BM,ANGMM)
      CALL AGSBLW(MOD2,NOEL,ANGMM,NBLW,WN,WA,
     >                                       BM,I3)
      DEV = ANGMM
      DO I = 1, I3
        DB(I) = A(NOEL,12+I)
        BM(I) = SCAL * BM(I) * (1.D0 + DB(I))
      ENDDO

C            write(*,*) ' agsmm scal ',scal

        XE =A(NOEL,30)
        DO IM=1,3
          DLE(IM) =A(NOEL,30+IM)
        ENDDO
        CALL RAZ(CE,MCOEF)
        NCE = NINT(A(NOEL,40))
        DO I=1, NCE
          CE(I) =A(NOEL,40+I)
        ENDDO

        XLS =A(NOEL,50)
        DO IM = 1,3
          DLS(IM) =A(NOEL,50+IM)
        ENDDO
        CALL RAZ(CS,MCOEF)
        NCS = NINT(A(NOEL,60))
         DO I=1,NCS
          CS(I) =A(NOEL,60+I)
        ENDDO
C------- FRINGE FIELD NOT INSTALLED FOR
C        DECA, DODECA, ... 18-POLE
        DLE(4)=ZERO
        DLS(4)=ZERO
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
  
C----- Pole rotation
        DO 35 IM=1,3
          RT(IM)=A(NOEL,70+IM-1)
          SKEW=SKEW .OR. RT(IM) .NE. ZERO
 35     CONTINUE
 
        NM0 = 1
        NM = NPOL
 
C------- 
        DO IM = NM0+1,NM
          DLE(IM)  = DLE(NM0)*DLE(IM)
          DLS(IM)  = DLS(NM0)*DLS(IM)
        ENDDO

      IF(NRES.GT.0) THEN
        WRITE(NRES,100) 'AGS Dipole',XL,DEV,DEV*DEG,RO
 100    FORMAT(/,5X,' -----  ',A10,'  : ', 1P
     >  ,/,15X,' Length  of  element         = ',G17.8,'  cm'
     >  ,/,15X,' Deviation                   = ',G17.8,'  rad'
     >  ,5X,                                ' (',G17.8,' deg)'
     >  ,/,15X,' Reference  pole  radius RO  = ',G14.5,'  cm')
        WRITE(NRES,103) (BE(KFL),LMNT(IM),BM(IM),DIM(KFL),IM=NM0,NM)
 103    FORMAT(15X,2A,'  =',1P,G17.8,1X,A)

        IF    (MOD.EQ.1) THEN 
          TXT30=' centered multipol model'
        ELSEIF(MOD.EQ.2) THEN 
          TXT30=' long-shifted dipole model'
        ELSE
          TXT30=' short-shifted dipole model'
        ENDIF
        WRITE(NRES,FMT='(/,15X,'' Mode  = '',I2,'',   '',A)') MOD,TXT30

        IF(SKEW) WRITE(NRES,101) (LMNT(IM),RT(IM),IM=NM0,NM)
 101    FORMAT(15X,A,'  Skew  angle =',1P,G17.8,' RAD')
        IF(XL .NE. 0.D0) THEN
          IF( (XL-DLE(NM)-DLS(NM)) .LT. 0.D0) WRITE(NRES,102)
 102      FORMAT(/,10X,'Entrance  &  exit  fringe  fields  overlap, ',
     >    /,10X,'  =>  computed  gradient  is ',' G = GE + GS - 1 ')

        WRITE(NRES,108) NBLW
 108    FORMAT(/,15X,' Nbr of backleg windings : ',I1,/)
        IF(NBLW .GE.1) THEN
          DO IBLW = 1, NBLW
            WRITE(NRES,109) IBLW, NINT(WN(IBLW)), WA(IBLW)
 109        FORMAT(15X,' Backleg winding # ',I1',  Nbr of windings : '
     >      ,I1,',  intensity in that winding : ',1p,e15.6,' A')      
          ENDDO
        ENDIF
        ELSE
          GOTO 98
        ENDIF
      ENDIF
 
      DL0=0.D0
      SUM=0.D0

      DO IM=NM0,NM
        DL0=DL0+DLE(IM)+DLS(IM)
        SUM=SUM+BM(IM)*BM(IM)
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
        CALL INTEG1(ZERO,FINTE,GAP)
        CALL INTEG2(ZERO,FINTS,GAP)
        
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
          CALL INTEG1(ZERO,FINTE,GAP)
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
          CALL INTEG2(ZERO,FINTS,GAP)
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
        IF(XE .EQ. 0.D0) THEN
C------- Entrance sharp edge field model
C          IF(NM .EQ. 1 .AND. BM(1) .NE. 0.D0) THEN
C          IF(BM(1) .NE. 0.D0) THEN
            IF(NRES.GT.0) 
     >      WRITE(NRES,FMT='(/,''  ***  Warning : sharp edge model,'',
     >      '' vertical wedge focusing approximated with '',
     >      ''first order kick, FINT at entrance ='',1P,2G12.4)') 
     >      FINTE
C          ENDIF
        ENDIF
        IF(XLS .EQ. 0.D0) THEN
C------- Exit sharp edge field model
C          IF(NM .EQ. 1 .AND. BM(1) .NE. 0.D0) THEN
C          IF(BM(1) .NE. 0.D0) THEN
            IF(NRES.GT.0) 
     >      WRITE(NRES,FMT='(/,''  ***  Warning : sharp edge model,'',
     >      '' vertical wedge focusing approximated with '',
     >      ''first order kick, FINT at exit ='',1P,2G12.4)') 
     >      FINTS
C          ENDIF
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

 98   RETURN
      END