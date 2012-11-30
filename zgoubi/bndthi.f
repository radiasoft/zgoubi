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
      SUBROUTINE BNDTHI(ND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/AIM/ BO,RO,FG,GF,XI,XF,EN,EB1,EB2,EG1,EG2
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      PARAMETER(MCOEF=6)
      COMMON/CHAFUI/ XE,XS,CE(MCOEF),CS(MCOEF),QCE(MCOEF),QCS(MCOEF)
      INCLUDE "MAXTRA.H"
      COMMON/CHAMBR/ LIMIT,IFORM,YLIM2,ZLIM2,SORT(MXT),FMAG,BMAX
     > ,YCH,ZCH
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      COMMON/CONST2/ ZERO, UN
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      COMMON/DROITE/ CA(9),SA(9),CM(9),IDRT
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      CHARACTER LET
      COMMON/FAISCT/ LET(MXT)
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      COMMON/MARK/ KART,KALC,KERK,KUASEX
      PARAMETER(MPOL=10)
      COMMON/MULTPL/ BBM(MPOL),DLE(MPOL),DLS(MPOL)
     >,DE(MPOL,MCOEF),DS(MPOL,MCOEF),RTB(MPOL)
      LOGICAL ZSYM
      COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      INCLUDE 'MXFS.H'
      COMMON/SCAL/ SCL(MXF,MXS),TIM(MXF,MXS),NTIM(MXF),KSCL
C      COMMON/SCAL/SCL(MXF,MXS),TIM(MXF,MXS),NTIM(MXF),JPA(MXF,MXP),KSCL
      COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)
      COMMON/SYNRA/ KSYN
      COMMON/TRAJ/ Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT
      COMMON/UNITS/ UNIT(MXJ)
 
      EQUIVALENCE (RTB(1),CTE),(RTB(2),STE),(RTB(4),CTS),(RTB(5),STS)

      LOGICAL DRT

      LOGICAL CHGRFE,CHGRFS,EVNT
      LOGICAL BONUL, LNUL
      SAVE NND
      DATA NND / 10 / 

C----- Magnet length, straight axis
      XL =A(NOEL,10)
C----- Field
      BFLD  =A(NOEL,12)
      KFLD = MG
      IF(BFLD .EQ. 0.D0) KFLD = 0
      
C FIELDS ARE DEFINED IN CARTESIAN COORDINATES
      KART = 1
      IF ( LNUL(XL) ) RETURN
      IF ( BONUL(XL,XL) ) THEN
        CALL SCUMW(XL)
        RETURN
      ENDIF

C     ... FACTEUR D'ECHELLE DES ChampS. UTILISE PAR 'SCALING'
      SCAL = SCAL0()
      IF(KSCL .EQ. 1) SCAL = SCAL0()*SCALER(IPASS,NOEL,
     >                                                 DTA1)
      BFLD = BFLD * SCAL

      IF(BFLD .EQ. 0.D0) KFLD = 0
      DEV = 2.D0 * ASIN(XL/2.D0/(BORO/BFLD))

          KP = NINT(A(NOEL,ND+NND))
      IF( KP .EQ. 3 ) THEN
        IF(A(NOEL,73) .NE. 0.D0) DEV = -A(NOEL,73) * 2.D0
      ENDIF

      XE = A(NOEL,20)
      DLE(1) = A(NOEL,21)
      WE = A(NOEL,22)
      TE=.5D0*DEV-WE
      XLS = A(NOEL,40)
      DLS(1) = A(NOEL,41)
      WS =  A(NOEL,42)
      TS=.5D0*DEV-WS

C----------- No fringe field but for 1st order V-kick
        FINTE = XE
        XE=0.D0
c        if (xe.lt.1d-10) xe = 5.d0*tan(abs(te))
        GAPE = -DLE(1)
        FINTS = XLS
        XLS=0.D0
C Causes appropriate change-ref once EFB has been reached
C        if (xls.lt.1d-10) xls = 5.d0*tan(abs(ts))
        GAPS = -DLS(1)
      XI = 0.D0
      XLIM = XL + XE + XLS
      XF = XLIM
      XS = XL + XE
      CTE=COS(TE)
      STE=SIN(TE)
      CTS=COS(TS)
      STS=SIN(TS)
 
C----- SHARP EDGE => INTEGR STOPPE SUR DR. DE COUPURE
        IDRT = 2
        CA(1)=CTE
        SA(1)=STE
        CM(1)=-XE*CA(1)
        CA(2)= CTS
        SA(2)=-STS
        CM(2)=-XS*CA(2)

      IF(NRES.GT.0) THEN
        WRITE(NRES,107) ' Theoretical BEND',
     >          XL,BORO/BFLD*DEV,DEV*DEG,DEV
 107    FORMAT(1P, /,5X,' +++++  ',A10,'  : ',
     >       //,15X, ' Length    = ',E14.6,' cm'
     >       ,/,15X, ' Arc length    = ',E14.6,' cm'
     >       ,/,15X, ' Deviation    = ',E14.6,' deg.,  ',E14.6,' rad',/)
        WRITE(NRES,103) BFLD,BORO/BFLD
 103    FORMAT(1P,15X,' Field  =',E14.6,'  kG ',
     >        /,15X, ' Reference  radius  (BRo/B)  = ',E14.6,'  cm')
 
C        WRITE(NRES,104) 'D''ENTREE'
C 104    FORMAT(/,15X,' FACE  ',A)
        WRITE(NRES,104) 'Entrance ' 
 104    FORMAT(/,15X,A9,' face  ')
        WRITE(NRES,101) XE,DLE(1),WE
 101    FORMAT(15X,' DX = ',F10.3,'    LAMBDA = ',F10.3
     >        ,/,15X,' Wedge  angle  =',F10.6,' RD')

        WRITE(NRES,104) 'Exit     '
        WRITE(NRES,101) XLS,DLS(1),WS

        WRITE(NRES,FMT='(/,
     >    ''***************************************'',
     >    ''Vertical wedge focusing still needs be approximated with'',
     >    '' first order kick, FINT values entr/exit : '',1P,2G12.4)') 
     >         FINTE, FINTS

      ENDIF

          IF( KP .EQ. 3 ) THEN
C            Stop if XCE .ne. 0. To be provisionned...
            IF(A(NOEL,ND+NND+1) .NE. 0.D0) 
     >                  STOP ' KPOS=3 does not support XCE.ne.0'
C             Calculate ALE as half deviation. 
            IF(A(NOEL,ND+NND+3).EQ.0.D0) 
     >                      A(NOEL,ND+NND+3)=-DEV/2.D0
C Modified, FM, Dec 05 :
C             Calculate XCE, YCE for entrance change of referential    
             YSHFT = A(NOEL,ND+NND+2)
             A(NOEL,ND+NND+1) = - YSHFT * SIN(A(NOEL,ND+NND+3))
             A(NOEL,ND+NND+2) =   YSHFT * COS(A(NOEL,ND+NND+3))
          ENDIF
          DSREF = ABS(DEV * (XL/(2.D0 * SIN(DEV/2.D0))))
          
          CALL SCUMW(DSREF)

C-------------------------------------------------------------------
      XCE  = A(NOEL,ND+NND+1)
      YCE  = A(NOEL,ND+NND+2)
      ALE  = A(NOEL,ND+NND+3)

      IF(KSCL .EQ. 1
C------ Cavity
     >  .OR. KSYN.EQ.1) THEN
C------------SR Loss

        IF(NRES .GT. 0)  WRITE(NRES,199) SCAL
 199    FORMAT(/,20X,'Field has been * by scaling factor ',1P,G16.8)

      ENDIF

      IF( KP .EQ. 3 )  THEN
        IF(NRES.GT.0) WRITE(NRES,FMT='(/,15X,
     >     ''Automatic positioning of element, XCE, YCE, ALE  ='',
     >                 1P,3G16.8,'' cm/cm/rad'' : )') XCE, YCE, ALE
      ENDIF

      XFE=-XE
      XFS=XS-XLIM
      IF(KP .EQ. 1)  THEN
C------- Optical elmnt undergoes new positionning. Its frame becomes the new 
C        reference frame. No exit CHANGREF needed
        XCS = 0.D0
        YCS = 0.D0
        ALS = 0.D0
      ELSEIF(KP .EQ. 2)  THEN
C------- KP=2, element is misaligned. Hence compute XCS,YCS,ALS for automatic 
C        re-positionning of the exit frame to old position:
        XLM = XS-XE
        COL=COS(ALE)
        SIL=SIN(ALE)
        XTEMP=XCE-XLM*(1.D0-COL)
        YTEMP=YCE+XLM*SIL
        XCS=-XTEMP*COL-YTEMP*SIL
        YCS=XTEMP*SIL-YTEMP*COL
        ALS=-ALE
        IF(NRES.GT.0) WRITE(NRES,100) XCE,YCE,ALE
      ELSEIF(KP .EQ. 3)  THEN
C------- Optical elmnt is Z-rotated. Entrance and exit frames are 
C        tilted by (normally) half the deviation. 
        XCS = 0.D0
        YCS = -YCE/COS(ALE)
        ALS=ALE
      ENDIF

      X=XI
      XLIM=XF
 
C--------------------------------
C------- Now equivalent to TRANSF

      CHGRFE= ( KP .GE. 2 .AND. XCE*XCE + YCE*YCE + ALE*ALE .GT. ZERO )
      CHGRFS= ( KP .GE. 2 .AND. XCS*XCS + YCS*YCS + ALS*ALS .GT. ZERO )

C----- Events, such as spin tracking, in-flight decay, etc...
      EVNT = KSPN .EQ. 1 

C----- Droite de coupure entree
      DRT = IDRT .EQ. -1 .OR. IDRT .GE. 2

      XO=X
      DO 1 IT=1,IMAX
 
C-------- IEX<-1 <=> Particle stopped
        IF(IEX(IT) .LT. -1) GOTO 1
 
          CALL INITRA(IT)
          X=XO

          IF( CHGRFE ) CALL CHAREF(EVNT,XCE,YCE,ALE)
C          CALL CHAREF(EVNT,ZERO,ZERO,TE)

          IF( DRT  ) CALL DRTENT
 
          CALL BNDTH(BFLD,XL,-TE,FINTE,GAPE,-TS,FINTS,GAPS)
 
C          CALL CHAREF(EVNT,XFS,ZERO,ZERO)
          IF( CHGRFS ) THEN
            IF(NRES.GT.0)
     >        WRITE(NRES,108)LET(IT),KEX,(FO(J,IT),J=1,5),X ,Y,T,Z,P,IT
 108        FORMAT(2X,A1,2X,I2,F8.4,4F10.3,8X,F9.3,4F10.3,8X,I4)
            CALL CHAREF(EVNT,XCS,YCS,ALS)
          ENDIF

          CALL MAJTRA(IT)

 1    CONTINUE
 
      IF(NRES .GT. 0) THEN
        CALL SCUMR(
     >             DUM,SCUM,TCUM) 
        WRITE(NRES,FMT='(/,'' Cumulative length of optical axis = '',
     >  1P,G17.9,
     >'' m ;  Time  (for ref. rigidity & particle) = '', 
     >  1P,G14.6,'' s '')')  SCUM*UNIT(5), TCUM
      ENDIF

      RETURN
C100   FORMAT(/,5X,'ELEMENT  DECENTRE  PAR  RAPPORT  A  L''AXE  OPTIQUE'
C     1,/,10X,'CENTRE  DE  LA  FACE  D''ENTREE  EN  X =',1P,G10.2
C     2,' CM   Y =',G10.2,' CM   INCLINAISON =',G12.4,' RAD',/)
100   FORMAT(/,5X,'Element  is  mis-aligned  wrt.  the  optical  axis'
     1,/,10X,'Center  of  entrance  EFB  is  at    X =',1P,G12.4
     2,' CM   Y =',G12.4,' cm,  tilt  angle =',G14.6,' RAD',/)
 
      RETURN
      END
