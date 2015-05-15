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
      SUBROUTINE CHXP(ND,KALC,KUASEX,
     >                               XL,DSREF,NDD)
      USE DYNHC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.AIM_2.H"     ! COMMON/AIM/ AE,AT,AS,RM,XI,XF,EN,EB1,EB2,EG1,EG2
      INCLUDE 'PARIZ.H'
      INCLUDE "XYZHC.H"
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.CONST2.H"     ! COMMON/CONST2/ ZERO, UN
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "C.DROITE.H"     ! COMMON/DROITE/ CA(9),SA(9),CM(9),IDRT
      INCLUDE "C.INTEG.H"     ! COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
C      LOGICAL ZSYM
      INCLUDE "C.TYPFLD.H"     ! COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      INCLUDE "C.ORDRES.H"     ! COMMON/ORDRES/ KORD,IRD,IDS,IDB,IDE,IDZ
      INCLUDE "C.PTICUL.H"     ! COMMON/PTICUL/ AM,Q,G,TO
      INCLUDE "C.REBELO.H"   ! COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      INCLUDE 'MXFS.H'
      INCLUDE 'MXSCL.H'
      INCLUDE "C.SCAL.H"     ! COMMON/SCAL/ SCL(MXF,MXS,MXSCL),TIM(MXF,MXS),NTIM(MXF),KSCL
C      COMMON/SCAL/SCL(MXF,MXS),TIM(MXF,MXS),NTIM(MXF),JPA(MXF,MXP),KSCL
      INCLUDE "C.STEP.H"     ! COMMON/STEP/ TPAS(3), KPAS
C      COMMON/STEP/ KPAS, TPAS(3) 
      INCLUDE "C.SYNRA.H"     ! COMMON/SYNRA/ KSYN
      INCLUDE "C.VITES.H"     ! COMMON/VITES/ U(6,3),DQBR(6),DDT(6)
  
      CHARACTER(39) TXTT, TXTS, TXTA
      CHARACTER(11) TXTEMP
      SAVE IPREC

      LOGICAL FITING
      LOGICAL NEWFIC

      INCLUDE 'FILPLT.H'
      DIMENSION DBDX(3)

      PARAMETER (PI = 4.D0 * ATAN(1.D0))
      PARAMETER (RAD = PI/180.D0)

      DATA DTA1 / 0.D0 /
      DATA NEWFIC / .TRUE. /
      data scal / 1.d0 /

c         open(unit=88,file='fort.88')

      ZSYM=.TRUE.
      CALL CHAMK2(UN)
      CALL RAZ(DBDX,3)
      CALL CHAMK4(DBDX,3)

C- KALC = TYPE CALCUL : ANALYTIQUE + SYM PLAN MEDIAN (1) , ANALYTIQUE 3D (3)
C   &  CARTE (2)

      IF    (KALC.EQ.2 ) THEN
C------- Field is defined by maps

        IF(KUASEX .EQ. 22) THEN
C--------- POLARMES keyword
          KP = NINT(A(NOEL,ND+10))
          NDD = ND+20
        ELSEIF(KUASEX .EQ. 20 .OR. KUASEX .EQ. 21) THEN
C--------- AIMANT, DIPOLE-M keywords
C Modif, FM, Dec. 05
C          KP = NINT(A(NOEL,ND+1))
C          NDD = ND+2
          KP = NINT(A(NOEL,ND+3))
          NDD = ND+4
        ELSEIF(KUASEX .EQ. 7 .OR. KUASEX .EQ. 2) THEN
C--------- TOSCA keyword using cylindrical mesh (MOD.ge.20)
          KP = NINT(A(NOEL,ND+10))
          NDD = ND+20
        ELSEIF(KUASEX.EQ.34 .OR. KUASEX.EQ.35) THEN
C--------- EMMA keyword using cylindrical mesh (MOD.ge.20)
          KP = NINT(A(NOEL,ND+10))
          NDD = ND+20
        ENDIF

      ELSE
C        KALC = 1 or 3
C        Field is defined by analytical models, mid plane (1) or 3D (3). 

        IF(KUASEX .EQ. 30) THEN
C--------- EMPTY

        ELSE

          KP = NINT(A(NOEL,ND+10))
          NDD = ND+11
        ENDIF

      ENDIF

      CALL FITSTA(5,FITING)
      IF(.NOT. FITING) THEN        
        IF    (KALC.EQ.2 ) THEN
C--------- Field is defined by maps
          IF(KUASEX .EQ. 22) THEN
C----------- POLARMES keyword
            LF = NINT(A(NOEL,1))
C  LST=1(2) : PRINT step by step coord. and field in zgoubi.res (zgoubi.plt)
            LST =LSTSET(NINT(A(NOEL,2)))
          ELSEIF(KUASEX .EQ. 20 .OR. KUASEX .EQ. 21) THEN
C----------- AIMANT, DIPOLE-M keywords
            LF =  NINT(A(NOEL,2))
            LST = LSTSET(NINT(A(NOEL,3)))
C Modif, FM, Dec. 05
C          KP = NINT(A(NOEL,ND+1))
C          NDD = ND+2
          ELSEIF(KUASEX .EQ. 7 .OR. KUASEX .EQ. 2) THEN
C----------- TOSCA keyword using cylindrical mesh (MOD.ge.20)
            LF =  NINT(A(NOEL,1))
            LST =LSTSET(NINT(A(NOEL,2)))
          ELSEIF(KUASEX.EQ.34 .OR. KUASEX.EQ.35) THEN
C----------- EMMA keyword using cylindrical mesh (MOD.ge.20)
            LF =  NINT(A(NOEL,1))
            LST =LSTSET(NINT(A(NOEL,2)))
          ENDIF
        ELSE
C          KALC = 1 or 3
C          Field is defined by analytical models, mid plane (1) or 3D (3). 
          LF = 0
          LST =LSTSET(NINT(A(NOEL,1)))
        ENDIF
      ELSE
        LF = 0
        LST = 0
      ENDIF

      IF(LST.EQ.2 .OR. LST.GE.4) CALL OPEN2('CHXP',NPLT,FILPLT)

C----- FACTEUR D'ECHELLE DES Champs. UTILISE PAR 'SCALING'
C Field scale factor. Used by  'SCALING'

C--------------------------------------------------------------
C Problem here with my laptop (ok w owl !) : this write(89 is necessary for the FIT problem 
C /home/meot/zgoubi/struct/folks/thomasPlanche/FITBugWithPARTICLE/FitWorks.res
C to work (otherwise does run but won't fit to the expected values). 
C Otherwise scal=NaN is passed to dipi
          CALL FITSTA(5,FITING)
            if(FITING) then
       write(89,*) ' chxp  SCAL , SCAL0() ',SCAL , SCAL0()
        rewind(89)
        endif
C--------------------------------------------------------------

      SCAL = SCAL0()
      IF(KSCL .EQ. 1) SCAL = SCAL0()*SCALER(IPASS,NOEL,
     >                                                 DTA1)
      AE = 0.D0
      AS = 0.D0
      IDRT = 0

        IDE=4
        IDB=2
        IDZ=2

      RFR = 0.D0

      GOTO (2001, 2002, 2003) KALC
 
 
 2001 CONTINUE
C----------- KALC = 1: Define field in the median plane 
C                  with median plane symetry

C Default IDB value is set here : 
      IRD = KORD
      IF(IRD.EQ.4) IDB=4

      IF(KUASEX .EQ. 27 )   THEN
C-------- FFAG                ffag radial

        CALL FFAGI(SCAL,
     >                  DSREF,IRD,IDB)
        IDZ=3
C Modif, FM, Dec. 05
C        KP = NINT(A(NOEL,ND+1))
C        NDD = ND+2
C Modif, FM, Sept 2014
C        KP = NINT(A(NOEL,ND+3))

        KP = NINT(A(NOEL,ND+1))
        NDD = ND+2
C        DSREF = ABS(DEV * (XL/(2.D0 * SIN(DEV/2.D0))))

      ELSEIF(KUASEX .EQ. 30) THEN
C------- Installiing SBEND, not yet complete 
        
          CALL SBENDI(SCAL,
     >                    DSREF)
          KP = NINT(A(NOEL,ND+NND))
C          DSREF = ABS(DEV * (XL/(2.D0 * SIN(DEV/2.D0))))

      ELSEIF(KUASEX .EQ. 31 )   THEN
C-------- DIPOLE

        CALL DIPI(SCAL,
     >                 DSREF)
C Modif, FM, Dec. 05
C        KP = NINT(A(NOEL,ND+1))
C        NDD = ND+2
        KP = NINT(A(NOEL,ND+3))
        NDD = ND+4
C        DSREF = ABS(DEV * (XL/(2.D0 * SIN(DEV/2.D0))))

      ELSEIF(KUASEX .EQ. 32 )   THEN
C-------- DIPOLES

        CALL DIPSI(SCAL,
     >                  DSREF,IRD,IDB)

C Modif, FM, Dec. 05
C        KP = NINT(A(NOEL,ND+1))
C        NDD = ND+2
        KP = NINT(A(NOEL,ND+3))
        NDD = ND+4
C        DSREF = ABS(DEV * (XL/(2.D0 * SIN(DEV/2.D0))))

      ELSEIF(KUASEX .EQ. 33 )   THEN
C-------- FFAG-SPI     spiral ffag

        CALL FFGSPI(SCAL,
     >                  DSREF,IRD,IDB)
C     >                  XL,DEV)
        IDZ=3

        NP = ND + 1
        KP = NINT(A(NOEL,NP))
        NDD = NP+1
C        DSREF = ABS(DEV * (XL/(2.D0 * SIN(DEV/2.D0))))

      ELSEIF(KUASEX .EQ. 40 )   THEN
C-------- CYCLOTRON

        CALL CYCLO(SCAL,
     >                  DSREF,IRD,IDB)

        NP = ND + 3
        KP = NINT(A(NOEL,NP))
        NDD = NP+1
C        DSREF = ABS(DEV * (XL/(2.D0 * SIN(DEV/2.D0))))

      ENDIF

      GOTO 99
 
C---------------------------------------------------------------
 2002 CONTINUE
C----- KALC = 2 : either read (TOSCAP, POLMES), or generate (CARLA, DIPOLM) a field map

      IF(KUASEX.EQ.2 .OR. KUASEX.EQ.7) THEN
C TOSCA keyword with MOD.ge.20. 

        IF(KUASEX .EQ. 2) THEN
          NDIM = 2
          CALL TOSCAP(SCAL,NDIM,
     >                          BMIN,BMAX,BNORM,XNORM,YNORM,ZNORM,
C     >                          BMIN,BMAX,BNORM,
     >                          ABMI,RBMI,ZBMI,ABMA,RBMA,ZBMA,NEWFIC)

        ELSEIF(KUASEX .EQ. 7) THEN
          NDIM = 3
          CALL TOSCAP(SCAL,NDIM,
     >                          BMIN,BMAX,BNORM,XNORM,YNORM,ZNORM,
C     >                          BMIN,BMAX,BNORM,
     >                          ABMI,RBMI,ZBMI,ABMA,RBMA,ZBMA,NEWFIC)

        ENDIF

        RFR = RM

      ELSEIF(KUASEX.EQ.34 .OR. KUASEX.EQ.35) THEN
C EMMA. Polar -> MOD.ge.20. 

        IF(KUASEX .EQ. 34) THEN

          NDIM = 2
          CALL EMMAP(SCAL,NDIM,
     >                          BMIN,BMAX,BNORM,
     >                          ABMI,RBMI,ZBMI,ABMA,RBMA,ZBMA)

        ELSEIF(KUASEX .EQ. 35) THEN
          NDIM = 3
          CALL EMMAP(SCAL,NDIM,
     >                          BMIN,BMAX,BNORM,
     >                          ABMI,RBMI,ZBMI,ABMA,RBMA,ZBMA)

        ENDIF

        RFR = RM

      ELSEIF(KUASEX .EQ. 20) THEN
C AIMANT keyword

        CALL CARLA(SCAL,
     >                          BMIN,BMAX,BNORM,
     >                          ABMI,RBMI,ZBMI,ABMA,RBMA,ZBMA)
        NDIM = 2

      ELSEIF(KUASEX .EQ. 21) THEN
C DIPOLE-M keyword

        CALL DIPOLM(SCAL,
     >                          BMIN,BMAX,BNORM,
     >                          ABMI,RBMI,ZBMI,ABMA,RBMA,ZBMA)
        NDIM = 2
        IF(NINT(A(NOEL,ND+3)) .EQ. 2) RFR = A(NOEL,ND+4)

      ELSEIF(KUASEX .EQ. 22) THEN
C POLARMES keyword

        CALL POLMES(SCAL,KUASEX,
     >                          BMIN,BMAX,BNORM,
     >                          ABMI,RBMI,ZBMI,ABMA,RBMA,ZBMA,NEWFIC)
        NDIM = 2
        RFR = RM

      ENDIF

      IF(NRES.GT.0) THEN
C        I2=2 introduced to avoid compiler complainig when IZ=1...
        I2 = 2
        WRITE(NRES,203) 
     >      XH(1),XH(IXMA),XH(IXMA)-XH(1),
     >      YH(1),YH(JYMA),YH(JYMA)-YH(1),
     >      BMIN/BNORM,BMAX/BNORM,
     >      ABMI,RBMI,ZBMI,ABMA,RBMA,ZBMA, 
     >      BNORM,BMIN,BMAX,IXMA,JYMA,KZMA,
     >      XH(2)-XH(1),YH(2)-YH(1),ZH(I2)-ZH(1)
  203   FORMAT(
     >   //,5X,'Field map limits, angle, min, max, max-min (rad) :',
     >             3E15.6,
     >    /,5X,'Field map limits, radius, min, max, max-min (cm) :',
     >             3E15.6,
     >    /,5X,' Min / max  fields  drawn  from  map  data : ', 
     >                           1P,G11.3,T80,' / ',G11.3,
     >    /,5X,'  @  x-node, y-node, z-node :              ', 
     >                             3G10.3,T80,' / ',3G10.3,
     >    /,5X,'Normalisation  coeff.  BNORM   :', G12.4,
     >    /,5X,'Field  min/max  normalised  :             ', 
     >                           1P,G11.3,T64,' / ',G11.3,
     >    /,5X,'Nbre  of  nodes  in  x/y/z : ',I4,'/',I4,'/',I4,
     >    /,5X,'Node  distance  in   x/y/z : ',G12.4,'/',
     >                  G12.4,'/',G12.4)
CCCCCCCC     >                  G12.4,'/',G12.4,'  (rad/cm/cm)')

        IF    (NDIM .EQ. 1) THEN
          WRITE(NRES,111) IRD
        ELSEIF(NDIM .EQ. 2) THEN
          WRITE(NRES,111) IRD
 111      FORMAT(/,20X,' Option  for  interpolation :',I2)
          IF(IRD .EQ. 2) THEN
            WRITE(NRES,113)
 113        FORMAT(20X,' Smoothing  using  9  points ')
          ELSE
C           .... IRD=4 OU 25
            WRITE(NRES,115)
 115        FORMAT(20X,' Smoothing  using  25  point ')
          ENDIF
        ELSEIF(NDIM .EQ. 3) THEN
          WRITE(NRES,119) IRD
 119      FORMAT(/,20X,' Option  for  interpolation :  3-D  grid',2X,
     >    '  3*3*3  points,   interpolation  at  ordre ',I2)
        ENDIF
 
        IF(LF .NE. 0) CALL FMAPW(ZERO,RFR,2)
 
      ENDIF
C     ... endif NRES>0
 
      GOTO 99
 
C---------------------------------------------------------------
 2003 CONTINUE
C------- KALC = 3 : Full 3D calculation from analytical model 
C ELMIR,  ELCYLDEF
        IRD = KORD
        IF(KUASEX .EQ. 24)   THEN
C--------- ELCYLDEF
          IDE = 2
          CALL ELCYLD(SCAL,
     >                     XL)
C           Motion in this lmnt has no z-symm. 
          ZSYM=.FALSE.

          KP = NINT(A(NOEL,ND+10))
          NDD = ND+20

        ELSEIF(KUASEX .EQ. 26 )   THEN
C--------- ELCMIR
          CALL ELCMII(SCAL,  
     >                     XL)
C           Motion in this lmnt has no z-symm. 
          ZSYM=.FALSE.

        ENDIF

      GOTO 99



C-----------------------------------------------------------------
 99   CONTINUE

      PAS = A(NOEL,ND)

      IF(KSCL .EQ. 1
C------ Cavity
     >               .OR. KSYN .EQ. 1) THEN
C--------------------- SR Loss
        IF(NRES .GT. 0)  WRITE(NRES,199) SCAL
 199    FORMAT(/,20X,'Field has been * by scaling factor ',1P,G16.8)
      ENDIF

      IF( KPAS .GT. 0 ) THEN
C------- Coded step size of form  #stp_1|stp_2|stp_3, stp_i= arbitrary integer
C                                    entr body exit
C      (max stp_i is MXSTEP as assigned in SBR INTEGR), 
C                      step size = length/stp_i
C------- KPAS=2 -> Variable step

C Modif, FM, Dec. 05
C        STP1 = A(NOEL,ND-2)
C        STP2 = A(NOEL,ND-1)
C        STP3 = A(NOEL,ND)
        STP1 = A(NOEL,ND)
        STP2 = A(NOEL,ND+1)
        STP3 = A(NOEL,ND+2)
        TPAS(1) = 2.D0 * AE * RM / STP1
        TPAS(2) =  (AT-2.D0 *(AE+AS))*RM / STP2
        TPAS(3) = 2.D0 *  AS * RM / STP3

        IF(NRES.GT.0) THEN
          WRITE(TXTA,FMT='(1P,'' / '',G11.3,'' / '')') TPAS(2)/RM
          WRITE(TXTS,FMT='(1P,'' / '',G11.3,'' / '')') TPAS(2)
          WRITE(TXTT,FMT='('' /   central   / '')') 
          IF(TPAS(1) .NE. 0.D0) THEN
            WRITE(TXTEMP,FMT='(1P,G11.3)') TPAS(1)/RM
            TXTA = TXTEMP//TXTA(1:28)
            WRITE(TXTEMP,FMT='(1P,G11.3)') TPAS(1)
            TXTS = TXTEMP//TXTS(1:28)
            TXTT = ' entrance  '//TXTT(1:28)
          ENDIF 
          IF(TPAS(3) .NE. 0.D0) THEN
            WRITE(TXTEMP,FMT='(1P,G11.3)') TPAS(3)/RM
            TXTA = TXTA(1:17)//TXTEMP
            WRITE(TXTEMP,FMT='(1P,G11.3)') TPAS(3)
            TXTS = TXTS(1:17)//TXTEMP
            TXTT = TXTT(1:28)//'   exit    '
          ENDIF 

          WRITE(NRES,FMT='(/,25X,'' Integration  step :'')') 
          WRITE(NRES,FMT='(30X,A,'' (rad) '')') TXTA
          WRITE(NRES,FMT='(30X,A,'' (cm,  at mean radius)'')') TXTS
          WRITE(NRES,FMT='(30X,A,/,35X,''region'')') TXTT

          IF(KPAS .EQ. 2) THEN
              IF(NRES.GT.0)
     >         WRITE(NRES,FMT='(/,25X,'' ++ Variable step ++'',/,
     >            25X,'' PRECISION ='',1P,G12.4,/)') 10.D0**(-IPREC)
          ENDIF


        ENDIF

      ELSE

        IF(NRES.GT.0) WRITE(NRES,FMT='(/,20X,''Integration step :'',
     >    1P,G12.4,'' cm   (i.e., '',G12.4,'' rad  at mean radius'',
     >    '' RM = '',G12.4,'')'')') 
     >      PAS, PAS/RM, RM

      ENDIF

      IF(KFLD.GE.LC) THEN
        IF(Q*AM.EQ.0D0) 
     >  CALL ENDJOB('Provide  mass  and  charge - using PARTICUL',-99)
      ENDIF

      RETURN

      ENTRY CHXP1(
     >            KPASO)
      KPASO = KPAS
      RETURN
      ENTRY CHXP2(KPASI,IPRECI)
      KPAS = KPASI
      IPREC = IPRECI
      IF(KPAS .EQ. 2) CALL DEPLAW(.TRUE.,IPREC)
      RETURN
 
      END
