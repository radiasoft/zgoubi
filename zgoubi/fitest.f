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
      SUBROUTINE FITEST(SAVFT,FNAME,
     >                             IER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SAVFT
      CHARACTER(*) FNAME
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     !  COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
                              ! >IREP(MXT),AMQLU,PABSLU
      PARAMETER (MXV=60) 
      INCLUDE "C.VARY.H"  !  COMMON/VARY/ NV,IR(MXV),NC,I1(MXV),I2(MXV),V(MXV),IS(MXV),W(MXV),
                          ! >IC(MXV),IC2(MXV),I3(MXV),XCOU(MXV),CPAR(MXV,27)
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      PARAMETER (INT1=1) 
      PARAMETER (I0=0) 

      IER=0
C----- CONTROLE VARIABLES
      IF(NV.LT.1) IER=1
      IF(NV.GT.MXV)IER=1
      IF(IER .EQ. 1) THEN
        IF(NRES.GT.0) WRITE(NRES,102) '''NV''',NV,'VARIABLES'
 102    FORMAT(/,5X,' SBR  FITEST - WARNING *  DATA  ',A,
     >  '=',I3,'  IN ',A,'  IS  OUT  OF  RANGE')
      ENDIF
      DO 1 I=1,NV
C------- LMNT VARIABLE
        IF(IR(I).LT.1 .OR. IR(I).GE.NOEL) THEN
          IF(NRES.GT.0) WRITE(NRES,101) '''IR''',IR(I),'VARIABLE',I
 101      FORMAT(/,5X,' SBR  FITEST - WARNING *  data  ',A,'=',
     >    I2,'  in ',A,' #',I3,'  is  out  of  range')
          IER = 1
        ELSE
          IF(NRES.GT.0) WRITE(NRES,*) '          variable # ',I,
     >    '      IR = ',IR(I),',   ok.'
        ENDIF
 
C------- PARAMTRE VARIABLE DANS L'LMNT
        IF(IS(I).LT.1 .OR. IS(I).GT.MXD) THEN
          IF(NRES.GT.0) WRITE(NRES,101) '''IP''',IS(I),'VARIABLE',I
          IER = 1
        ELSE
          IF(NRES.GT.0) WRITE(NRES,*) '          variable # ',I,
     >    '      IP = ',IS(I),',   ok.'
        ENDIF
 
C------- COUPLAGE AVEC LMNT #KL, PRMTR #KP
        KL=INT(XCOU(I))
        IF(KL .NE. 0) THEN
          KP=NINT((1D3*XCOU(I)-1D3*KL))

c                  write(*,*) ' fitest ',i,xcou(i),kl,kp
c                    read(*,*)

          CALL VRBLE(IER,I,KL,KP)
C--------- # LMNT COUPLE
          IF(KL.GE.NOEL) THEN
            IF(NRES.GT.0) WRITE(NRES,101) '''XC.''',KL,'VARIABLE',I
            IER = 1
          ELSE
            IF(NRES.GT.0) WRITE(NRES,*) '          variable # ',I,
     >      '      XC.= ',KL,',   ok.'
          ENDIF
C--------- # PARAMETR COUPLE
          IF(KP .LT. 1 .OR. KP .GT. MXD) THEN
            IF(NRES.GT.0) WRITE(NRES,101) '''.XC''',KP,'VARIABLE',I
            IER = 1
          ELSE
            IF(NRES.GT.0) WRITE(NRES,*) '          variable # ',I,
     >      '      .XC= ',KP,',   ok.'
          ENDIF
        ENDIF
 1    CONTINUE
 
C----- CONTROLE CONTRAINTES
      IF(NC.LT.1 .OR. NC.GT.MXV) THEN
        IF(NRES.GT.0) WRITE(NRES,102) '''NC''',NC,'CONSTRAINTS'
        IER = 1
      ENDIF
 
      CALL FITMM6(I0)

      DO 8 I=1,NC
        IF(I3(I).LT. 1 .OR. I3(I).GE.NOEL ) THEN
          IF(NRES.GT.0) WRITE(NRES,101) '''IR''',I3(I),'CONSTRAINT',I
          IER = 1
        ELSE
          IF(NRES.GT.0) WRITE(NRES,*) '          constraint # ',I,
     >    '      IR = ',I3(I),',   ok.'
        ENDIF
        IF(IC(I) .EQ. 3) THEN
          IF(I1(I) .GT. IMAX ) THEN
            IF(NRES.GT.0) WRITE(NRES,101) '''I''',I1(I),'CONSTRAINT',I
            IER = 1
          ELSE
            IF(NRES.GT.0) WRITE(NRES,*) '          constraint # ',I,
     >      '      I  = ',I1(I),',   ok.'
          ENDIF
        ELSEIF(IC(I) .EQ. 0) THEN
c          IF(IC2(I) .EQ. 0 ) THEN
c            IF(NRES.GT.0) WRITE(NRES,101)'''IC2''',IC2(I),'CONSTRAINT',I
c            IER = 1
c          ENDIF
        ENDIF
        IF(IC(I) .EQ. 7) THEN
C------------ Constraint on  coordinate or field in optical element, 
            CALL INTEG8(INT1)
            CALL FITMM4(I3(I))
        ENDIF
 8    CONTINUE
 
      IF(IER .EQ. 1) THEN
        IF(NRES.GT.0) WRITE(NRES,FMT=
     >  '(/,20X,''** NO  FIT  WILL  BE  PERFORMED **'')')
      ELSE
        IF(NRES.GT.0) WRITE(NRES,FMT=
     >  '(/,20X,''FIT  variables  in  good  order,'',  
     >              ''  FIT  will proceed. '')')
      ENDIF

      IF(SAVFT) THEN
        IF(NRES.GT.0) WRITE(NRES,FMT=
     >  '(/,20X,''Final FIT status will be saved in '',A,/)')
     >  FNAME
      ELSE
        IF(NRES.GT.0) WRITE(NRES,FMT=
     >  '(/,20X,''Final FIT status will NOT be saved. For so, use the'',
     >   '' ''''save [FileName]'''' command'')')
      ENDIF

      RETURN    
      END
