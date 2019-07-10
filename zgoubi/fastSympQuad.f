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
C  Upton, NY, 11973, USA
C  -------
      SUBROUTINE fastSympQuad(KUASEX, NRES, newIntegQ)
      IMPLICIT NONE
C---------------------------------------------------------
C     (1) Optical elements defined in cartesian frame
C     (2) Fast drift-kick-drift particle motion integrator
C---------------------------------------------------------
      INCLUDE "MAXTRA.H"   ! PARAMETER (MXJ=7)
      INCLUDE "MAXCOO.H"   ! PARAMETER (MXT=10000)
      INCLUDE "MXFS.H"     ! PARAMETRE (MXF=45, MXS=2000, MLF=9, MXP=11)
      INCLUDE "MXLD.H"     ! PARAMETER (MXL=15000 , MXD=1400)
      INCLUDE "MXSCL.H"    ! PARAMETER (MXSCL=10)
      INCLUDE "C.CONST.H"  ! COMMON/CONST/ CL9,CL,PI,RAD,DEG,QE,AMPROT,CM2M
      INCLUDE "C.DON.H"    ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE 'C.SCAL.H'   ! COMMON/SCAL/ SCL(MXF,MXS,MXSCL),TIM(MXF,MXS),NTIM(MXF),KSCL
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"  ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),IREP(MXT),AMQLU,PABSLU
      INCLUDE "C.PTICUL.H" ! COMMON/PTICUL/ AM,Q,G,TO
      INCLUDE 'C.INTEG.H'  ! COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      INCLUDE "C.OBJET.H"  ! COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT,KZOB
      INCLUDE "C.RIGID.H"  ! COMMON/RIGID/ BORO,DPREF,HDPRF,DP,QBR,BRI
      INCLUDE "SCALE_symp.H"  ! PARAMETER (scl_l0, scl_w0)
      INCLUDE "C.UNITS.H"    ! COMMON/UNITS/ UNIT(MXJ)

      DOUBLE PRECISION XL, XE, RO, GAP, SCAL, SCAL0, B_mag, b2_mag
      DOUBLE PRECISION DLE, XLS, DLS, DL0, TCUM, SCUM, X_axis(MXT)
      DOUBLE PRECISION FINTE, FINTS, PREF, EREF, BETA0, GAMMA0
      DOUBLE PRECISION eb2, gb0, PAS0, PAS_scl, RO_scl
      INTEGER KUASEX, NRES, newIntegQ, MAX_STEP, J, IT

C-------- Get the parameters for the magnetic qradrupole

      XL = A(NOEL,10)
      RO = A(NOEL,11)
      GAP = RO/KUASEX
      SCAL = SCAL0()
      B_mag = A(NOEL,12)*SCAL
      XE = A(NOEL,20)
      DLE = A(NOEL,21)
      XLS = A(NOEL,40)
      DLS =A(NOEL,41)
      PAS = A(NOEL,60)
      KP = NINT(A(NOEL,70))
      XCE = A(NOEL,71)
      YCE = A(NOEL,72)
      ALE = A(NOEL,73)

      b2_mag = B_mag/RO    ! The strength of the mag. quad.
      IF(NRES.GT.0) THEN
        WRITE(NRES, 100) XL, RO, B_mag, B_mag/SCAL, b2_mag
      ENDIF

      DL0 = DLE + DLS
      IF(DL0 .EQ. 0.D0) THEN
C-------- Sharp edge at entrance and exit
        FINTE = XE
        XE = 0.D0
        FINTS = XLS
        XLS = 0.D0
        IF(NRES.GT.0) THEN
          WRITE(NRES,102) FINTE, FINTS, GAP
        ENDIF
      ENDIF

      PAS0 = PAS           ! Integration step from input
      MAX_STEP = CEILING(XL/PAS0)
      PAS = XL/MAX_STEP    ! Adjust the step size to fit the length
      IF(NRES.GT.0) THEN
        WRITE(NRES,103) PAS0, PAS
      ENDIF

 100  FORMAT(/, 5X, ' -----  QUADRUPOLE  : ', 1P
     >  , /, 15X,' Length  of  element  = ',G16.8,'  cm'
     >  , /, 15X,' Bore  radius      RO = ',G13.5,'  cm'
     >  , /, 15X, 'B-QUADRUPOLE  =', 1P, E15.7, 1X, 'kG'
     >  , '  (i.e., ', E15.7, ' * SCAL)'
     >  , /, 15X, 'The strength of the magnetic quadrupole b2 = '
     >  , E15.7, ' kG/cm')
 102  FORMAT(/, 15X, 'Entrance/exit field models are sharp edge',
     >  /, 15X, 'FINTE, FINTS, gap : ', 1P, 3(1X, E12.4))
 103  FORMAT(
     >  /, 20X, 'The input  integration step :', 1P, G15.8, ' cm',
     >  /, 20X, 'The ACTUAL integration step :', 1P, G15.8, ' cm', /)

C---------Compute some constants for the sympletic integration
      IF(Q .EQ. 0.D0) Q = 1
      PREF = BORO*CL9*Q ! the reference momentum; BORO in kG.cm
                        ! CL9 = 0.29979246; Q=1 for proton;
      EREF = SQRT(PREF*PREF + AM*AM) ! the total energy
      BETA0 = PREF / EREF
      GAMMA0 = EREF / AM
      gb0 = 1.D0/(BETA0*GAMMA0)
      PAS_scl = PAS/scl_l0
      RO_scl = RO/scl_l0
      eb2 = ((B_mag*RO)/BORO)/(RO_scl*RO_scl)

      IF (newIntegQ .EQ. 1) THEN
C---------Drift-kick-drift motion integrator
        CALL DKD_INTEGR(MAX_STEP, PAS_scl,
     >           gb0, eb2, PREF, IMAX, AM, F, X_axis)
      ELSE IF (newIntegQ .EQ. 2) THEN
C---------Matrix-kick-Matrix motion integrator
        CALL MKM_INTEGR(MAX_STEP, PAS_scl,
     >           gb0, eb2, PREF, IMAX, AM, F, X_axis)
      ELSE
        WRITE (*,'(/, 10X, A, I0, A, /)') 'newIntegQ = ',
     >    newIntegQ, ', NOT implemented, stop the job.'
        IF (NRES .GT. 0) WRITE(NRES, '(/, 10X, A, I0, A, /)')
     >    'INTEGRATOR OPTION = ', newIntegQ,
     >    '. It is NOT implemented, stop the job.'
        CALL ENDJOB(' The choice of integrator has NOT been
     >implemented', -99)
      END IF

      DO IT = 1, IMAX
        IF(NRES.GT.0) THEN
          IF(IT.EQ.1) THEN
            WRITE(NRES,199)
          ENDIF
C----------NOTE: P and T in unit of mrad (as in F(J,I))
          WRITE(NRES,200) IEX(IT), (FO(J,IT), J=1,5),
     >      X_axis(IT), (F(J,IT), J=2,5), IT
        ENDIF

      ENDDO

      SCUM = F(6,1)
      TCUM = F(7,1)
      IF(NRES .GT. 0) THEN
        WRITE(NRES,201) KP, XCE, YCE, ALE, SCUM*UNIT(5), TCUM
      ENDIF

 199  FORMAT(2X, '  KPOS  DP         Y(cm)   T(mrad)  ',
     >  '   Z(cm)   P(mrad)  ', 8X, ' X(cm)         Y(cm)      ',
     >  '  T(mrdd)      Z(cm)        P(mrad)   ', 8X, '  I')
 200  FORMAT(2X, 'A', 2X, I3, F8.4, 4F10.3, 8X, F12.6, 4F13.6, 8X, I5)
 201  FORMAT(/, 5X, 'KPOS =  ', I0,
     >  '.  Change  of  frame  at  exit  of  element.', /, 10X,
     >  'X =', 1P, G12.4, ' CM   Y =', G12.4, 1P,
     >  'cm,  tilt  angle =', G14.6, ' RAD', /, /, /, 1P,
     >  'Cumulative length of optical axis = ', 1P, G17.9, ' m ;  ',
     >  'Time  (for ref. rigidity & particle) = ', 1P, G14.6, ' s')

      RETURN
      END
