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

C     written by Frédéric Desforges (2013)

      SUBROUTINE MATIMC(NRES)

      IMPLICIT DOUBLE PRECISION (A-H,M-Z)
      DIMENSION C(2,2),P(4,4)
      INTEGER NRES
      DIMENSION X4(2,2), X19(4,4)

      SAVE RPARAM,C,NU1,NU2,ALPHA1, ALPHA2, BETA1,    
     >BETA2,GAMMA1,GAMMA2,CMOINS,CPLUS,DELTA,DELTA2,NUX0,NUY0,P
      
      IF(NRES .GT. 0) THEN
      WRITE(NRES,FMT='(//,40X,
     >                     ''----------------------------------------''
     >)')
      WRITE(NRES,FMT='(40X,''-  EDWARDS AND TENG`S PARAMETRIZATION  -''
     >)')
      WRITE(NRES,FMT='(40X,''----------------------------------------''
     >,//)')

      WRITE(NRES,FMT='(6X,''COUPLING PARAMETERS:'')')
      
      WRITE(NRES,FMT='(27X,''- r:'',2X,F13.8)') rPARAM
      
      WRITE(NRES,FMT='(27X,''- C:'')')
      
      WRITE(NRES,200) ((C(I,J),J=1,2),I=1,2)
 200     FORMAT(2(30X,2(F15.8,1X),/))

      WRITE(NRES,FMT='(78X,''MODE 1'',15X,''MODE 2'',/)')
      
      WRITE(NRES,FMT='(6X,''FRACTIONAL PART OF THE BETATRON TUNES IN THE
     > DECOUPLED FRAME: '',4X,F13.8,8X,F13.8)') NU1,NU2
      
      WRITE(NRES,FMT='(/,6X,''EDWARDS-TENG`S PARAMETERS:'')')
      
      WRITE(NRES,FMT='(32X,''- ALPHA:'',33X,F13.8,8X,F13.8)') ALPHA1,
     >ALPHA2
      
      WRITE(NRES,FMT='(32X,''- BETA:'',34X,F13.8,8X,F13.8)') BETA1,BETA2
      
      WRITE(NRES,FMT='(32X,''- GAMMA:'',33X,F13.8,8X,F13.8)') GAMMA1,
     >GAMMA2

      WRITE(NRES,FMT='(/,6X,''HAMILTONIAN PERTURBATION PARAMETERS:'',/
     >)')

      WRITE(NRES,FMT='(27X,''- DISTANCE FROM THE NEAREST DIFFERENCE LINE
     >AR RESONANCE:'',4X,F11.8)') DELTA

      WRITE(NRES,FMT='(27X,''- COUPLING STRENGTH OF THE DIFFERENCE LINEA
     >R RESONANCE:'',3X,F13.8,/)') CMOINS

      WRITE(NRES,FMT='(27X,''- DISTANCE FROM THE NEAREST SUM LINEAR RESO
     >NANCE:'',11X,F11.8)') DELTA2

c      WRITE(NRES,FMT='(27X,''- COUPLING STRENGTH OF THE SUM LINEAR RESON
c     >ANCE:'',10X,F13.8,/)') CPLUS

      WRITE(NRES,FMT='(27X,''- UNPERTURBED HORIZONTAL TUNE:'',28X,F13.8)
     >') NUX0

      WRITE(NRES,FMT='(27X,''- UNPERTURBED VERTICAL TUNE:'',30X,F13.8,/)
     >') NUY0

      WRITE(NRES,FMT='(6X,''P MATRIX :'',/)')
      WRITE(NRES,188) ((P(I,J),J=1,4),I=1,4)
 188  FORMAT(1p,4(2X,E18.10))
      END IF

      RETURN

      ENTRY MATIC2(X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15,X16,
     >X17, X18, X19)    

      RPARAM=X3
      C=X4
      NU1=X5
      NU2=X6
      ALPHA1=X7
      ALPHA2=X8
      BETA1=X9
      BETA2=X10
      GAMMA1=X11
      GAMMA2=X12
      CMOINS=X13
      CPLUS=X14
      DELTA=X15
      DELTA2=X16
      NUX0=X17
      NUY0=X18
      P=X19      

      RETURN      
      END
