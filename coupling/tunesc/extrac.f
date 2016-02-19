c     extrac.f
c
c     written by Frédéric Desforges, frederic.desforges@grenoble-inp.org
c
c
      SUBROUTINE EXTRAC(NRES,ARCLEN,R,rPARAM,C,NU1,NU2,ALPHA1,ALPHA2,B
     >ETA1,BETA2,GAMMA1,GAMMA2,CMOINS,CPLUS,DELTA,DELTA2,NUX0,NUY0,P)

      IMPLICIT DOUBLE PRECISION (A-H,M-Z)
      DIMENSION R(4,4),C(2,2),P(4,4)
      INTEGER NRES

      dum = cplus

      WRITE(NRES,FMT='(/,''-----------------------------------------'',/
     >)')
      
      WRITE(NRES,FMT='(''POSITION (IN M): '',F13.8)') ARCLEN/100.D0
      
      WRITE(NRES,FMT='(/,40X,''########################################'
     >')')
      
      WRITE(NRES,FMT='(40X,''#  EDWARDS AND TENG`S PARAMETRIZATION  #'')
     >')
      
      WRITE(NRES,FMT='(40X,''########################################'',
     >/,/)')

      WRITE(NRES,FMT='(/,6X,''TRANSFERT MATRIX (Rij) IN THE COUPLED FRAM
     >E'',/)')
      WRITE(NRES,100) ((R(I,J),J=1,4),I=1,4)
 100     FORMAT(6X,4F13.8)

      WRITE(NRES,FMT='(/,/,6X,''COUPLING PARAMETERS:'')')
      
      WRITE(NRES,FMT='(27X,''- r:'',2X,F13.8)') rPARAM
      
      WRITE(NRES,FMT='(27X,''- C:'')')
      
      WRITE(NRES,200) ((C(I,J),J=1,2),I=1,2)
 200     FORMAT(2(30X,2(F15.8,1X),/))

      WRITE(NRES,FMT='(/,78X,''MODE 1'',15X,''MODE 2'',/)')
      
      WRITE(NRES,FMT='(6X,''FRACTIONAL PART OF THE BETATRON TUNES IN THE
     > DECOUPLED FRAME: '',4X,F13.8,8X,F13.8)') NU1,NU2
      
      WRITE(NRES,FMT='(/,6X,''EDWARDS-TENG`S PARAMETERS:'')')
      
      WRITE(NRES,FMT='(32X,''- ALPHA:'',33X,F13.8,8X,F13.8)') ALPHA1,ALP
     >HA2
      
      WRITE(NRES,FMT='(32X,''- BETA:'',34X,F13.8,8X,F13.8)') BETA1,BETA2
      
      WRITE(NRES,FMT='(32X,''- GAMMA:'',33X,F13.8,8X,F13.8)') GAMMA1,GAM
     >MA2

      WRITE(NRES,FMT='(/,/,6X,''HAMILTONIAN PERTURBATION PARAMETERS:'',/
     >)')

      WRITE(NRES,FMT='(27X,''- DISTANCE FROM THE NEAREST DIFFERENCE LINE
     >AR RESONANCE:'',4X,F11.8)') DELTA

      WRITE(NRES,FMT='(27X,''- COUPLING STRENGHT OF THE DIFFERENCE LINEA
     >R RESONANCE:'',3X,F13.8,/)') CMOINS

      WRITE(NRES,FMT='(27X,''- DISTANCE FROM THE NEAREST SUM LINEAR RESO
     >NANCE:'',11X,F11.8)') DELTA2

c      WRITE(NRES,FMT='(27X,''- COUPLING STRENGHT OF THE SUM LINEAR RESON
c     >ANCE:'',10X,F13.8,/)') CPLUS

      WRITE(NRES,FMT='(27X,''- UNPERTURBED HORIZONTAL TUNE:'',28X,F13.8)
     >') NUX0

      WRITE(NRES,FMT='(27X,''- UNPERTURBED VERTICAL TUNE:'',30X,F13.8,/,
     >/)') NUY0

      WRITE(NRES,FMT='(/,6X,''P MATRIX :'',/)')
      WRITE(NRES,188) ((P(I,J),J=1,4),I=1,4)
 188  FORMAT(1p,4(2X,E18.10))
      
      
      END
