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
C  Brookhaven National Laboratory                                               és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE MATCH(FNCT,DFNCT,NV,V,VMI,VMA,NC,CX,CY,MORE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(*),VMI(*),VMA(*), CX(*), CY(*)
      LOGICAL MORE

      PARAMETER(MXV=10)
      PARAMETER(MXC=400)
      DIMENSION RNU(MXV),RO(MXV),FCT(MXC)
      EXTERNAL FNCT,DFNCT

 1    CONTINUE

C------- INITIALISATION DE L'ALGORYTHME
        CALL FUNK(FNCT,V,CX(MXC-1),CX(MXC),NC,NV,CX,FCT)
        CALL GRAD(FNCT,DFNCT,V,CX,CY,NC,NV,FCT,RO)
        DO I=1,NV
          RNU(I)=-RO(I)
        ENDDO

C------- ALGORYTHME DE FLETCHER-REEVES
        DO 50 K=1,NV
          RLAMDA=AMINI(FNCT,CX,CY,V,VMI,VMA,NC,NV,RNU)      

          DO 52 I=1,NV
 52         V(I)=V(I)+RLAMDA*RNU(I) 

          RMOD1=0.D0 
          DO 55 I=1,NV
 55         RMOD1=RMOD1+RO(I)*RO(I)

          CALL FUNK(FNCT,V,CX(MXC-1),CX(MXC),NC,NV,CX,FCT)
          CALL GRAD(FNCT,DFNCT,V,CX,CY,NC,NV,FCT,RO)

          RMOD2=0.D0
          DO 56 I=1,NV
            RMOD2=RMOD2+RO(I)*RO(I)
 56       CONTINUE     

          DO 57 I=1,NV
 57         RNU(I)=-RO(I)+RNU(I)*RMOD2/RMOD1

 50     CONTINUE

        RNORMP=0.D0
        DO 31 I=1,NV
 31       RNORMP=RNORMP+V(I)*V(I)

        RLAMDA=AMINI(FNCT,CX,CY,V,VMI,VMA,NC,NV,RNU)      

        DO 32 I=1,NV
 32       V(I)=V(I)+RLAMDA*RNU(I) 

C------- TEST D'ARRET
        RNORM =0.D0
        DO 33 I=1,NV
 33       RNORM =RNORM + V(I)* V(I)

        TEST=1.D0-RNORM/RNORMP 
        TEST=TEST*TEST

C      IF(TEST .GT. 1.D-20) GOTO 1   
      IF(TEST .GT. 1.D-12) GOTO 1   

C------------ Feedback on the fit, for adjusting the 
C            coefficients so as to have the magnetic 
C            face for 1/(1+exp(P(x)) at x=0 (while the
C            initial value XFB is just an approximate value
C            computed from the data) :
       IF(MORE) CALL FMAG(FNCT,CX,NC,NV,V,
     >                                    *1)

      RETURN
      END
