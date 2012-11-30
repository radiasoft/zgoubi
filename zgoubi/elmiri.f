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
      SUBROUTINE ELMIRI(SCAL,
     >                       XL,DEV)
C Mirror, straight slits. Frame is Cartesian
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/AIM/ BO,RO,FG,GF,XI,XF,EN,EB1,EB2,EG1,EG2
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      PARAMETER(MCOEF=6)
      COMMON/CHAFUI/ XE,XS,CE(MCOEF),CS(MCOEF),QCE(MCOEF),QCS(MCOEF)
C      COMMON/CHAFUI/ XE,XS,CE(6),CS(6),QCE(6),QCS(6)
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QEL,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      COMMON/DROITE/ CA(9),SA(9),CM(9),IDRT
      COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      PARAMETER(MPOL=10)
      COMMON/MULTPE/ EM(MPOL),QLE(MPOL),QLS(MPOL)
     >,QE(MPOL,MCOEF),QS(MPOL,MCOEF),RTQ(MPOL)
      LOGICAL ZSYM
      COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
C      COMMON/PTICUL/ AM,Q,G,TO
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI

      DIMENSION AL(7), V(7)

      CHARACTER SMT(2)*6, HV(2,2)*10
      DATA SMT /          'mirror',                'lens' /
      DATA HV / 'horizontal', ' vertical', ' vertical', 'horizontal' /

      NP =A(NOEL,10)
      DO 10 I=1, NP
        AL(I) = A(NOEL,10+I)
        V(I) = A(NOEL,19+I)*SCAL
 10   CONTINUE

      D = A(NOEL,10+NP+1)
      MT = A(NOEL,10+NP+2)

      KP = NINT(A(NOEL,40))
      IF( KP .EQ. 3 ) THEN
        IF(A(NOEL,43) .EQ. 0.D0) THEN
C          DEV = 
        ELSE
          DEV = -A(NOEL,73) * 2.D0
        ENDIF
      ENDIF

C------ Normally L1 ~ 3D
C----- MT = 11, 12 for mirror H, V
C           21, 22 for lense V, H
      IF(MT/10.NE.1 .AND. MT/10.NE.2) 
     >        CALL ENDJOB('Wrong data MT',-99)
      IF(MT-10*(MT/10).NE.1 .AND. MT-10*(MT/10).NE.2) 
     >        CALL ENDJOB('Wrong data MT',-99)

      IF(NRES.GT.0) THEN
        ML=MT/10
        KPL=MT-10*(MT/10)
        WRITE(NRES,110) HV(ML,KPL), SMT(ML),NP
 110    FORMAT(/,40X,' --- Electrostatic mirror/lens ---'
     >        ,/,40X,'     * used as ',A,' ',A,' *', /,    1P
     >        ,/,30X,' Number of plates (N) = ',I2
     >        ,/,30X,' Length of plates  (Li) :')
        WRITE(NRES,111) (I,AL(I),I=1,NP)
 111    FORMAT(  30X,'         ',I2,' :',5X,1P,G12.5,' m')
        WRITE(NRES,FMT='(
     >           30X,'' Gap                  = '',1P,G12.5,'' m''
     >        ,/,30X,'' Potentials (Vi) :'')') D
        WRITE(NRES,112) (I,V(I),I=1,NP)
 112    FORMAT(  30X,'         ',I2,' :',5X,1P,G12.5,' V')

      ENDIF
 
C----- change units : AL,D -> cm  and V1-3 -> MeV
C----- V -> V/D, MeV/cm
      D=D*100.D0
      XL = 0.D0
      DO 11 I=1, NP
        AL(I) = AL(I)*100.D0
        XL = XL + AL(I)
 11     V(I) = V(I)*1.D-6/D

      XI = -AL(1)
      XLIM= XL + XI
      XF=XLIM
      XS = XF

      IF(MT-10*(MT/10).EQ.1) THEN
C----- MT=11 or 21. 
C----- plates are parallel to H-plane -> h-mirror or v-lens
        EM(6) = 0.D0
      ELSEIF(MT-10*(MT/10).EQ.2) THEN
C----- MT=12 or 22. 
C----- plates are parallel to V-plane -> v-mirror or h-lens
        EM(6) = 0.5D0*PI
      ENDIF

      IF(SMT(MT/10).EQ.'mirror') THEN
        IDRT = 2
        CA(1)=1.D0
        SA(1)=0.D0
        CM(1)= -XI
        CA(2)=1.D0
        SA(2)=0.D0
        CM(2)= -XI
        CALL TRANSW(.TRUE.,.TRUE.)
      ENDIF
      CALL ELMIRW(D,AL,V,NP)
      RETURN
      END
