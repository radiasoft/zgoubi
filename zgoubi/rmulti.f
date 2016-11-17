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
      SUBROUTINE RMULTI(NDAT,NOEL,MXL,A,MPOL,
     >                                       ND)
C     ------------------------
C     READS DATA FOR MULTIPOLE
C     ------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(MXL,*)
      PARAMETER (KSIZ=10)
      CHARACTER(KSIZ) KLE

C-- eRHIC, skew ffag dipoles
c          logical first
c             data first / .true. /
C-- eRHIC, skew ffag dipoles

      IA = 1
      LINE =1 
      READ(NDAT,*,END=90,ERR=90) IL
      A(NOEL,IA) = IL

      IA = IA+1                
      IB = IA+MPOL+1           
      LINE = LINE + 1 
      READ(NDAT,*,END=90,ERR=90) (A(NOEL,I),I=IA,IB)

      IA = IB+1                
C----- CHP FUITE ENTREE
C        XE, LambdaE, B1-Bmpol
      IB = IA + MPOL           
      LINE = LINE + 1 
      READ(NDAT,*,ERR=11) (A(NOEL,I),I=IA,IB)
 11   CONTINUE
      IA = IB+1                
      IB = IA + 6              
      LINE = LINE + 1 
      READ(NDAT,*,END=90,ERR=90) II,(A(NOEL,I),I=IA+1,IB)
      A(NOEL,IA) = II
C----- CHP FUITE SORTIE
      IA = IB + 1              
      IB = IA + MPOL           
      LINE = LINE + 1 
      READ(NDAT,*,ERR=12) (A(NOEL,I),I=IA,IB)
 12   CONTINUE
      IA = IB+1                
      IB = IA + 6              
      LINE = LINE + 1 
      READ(NDAT,*,END=90,ERR=90) II,(A(NOEL,I),I=IA+1,IB)
      A(NOEL,IA) = II
C----- Rotation of multipole components
      IA = IB +1               
      IB = IA + MPOL - 1       
      LINE = LINE + 1 
      READ(NDAT,*,END=90,ERR=90) (A(NOEL,I),I=IA,IB)

      ND = IB + 1              
      CALL STPSIZ(NDAT,NOEL,ND,
     >                         A)

c           write(*,*) ' multip nd ',nd

C Modif, FM, Dec. 05
C      READ(NDAT,*,END=90,ERR=90) II,(A(NOEL,I),I=ND+2,ND+4)
C      A(NOEL,ND+1) = II
      LINE = LINE + 1 
      READ(NDAT,*,END=90,ERR=90) II,(A(NOEL,I),I=ND+4,ND+6)
      A(NOEL,ND+3) = II 
 
c          write(*,*) ' multip  ',ND+3,A(NOEL,ND+3), II 
c          write(*,*) ' rfit '
c          write(*,fmt='(1p,10(e12.4,1x))') (a(noel,i),i=1,70)  
c                 stop

      RETURN

 90   CONTINUE
      CALL ZGKLEY( 
     >            KLE)
      CALL ENDJOB('*** Pgm rmulti, keyword '//KLE//' : '// 
     >'input data error, at line ',LINE)
      RETURN
      END
