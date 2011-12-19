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

C-- Sizing of the field map array HC(ID,MXX,MXY,IZ)
C-- Max mesh nodes in Z direction, # of field coordinates
C      PARAMETER (IZ=9 , ID=3)
C      PARAMETER (IZ=5 , ID=3)
C      PARAMETER (IZ=3 , ID=3)

C-- Max mesh nodes in X (longitudinal) and Y (transverse horiz) dirctions
C      PARAMETER (MXX=400, MXY=200)
C      PARAMETER (MXX=3202, MXY=200)
C      PARAMETER (MXX=601, MXY=200)
C      PARAMETER (MXX=2000, MXY=600)

C-- Sizing of the field map array HC(ID,MXX,MXY,IZ)
CC----- PRISM.  Yasutoshi Kuriyama
C      PARAMETER (IZ=35 , ID=3)
C      PARAMETER (IZ=41 , ID=3)
C      PARAMETER (MXX=361, MXY=171)

C-- Sizing of the field map array HC(ID,MXX,MXY,IZ)
CC----- RACCAM spiral field maps
C       PARAMETER (IZ=21 , ID=3, MMAP=1)
C       PARAMETER (MXX=671, MXY=53)

C-- Sizing of the field map array HC(ID,MXX,MXY,IZ)
C      PARAMETER (IZ=41 , ID=3)
C      PARAMETER (MXX=1000, MXY=300)

C-- Sizing of the field map array HC(ID,MXX,MXY,IZ)
CC---- EMMA 
C      PARAMETER (IZ=5 , ID=3, MMAP=1)
C      PARAMETER (MXX=401, MXY=300)

C-- Sizing of the field map array HC(ID,MXX,MXY,IZ,MMAP)
C-- IZ : # of mesh nodes in Z direction or 3 with MAP2. ID : # of field coordinates (3 if IZ>1)
C-- MMAP : maximum number of maps in a problem
C      PARAMETER (IZ=3 , ID=3, MMAP=2)
C      PARAMETER (MXX=801, MXY=29)

C-- Sizing of the field map array HC(ID,MXX,MXY,IZ)
CC----- RACCAM spiral field maps
c       PARAMETER (IZ=1 , ID=3, MMAP=1)
c       PARAMETER (MXX=361, MXY=221)

C-- Sizing of the field map array HC(ID,MXX,MXY,IZ)
CC----- RACCAM spiral field maps
C       PARAMETER (IZ=21 , ID=3, MMAP=1)
C       PARAMETER (MXX=671, MXY=53)

C-- Sizing of the field map array HC(ID,MXX,MXY,IZ)
C----- KEK/TOSCA.  Masamitsu Aiba
c      PARAMETER (IZ=41 , ID=3, MMAP=1)
c      PARAMETER (MXX=301, MXY=121)

C-- Sizing of the field map array HC(ID,MXX,MXY,IZ)
C----- Jaroslaw problem, Oct 2010
c      PARAMETER (IZ=1 , ID=3, MMAP=1)
c      PARAMETER (MXX=1271, MXY=97)

C AGS with 2 (=mmap) snake maps
C Using MAP2D
C      PARAMETER (IZ=1 , ID=3, MMAP=2)
C      PARAMETER (MXX=801, MXY=29)
      

C AGS Magnet maps
C Using TOSCA
C      PARAMETER (IZ=1 , ID=3, MMAP=2)
C      PARAMETER (MXX=498, MXY=161)
      
C AGS with 2 (mmap -> +2) 2D snake maps (using MAP2D)  
C and using AGS dipole 2D maps (mmap -> +6) (using TOSCA)  
C      PARAMETER (IZ=1 , ID=3, MMAP=8)
C      PARAMETER (MXX=1000, MXY=201)

C AGS with 2 (mmap -> +2) 2D snake maps (using MAP2D)  
C and using AGS dipole 2D maps (mmap -> +6) (using TOSCA)  
      PARAMETER (IZ=30 , ID=3, MMAP=8)
      PARAMETER (MXX=1000, MXY=201)

C-- Sizing of the field map array HC(ID,MXX,MXY,IZ)
CC---- GSI KAOS
C      PARAMETER (IZ=5 , ID=3, MMAP=2)
C      PARAMETER (MXX=401, MXY=300)

