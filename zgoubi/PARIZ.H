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

C MMAP : The total number of field maps declared in zgoubi.dat, whether repeatedly for some 
C of them or not.     
C IZ :   
C ID :   

C-- Sizing of the field map array HC(ID,MXX,MXY,IZ,MMAP)
C-- Max mesh nodes in Z direction, # of field coordinates
C      PARAMETER (IZ=9 , ID=3)
C      PARAMETER (IZ=5 , ID=3)
C      PARAMETER (IZ=3 , ID=3)

C-- Max mesh nodes in X (longitudinal) and Y (transverse horiz) dirctions
C      PARAMETER (MXX=400, MXY=200)
C      PARAMETER (MXX=3202, MXY=200)
C      PARAMETER (MXX=601, MXY=200)
C      PARAMETER (MXX=2000, MXY=600)

C-- Sizing of the field map array HC(ID,MXX,MXY,IZ,MMAP)
CC----- PRISM.  Yasutoshi Kuriyama
C      PARAMETER (IZ=35 , ID=3)
C      PARAMETER (IZ=41 , ID=3)
C      PARAMETER (MXX=361, MXY=171)

C-- Sizing of the field map array HC(ID,MXX,MXY,IZ,MMAP)
CC----- RACCAM spiral field maps
C       PARAMETER (IZ=21 , ID=3, MMAP=1)
C       PARAMETER (MXX=671, MXY=53)

C-- Sizing of the field map array HC(ID,MXX,MXY,IZ,MMAP)
C      PARAMETER (IZ=41 , ID=3, MMAP=1)
C      PARAMETER (MXX=1000, MXY=300)

C-- Sizing of the field map array HC(ID,MXX,MXY,IZ,MMAP)
CC---- EMMA 
C      PARAMETER (IZ=5 , ID=3, MMAP=1)
C      PARAMETER (MXX=401, MXY=300)

C-- Sizing of the field map array HC(ID,MXX,MXY,IZ,MMAP)
CC---- EMMA. Two 2D maps, Ben, Feb 2013 
c      PARAMETER (IZ=3 , ID=3, MMAP=2)
c      PARAMETER (MXX=961, MXY=161)
cccccccccccc      PARAMETER (IZ=3 , ID=3, MMAP=1)
cccccccccccc      PARAMETER (MXX=961, MXY=161)

C-- Sizing of the field map array HC(ID,MXX,MXY,IZ,MMAP)
C-- IZ : # of mesh nodes in Z direction or 3 with MAP2. ID : # of field coordinates (3 if IZ>1)
C-- MMAP : maximum number of maps in a problem
C      PARAMETER (IZ=3 , ID=3, MMAP=2)
C      PARAMETER (MXX=801, MXY=29)

C-- Sizing of the field map array HC(ID,MXX,MXY,IZ,MMAP)
CC----- RACCAM spiral field maps
c       PARAMETER (IZ=3 , ID=3, MMAP=1)
c       PARAMETER (MXX=361, MXY=221)

C-- Sizing of the field map array HC(ID,MXX,MXY,IZ,MMAP)
CC----- RACCAM spiral field maps
C       PARAMETER (IZ=21 , ID=3, MMAP=1)
C       PARAMETER (MXX=671, MXY=53)

C-- Sizing of the field map array HC(ID,MXX,MXY,IZ,MMAP)
C----- KEK/TOSCA.  Masamitsu Aiba
c      PARAMETER (IZ=41 , ID=3, MMAP=1)
c      PARAMETER (MXX=301, MXY=121)

C-- Sizing of the field map array HC(ID,MXX,MXY,IZ,MMAP)
C----- Jaroslaw problem, Oct 2010
c      PARAMETER (IZ=3 , ID=3, MMAP=1)
c      PARAMETER (MXX=1271, MXY=97)

C AGS ----------------------------------------------------
CC AGS with 2 (mmap -> +2) 2D snake maps (using MAP2D)  
CC and using AGS dipole 2D maps (mmap -> +6) (using TOSCA)  
C      PARAMETER (IZ=3 , ID=3, MMAP=8)
C      PARAMETER (MXX=801, MXY=161)

CC D5 IPM magnet in bare AGS.
c      PARAMETER (IZ=41 , ID=3, MMAP=2)
c      PARAMETER (MXX=1801, MXY=41)

C-- Sizing of the field map array HC(ID,MXX,MXY,IZ,MMAP)
C AGS warm snake 
c       PARAMETER (IZ=29 , ID=3, MMAP=1)
c       PARAMETER (MXX=801, MXY=29)
cccccccccccccccc      PARAMETER (IZ=61 , ID=3, MMAP=4)
c      PARAMETER (MXX=1001, MXY=150)

C-- Sizing of the field map array HC(ID,MXX,MXY,IZ,MMAP)
C AGS warm snake 
c       PARAMETER (IZ=29 , ID=3, MMAP=1)
c       PARAMETER (MXX=801, MXY=29)

CC AGS with (mmap -> +2) 3D snake maps (using TOSCA)  
CC and using AGS dipole 2D maps (mmap -> +6) (using TOSCA)  
c        PARAMETER (IZ=61 , ID=3, MMAP=5)
c      PARAMETER (MXX=1001, MXY=61)
c      PARAMETER (MXX=801, MXY=201)
c      PARAMETER (MXX=321, MXY=29)
c       compatible eRHIC Halbach magnet
c      PARAMETER (MXX=1601, MXY=200)
C      PARAMETER (MXX=991, MXY=161)   ! compatible EMMA
C AGS ----------------------------------------------------

CC Cyl. condenser end, pEDM, MAP2D-E
c      PARAMETER (IZ=3 , ID=3, MMAP=1)
c      PARAMETER (MXX=71, MXY=151)

CC Daedalus cyclo
c      PARAMETER (IZ=3 , ID=3, MMAP=6)
c      PARAMETER (MXX=596, MXY=596)

CC-- Sizing of the field map array HC(ID,MXX,MXY,IZ,MMAP)
CC---- GSI KAOS
c      PARAMETER (IZ=3, ID=3, MMAP=2)
c      PARAMETER (MXX=47, MXY=37) 

CC Laurent, DESIR, tmag field maps
C Using TOSCA
C      PARAMETER (IZ=3 , ID=3, MMAP=8)
C      PARAMETER (MXX=651, MXY=361)

CC FFAG 14 tutorial - POLARMES with PSI field map
c      PARAMETER (IZ=41 , ID=3, MMAP=2)
c      PARAMETER (MXX=351, MXY=200)

CC ~/zgoubi/SVN/current/exemples/tools
c      PARAMETER (IZ=61 , ID=3, MMAP=3)
c      PARAMETER (MXX=801, MXY=200)

CC Malek FFAG 150MeV
C       PARAMETER (IZ=21 , ID=3, MMAP=2)
C       PARAMETER (MXX=301, MXY=319)

C-- Sizing of the field map array HC(ID,MXX,MXY,IZ,MMAP)
C Carols 1GeV FFAG
c       PARAMETER (IZ=3 , ID=3, MMAP=4)
c       PARAMETER (MXX=370, MXY=200)
C       PARAMETER (IZ=3 , ID=3, MMAP=4)
C       PARAMETER (MXX=370, MXY=200)

C-- Sizing of the field map array HC(ID,MXX,MXY,IZ,MMAP)
C For ~/zgoubi/SVN/zgoubi-code/exemples/tools/compare
C       PARAMETER (IZ=29 , ID=3, MMAP=6)      ! strictu sensu
       PARAMETER (IZ=31 , ID=3, MMAP=6)    ! compatible Cbeta
       PARAMETER (MXX=801, MXY=200)

C-- Sizing of the field map array HC(ID,MXX,MXY,IZ,MMAP)
C For ~/zgoubi/SVN/zgoubi-code/exemples/tools/compare
       PARAMETER (IZ=31 , ID=3, MMAP=6)
       PARAMETER (MXX=801, MXY=200)

C-- Sizing of the field map array HC(ID,MXX,MXY,IZ,MMAP)
C For Scott/Holger Cbeta iron quad map
c       PARAMETER (IZ=9 , ID=3, MMAP=2)
c       PARAMETER (MXX=447, MXY=87)

C-- Sizing of the field map array HC(ID,MXX,MXY,IZ,MMAP)
C For CBeta - Nick
c       PARAMETER (IZ=31 , ID=3, MMAP=2)
c       PARAMETER (MXX=801, MXY=31)


C-- Sizing of the field map array HC(ID,MXX,MXY,IZ,MMAP)
C For Scott/Holger Cbeta iron quad map
c       PARAMETER (IZ=9 , ID=3, MMAP=2)
c       PARAMETER (MXX=447, MXY=87)

C-- Sizing of the field map array HC(ID,MXX,MXY,IZ,MMAP)
C For CBeta - Nick
c       PARAMETER (IZ=31 , ID=3, MMAP=2)
c       PARAMETER (MXX=801, MXY=31)

 
