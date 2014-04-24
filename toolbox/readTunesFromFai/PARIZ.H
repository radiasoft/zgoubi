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
C  François Méot <meot@lpsc.in2p3.fr>
C  Service Accélerateurs
C  LPSC Grenoble
C  53 Avenue des Martyrs
C  38026 Grenoble Cedex
C  France
C-- Max mesh nodes in Z direction, # of field coordinates
C      PARAMETER (IZ=9 , ID=3, IZ2=2, ID2=2, ID3=3)
C      PARAMETER (IZ=5 , ID=3, IZ2=2, ID2=2, ID3=3)
C      PARAMETER (IZ=3 , ID=3, IZ2=2, ID2=2, ID3=3)
C      PARAMETER (IZ=1 , ID=1, IZ2=1, ID2=1, ID3=1)

C-- Max mesh nodes in X (longitudinal) and Y (transverse horiz) dirctions
C      PARAMETER (MXX=400, MXY=200)
C      PARAMETER (MXX=3202, MXY=200)
C      PARAMETER (MXX=601, MXY=200)
C      PARAMETER (MXX=2000, MXY=600)

C----- KEK/TOSCA.  Masamitsu Aiba
C      PARAMETER (IZ=41 , ID=3, IZ2=2, ID2=2, ID3=3)
C      PARAMETER (MXX=301, MXY=121)

CC----- PRISM.  Yasutoshi Kuriyama
C      PARAMETER (IZ=35 , ID=3, IZ2=2, ID2=2, ID3=3)
C      PARAMETER (IZ=41 , ID=3, IZ2=2, ID2=2, ID3=3)
C      PARAMETER (MXX=361, MXY=171)

CC----- RACCAM spiral field maps
C       PARAMETER (IZ=21 , ID=3, IZ2=2, ID2=2, ID3=3)
C       PARAMETER (MXX=671, MXY=53)

      PARAMETER (IZ=41 , ID=3, IZ2=2, ID2=2, ID3=3)
      PARAMETER (MXX=1000, MXY=500)
