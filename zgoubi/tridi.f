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
      SUBROUTINE TRIDI(L,D,U,COE,nd,n)
C TRIDIAGONAL MATRIX
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer m,i
      doubleprecision L(nd),D(nd),U(nd),COE(nd)
      m = n - 1 
      do 10 i = 1,m
         L(i+1) = L(i+1) / D(i)
         D(i+1) = D(i+1) - L(i+1) *U(i)
         COE(i+1) = COE(i+1) - L(i+1) *COE(i)
10    continue
C  THE COEFFICIENT VECTOR WILL TRANSFORM TO SOLUTION VECTOR*
      COE(n) = COE(n)/D(n)
      do 20 i = m,1,-1
         COE(i) = (COE(i) - U(i) *COE(i+1)) / D(i)
20    continue
      return
      end
