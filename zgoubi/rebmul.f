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
      SUBROUTINE REBMUL(ND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ----------------------------
C     READS DATA FOR E-B MULTIPOLE
C     ----------------------------
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      COMMON/MARK/ KART,KALC,KERK,KUASEX
      PARAMETER(MPOL=10)
 
      IB = 1
      READ(NDAT,*) A(NOEL,IB)
 
C----- E FIELD
C----- X,RO,E1,E2,E3,E4,E5,E6,...E10
      IA = IB + 1             
      IB = IA+MPOL+1          
      READ(NDAT,*) (A(NOEL,I),I=IA,IB)
C----- Entrance fringe field
      IA = IB + 1             
      IB = IA+MPOL            
      READ(NDAT,*) (A(NOEL,I),I=IA,IB)
      IA = IB+1               
      IB = IA + 6             
      READ(NDAT,*) II,(A(NOEL,I),I=IA+1,IB)
      A(NOEL,IA) = II
C----- Exit fringe field
      IA = IB + 1             
      IB = IA + MPOL          
       READ(NDAT,*) (A(NOEL,I),I=IA,IB)
      IA = IB+1               
      IB = IA + 6             
       READ(NDAT,*) II,(A(NOEL,I),I=IA+1,IB)
      A(NOEL,IA) = II
C------- Rotation of multipole components
      IA = IB +1              
      IB = IA + MPOL - 1      
      READ(NDAT,*) (A(NOEL,I),I=IA,IB)
 
C----- B FIELD
C----- X,RO,B1,B2,B3,B4,B5,B6,...B10
      IA = IB + 1             
      IB = IA+MPOL+1          
      READ(NDAT,*) (A(NOEL,I),I=IA,IB)
C----- Entrance fringe field
      IA = IB + 1             
      IB = IA+MPOL            
      READ(NDAT,*) (A(NOEL,I),I=IA,IB)
      IA = IB+1               
      IB = IA + 6             
      READ(NDAT,*) II,(A(NOEL,I),I=IA+1,IB)
      A(NOEL,IA) = II
C----- Exit fringe field
      IA = IB + 1             
      IB = IA+MPOL            
      READ(NDAT,*) (A(NOEL,I),I=IA,IB)
      IA = IB+1               
      IB = IA + 6             
      READ(NDAT,*) II,(A(NOEL,I),I=IA+1,IB)
      A(NOEL,IA) = II
C----- Rotation of multipole components
      IA = IB +1              
      IB = IA + MPOL - 1      
      READ(NDAT,*) (A(NOEL,I),I=IA,IB)
 
      ND= IB + 1              
C----- XPAS
      READ(NDAT,*) A(NOEL,ND)
C----- KP,XCE,YCE,ALE
C Modif, FM, Dec. 05
C      READ(NDAT,*) II,(A(NOEL,I),I=ND+2,ND+4)
C      A(NOEL,ND+1) = II
      READ(NDAT,*) II,(A(NOEL,I),I=ND+4,ND+6)
      A(NOEL,ND+3) = II
 
      RETURN
      END
