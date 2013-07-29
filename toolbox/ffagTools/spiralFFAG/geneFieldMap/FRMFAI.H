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
C Storage format advocated in
C      SUBROUTINE IMPFAI
C Add space(1X) in the line of C5 to ensure a blank  between the 
C data (  Takeichiro Yokoi)16/11/2006 

 110     FORMAT(1X,

C1   LET(IT),KEX,      XXXO,(FO(J,IT),J=2,MXJ),
     1           A1,I2,1P,7E16.8,

C2        XXX,Y,T*1.D3,
     2         /,3E24.16,

C3       Z,P*1.D3,SAR,     TAR,  ENEKI,
     3         /,4E24.16,E16.8,

C4       IT,IREP(IT),SORT(IT), (AMQ(J,I),J=1,5), RET(IT), DPR(IT),
     4         /,2I6,8E16.8,

C5        BORO,  IPASS,  KLEY,  (LABEL(NOEL,I),I=1,2),NOEL
C     5  /,E16.8,   I6,  A8,  2A8,                I5)
     5  /,E16.8,   I6,1X,  A8, 1X , 2A8,                I5)

