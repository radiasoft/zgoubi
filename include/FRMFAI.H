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
C Storage format advocated in
C      SUBROUTINE IMPFAI
C Add space(1X) in the line of C5 to ensure a blank  between the 
C data (  Takeichiro Yokoi)16/11/2006 

 110     FORMAT(1X

C1               KEX,      XXXO,(FO(J,IT),J=2,MXJ)
     >  ,1P,   1X, I3,    7(1X,E16.8)

C2              XXX,Y,T*1.D3
     >         ,3(1X,E24.16)

C3       Z,P*1.D3,SAR,TAR 
     >   ,4(1X,E24.16)

C4       SXo, SYo, SZo, So, SX, SY, SZ, S
     >   ,8(1X,E15.7)

C5        ENEKI,ENERG
     >   ,2(1X,E16.8)

C6       IT,IREP(IT),   SORT(IT), (AMQ(J,I),J=1,5), RET(IT), DPR(IT), PS
     >   ,2(1X,I6),   9(1X,E16.8)

C7         BORO,    IPASS, NOEL,    
     >   ,1X,E16.8, 1X,I8, 1X,I6

C8            'KLEY',    ('LABEL(NOEL,I)',I=1,2),           'LET(IT)'
     >    ,1X,A1,A10,A1,   2(1X,A1,A ,A1),                   1X,3A1)
