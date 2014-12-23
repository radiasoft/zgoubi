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
C  USA
C  -------
C
C     PHYSICAL CONSTANTS
C
C     c        (m/s)
      PARAMETER (CL=2.99792458E8)
C     Pi       (rad)
      PARAMETER (PI=3.1415926535897932384)
C     2*Pi     (rad)
      PARAMETER (DEUXPI=2.D+0*PI)
C     4*Pi     (rad)
      PARAMETER (QUATPI=4.D+0*PI)
C     Pi/2     (rad)
      PARAMETER (PISUR2=PI/2.D+0)
C     mu0      (H/m)
      DOUBLE PRECISION MU0
      PARAMETER (MU0=QUATPI*1.D-07)
C     epsilon0=1/(muo*c^2) (F/m)
      PARAMETER (EPS0=1.D+00/MU0/CL/CL)
C     h*c      (eV*m)
      PARAMETER (CHPLANCK=1.23984244D-6)
C     h/(2*Pi) (J.s)
      PARAMETER (HBARRE=1.05457266E-34)
C     Qelectron(C)
      PARAMETER (ECHARGE=1.60217733E-19)
C     Melectron(eV)
      PARAMETER (EMASSE0=0.511003E+6) 
      PARAMETER (EMAS0G=EMASSE0/1.E+9)
C     Relectron(m) (q^2/(4Pi*eps0*(mc^2*q))
      PARAMETER (ERAYON=2.81794092E-15)
C     Mproton  (eV)
      PARAMETER (PMASSE0=938.2796E+6) 
C     Rproton (m)  (q^2/(4Pi*eps0*(mc^2*q))
      PARAMETER (PRAYON=ECHARGE/QUATPI/EPS0/PMASSE0)
C     Cq       (m)
C     Cq = 55/(32 Sqrt(3)) (h/2Pi) /(mc)
C     Cq = 3.83191E-13
      PARAMETER (SQRT3=1.732050808)
      PARAMETER (CQ = 55.*HBARRE*CL/32./SQRT3/EMASSE0/ECHARGE)
C     Cgamma   (m/GeV^3)
C     Ggamma = 4 Pi Re /(3 (m c^2)^3) 
C     Ggamma = 8.84607E-5
      PARAMETER (CGAM=QUATPI*ERAYON/EMAS0G/EMAS0G/EMAS0G/3.)      
C     Square Root (2)
      PARAMETER (SQRT2=1.414213562)      
