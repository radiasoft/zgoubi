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
C  Brookhaven National Laboratory                    és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      FUNCTION POIDEV(XM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (PI=3.141592653589793238D0)
      DATA OLDM / -1.D0 /
C Going through GAMMLN calculation (namely, whenever XM.ge.12.) 
C increase the CPU time for POIDEV by a factor of about 2.5 (HP 
C station, 1999) 
      IF(XM.LT.12.D0) THEN
        IF(XM.NE.OLDM) THEN
          OLDM=XM
          G=EXP(-XM)
        ENDIF
        EM=-1
        T=1.D0
 2      EM=EM+1.D0
        R1=RNDM()
        T=T*R1
        IF(T.GT.G) GOTO 2
      ELSE
        IF(XM.NE.OLDM) THEN
          OLDM=XM
          SQ=SQRT(2.D0*XM)
          ALXM=LOG(XM)
          G=XM*ALXM-GAMMLN(XM+1D0)
        ENDIF
 1      Y=TAN(PI*RNDM())
        EM=SQ*Y+XM
        IF(EM.LT.0.D0) GOTO 1
        EM=INT(EM)
        GG=GAMMLN(EM+1.D0)
        T=0.9D0*(1.D0+Y*Y)*EXP(EM*ALXM-GG-G)
        IF(RNDM().GT.T) GOTO 1
      ENDIF
      POIDEV=EM
      RETURN
      END
