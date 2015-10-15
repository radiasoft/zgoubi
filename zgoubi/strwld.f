C  ZGOUBI, a program for computing the trajectories of charged particles
C  in electric and magnetic fields
C  Copyright (C) 1988-2007  Fran�ois M�ot
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
C  Fran�ois M�ot <fmeot@bnl.gov>
C  Brookhaven National Laboratory     
C  C-AD, Bldg 911
C  Upton, NY, 11973, USA
C  -------
      FUNCTION STRWLD(NLB,LBLST,LABEL,
     >                                IS)
C A string LBLST(i=1,NLB) may have a wild card '*' at start eor at end
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL STRWLD
      CHARACTER(*) LABEL
      CHARACTER(*) LBLST(*)

      INTEGER DEBSTR, FINSTR

      STRWLD = .FALSE.
      
      DO I = 1, NLB
        IA = DEBSTR(LBLST(I))
        IB = FINSTR(LBLST(I))
c             write(*,*) ' strwld ',ia,ib,LBLST(I)(IA:IA)
c             write(*,*) ' strwld ',ia,ib,LBLST(I)(Ib:Ib)
        IF    (LBLST(I)(IA:IA) .EQ. '*') THEN
          IC = FINSTR(LABEL)
          IF(LBLST(I)(IA+1:IB) .EQ. LABEL(IC-(IB-IA)+1:IC)) THEN
c             write(*,*) LBLST(I)(IA+1:IB),' ', LABEL(IC-(IB-IA)+1:IC)
            STRWLD = .TRUE.
            IS = I
            RETURN
          ELSE
            STRWLD = .FALSE.
          ENDIF
        ELSEIF(LBLST(I)(IB:IB) .EQ. '*') THEN
          IC = DEBSTR(LABEL)
          IF(LBLST(I)(IA:IB-1) .EQ. LABEL(IC:IC+IB-IA-1)) THEN
c             write(*,*) LBLST(I)(IA:IB-1),' ', LABEL(IC:IC+IB-IA-1)
            STRWLD = .TRUE.
            IS = I
            RETURN
          ELSE
            STRWLD = .FALSE.
          ENDIF
        ENDIF
      ENDDO

      RETURN
      END
