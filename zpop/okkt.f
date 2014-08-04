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
      FUNCTION OKKT(KT1,KT2,KT3,IT,KEX,LET,
     >                                 IEND)
      LOGICAL OKKT
      CHARACTER*1 LET,KLETO
      LOGICAL OKKT5

      INCLUDE 'MAXNTR.H'          
      LOGICAL LOST(NTRMAX)
      SAVE LOST
      DATA LOST / NTRMAX*.FALSE. / 

      CALL FBGTXT

      INCLUDE "OKKT.H"

      CALL READC9(KEXO,KLETO)

C------ Plot or store only as long as KEX is correct 
      IF(KEXO.NE.99) OKKT = OKKT .AND. (KEX.EQ.KEXO)
      IF(KLETO.NE.'*') THEN
        IF(KLETO.EQ.'S') THEN
C------ Only secondary particles (from decay process using MCDESINT) can be plotted
          OKKT = OKKT .AND. (LET.EQ.'S')
        ELSEIF(KLETO.EQ.'P') THEN
C------ Only parent particles (if decay process using MCDESINT) can be plotted
          OKKT = OKKT .AND. (LET.NE.'S')
        ENDIF
      ENDIF

      IF(.NOT.LOST(IT)) LOST(IT) = KEX.LE.0
      RETURN

      ENTRY OKKT5(KT)      
      OKKT5 = .NOT. LOST(KT)
c         write(89,*) 'oeiuvc',okkt, okkt5,it
c         pause
      RETURN
      END
