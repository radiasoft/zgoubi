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
      SUBROUTINE FITWDA
C Will cause save of zgoubi.dat list with updated variables as following from FIT[2].
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      
      CHARACTER(2000) TXT132

      PARAMETER (KSIZ=10)
      CHARACTER(KSIZ)  KLEY
      INTEGER DEBSTR, FINSTR
      LOGICAL OK, IDLUNI

      OK = IDLUNI(
     >            LWDAT)
      CALL SYSTEM('\cp zgoubi.dat zgoubi.FIT.out.dat')
      OPEN(UNIT=LWDAT,FILE='zgoubi.FIT.out.dat')
      OPEN(UNIT=LTEMP,FILE='zgoubi.temp.dat')
      IF(NRES.GT.0) WRITE(NRES,FMT='(/,20X,
     >''Saved new version of zgoubi.dat with variables updated.'')')

      TXT132 = '***' 
      DOWHILE(TXT132(1:5) .NE.'''END''')
        READ(LWDAT,FMT='(A)',ERR=10,END=10) TXT132
        TXT132 = TXT132(DEBSTR(TXT132):FINSTR(TXT132))
        WRITE(LTEMP,FMT='(A)') 
     >       TXT132(DEBSTR(TXT132):FINSTR(TXT132))
        IF(TXT132(1:1) .EQ. '''') THEN
          READ(TXT132(105:132),*,err=11,end=11) NUEL      ! Position follows from prdata
          CALL ZGKLE(IQ(NUEL), 
     >                        KLEY)
          IF    (KLEY(1:8) .EQ. 'MULTIPOL') THEN 
            READ(LWDAT,FMT='(A)',err=10,end=10) TXT132
            WRITE(LTEMP,FMT='(A)') 
     >                    TXT132(DEBSTR(TXT132):FINSTR(TXT132))
            READ(LWDAT,FMT='(A)',err=10,end=10) TXT132
            WRITE(LTEMP,FMT='(F11.6,F7.2,3F15.10,7F4.1)')
     >                                  (A(NUEL,J),J=2,13)
          ELSEIF(KLEY(1:8) .EQ. 'CHANGREF') THEN 
            READ(LWDAT,FMT='(A)',err=10,end=10) TXT132
            WRITE(LTEMP,FMT='(3F14.8)') (A(NUEL,J),J=1,3)
          ENDIF
        ENDIF
      ENDDO


 10   CONTINUE

      IF(NRES.GT.0) WRITE(NRES,FMT='(/,20X,
     >''Updated version of zgoubi.dat saved in  '',a)')
     >'zgoubi.FIT.out.dat'
      WRITE(*,FMT='(/,20X,
     >''Updated version of zgoubi.dat saved in  '',a)')
     >'zgoubi.FIT.out.dat'

      CLOSE(LWDAT)
      CLOSE(LTEMP)

        CALL SYSTEM('\cp zgoubi.temp.dat zgoubi.FIT.out.dat')
        CALL SYSTEM('rm -f zgoubi.temp.dat')

      RETURN

 11   CONTINUE
      CALL ENDJOB('Pgm fitwda. Need number at each element.'
     >//' This element may have no number in zgoubi.dat ?',-99)
      RETURN
      END
