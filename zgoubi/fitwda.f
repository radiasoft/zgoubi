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
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER(2000) TXT132
      PARAMETER (KSIZ=10)
      CHARACTER(KSIZ)  KLEY
      INTEGER DEBSTR, FINSTR
      LOGICAL OK, IDLUNI, STRCON, ISNUM
      PARAMETER (MSR=8,MSR2=2*MSR)
      CHARACTER(30) STRA(MSR2)
      CHARACTER(40) FRMT
      CHARACTER(1) LET

      OK = IDLUNI(
     >            LWDAT)
      CALL SYSTEM('\cp zgoubi.dat zgoubi.FIT.out.dat')
      OPEN(UNIT=LWDAT,FILE='zgoubi.FIT.out.dat')
      OK = IDLUNI(
     >            LTEMP)
      OPEN(UNIT=LTEMP,FILE='zgoubi.temp.dat')
      IF(NRES.GT.0) WRITE(NRES,FMT='(/,20X,''Saving new version of '',
     >''zgoubi.dat to zgoubi.FIT.out.dat, with variables updated.'')')

      TXT132 = 'not END !' 
      DOWHILE(TXT132(1:5) .NE.'''END''')
        READ(LWDAT,FMT='(A)',ERR=10,END=10) TXT132
        TXT132 = TXT132(DEBSTR(TXT132):FINSTR(TXT132))
        WRITE(LTEMP,FMT='(A)') 
     >       TXT132(DEBSTR(TXT132):FINSTR(TXT132))
        IF(TXT132(1:1) .EQ. '''') THEN
          READ(TXT132(105:132),*,ERR=11,END=11) NUEL      ! Position follows from prdata
          IF(NUEL .LT. 0 .OR. NUEL .GT. MXL) GOTO 11
          CALL ZGKLE(IQ(NUEL), 
     >                        KLEY)
          IF    (KLEY(1:8) .EQ. 'OBJET') THEN 
            READ(LWDAT,FMT='(A)',ERR=10,END=10) TXT132              ! BORO
            WRITE(LTEMP,FMT='(A)')                                  
     >                    TXT132(DEBSTR(TXT132):FINSTR(TXT132))
            READ(LWDAT,FMT='(A)',ERR=10,END=10) TXT132              ! KOBJ[.KOBJ2]
            WRITE(LTEMP,FMT='(A)') 
     >                    TXT132(DEBSTR(TXT132):FINSTR(TXT132))
             IF(STRCON(TXT132,'.',
     >                           IS)) THEN
              READ(TXT132(1:IS-1),*) KOBJ 
c                   write(*,*) ' fitwda ',TXT132(IS+1:FINSTR(TXT132))
               READ(TXT132(IS+1:FINSTR(TXT132)),*) KOBJ2 
c                   write(*,*) ' fitwda kobj, kobj2 : ',kobj, kobj2
c                       read(*,*)
            ELSE
              READ(TXT132,*) KOBJ 
              KOBJ2 = 0
            ENDIF
            IF    (KOBJ.EQ.2) THEN
              READ(LWDAT,FMT='(A)',ERR=10,END=10) TXT132
              WRITE(LTEMP,FMT='(A)') 
     >                    TXT132(DEBSTR(TXT132):FINSTR(TXT132))
              READ(LWDAT,FMT='(A)',ERR=10,END=10) TXT132
              READ(TXT132,*,ERR=10,END=10) Y,T,Z,P,X,D,LET
              WRITE(LTEMP,FMT='(1P,4(1X,E14.6),F7.2,1X,E13.6,3A)')
     >              (A(NUEL,J),J=30,35),' ''',LET,''' '
            ELSEIF(KOBJ.EQ.5) THEN
              READ(LWDAT,FMT='(A)',ERR=10,END=10) TXT132                 ! Sampling
              WRITE(LTEMP,FMT='(1P,4(1X,E14.6),F7.2,1X,E13.6)')
     >              (A(NUEL,J),J=20,25)
              READ(LWDAT,FMT='(A)',ERR=10,END=10) TXT132                 ! (1st) reference
              WRITE(LTEMP,FMT='(1P,4(1X,E14.6),F7.2,1X,E13.6)')
     >              (A(NUEL,J),J=30,35)
              IF(KOBJ2 .GE. 1) THEN
                READ(LWDAT,FMT='(A)',ERR=10,END=10) TXT132
                IF(KOBJ2 .EQ. 1) THEN
                  WRITE(LTEMP,FMT='(1P,10(E10.3,1X))')         ! Initial beta 
     >            (A(NUEL,J),J=40,49)
                ELSE
                  DO I = 1, KOBJ2-1
                    READ(LWDAT,FMT='(A)',ERR=10,END=10) TXT132           ! Next references in a total of KOBJ2
                    WRITE(LTEMP,FMT='(1P,4(1X,E14.6),F7.2,1X,E13.6)')
     >              (A(NUEL,J),J=(I+4)*10,(I+4)*10+5)
                  ENDDO
                ENDIF
              ENDIF
            ENDIF

          ELSEIF(KLEY(1:8) .EQ. 'AGSMM') THEN 

            READ(LWDAT,FMT='(A)',ERR=10,END=10) TXT132
            WRITE(LTEMP,FMT='(A)') 
     >                    TXT132(DEBSTR(TXT132):FINSTR(TXT132))
            READ(LWDAT,FMT='(A)',ERR=10,END=10) TXT132
            WRITE(LTEMP,FMT='(F11.6,2F7.2,3E16.8)')
     >                                  (A(NUEL,J),J=10,15)

          ELSEIF(KLEY(1:8) .EQ. 'CHANGREF') THEN 

            READ(LWDAT,FMT='(A)',ERR=10,END=10) TXT132
            IF(STRCON(TXT132,'!',
     >                           IS))
     >      TXT132 = TXT132(DEBSTR(TXT132):IS-1)
            CALL STRGET(TXT132,MSR2,
     >                             NSR2,STRA)
            IF(.NOT. ISNUM(STRA(1))) THEN
C New style CHANGREF
              WRITE(FRMT,FMT='(A,I0,A)') '(,',NSR2/2,'(A2,1X,F12.8,1X))'
              WRITE(LTEMP,FRMT) 
     >        (STRA(J)(DEBSTR(STRA(J)):FINSTR(STRA(J))),
     >        A(NUEL,J/2+1),J=1,NSR2,2)
            ELSE
C Old style CHANGREF
              WRITE(LTEMP,FMT='(3(E14.6,1X))') (A(NUEL,J),J=1,3)
            ENDIF

          ELSEIF(KLEY(1:6) .EQ. 'SPNTRK') THEN 

            READ(LWDAT,FMT='(A)',ERR=10,END=10) TXT132
            WRITE(LTEMP,FMT='(A)') 
     >                    TXT132(DEBSTR(TXT132):FINSTR(TXT132))
            READ(TXT132,*,ERR=10,END=10) XSO
            IF(XSO .GE. 4) THEN
              IF    (NINT(10*XSO) .EQ. 40) THEN
                CALL ZGIMAX(
     >                      IMAX) 
                DO IT = 1, IMAX
                  READ(LWDAT,FMT='(A)',ERR=10,END=10) TXT132
                  WRITE(LTEMP,FMT='(F11.6,F7.2,3E16.8,7F4.1)')
     >                                  (A(NUEL,J),J=2,13)
                ENDDO
              ELSEIF(NINT(10*XSO) .EQ. 41) THEN
                READ(LWDAT,FMT='(A)',ERR=10,END=10) TXT132
                WRITE(LTEMP,FMT='(3(E16.8,2X))')
     >                                  (A(NUEL,J),J=10,12)
              ELSEIF(NINT(10*XSO) .EQ. 50) THEN
                READ(LWDAT,FMT='(A)',ERR=10,END=10) TXT132
                WRITE(LTEMP,FMT='(F11.6,F7.2,3E16.8,7F4.1)')
     >                                  (A(NUEL,J),J=2,13)
                READ(LWDAT,FMT='(A)',ERR=10,END=10) TXT132
                WRITE(LTEMP,FMT='(F11.6,F7.2,3E16.8,7F4.1)')
     >                                  (A(NUEL,J),J=2,13)
              ELSE
                CALL ENDJOB('Pgm fitwda. Present KSO option in SPNTRK '
     >          //' needs be installed in fitwda. ',-99)
              ENDIF
            ENDIF

          ELSEIF(KLEY(1:8) .EQ. 'MULTIPOL') THEN 

            READ(LWDAT,FMT='(A)',ERR=10,END=10) TXT132
            WRITE(LTEMP,FMT='(A)') 
     >                    TXT132(DEBSTR(TXT132):FINSTR(TXT132))
            READ(LWDAT,FMT='(A)',ERR=10,END=10) TXT132
            WRITE(LTEMP,FMT='(F11.6,F7.2,1P,3E16.8,0P,7F4.1)')
     >                                  (A(NUEL,J),J=2,13)
            DO I = 1, 6
              READ(LWDAT,FMT='(A)',ERR=10,END=10) TXT132
              WRITE(LTEMP,FMT='(A)') 
     >                    TXT132(DEBSTR(TXT132):FINSTR(TXT132))
            ENDDO
            READ(LWDAT,FMT='(A)',ERR=10,END=10) TXT132
            WRITE(LTEMP,FMT='(I0,1X,1P,3E16.8)')
     >                           NINT(A(NUEL,63)),(A(NUEL,J),J=64,66)  ! KPOS, XCE, YC, ALE

          ELSEIF(KLEY(1:6) .EQ. 'DIPOLE') THEN 
 
            READ(LWDAT,FMT='(A)',ERR=10,END=10) TXT132                 ! IL
            WRITE(*,FMT='(A)')  
     >                    TXT132(DEBSTR(TXT132):FINSTR(TXT132))
            WRITE(LTEMP,FMT='(A)')  
     >                    TXT132(DEBSTR(TXT132):FINSTR(TXT132))
            READ(LWDAT,FMT='(A)',ERR=10,END=10) TXT132                 ! AT, RM
            WRITE(*,FMT='(A)')  
     >                    TXT132(DEBSTR(TXT132):FINSTR(TXT132))
            WRITE(LTEMP,FMT='(A)') 
     >                    TXT132(DEBSTR(TXT132):FINSTR(TXT132))
            READ(LWDAT,FMT='(A)',ERR=10,END=10) TXT132
            WRITE(*,FMT='(A)')  
     >                    TXT132(DEBSTR(TXT132):FINSTR(TXT132))
            WRITE(*,FMT='(5F13.7,T70,A)')
     >             (A(NUEL,J),J=4,8), ' ! ACNT,  HNORM, indices'
            WRITE(LTEMP,FMT='(5F13.7,T70,A)')
     >             (A(NUEL,J),J=4,8), ' ! ACNT,  HNORM, indices'
            DO I = 1, 11
              READ(LWDAT,FMT='(A)',ERR=10,END=10) TXT132
              WRITE(*,FMT='(A)')  
     >                    TXT132(DEBSTR(TXT132):FINSTR(TXT132))
              WRITE(LTEMP,FMT='(A)') 
     >                    TXT132(DEBSTR(TXT132):FINSTR(TXT132))
            ENDDO
            READ(LWDAT,FMT='(A)',ERR=10,END=10) TXT132
            WRITE(*,FMT='(A)')  
     >                    TXT132(DEBSTR(TXT132):FINSTR(TXT132))
            WRITE(LTEMP,FMT='(I0,4(1X,F13.7),T70,A)')
     >                      NINT(A(NUEL,63)),(A(NUEL,J),J=64,67),
     >                                  ' ! KPOS, RE, TE, RS, TS' 
C                read(*,*)
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
      WRITE(*,fmt='(/,A)') 
     >'Pgm fitwda. Need number at each element.'
     >//' Some element may have wrong/no number in zgoubi.dat ?'
      WRITE(*,fmt='(A,/)') 'Write to zgoubi.FIT.out.dat skipped.'
      NRESA = ABS(NRES)
      WRITE(NRESA,fmt='(/,10X,A)') 
     >'Pgm fitwda. Need number at each element.'
     >//' Some element may have wrong/no number in zgoubi.dat ?'
      WRITE(NRESA,fmt='(10X,A,/)') 
     >'Write to zgoubi.FIT.out.dat skipped.'


      RETURN
      END
