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
      SUBROUTINE GO2KEY(NUMKEY,IQCNT,mxkle,KLE,
     >                                   KLEY, LBL1, LBL2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     -------------------------------------------------------
C     Places pointer right before KLEY # numkey in zgoubi.dat,
C     ready to carry on reading in sequence. 
C     Returns keyword, label1, label2
C     -------------------------------------------------------
      CHARACTER(*)  KLE(*)
      PARAMETER (KSIZ=10)
      CHARACTER(KSIZ) KLEY, LBL1, LBL2

      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL

      CHARACTER(40) TXT40
      INTEGER DEBSTR, FINSTR
      LOGICAL STRCON
      CHARACTER(KSIZ) STRA(3)
      SAVE NKMX
      DATA NKMX / -1 /

      REWIND(NDAT)

      NK = 0
 1    CONTINUE
        READ(NDAT,FMT='(A)',END=98,ERR=99) TXT40
        IF(STRCON(TXT40,'!',
     >                      IS)) THEN
           IF(IS.GE.2) THEN
             TXT40 = TXT40(1:IS-1)
           ELSE
             TXT40 = ' '
           ENDIF
        ENDIF
        IF(DEBSTR(TXT40) .GE. 1) THEN
          TXT40 = TXT40(DEBSTR(TXT40):FINSTR(TXT40))
        ELSE
          DO I = 1, 40
            TXT40(I:I) = ' '
          ENDDO
        ENDIF
        IF(TXT40(1:1).NE.'''') GOTO 1
        NK = NK+1

        CALL STRGET(TXT40,3,
     >                    NST,STRA)
        KLEY = STRA(1)(debstr(STRA(1))+1:finstr(STRA(1))-1)
        IF(NST.GE.2) THEN
          LBL1 = STRA(2) 
          IF(NST.EQ.3) LBL2 = STRA(3)
        ENDIF

        IF(IQCNT .EQ. 1) THEN
          IF(NKMX.LT.NK) THEN
            DO IKLE=1,MXKLE

c       write(*,*) 'go2key iqcnt  *',kley,'*',kle(ikle),'*'
c     > ,kley .eq. kle(ikle)

              IF(KLEY .EQ. KLE(IKLE)) then 
                IQ(NK) =  IKLE

c                write(*,*) 'go2key iqcnt ',iqcnt,nk,iq(nk)
c                write(*,*) 'go2key iqcnt ',iqcnt,nk,iq(nk)
c                     read(*,*)
              endif
            ENDDO          
            NKMX = NK
          ENDIF
        ENDIF

      IF(NK .LT. NUMKEY) GOTO 1

      BACKSPACE(NDAT)

C Test
c        READ(NDAT,FMT='(A)',END=98,ERR=99) TXT40
c      backspace(ndat)


c      write(*,*) ' go2key txt40 : ',TXT40(DEBSTR(TXT40):FINSTR(TXT40))
c      write(*,*) ' go2key ::::::: ',nst,numkey,(stra(i),i=1,nst)
c               read(*,*)

 98   RETURN
 99   STOP ' *** ERROR IN SBR GO2KEY'
      END 
