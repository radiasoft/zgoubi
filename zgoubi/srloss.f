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
      SUBROUTINE SRLOSS(IMAX,IPASS,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER*80 TA
      COMMON/DONT/ TA(MXL,40)
      COMMON/PTICUL/ AM,Q,G,TO
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      COMMON/SYNRA/ KSYN
      LOGICAL SCALE, SCALO
      SAVE SCALE
      CHARACTER*80 TXT
      INTEGER DEBSTR, FINSTR
      LOGICAL OKSR, OKSRO

      LOGICAL OKOPEN, OKIMP, IDLUNI
      SAVE LUN, OKOPEN, OKIMP

      DATA OKSR / .FALSE. /
      DATA OKOPEN, OKIMP /.FALSE., .FALSE. /

      KSYN= A(NOEL,1)
      IF(KSYN.EQ.0) THEN
        IF(NRES.GT.0) WRITE(NRES,FMT='(/,15X,''S.R. LOSS IS OFF'')')
        OKSR = .FALSE.
        RETURN
      ENDIF

      IF(IPASS .EQ. 1) THEN 
       OKIMP = (NINT(10.D0*A(NOEL,1)) - 10*KSYN) .EQ. 1
       IF(OKIMP) THEN 
        IF(.NOT.OKOPEN) THEN
          IF(IDLUNI(
     >              LUN)) THEN
            OPEN(UNIT=LUN,FILE='zgoubi.SRLOSS.Out',
     >                     FORM='FORMATTED',ERR=99, IOSTAT=IOS)
          ELSE
            OKIMP = .FALSE.
            GOTO 99
          ENDIF
          IF(IOS.NE.0) GOTO 99
          OKOPEN = .TRUE.
        ENDIF
       ENDIF
      ENDIF

      TXT = TA(NOEL,2)
      TXT = TXT(DEBSTR(TXT):FINSTR(TXT))
      SCALE=TXT.EQ.'SCALE' .OR. TXT.EQ.'scale'

C----- Set SR loss tracking
      IF(NRES.GT.0) THEN 
          WRITE(NRES,FMT='(/,15X,'' S.R.  TRACKING  REQUESTED'')')
          IF(TA(NOEL,1).NE.'ALL' .AND. TA(NOEL,1).NE.'all') 
     >                                    WRITE(NRES,FMT='(20X,
     >           '' Accounted for only in '',A)') TA(NOEL,1)
          IF(TA(NOEL,2).EQ.'SCALE') THEN
            WRITE(NRES,FMT='(20X,'' Magnetic strengths will '',
     >      ''scale with energy lost in dipole fields'')')
          ELSE
            WRITE(NRES,FMT='(20X,'' Magnetic strengths will NOT '',
     >      ''scale with energy lost in dipole fields'')')
          ENDIF 
          IF(A(NOEL,10).NE.1) WRITE(NRES,FMT='(20X,
     >      '' Loss entails dp only, no angle kick installed!! '')') 
          P=BORO*CL9*Q
          E=SQRT(P*P + AM*AM)
          BTA = P/E
          GAMMA=E/AM
          WRITE(NRES,102) BORO,BTA,GAMMA,E-AM
 102      FORMAT(
     >    /,25X,' Reference dynamical parameters :', 1P,
     >    /,30X,' B.rho           =  ',G15.7,' kG*cm',
     >    /,30X,' beta=v/c        =  ',G15.7,
     >    /,30X,' gamma           =  ',G15.7,
     >    /,30X,' Kinetic  energy =  ',G15.7,' MeV') 
      ENDIF
 
      IF(Q*AM.EQ.0D0) THEN
        WRITE(NRES,106)
 106    FORMAT(//,15X,' Please provide mass and charge of particles !'
     >         ,/,15X,' - use keyword ''PARTICUL''',/)
        RETURN 1
      ENDIF

      OKSR = .TRUE.

      IRA=1+(NINT(A(NOEL,11))/2)*2    
      IF(IPASS.EQ.1) THEN 
        CALL RAYSY1(IMAX,IRA)
        CALL RAYSY3(TA(NOEL,1))
      ENDIF

      IF(OKIMP) CALL SRPRN(I0,LUN,IMAX)

      RETURN

      ENTRY SRLOSR(
     >             SCALO)
      SCALO=SCALE
      RETURN

      ENTRY SRLOS3(
     >             OKSRO)
      OKSRO = OKSR
      RETURN

 99   IF(NRES .GT. 0) WRITE(NRES,101) 'SRLOSS',
     >                          ' OPEN zgoubi.SRLOSS.Out'
 101  FORMAT(/,' ******  SBR ',A,' : ERROR',A,/)
      RETURN
      END
