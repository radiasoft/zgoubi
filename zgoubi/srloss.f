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
      SUBROUTINE SRLOSS(IMAX,IPASS,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER(80) TA
      PARAMETER (MXTA=45)
      INCLUDE "C.DONT.H"     ! COMMON/DONT/ TA(MXL,MXTA)
      INCLUDE "C.PTICUL.H"     ! COMMON/PTICUL/ AM,Q,G,TO
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      INCLUDE "C.SYNRA.H"     ! COMMON/SYNRA/ KSYN
      LOGICAL SCLFLD, SCLFLO
      SAVE SCLFLD
      CHARACTER(80) TSCAL, LIST
      INTEGER DEBSTR, FINSTR
      LOGICAL OKSR, OKSRO

      LOGICAL OKOPEN, OKPRSR, OKPRSO, IDLUNI
      LOGICAL EMPTY

      SAVE LNSR, OKOPEN, OKSR, OKPRSR

      DATA OKSR / .FALSE. /
      DATA OKOPEN, OKPRSR /.FALSE., .FALSE. /
      DATA TSCAL / ' ' /
      DATA LIST / ' ' /
      DATA SCLFLD / .FALSE. /

      KSYN= A(NOEL,1)
      IF(KSYN.EQ.0) THEN
        IF(NRES.GT.0) WRITE(NRES,FMT='(/,15X,''S.R. LOSS IS OFF'')')
        OKSR = .FALSE.
        RETURN
      ENDIF

      IF(IPASS .EQ. 1) THEN 
       OKPRSR = (NINT(10.D0*A(NOEL,1)) - 10*KSYN) .EQ. 1
       IF(OKPRSR) THEN 
        IF(.NOT.OKOPEN) THEN
          IF(IDLUNI(
     >              LNSR)) THEN
            OPEN(UNIT=LNSR,FILE='zgoubi.SRLOSS.Out',
     >                     FORM='FORMATTED',ERR=99, IOSTAT=IOS)
            WRITE(LNSR,FMT='(A)') 
     >      '# File created by srloss, print out by prsr'
            WRITE(LNSR,FMT='(A)') '# '
            WRITE(LNSR,FMT='(A)') '# 1 2 3 4 5 6 7 8 ...' 
            WRITE(LNSR,FMT='(A)') 
     >      '# KLE, NOEL, IPASS, BORO, DPREF, AM, Q, G, IMAX, '
     >      //'PI*EMIT(1), ALP(1), BET(1), XM(1), XPM(1), '
     >      //'NLIV(1), NINL(1), RATIN(1), '
     >      //'PI*EMIT(2), ALP(2), BET(2), XM(2), XPM(2), '
     >      //'NLIV(2), NINL(2), RATIN(2), '
     >      //'PI*EMIT(3), ALP(3), BET(3), XM(3), XPM(3), '
     >      //'NLIV(3), NINL(3), RATIN(3), '
     >      //'dE local, sigE local ; theoretical '
     >      //'dEav (kEv/prtcl), Eav_phot (keV), Erms_phot (keV)'
          ELSE
            OKPRSR = .FALSE.
            GOTO 99
          ENDIF
          IF(IOS.NE.0) GOTO 99
          OKOPEN = .TRUE.
        ENDIF
       ENDIF
      ENDIF

      TSCAL = TA(NOEL,2)
      LIST = TA(NOEL,3)
      TSCAL = TSCAL(DEBSTR(TSCAL):FINSTR(TSCAL))
      IF(TSCAL(1:5) .EQ. 'scale' .OR. TSCAL(1:5) .EQ. 'SCALE') THEN
        SCLFLD = .TRUE. 
      ELSE
        SCLFLD = .FALSE. 
        TSCAL = ' ' 
      ENDIF
      IF(EMPTY(LIST)) THEN
        LIST = ' '
      ELSE
        LIST = LIST(DEBSTR(LIST):FINSTR(LIST))
      ENDIF

C----- Set SR loss tracking
      IF(NRES.GT.0) THEN 
          WRITE(NRES,FMT='(/,15X,'' S.R.  TRACKING  REQUESTED'')')
          IF(TA(NOEL,1).NE.'ALL' .AND. TA(NOEL,1).NE.'all') 
     >                                    WRITE(NRES,FMT='(20X,
     >           '' Accounted for only in '',A)') TA(NOEL,1)
          IF(SCLFLD) THEN
            WRITE(NRES,FMT='(20X,'' Magnetic strengths will scale '',
     >      ''with energy lost in the following list of elements :'',
     >    /,25X,''{'',A,''}'')') LIST(DEBSTR(LIST):FINSTR(LIST))
          ELSE
            WRITE(NRES,FMT='(20X,'' Magnetic strengths will NOT '',
     >      ''scale with energy lost in dipole fields.'')')
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

      RETURN

      ENTRY SRLOSR(
     >             SCLFLO)
      SCLFLO=SCLFLD
      RETURN

      ENTRY SRLOS1(
     >             OKPRSO,LNSRO)
      OKPRSO = OKPRSR
      LNSRO = LNSR
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
