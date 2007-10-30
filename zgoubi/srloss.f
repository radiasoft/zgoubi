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
C  François Méot <meot@lpsc.in2p3.fr>
C  Service Accélerateurs
C  LPSC Grenoble
C  53 Avenue des Martyrs
C  38026 Grenoble Cedex
C  France
      SUBROUTINE SRLOSS(IMAX,IPASS,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER*80 TA
      COMMON/DONT/ TA(MXL,20)
      COMMON/PTICUL/ AM,Q,G,TO
      COMMON/RIGID/ BORO,DPREF,DP,BR
      COMMON/SYNRA/ KSYN
      LOGICAL SCALE, SCALO
      SAVE SCALE
      CHARACTER*80 TXT
      INTEGER DEBSTR, FINSTR
      KSYN= A(NOEL,1)
      IF(KSYN.EQ.0) THEN
        IF(NRES.GT.0) WRITE(NRES,FMT='(/,15X,''S.R. TRACKING IS OFF'')')
        RETURN
      ENDIF

      TXT = TA(NOEL,2)
      TXT = TXT(DEBSTR(TXT):FINSTR(TXT))
      SCALE=TXT.EQ.'SCALE' .OR. TXT.EQ.'scale'

C----- Set SR loss tracking
 1    CONTINUE
      IF(NRES.GT.0) THEN 
          WRITE(NRES,FMT='(/,15X,'' S.R.  TRACKING  REQUESTED'')')
          IF(TA(NOEL,1).NE.'ALL') WRITE(NRES,FMT='(20X,
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
C          P=BORO*CL*1.D-9*Q/QE
C          P=BORO*CL*1.D-9*Q
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
 
C      IF(AM*Q .EQ. 0.D0) THEN
      IF(Q*AM.EQ.0D0) THEN
        WRITE(NRES,106)
 106    FORMAT(//,15X,' PLEASE GIVE MASS AND CHARGE OF PARTICLES !'
     >         ,/,15X,' - USE KEWORD ''PARTICUL''',/)
        RETURN 1
      ENDIF

C      TXT = TA(NOEL,2)
C      TXT = TXT(DEBSTR(TXT):FINSTR(TXT))
C      SCALE=TXT.EQ.'SCALE' .OR. TXT.EQ.'scale'
      IRA=1+(NINT(A(NOEL,11))/2)*2    
      IF(IPASS.EQ.1) THEN 
        CALL RAYSY1(IMAX,IRA)
        CALL RAYSY3(TA(NOEL,1))
      ENDIF
      RETURN
      ENTRY SRLOSR(SCALO)
      SCALO=SCALE
      RETURN
      END
