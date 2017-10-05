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
      SUBROUTINE SYNPAR(B,XTLT,XL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "C.PTICUL.H"     ! COMMON/PTICUL/ AM,Q,G,TO
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      LOGICAL SCALE
      SAVE SMELPP,E
      DATA UNITB / 1.D-1/
      DATA SMELPP, SNMPP, SRMSE2
     >                      / 0.D0, 0.D0, 0.D0 /

      E0=SQRT((BORO*CL*1.D-9*Q/Q)**2+AM*AM)
      BRO=BORO * SCAL0()
      P = BRO*CL9*Q 
      E = SQRT(P*P+AM*AM)
      BTA = P/E
      GAMMA=E/AM
      G3=GAMMA**3
      R0=QE*CL*CL*1.D-7/(AM*1.D6)
      EKEV = (4.D0*PI/3.D0)*R0*CL*BTA*BTA*G3*ABS(B*UNITB)*1.D-3
C Fanglei, FM - 16-02-19
      IF(ABS(B).GT. 1.D-10) THEN
        RHO=BRO/B *1.D-2
      ELSE
        RHO=1.D20
      ENDIF
      ANG=2.D0*ASIN(0.5D0*(XL*1.D-2)/RHO)
      HBAR=6.6260755D-34/(2.D0*PI)
      EC=3.D0*G3*CL/(2.D0*ABS(RHO))*HBAR/QE*1.D-3
      EPHOT=8.D0/(15.D0*SQRT(3.D0))*EC
      SMELPP=SMELPP+EKEV *ABS(ANG)/(2.D0*PI)
      CALL SRLOSR(
     >            SCALE)
      IF(SCALE) SCL=SCAL0W(1.D0-SMELPP*1.D-3/E0)
      SNMPP=SNMPP+EKEV/EPHOT *ABS(ANG)/(2.D0*PI)
C Bug. FM. Sept 2015
C      SRMSE2=SRMSE2+11.d0/27.d0* EC**2  *ABS(ANG)/(2.D0*PI)
      SRMSE2=SRMSE2+11.d0/27.d0* EC**2  

      IF(NRES.LE.0) RETURN

      WRITE(NRES,FMT='(/,2X,
     >'' * Theoretical S.R. parameters in local *dipole* field :'')')
      IF(XTLT .GT. 1D-10) WRITE(NRES,FMT='(5X,''Bend plane is X-tilted''
     >,'' tilt angle is '',1P,E14.6,'' rad'')') XTLT 
      WRITE(NRES,FMT='(5X,''Bending radius (Brho/B) :'',1P,G16.8,
     > ''m,   deviation angle :'',G16.8,'' rad'')') RHO, ANG
      WRITE(NRES,FMT='(5X,''Average energy loss per particle :'',
     >'' Eloss = (2/3).r0.c.gamma^3.B.Ang/1000  ='',1P,T80,G16.8,
     >'' keV'',/,30X,
     >''(elctrn with bta~1 : 88.463*E[GeV]^4/rho[m]*(Ang/2pi))'')')
     >EKEV *ABS(ANG)/(2.D0*PI)
      WRITE(NRES,FMT='(5X,''Critical photon energy :'',
     >'' Ec = 3.gamma^3.c/(2.rho)*(Hbar/e)/1000 ='',1P,T80,G16.8,
     >'' keV'')') EC
      WRITE(NRES,FMT='(5X,''Average photon energy :'',
     >'' <Eph> = 8/(15.sqrt(3)) *Ec ='',1P,T80,G16.8,'' keV'')') EPHOT
      WRITE(NRES,FMT='(5X,''rms photon energy :'',
     >'' Eph_rms = 0.6383.Ec ='',1P,T80,G16.8,'' keV'')') 0.6383D0*EC
      WRITE(NRES,FMT='(5X,''Number of average photons per particle'',
     >'' inside dipole :'','' N = Eloss/<Eph> ='',1P,T80,G16.8)') 
     > EKEV/EPHOT *ABS(ANG)/(2.D0*PI)

      WRITE(NRES,FMT='(/,2X,'' * Theoretical S.R. parameters, '',
     >''summed over magnets up to this point : '')')
      WRITE(NRES,FMT='(5X,1P,
     >''Average energy loss : '',T80,G16.8,'' keV/particle'',6X,
     >/,15X,''- relative to initial energy :'',T80,G16.8)') 
     >SMELPP, SMELPP*1.D-3/E
      WRITE(NRES,FMT='(5X,''rms photon energy : '',
     >T80,1P,G16.8,'' keV'')') SQRT(SRMSE2)
      WRITE(NRES,FMT='(5X,''Number of average photons : '',
     >T80,1P,G16.8,'' /particle'')') SNMPP

      CALL PRSR2(SMELPP,EPHOT,SQRT(SRMSE2))

      RETURN

      ENTRY SYNPA0
        SMELPP=0.D0
        SNMPP=0.D0
        SRMSE2=0.D0
      RETURN

      ENTRY SYNPA1(
     >             SMELPO)
      SMELPO = SMELPP*1.D-3
      RETURN

      ENTRY SYNPA3(LUN,
     >                SMELPO,EO)

      SMELPO = SMELPP*1.D-3
      EO = E

      IF(LUN.GT.0) 
     >WRITE(LUN,FMT='(/,5X,''Average energy loss per particle, summed''
     >,'' UP TO THIS POINT :'', 1P,G16.8,'' keV.'',3X,
     >''Relative to initial energy :'',G16.8,/)') SMELPP, SMELPP*1.D-3/E

      RETURN

      END
