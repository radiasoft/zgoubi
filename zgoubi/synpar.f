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
      SUBROUTINE SYNPAR(B,XL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      COMMON/PTICUL/ AM,Q,G,TO
      COMMON/RIGID/ BORO,DPREF,DP,BR
      LOGICAL SCALE
      DATA UNITB / 1.D-1/
      DATA OCCUR, SMELPP, SNMPP, SRMSE2
     >                      / 0.D0, 0.D0, 0.D0 , 0.D0/
      E0=SQRT((BORO*CL*1.D-9*Q/Q)**2+AM*AM)
      BRO=BORO * SCAL0()
C      P = BRO*CL*1.D-9*Q/QE  
C      P = BRO*CL*1.D-9*Q 
      P = BRO*CL9*Q 
      E = SQRT(P*P+AM*AM)
      BTA = P/E
      GAMMA=E/AM
      G3=GAMMA**3
      R0=QE*CL*CL*1.D-7/(AM*1.D6)
      EKEV = (4.D0*PI/3.D0)*R0*CL*BTA*BTA*G3*ABS(B*UNITB)*1.D-3
      RHO=BRO/B *1.D-2
      ANG=2.D0*ASIN(0.5D0*(XL*1.D-2)/RHO)
      HBAR=6.6260755D-34/(2.D0*PI)
      EC=3.D0*G3*CL/(2.D0*ABS(RHO))*HBAR/QE*1.D-3
      EPHOT=8.D0/(15.D0*SQRT(3.D0))*EC
      SMELPP=SMELPP+EKEV *ABS(ANG)/(2.D0*PI)
      CALL SRLOSR(SCALE)
      IF(SCALE) SCL=SCAL0W(1.D0-SMELPP*1.D-3/E0)
      SNMPP=SNMPP+EKEV/EPHOT *ABS(ANG)/(2.D0*PI)
C      SRMSE2=SRMSE2 + 11.d0/27.d0* EC**2  *ABS(ANG)/(2.D0*PI)
      SRMSE2=SRMSE2+11.d0/27.d0* EC**2  *ABS(ANG)/(2.D0*PI)
C      SRMSE2=SRMSE2 + (6.72d-14*GAMMA**2.5/RHO)**2
      OCCUR=OCCUR+1
      IF(NRES.LE.0) RETURN
      WRITE(NRES,FMT='(/,2X,
     >'' * Theoretical S.R. parameters in the *dipole* field :'',/)')
      WRITE(NRES,FMT='(5X,'' Deviation Ang. :'',1P,G12.4,
     > ''rad.    Bending radius :'',G12.4,''m'')') ANG,RHO
      WRITE(NRES,FMT='(5X,'' Mean energy loss per particle :'',
     >'' Eloss = (2/3).r0.c.gamma^3.B/1000 .Ang ='',1P,T80,G12.4,
     >'' keV'')') EKEV *ABS(ANG)/(2.D0*PI)
      WRITE(NRES,FMT='(5X,'' Critical energy :'',
     >'' Ec = 3.gamma^3.c/(2.rho)*(Hbar/e)/1000 ='',1P,T80,G12.4,
     >'' keV'')') EC
      WRITE(NRES,FMT='(5X,'' Mean energy of radiated photons :'',
     >'' <Eph> = 8/(15.sqrt(3)).Ec ='',1P,T80,G12.4,'' keV'')') EPHOT
      WRITE(NRES,FMT='(5X,'' rms energy of radiated photons :'',
C     >'' Eph_rms = 0.5591.Ec ='',1P,T80,G12.4,'' keV'')') .5591D0*EC
     >'' Eph_rms = 0.6383.Ec ='',1P,T80,G12.4,'' keV'')') .6383D0*EC
      WRITE(NRES,FMT='(5X,'' Number of mean photons per particle'',
     >'' inside dipole :'','' N = Eloss/<Eph> ='',1P,T80,G12.4)') 
     > EKEV/EPHOT *ABS(ANG)/(2.D0*PI)

      WRITE(NRES,FMT='(//,5X,'' Mean energy loss per particle, summed'',
     >'' UP TO THIS MAGNET :'',1P,G12.4,'' keV'',6X,
     >''Relative to initial energy :'',G12.4)') SMELPP, SMELPP*1.D-3/E
      WRITE(NRES,FMT='(5X,'' # of mean photons per particle, summed'',
     >'' UP TO THIS MAGNET :'',1P,G12.4)') SNMPP
      WRITE(NRES,FMT='(5X,'' rms energy of radiated photons'',
     > 1P,G12.4)') sqrt(SRMSE2)

C      write(30,fmt='(2I5,1P,3G12.4)') 
C     >  noel,occur,SMELPP,SNMPP,sqrt(SRMSE2)

      RETURN

      ENTRY SYNPA0
        SMELPP=0.D0
        SNMPP=0.D0
        SRMSE2=0.D0
        OCCUR=0
      RETURN

C      ENTRY SYNPAI(SCLIN)
C      SCL=SCLIN
C      RETURN

      END
