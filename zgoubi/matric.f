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
      SUBROUTINE MATRIC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ------------------------------------
C     Compute transfer matrix coefficients
C     ------------------------------------
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
C      COMMON/DON/ A(09876,99),IQ(09876),IP(09876),NB,NOEL
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
 
      LOGICAL PRDIC

C------         R_ref    +dp/p     -dp/p
      DIMENSION R(6,6), RPD(6,6), RMD(6,6) 
      DIMENSION T(6,6,6)
      DIMENSION T3(5,6) , T4(5,6)
      SAVE R,T, T3,       T4

      LOGICAL PRBEAM
      SAVE PRBEAM

C------        Beam_ref    +dp/p     -dp/p
      DIMENSION F0(6,6), F0PD(6,6), F0MD(6,6) 

      LOGICAL KWRMAT
      DIMENSION F0P(6,6)

      DATA KWRMAT / .FALSE. /

      DATA PRBEAM / .FALSE. /

      IF(NRES.LE.0) RETURN

      IF(.NOT. (KOBJ.EQ.5 .OR. KOBJ.EQ.6)) THEN
        WRITE(NRES,FMT='('' Matrix  cannot  be  computed :  need "OBJET" 
     >  with  KOBJ=5 or 6'')')
        RETURN
      ENDIF

      IORD = A(NOEL,1)
      IF(IORD .EQ. 0) THEN
        CALL IMPTRA(1,IMAX,NRES)
        RETURN
      ENDIF
      IF(KOBJ .EQ. 5) THEN
        IORD=1
      ELSEIF(KOBJ .EQ. 6) THEN
        IORD=2
      ENDIF 

      IFOC = A(NOEL,2) 
      PRDIC = IFOC .GT. 10
 
      KWRMAT = NINT(A(NOEL,3)) .EQ. 1
      IF(KWRMAT) CALL MATIM6(KWRMAT)

      IF    (IORD .EQ. 1) THEN
        CALL OBJ51(NBREF)
        IREF = 0
 1      CONTINUE
          IREF = IREF + 1

          IT1 = 1 + 11 * (IREF-1)
          IT2 = IT1+3
          IT3 = IT1+4

C FM, Nov. 2008
C          CALL REFER(1,IORD,IFOC,IT1,IT2,IT3)
          IFC = IFOC
          IF(PRDIC) IFC = 0
          CALL REFER(1,IORD,IFC,IT1,IT2,IT3)
          CALL MAT1(R,T,IT1)
          CALL MKSA(IORD,R,T,T3,T4)
          CALL MATIMP(R)
          IF(PRDIC) CALL TUNES(R,F0,IFOC-10,IERY,IERZ,.TRUE.,
     >                                                YNU,ZNU,CMUY,CMUZ)
C FM, Nov. 2008
C          CALL REFER(2,IORD,IFOC,IT1,IT2,IT3)
          CALL REFER(2,IORD,IFC,IT1,IT2,IT3)
          IF(IREF.LT.NBREF) GOTO 1
          
      ELSEIF(IORD .EQ. 2) THEN

C FM, Nov. 2008
C        CALL REFER(1,IORD,IFOC,1,6,7)
        IFC = IFOC
        IF(PRDIC) IFC = 0
        CALL REFER(1,IORD,IFC,1,6,7)
        CALL MAT2(R,T,T3,T4)
        CALL MKSA(IORD,R,T,T3,T4)
        CALL MATIMP(R)
        IF(PRDIC) CALL TUNES(R,F0,IFOC-10,IERY,IERZ,.TRUE.,
     >                                                YNU,ZNU,CMUY,CMUZ)
        CALL MATIM2(R,T,T3)
        IF(PRDIC) THEN 
          CALL MAT2P(RPD,DP)
          CALL MKSA(IORD,RPD,T,T3,T4)
          CALL MATIMP(RPD)
          CALL TUNES(RPD,F0PD,IFOC-10,IERY,IERZ,.TRUE.,
     >                                              YNUP,ZNUP,CMUY,CMUZ)
          CALL MAT2M(RMD,DP)
          CALL MKSA(IORD,RMD,T,T3,T4)
          CALL MATIMP(RMD)
          CALL TUNES(RMD,F0MD,IFOC-10,IERY,IERZ,.TRUE.,
     >                                              YNUM,ZNUM,CMUY,CMUZ)
C Momentum detuning
          NUML = 1
C          DNUYDP = (YNUP-YNUM)/2.D0/A(NUML,25)
C          DNUZDP = (ZNUP-ZNUM)/2.D0/A(NUML,25)
          DNUYDP = (YNUP-YNUM)/2.D0/DP
          DNUZDP = (ZNUP-ZNUM)/2.D0/DP
          IF(NRES .GT. 0) WRITE(NRES,FMT='(/,34X,'' Chromaticities : '',
     >      //,30X,''dNu_y / dp/p = '',G14.8,/, 
     >         30X,''dNu_z / dp/p = '',G14.8)') DNUYDP, DNUZDP
C             write(nres,*) dp, a(numl,25)
        ENDIF
C        CALL REFER(2,IORD,IFOC,1,6,7)
        CALL REFER(2,IORD,IFC,1,6,7)
      ENDIF

      IF(.NOT. PRDIC) THEN
        IF(PRBEAM) THEN
C------- Transported beam matrix. Initial beam input with OBJET, Kobj=5.1 or 6.1
          CALL BEAMAT(R, 
     >                  F0P)
          SUM=0.D0
          DO 2 J=1,6
            DO 2 I=1,6
 2            SUM = SUM + F0P(I,J)
          IF(SUM.NE.0.D0) CALL BEAIMP(F0P)
        ENDIF
      ENDIF

      RETURN

      ENTRY MATRI2(KOBJ2)
      PRBEAM = KOBJ2.EQ.1
      RETURN

      END
