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
 
C------         R_ref    +dp/p     -dp/p
      DIMENSION R(6,6), RPD(6,6), RMD(6,6) 
      DIMENSION T(6,6,6)
      DIMENSION T3(5,6) , T4(5,6)
      SAVE R,T, T3,       T4

      LOGICAL PRBEAM
      SAVE PRBEAM

C------        Beam_ref    +dp/p     -dp/p
      DIMENSION F0(6,6), F0PD(6,6), F0MD(6,6) 

      DIMENSION F0P(6,6)

      DATA PRBEAM / .FALSE. /

      IF(.NOT. (KOBJ.EQ.5 .OR. KOBJ.EQ.6)) THEN
        WRITE(NRES,FMT='('' Matrix  cannot  be  computed :  need "OBJET" 
     >  with  KOBJ=5 or 6'')')
        RETURN
      ENDIF

      IF(KOBJ .EQ. 5) THEN
        IORD=1
      ELSEIF(KOBJ .EQ. 6) THEN
        IORD=2
      ENDIF 

      IFOC = 0
 
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
          CALL REFER(1,IORD,IFC,IT1,IT2,IT3)
          CALL MAT1(R,T,IT1)
          CALL MKSA(IORD,R,T,T3,T4)
          CALL MATIMP(R)
          CALL REFER(2,IORD,IFC,IT1,IT2,IT3)
          IF(IREF.LT.NBREF) GOTO 1
          
      ELSEIF(IORD .EQ. 2) THEN

        IFC = IFOC
        CALL REFER(1,IORD,IFC,1,6,7)
        CALL MAT2(R,T,T3,T4)
        CALL MKSA(IORD,R,T,T3,T4)
        CALL MATIMP(R)
        CALL MATIM2(R,T,T3)
        CALL REFER(2,IORD,IFC,1,6,7)
      ENDIF

        IF(PRBEAM) THEN
C------- Transported beam matrix. Initial beam input with OBJET, Kobj=5.ij or 6.ij
          CALL BEAMAT(R, 
     >                  F0P)
          SUM=0.D0
          DO 2 J=1,6
            DO 2 I=1,6
 2            SUM = SUM + F0P(I,J)
          IF(SUM.NE.0.D0) CALL BEAIMP(F0P)
        ENDIF

      RETURN

      ENTRY MATRI2(KOBJ2)
      PRBEAM = KOBJ2.EQ.1
      RETURN

      END
