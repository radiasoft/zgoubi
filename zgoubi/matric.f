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
C  François Meot <fmeot@bnl.gov>
C  Brookhaven National Laboratory
C  C-AD, Bldg 911
C  Upton, NY, 11973, USA
C  -------
      SUBROUTINE MATRIC(JORD,JFOC,KWR,SCPLD,
     >                                      IER)
C      SUBROUTINE MATRIC(JORD,JFOC,KWR,KCPLD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ------------------------------------
C     Compute transfer matrix coefficients
C     ------------------------------------
      CHARACTER(*) SCPLD
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
C      INCLUDE 'MXLD.H'
C      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
C      COMMON/DON/ A(09876,99),IQ(09876),IP(09876),NB,NOEL
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
C      INCLUDE "C.OBJET.H"     ! COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT,KZOB
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),IREP(MXT),AMQLU,PABSLU

C------         R_ref    +dp/p     -dp/p
      DIMENSION R(6,6), RPD(6,6), RMD(6,6)
      DIMENSION T(6,6,6)
      DIMENSION T3(5,6) , T4(5,6)
      SAVE R,T, T3,       T4

C------        Beam_ref    +dp/p     -dp/p
      DIMENSION F0(6,6), F0PD(6,6), F0MD(6,6)

      LOGICAL OKCPLD
      LOGICAL KWRMAT
      LOGICAL PRDIC

      DIMENSION RO(6,6)

      DATA KWRMAT / .FALSE. /

      IER = 0
      CALL OBJET1(
     >            KOBJ,KOBJ2)
      IF(.NOT. (KOBJ.EQ.5 .OR. KOBJ.EQ.6)) THEN
        IF(NRES.GT.0) THEN
          WRITE(NRES,FMT='('' Matrix  cannot  be  computed :  '',
     >    ''need "OBJET" with  KOBJ=5 or 6'')')
          CALL IMPTRA(1,IMAX,NRES)
          IER = 1
          RETURN
        ENDIF
      ENDIF

C      IORD = A(NOEL,1)
      IORD = JORD
      IF(IORD .EQ. 0) THEN
        IF(NRES.GT.0) THEN
          WRITE(NRES,FMT='(/,9X,'' Matrix not computed : IORD = 0'',/)')
          CALL IMPTRA(1,IMAX,NRES)
        ENDIF
        RETURN
      ENDIF
      IF(KOBJ .EQ. 5) THEN
        IORD=1
      ELSEIF(KOBJ .EQ. 6) THEN
        IF(KOBJ2 .EQ. 0) THEN
          IORD=2
        ELSEIF(KOBJ2 .EQ. 1) THEN
          IORD=3
        ENDIF
      ENDIF

      IFOC = JFOC
      PRDIC = IFOC .GT. 10
      NMAIL = IFOC-10

      KWRMAT = KWR .EQ. 1
      CALL MATIM6(KWRMAT)
      IF(KWRMAT .AND. NRES .GT. 0) THEN
        WRITE(NRES,*)
        WRITE(NRES,*)  ' Matrix coefficients are printed in '
     >  // ' zgoubi.MATRIX.out.'
      ENDIF

      OKCPLD = SCPLD .EQ. 'coupled'
      CALL MATIM4(OKCPLD)

      IF    (IORD .EQ. 1) THEN
        CALL OBJ51(
     >             NBREF)
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
          CALL MAT1(IT1,F,IMAX,
     >                  R,T)
          CALL MKSA(IORD,R,T,T3,T4)
          IF(PRDIC) THEN
            IF(OKCPLD) THEN
              CALL TUNESC(R,
     >                      F0,YNU,ZNU,CMUY,CMUZ,IERY,IERZ)
            ELSE
              CALL TUNES(R,F0,NMAIL,IERY,IERZ,.TRUE.,
     >                                                YNU,ZNU,CMUY,CMUZ)
            ENDIF
          ENDIF
          CALL MATIMP(R,F0,YNU,ZNU,CMUY,CMUZ,NMAIL,PRDIC,IT1)

C FM, Nov. 2008
          CALL REFER(2,IORD,IFC,IT1,IT2,IT3)

          IF(IREF.LT.NBREF) GOTO 1

      ELSEIF(IORD .EQ. 2) THEN

C FM, Nov. 2008
C        CALL REFER(1,IORD,IFOC,1,6,7)
        IREF = 1
        IFC = IFOC
        IF(PRDIC) IFC = 0
        CALL REFER(1,IORD,IFC,1,6,7)
        CALL MAT2(R,T,T3,T4)
        CALL MKSA(IORD,R,T,T3,T4)
        IF(PRDIC) THEN
          IF(OKCPLD) THEN
            CALL TUNESC(R,
     >                    F0,YNU,ZNU,CMUY,CMUZ,IERY,IERZ)
          ELSE
            CALL TUNES(R,F0,NMAIL,IERY,IERZ,.TRUE.,
     >                                            YNU,ZNU,CMUY,CMUZ)
          ENDIF
        ENDIF
        CALL MATIMP(R,F0,YNU,ZNU,CMUY,CMUZ,NMAIL,PRDIC,iref)
        CALL MATIM2(R,T,T3)
        IF(PRDIC) THEN
          CALL MAT2P(RPD,DP)
          CALL MKSA(IORD,RPD,T,T3,T4)
          IF(OKCPLD) THEN
            CALL TUNESC(R,
     >                    F0,YNU,ZNU,CMUY,CMUZ,IERY,IERZ)
          ELSE
            CALL TUNES(RPD,F0PD,NMAIL,IERY,IERZ,.TRUE.,
     >                                              YNUP,ZNUP,CMUY,CMUZ)
          ENDIF
          CALL MATIMP(RPD,F0PD,YNUP,ZNUP,CMUY,CMUZ,NMAIL,PRDIC,iref)
          CALL MAT2M(RMD,DP)
          CALL MKSA(IORD,RMD,T,T3,T4)
          IF(OKCPLD) THEN
            CALL TUNESC(R,
     >                    F0,YNU,ZNU,CMUY,CMUZ,IERY,IERZ)
          ELSE
            CALL TUNES(RMD,F0MD,NMAIL,IERY,IERZ,.TRUE.,
     >                                              YNUM,ZNUM,CMUY,CMUZ)
          ENDIF
          CALL MATIMP(RMD,F0MD,YNUM,ZNUM,CMUY,CMUZ,NMAIL,PRDIC,iref)
C Momentum detuning
          NUML = 1
C          DNUYDP = (YNUP-YNUM)/2.D0/A(NUML,25)
C          DNUZDP = (ZNUP-ZNUM)/2.D0/A(NUML,25)
          DNUYDP = (YNUP-YNUM)/2.D0/DP
          DNUZDP = (ZNUP-ZNUM)/2.D0/DP
          IF(NRES .GT. 0) WRITE(NRES,FMT='(/,34X,'' Chromaticities : '',
     >      //,30X,''dNu_y / dp/p = '',G15.8,/,
     >         30X,''dNu_z / dp/p = '',G15.8)') DNUYDP, DNUZDP
        ENDIF
        CALL REFER(2,IORD,IFC,1,6,7)

      ELSEIF(IORD .EQ. 3) THEN
C From Nick Tsoupas, RAYTRACE, 3rd order transport coeffs

        IREF = 1
        IFC = IFOC
        CALL REFER(1,IORD,IFC,1,6,7)
        CALL MAT3RD(NRES)
        CALL REFER(2,IORD,IFC,1,6,7)

      ENDIF

      RETURN

      ENTRY MATRI1(
     >             RO)
      DO IB = 1, 6
        DO IA = 1, 6
          RO(IA,IB)  = R(IA,IB)
        ENDDO
      ENDDO
      RETURN
      END
