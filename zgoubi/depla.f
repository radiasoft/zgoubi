C  ZGOUBI, a program for computing the trajectories of charged particles
C  in electric and magnetic fields
C  Copyright (C) 1988-2007  Fran�ois M�ot
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
C  Fran�ois M�ot <fmeot@bnl.gov>
C  Brookhaven National Laboratory   
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  -------
      SUBROUTINE DEPLA(DS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ----------------------------
C     CALCUL D'UN PAS EN CARTESIEN
C     ----------------------------
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      COMMON/DEPL/ XF(3),DXF(3),DQBRO,DTAR
      LOGICAL ZSYM
      COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      COMMON/ORDRES/ KORD,IRD,IDS,IDB,IDE,IDZ
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      COMMON/TRAJ/ Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT
      COMMON/VITES/ U(18),DQBR(6),DDT(6)

      PARAMETER (IMX=6)

      LOGICAL VARSTP,VARSTI
      SAVE VARSTP, PREC

      DATA VARSTP, PREC / .FALSE., 1D-2 /

 11   CONTINUE

CALCUL VECTEURS dR ET dU RESULTANT DE DS
 
      DO K=1,3
        XF(K)=0.D0
        DXF(K)=0.D0
      ENDDO

      DV=1.D0
      DO 3 I=1,IDS
         V=DV*DS/DBLE(I)
         IK = I - IMX
         DO 1 K=1,3
            IK = IK + IMX
            DXF(K)=DXF(K)+U(IK)*DV
            XF(K)=XF(K)+U(IK)*V
 1       CONTINUE
         DV=V
 3    CONTINUE

      IF(KFLD .GE. LC) THEN
C------ Electric or elctro-mag lmnt

        DQBRO=0.D0
        DTAR = 0.D0
        V=1.D0
        DO 31 I=1,IDS
           V=V*DS/DBLE(I)
           DQBRO = DQBRO+DQBR(I)*V
           DTAR = DTAR+DDT(I)*V
 31     CONTINUE

        IF (VARSTP) THEN

C----- Tests implementation EL[T]MIR
C          CALL VTHET(U(1))
C          CALL ENRGY(Error)
          DANG2 = 
     >    ( (DXF(1)-U(1))*(DXF(1)-U(1)) +
     >    (DXF(2)-U(7))*(DXF(2)-U(7)) +
     >    (DXF(3)-U(13))*(DXF(3)-U(13)) )

          IF(DQBRO*DQBRO .GT. PREC
     >      .OR.  DANG2 .GT. PREC
     >                             ) THEN
            DS = DS * 0.5D0
            GOTO 11
          ENDIF

        ENDIF

      ELSE
         QBRO = QBR*CL9
         DTAR = DS / (QBRO/SQRT(QBRO*QBRO+AMT*AMT)*CL9) 
      ENDIF

      RETURN

      ENTRY DEPLAW(VARSTI,IPREC)
      VARSTP = VARSTI
      PREC = 10.D0**(-IPREC)
      RETURN

      END
