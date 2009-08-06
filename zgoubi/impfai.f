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
      SUBROUTINE IMPFAI(KPR,NOEL,KLEY,LBL1,LBL2)
C------- Called by keyword FAISTORE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER KLEY*10
      CHARACTER*10 LBL1,LBL2
      PARAMETER(MLB=10)
      CHARACTER*10 LBL(MLB)
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      COMMON/CHAMBR/ LIMIT,IFORM,YLIM2,ZLIM2,SORT(MXT),FMAG,BMAX
     > ,YC,ZC
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      COMMON/DESIN/ FDES(7,MXT),IFDES,KINFO,IRSAR,IRTET,IRPHI,NDES
     >,AMS,AMP,AM3,TDVM,TETPHI(2,MXT)
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      CHARACTER LET
      COMMON/FAISCT/ LET(MXT)
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      COMMON/REBELO/ NPASS,IPASS,KWRT,NNDES,STDVM
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      COMMON/SYNCH/ RET(MXT), DPR(MXT),PS
 
      CHARACTER FNAME*80
      LOGICAL BINARY, FITING
      SAVE BINARY

      CALL FITSTA(5,FITING)
      IF(FITING) RETURN

      IF    (KPR .GE. 1) THEN
C------- Keyword FAISTORE: store at ipass=1 & every ipass=multiple of KPR
        IF(IPASS .GT. 1) THEN
          IF(IPASS .LT. NPASS) THEN
            IF(KPR*(IPASS/KPR) .NE. IPASS) RETURN
          ENDIF
        ENDIF
      ENDIF
 
C------- Dummies
        D = 0.D0
        ID = 0

      IF(BINARY) THEN
        DO 2 I=1,IMAX
            P = BORO*CL9 *F(1,I) * AMQ(2,I)
            ENERG = SQRT(P*P + AMQ(1,I)*AMQ(1,I))
            ENEKI = ENERG - AMQ(1,I)
            WRITE(NFAI)
     1      LET(I),IEX(I),-1.D0+FO(1,I),(FO(J,I),J=2,MXJ),
     2      -1.D0+F(1,I),F(2,I),F(3,I),
     3      (F(J,I),J=4,MXJ),ENEKI,
     4      I,IREP(I), SORT(I),(AMQ(J,I),J=1,5),RET(I),DPR(I),
     5                          BORO, IPASS, KLEY,LBL1,LBL2,NOEL
 2      CONTINUE
      ELSE
        DO 1 I=1,IMAX
            P = BORO*CL9 *F(1,I) * AMQ(2,I)
            ENERG = SQRT(P*P + AMQ(1,I)*AMQ(1,I))
            ENEKI = ENERG - AMQ(1,I)
            WRITE(NFAI,110)
     1      LET(I),IEX(I),-1.D0+FO(1,I),(FO(J,I),J=2,MXJ),
     2      -1.D0+F(1,I),F(2,I),F(3,I),
     3      (F(J,I),J=4,MXJ),ENEKI,
     4      I,IREP(I), SORT(I),(AMQ(J,I),J=1,5),RET(I),DPR(I),
     5                          BORO, IPASS, KLEY,LBL1,LBL2,NOEL
          INCLUDE "FRMFAI.H"
 1      CONTINUE
      ENDIF

CC--------- NuFact, jaroslaw ------------------------------------------
C      itot = 0
C      DO 3 I=1,IMAX
C        if(iex(i) .LE. -1) goto 3
C        itot = itot + 1
C 3    CONTINUE
CC---------------------------------------------------------------------

      CALL FLUSH2(NFAI,BINARY)

      RETURN

      ENTRY IMPFAW(FNAME,LBL,NLB)

      IF(IPASS .EQ. 1) CALL OPEN2('FAISCN',NFAI,FNAME)
C            write(*,*) ' impfaw ',NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      BINARY=FNAME(1:2).EQ.'B_' .OR. FNAME(1:2).EQ. 'b_'
      IF(NRES .GT. 0) THEN
        WRITE(NRES,FMT='(15X,
     >    ''Print will occur at element[s] labeled : '')') 
        WRITE(NRES,FMT='(20X,A)') (LBL(I),I=1,NLB)
      ENDIF

      RETURN
      END
