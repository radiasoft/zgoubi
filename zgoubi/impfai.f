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
      SUBROUTINE IMPFAI(KPR,NOEL,KLEY,LBL1,LBL2)
C------- Called by keyword FAISTORE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER(*) KLEY
      CHARACTER(*) LBL1,LBL2
      PARAMETER(MLB=10)
      PARAMETER (LBLSIZ=20)
      CHARACTER(LBLSIZ) LBLST(MLB)
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      INCLUDE "C.CHAMBR.H"     ! COMMON/CHAMBR/ LIMIT,IFORM,YLIM2,ZLIM2,SORT(MXT),FMAG,YCH,ZCH

      INCLUDE "C.CONST_2.H"     ! COMMON/CONST/ CL9,CL,PI,RAD,DEG,QEL,AMPROT,CM2M
      INCLUDE "C.DESIN.H"     ! COMMON/DESIN/ FDES(7,MXT),IFDES,KINFO,IRSAR,IRTET,IRPHI,NDES
C     >,AMS,AMP,AM3,TDVM,TETPHI(2,MXT)
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      CHARACTER(1) LET
      INCLUDE "C.FAISCT.H"     ! COMMON/FAISCT/ LET(MXT)
      INCLUDE "C.OBJET.H"     ! COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT,KZOB
      INCLUDE "C.REBELO.H"   ! COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,HDPRF,DP,QBR,BRI
      INCLUDE "C.SPIN.H"     ! COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)
      INCLUDE "C.SYNCH.H"     ! COMMON/SYNCH/ PH(MXT), DPR(MXT), PS

      CHARACTER(80) FNAME
      LOGICAL BINARY, FITING

      CHARACTER(1) TX1
      PARAMETER (TX1='''')
      SAVE BINARY
      LOGICAL OPN
      SAVE OPN

      LOGICAL SRLOSS
      DIMENSION SRLT(MXT)

      DATA OPN / .FALSE. /

      CALL FITSTA(5,
     >              FITING)
      IF(FITING) RETURN

      IF    (KPR .GE. 1) THEN
C------- Keyword FAISTORE: store at ipass=1 & every ipass=multiple of KPR
        IF(IPASS .GT. 1) THEN
          IF(IPASS .LT. NRBLT) THEN
            IF(KPR*(IPASS/KPR) .NE. IPASS) RETURN
          ENDIF
        ENDIF
      ENDIF

C--- Case SR loss
      CALL SRLOS3(
     >            SRLOSS)
      IF(SRLOSS) THEN
        CALL RAYSY7(
     >             SRLT)
c             write(*,*) ' impfai SRLT ',SRLT(1)
c             write(*,*)
      ENDIF


C------- Dummies
      D = 0.D0
      ID = 0

C       write(*,*) ' impfai f ',(f(7,i),i=1,imax)
C           read(*,*)

      IF(BINARY) THEN
        DO 2 I=1,IMAX
            P = BORO*CL9 *F(1,I) * AMQ(2,I)
            ENERG = SQRT(P*P + AMQ(1,I)*AMQ(1,I))
            ENEKI = ENERG - AMQ(1,I)
            WRITE(NFAI)
     1      IEX(I),-1.D0+FO(1,I),(FO(J,I),J=2,MXJ),
     2      -1.D0+F(1,I),F(2,I),F(3,I),(F(J,I),J=4,MXJ),
     >      (SI(J,I),J=1,4),(SF(J,I),J=1,4),
     >      ENEKI,ENERG,
     4      I,IREP(I), SORT(I),(AMQ(J,I),J=1,5),PH(I),DPR(I),PS,
     5      BORO, IPASS, NOEL, KLEY,LBL1,LBL2,LET(I),SRLT(I),
     6      DPREF,HDPRF
 2      CONTINUE
      ELSE
        DO 1 I=1,IMAX
          P = BORO*CL9 *F(1,I) * AMQ(2,I)
          ENERG = SQRT(P*P + AMQ(1,I)*AMQ(1,I))
          ENEKI = ENERG - AMQ(1,I)
          WRITE(NFAI,110)
     1    IEX(I),-1.D0+FO(1,I),(FO(J,I),J=2,MXJ),
     2    -1.D0+F(1,I),F(2,I),F(3,I),
     3    (F(J,I),J=4,MXJ),
     4    (SI(J,I),J=1,4),(SF(J,I),J=1,4),
     5    ENEKI,ENERG,
     6    I,IREP(I), SORT(I),(AMQ(J,I),J=1,5),PH(I),DPR(I),PS,
     7    BORO, IPASS, NOEL,
     8    TX1,KLEY,TX1,TX1,LBL1,TX1,TX1,LBL2,TX1,TX1,LET(I),TX1,
     9    SRLT(I),
     X    DPREF,HDPRF
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

      ENTRY IMPFAW(FNAME,LBLST,NLB)

C Bug found by Sam T, July 2015 : FAISTORE may come after REBELOTE and ipass>1
C      IF(.NOT. OPN) INQUIRE(FILE=FNAME,OPENED=OPN)
C      IF(IPASS .EQ. 1) CALL OPEN2('FAISCN',NFAI,FNAME)
C      IF(.NOT. OPN) CALL OPEN2('FAISCN',NFAI,FNAME)
      CALL OPEN2('FAISCN',NFAI,FNAME)
      BINARY=FNAME(1:2).EQ.'B_' .OR. FNAME(1:2).EQ. 'b_'
      IF(NRES .GT. 0) THEN
        WRITE(NRES,FMT='(15X,
     >    ''Print will occur at element[s] labeled : '')')
        WRITE(NRES,FMT='(20X,A)') (LBLST(I),I=1,NLB)
      ENDIF

      RETURN
      END
