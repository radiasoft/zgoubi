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
      SUBROUTINE AUTORF
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ------------------------------------------------
C     CHANGEMENT AUTOMTIQ DE REFERENTIEL
C     ------------------------------------------------
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.CONST_2.H"   ! COMMON/CONST/ CL9,CL,PI,RAD,DEG,QEL,AMPROT,CM2M
      INCLUDE "C.CONST2.H"  ! COMMON/CONST2/ ZERO, UN
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"   ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      INCLUDE "C.OBJET.H"   ! COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT,KZOB
      INCLUDE "C.SPIN.H"    ! COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)
C      LOGICAL ZSYM
      INCLUDE "C.TYPFLD.H"  ! COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM

      LOGICAL EVNT
      SAVE EVNT
      PARAMETER (MSR=8)
      CHARACTER(2) QSHRO(MSR)
      DIMENSION VSHRO(MSR)

      DATA EVNT / .FALSE. /

      IOP = NINT(A(NOEL,1))
      IOP2 = NINT(A(NOEL,2))

      IF    (IOP .LE. 3) THEN

        IF(IOP .EQ. 0) GOTO 99

        IF    (IOP .EQ. 1) THEN
C--------- POSITIONNEMENT REFERENCE = TRAJ. 1, ABSCISSE ACTUELLE
          IRF=1
          XC =ZERO
          YC =F(2,IRF)
          AA  =F(3,IRF) * 0.001D0
          IF(IOP2 .EQ. 5) THEN
            ZC =F(4,IRF)
            BB = F(5,IRF) * 0.001D0
          ENDIF
        ELSEIF(IOP .EQ. 2) THEN
C--------- POSITIONNEMENT REFERENCE = WAIST DES TRAJ. 1 ET 4-5
          IRF=1
          CALL FOCAL1(IRF,4,5,
     >                        XC,YC )
          AA  =F(3,IRF) * 0.001D0
        ELSEIF(IOP .EQ. 3) THEN
C--------- POSITIONNEMENT REFERENCE = WAIST DES TRAJ. IRF ET MX1-MX2
          IRF=NINT(A(NOEL,10))
          MX1=NINT(A(NOEL,11))
          MX2=NINT(A(NOEL,12))
          CALL FOCAL1(IRF,MX1,MX2,
     >                            XC,YC )
          AA  =F(3,IRF) * 0.001D0
        ENDIF

      ELSEIF(IOP .EQ. 4) THEN

        XC = ZERO
        YC = ZERO
        AA = ZERO
        II = 0
C First center the beam on Y=0, T=0, D
        DO I = 1, IMAX
          IF( IEX(I) .GT. 0) THEN
            II = II + 1
            YC = YC + F(2,I)
            AA = AA + F(3,I)
          ENDIF
        ENDDO

        YC = YC / DBLE(II)
        AA = AA / DBLE(II) * 1.D-3

        DO I=1,IMAX
C         +++ IEX<-1 <=> Particule stoppee
          IF( IEX(I) .GE. -1) THEN
            IF(I .EQ. IREP(I) .OR. .NOT. ZSYM) THEN
              CALL INITRA(I)
              CALL CHAREF(.FALSE.,XC,YC,AA)
              CALL MAJTRA(I)
            ELSE
              CALL DEJACA(I)
            ENDIF

          ENDIF
        ENDDO
        IF(NRES .GT. 0) THEN
          WRITE(NRES,FMT='(1X,A)') ' First,  center  beam  on'
     >    // ' Y=0,  T=0 : '
          WRITE(NRES,100) XC,YC,AA*DEG,AA
          WRITE(NRES,101) IEX(1),(F(J,1),J=1,7)
          WRITE(NRES,FMT='(/,1X,A)') ' Then,  move  beam  to '
     >    // ' new  centroid :'
        ENDIF

C Then update to requested beam centering coordinates
        XC =  -A(NOEL,10)
        YC =  -A(NOEL,11)
        AA =  -A(NOEL,12) * 1.D-3

      ELSEIF(IOP .EQ. 5) THEN

        ZC = ZERO
        BB = ZERO
        II = 0
C First center the beam on Z=0, P=0
        DO I = 1, IMAX
          IF( IEX(I) .GT. 0) THEN
            II = II + 1
            ZC = ZC + F(4,I)
            BB = BB + F(5,I)
          ENDIF
        ENDDO
        ZC = ZC / DBLE(II)
        BB = BB / DBLE(II)
        DO I=1,IMAX
C         +++ IEX<-1 <=> Particule stoppee
          IF( IEX(I) .GE. -1) THEN
            F(4,I) = F(4,I) - ZC
            F(5,I) = F(5,I) - BB
          ELSE
          ENDIF
        ENDDO
        IF(NRES .GT. 0) THEN
          BB = BB *1.D-3
          WRITE(NRES,107) ZC,BB*DEG,BB
          WRITE(NRES,101) IEX(1),(F(J,1),J=1,7)
        ENDIF

C Then update to requested beam centering coordinates
        ZC =  -A(NOEL,10)
        BB =  -A(NOEL,11) * 1.D-3

      ELSE
        CALL ENDJOB('Sbr autorf. No such option I = ',IOP)
      ENDIF

      IF    (IOP .LE. 4) THEN
        DO I=1,IMAX
C         +++ IEX<-1 <=> Particule stoppee
          IF( IEX(I) .GE. -1) THEN
            IF(I .EQ. IREP(I) .OR. .NOT.ZSYM) THEN
              CALL INITRA(I)
              CALL CHAREF(.FALSE.,XC,YC,AA)
              CALL MAJTRA(I)
            ELSE
              CALL DEJACA(I)
            ENDIF
C            IF(IOP2.EQ.1) F(1,I) = F(1,I) - DD

            IF(IOP2.EQ.5) THEN

              VSHRO(MSR) = 2
              VSHRO(1) = ZC
              VSHRO(2) = BB
              QSHRO(1) = 'ZS'
              QSHRO(2) = 'YR'
              CALL INITRA(I)
              CALL CHANRF(EVNT,QSHRO,VSHRO)
              CALL MAJTRA(I)
              IF(KSPN .EQ. 1 ) then
                if(nres .gt. 0) write(nres,*)
     >          ' WARNING : Y-rotation in autorf.f. Spin '
     >          //'rotation is to be implemented. '
              ENDIF
            ENDIF

          ENDIF
        ENDDO

            IF(IOP2.EQ.1) THEN
              II = 0
              DD = 0.D0
              TT = 0.D0
              DO I = 1, IMAX
               IF( IEX(I) .GT. 0) THEN
                II = II + 1
                DD = DD + F(1,I)
                TT = TT + F(7,I)
               ENDIF
              ENDDO
              DD = DD / DBLE(II)
              TT = TT / DBLE(II)
              DO I = 1, IMAX
                IF( IEX(I) .GT. 0) THEN
                  F(1,I) = F(1,I) - DD + A(NOEL,13)
                  F(7,I) = F(7,I) - TT + A(NOEL,14)    ! case eRHIC linac
                ENDIF
              ENDDO
            ELSEIF(IOP2.EQ.2) THEN
              II = 0
              DD = 0.D0
              DO I = 1, IMAX
               IF( IEX(I) .GT. 0) THEN
                II = II + 1
                DD = DD + F(1,I)
               ENDIF
              ENDDO
              DD = DD / DBLE(II)
              DO I = 1, IMAX
                IF( IEX(I) .GT. 0) THEN
                  F(1,I) = F(1,I) - DD + A(NOEL,13)
                  F(7,I) = A(NOEL,14)    ! case eRHIC linac
                ENDIF
              ENDDO
            ENDIF

      ELSEIF(IOP .EQ. 5) THEN
        VSHRO(MSR) = 2
        VSHRO(1) = ZC
        VSHRO(2) = BB
        QSHRO(1) = 'ZS'
        QSHRO(2) = 'YR'
        DO I=1,IMAX
C         +++ IEX<-1 <=> Particule stoppee
          IF( IEX(I) .GE. -1) THEN
            CALL INITRA(I)
            CALL CHANRF(EVNT,QSHRO,VSHRO)
            CALL MAJTRA(I)
          ENDIF
        ENDDO
      ENDIF

 99   CONTINUE

      IF(NRES .GT. 0) THEN

        IF    (IOP .EQ. 0) THEN
          WRITE(NRES,102)
 102      FORMAT(/,20X,' AUTOREF  is  off')
        ELSEIF(IOP .LE. 5) THEN
          IF(IOP .LE.4) THEN
            WRITE(NRES,100) XC,YC,AA*DEG,AA
 100        FORMAT(/,' Change  of  reference,  horizontal,   XC =',
     >      F15.8,' cm , YC =',
     >      F15.8,' cm ,   A =',F16.9,' deg  (i.e., ',F14.10,' rad)',/)
C 100      FORMAT(/,' CHANGEMENT  DE  REFERENCE  XC =',F9.3,' cm , YC =',
C     >     F10.3,' cm ,   A =',F12.5,' deg  (i.e., ',F10.6,' rad)',/)
            IF(IOP2.EQ.1) THEN
              WRITE(NRES,FMT='(/,'' Beam centerd on momentum p/p_Ref ''
     >        ,1P,E16.8,''  and tine '',E16.8,/)') A(NOEL,13),A(NOEL,14)
            ENDIF
          ELSEIF(IOP .EQ. 5) THEN
            WRITE(NRES,107) ZC,BB*DEG,BB
 107        FORMAT(/,' Change  of  reference,  vertical,   ZC =',
     >      F15.8,' cm,  A =',F16.9,' deg  (i.e., ',F14.10,' rad)',/)
          ENDIF
          WRITE(NRES,101) IEX(1),(F(J,1),J=1,7)
  101     FORMAT(' TRAJ 1 IEX,D,Y,T,Z,P,S,time :',I3,1P,5G12.4,2G17.5)
        ENDIF

      ENDIF

      RETURN
      END

