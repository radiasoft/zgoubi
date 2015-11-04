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
      SUBROUTINE OBJ102
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     **************************************
C     CONSTITUTION DE L'OBJET INITIAL KOBJ=6
C     **************************************
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON_2.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IIP(MXL),NB,NOEL
      PARAMETER (LNTA=132) ; CHARACTER(LNTA) TA
      PARAMETER (MXTA=45)
      INCLUDE "C.DONT.H"     ! COMMON/DONT/ TA(MXL,MXTA)
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      CHARACTER(1) LET
      INCLUDE "C.FAISCT.H"     ! COMMON/FAISCT/ LET(MXT)
C----- KAR: LETTRES AFFECTEES AUX TRAJECTOIRES ( 'S'  EST RESERVEE
C      POUR ETIQUETER LES PARTICULES SECONDAIRES -OPTION 'MCDESINT')
      CHARACTER(1) KAR(41)
      INCLUDE "C.KAR.H"     ! COMMON/KAR/ KAR
      INCLUDE "C.OBJET.H"     ! COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      INCLUDE "C.REBELO.H"   ! COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
 
      PARAMETER(MXJ1=MXJ-1)

      DIMENSION P(MXJ)
 
      DIMENSION REF(MXJ)

      IMAX = 102 
      IDMAX= 1

      P(2) = A(NOEL,20)
      P(3) = A(NOEL,21)
      P(4) = A(NOEL,22)
      P(5) = A(NOEL,23)
      P(6) = A(NOEL,24)
      P(1) = A(NOEL,25)

      K = 30
      DO 62 J = 2,MXJ1
        REF(J) = A(NOEL,K)
        K = K + 1
 62   CONTINUE
      REF(1) = A(NOEL,35)

      call ray102(p)

      DO 66 I=1,IMAX
         IEX(I) =1
         IREP(I)=I
         IKAR = IKAR+1
         IF(IKAR.GT.41)  IKAR=1
         LET(I)=KAR(IKAR)
         DO 66 J=1,6
            FO(J,I) = FO(J,I) + REF(J)
            F(J,I)=FO(J,I)
 66   CONTINUE
      IF(NRES.GT.0) THEN
        WRITE(NRES,100) KOBJ,IMAX
        WRITE(NRES,102) (P(J),J=1,MXJ1)
        WRITE(NRES,FMT='(/,19X,'' Reference trajectory : '',
     >    4X,F9.6,6(4X,F9.6),/)') (REF(I), I=1, MXJ1)
      ENDIF
  100 FORMAT(/,41X,'CALCUL  DES  TRAJECTOIRES',//,30X,'OBJET  (',I1,
     1')  FORME  DE ',I6,' POINTS ')
  102 FORMAT(/,19X,' Sampling : ',16X,F9.6,5(4X,F9.6))
      RETURN
      END
