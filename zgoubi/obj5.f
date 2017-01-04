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
      SUBROUTINE OBJ5
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     **************************************
C     CONSTITUTION DE L'OBJET INITIAL KOBJ=5
C     **************************************
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON_2.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IIP(MXL),NB,NOEL
C      PARAMETER (LNTA=132) ; CHARACTER(LNTA) TA
C      PARAMETER (MXTA=45)
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
C      DIMENSION IDE(MXJ),JDE(MXJ),P(MXJ)
C      EQUIVALENCE (IDE(2),IYMAX),(IDE(3),ITMAX),(IDE(4),IZMAX),
C     > (IDE(5),IPMAX),(IDE(1),IMAXD)
C      EQUIVALENCE (JDE(2),IY   ),(JDE(3),IT   ),(JDE(4),IZ   ),
C     > (JDE(5),IP   ),(JDE(1),ID)
 
      PARAMETER(MXREF=MIN(999,11*MXT))
      COMMON/OBJ5RF/ REF(MXJ,MXREF)

      DIMENSION FI(6,6)

      SAVE KOBJ2
      SAVE NBREF

      DATA NBREF / 1 /

      IMAX=11 * NBREF
      IDMAX=1
      IMAXT=IMAX/IDMAX
      P(2) = A(NOEL,20)
      P(3) = A(NOEL,21)
      P(4) = A(NOEL,22)
      P(5) = A(NOEL,23)
      P(6) = A(NOEL,24)
      P(1) = A(NOEL,25)
           
      IREF = 0

 1    CONTINUE
        IREF = IREF + 1
        IREF1 = IREF-1
        K = 30 + 10 * IREF1

        IF(IREF .LE. 7) THEN
C----------- For allowing possible use of the first 7 reference trajectories with FIT
          DO J = 2,MXJ1
            REF(J,IREF) = A(NOEL,K)
            K = K + 1
          ENDDO 
          REF(1,IREF) = A(NOEL,K)
        ENDIF

        I = 11 * IREF1 
        DO 53 J=2,5
          I=I+2
          DX = P(J)
          FO(J,I   ) = DX
          FO(J,I+1 ) = - DX
 53     CONTINUE
        FO(1,11*IREF-1) = P(1)
        FO(1,11*IREF) = - P(1)   

        IKAR = 1
        DO 51 I=11*IREF1+1,11*IREF
          IEX(I) = 1
          IREP(I) = I
          LET(I) = KAR(IKAR)
          IKAR = IKAR+1
          IF(IKAR .GT. 41) IKAR = 1
          DO 51 J = 1, 6
            FO(J,I) = FO(J,I) + REF(J,IREF)
            F(J,I) = FO(J,I)
 51     CONTINUE

      IF(IREF.LT.NBREF) GOTO 1

      IF(KOBJ2 .EQ. 1) THEN
C        if(nres.gt.0) write(nres,fmt='(a)')
C     >  'beam line initial alpha_y, beta_y, *_z, *_d are as follows :'
C In case of periodic structure, FI is filled up by tunes.f (uncoupled option) or tunesc.f (coupled option)
        FI(1,1) = A(NOEL,41)  
        IF(FI(1,1) .EQ. 0.D0) FI(1,1) = 1.D0    
        FI(2,1) = A(NOEL,40)      ! +alpha
        FI(1,2) = FI(2,1)
        FI(2,2) = (1.D0+FI(2,1)*FI(2,1))/FI(1,1)
        FI(3,3) = A(NOEL,43)  
        IF(FI(3,3) .EQ. 0.D0) FI(3,3) = 1.D0    
        FI(4,3) = A(NOEL,42)      ! +alpha
        FI(3,4) = FI(4,3)
        FI(4,4) = (1.D0+FI(4,3)*FI(4,3))/FI(3,3)
        FI(5,5) = A(NOEL,45)  
        IF(FI(5,5) .EQ. 0.D0) FI(5,5) = 1.D0    
        FI(6,5) = A(NOEL,44)      
        FI(5,6) = FI(6,5)
        FI(6,6) = (1.D0+FI(6,5)*FI(6,5))/FI(5,5)        
C Dy, Dy', Dz, Dz'
        FI(1,6) = A(NOEL,46)      
        FI(6,1) = FI(1,6)
        FI(2,6) = A(NOEL,47)      
        FI(6,2) = FI(2,6)
        FI(3,6) = A(NOEL,48)      
        FI(6,3) = FI(3,6)
        FI(4,6) = A(NOEL,49)      
        FI(6,4) = FI(4,6)
        SIGN = +1.D0
        CALL BEAMA2(FI,SIGN)
      ENDIF

      IF(NRES.GT.0) THEN
        WRITE(NRES,100) KOBJ,IMAX
  100   FORMAT(/,41X,'CALCUL  DES  TRAJECTOIRES',//,30X,'OBJET  (',I1,
     >  ')  FORME  DE ',I6,' POINTS ',//)
        WRITE(NRES,FMT='(/,T33,''Y (cm)'',T48,''T (mrd)'',T62,
     >  ''Z (cm)'',T76,''P (mrd)'',T90,''S (cm)'',T103,'' dp/p '')')
       WRITE(NRES,FMT='(14X,'' Sampling : '',T30, 
     >  5(4X,G10.2),4X,G12.4)') (P(J), J=2,6), P(1)
        IREF = 1
        WRITE(NRES,FMT='(2X,''Reference trajectry # '',I6,'' : '',T30,
     >                         5(4X,G10.2),4X,G12.4)') 
     >  IREF,(REF(J,IREF), J=2,6), REF(1,IREF)
        DO 20 IREF=2, NBREF
          WRITE(NRES,FMT='(
     >    20X,''# '',I6,'' : '',T30,5(4X,G10.2),4X,G12.4)') 
     >    IREF,(REF(J,IREF), J=2,6), REF(1,IREF)
 20     CONTINUE
      ENDIF

      RETURN

      ENTRY OBJ51(
     >            NBREFO)
      NBREFO = NBREF
      RETURN

      ENTRY OBJ52(KOBJ2I)
      KOBJ2 = KOBJ2I
      IF(KOBJ2 .EQ. 0) THEN
        NBREF = 1
      ELSEIF(KOBJ2 .EQ. 1) THEN
        NBREF = 1
      ELSEIF(KOBJ2 .GE. 2 .AND. KOBJ2 .LE.MXREF) THEN
        NBREF = KOBJ2
      ELSE
        CALL ENDJOB(
     >  'Pgm obj5, wrong value KOBJ2 in OBJET. Max is MXREF=',MXREF)
      ENDIF 
      RETURN
      END
