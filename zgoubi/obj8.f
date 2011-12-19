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
C  Brookhaven National Laboratory               és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE OBJ8(KREB31)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     -----------------------------------------
C      CONSTITUTION DE L'OBJET SUR UN ELLIPSOID
C     -----------------------------------------
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      CHARACTER LET
      COMMON/FAISCT/ LET(MXT)
      COMMON/FITEX/ DN0,C0,C1,C2,C3,DL
      CHARACTER  KAR(41)
      COMMON/KAR/ KAR
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IMAXD,IMAXT
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      COMMON/UNITS/ UNIT(MXJ)

      DIMENSION ALP(MXJ),BET(MXJ),EPS(MXJ)
      DIMENSION CENTRE(MXJ)
      PARAMETER(MXJ1=MXJ-1)
      CHARACTER*80 TXT
      PARAMETER(NDIM = 3)
      DIMENSION ISAM(NDIM)

      DIMENSION NVRNT(NDIM) 
      LOGICAL STRCON
      CHARACTER TXT80*80, STRA(2)*40


C----- IX, IY, IZ
      DO I= 1, NDIM
        ISAM(I) = INT(A(NOEL,20+I-1))
        NVRNT(I) = 1
        WRITE(TXT80,*) A(NOEL,20+I-1)
        IF(STRCON(TXT80,'.',
     >                      IS)) THEN
          IF(IS.LE.80-2) THEN
            ln = nint(a(noel,22+i))
            if(ln.gt.0) then
              READ(TXT80(IS+1:is+ln),*) ITMP
            endif
            IF(ITMP .GT. 0) NVRNT(I) = ITMP
          ENDIF
        ENDIF
      ENDDO

C----- CENTRE ELLIPSES
C       2-Y, 3-T, 4-Z, 5-Z, 6-X, 1-D
      CENTRE(2) = A(NOEL,30)
      CENTRE(3) = A(NOEL,31)
      CENTRE(4) = A(NOEL,32)
      CENTRE(5) = A(NOEL,33)
      CENTRE(6) = A(NOEL,34)
      CENTRE(1) = A(NOEL,35)
C----- PARAMETRE ELLIPSES
      ALP(2)=A(NOEL,40)
      ALP(4)=A(NOEL,50)
      ALP(6)=A(NOEL,60)
      BET(2)=A(NOEL,41)
      BET(4)=A(NOEL,51)
      BET(6)=A(NOEL,61)
      EPS(2)=A(NOEL,42)
      EPS(4)=A(NOEL,52)
      EPS(6)=A(NOEL,62)
      DO I= 1, NDIM
        IF(EPS(2*I) .EQ. 0) NVRNT(I) = 1
      ENDDO
 
      IMAX = 1
      DO  I= 1, 3
        IF(ISAM(I).NE.0) IMAX = IMAX * ISAM(I)*NVRNT(I)
      ENDDO
      IF(IMAX .GT. MXT) GOTO 98

      IF(KREB31 .EQ. 0) THEN
C------- Normal case
        IMI  = 1
        IMA = IMAX
        IF(IMAX .GT. MXT) GOTO 98
      ELSE
C------- Multiturn injection
        IF(IMAX*(KREB31+1) .GT. MXT) GOTO 98
        IMI  = 1 + IMAX*(IPASS-1)
        IMA = IMAX*IPASS
        IMAX=IMA
      ENDIF

      IF(NRES.GT.0) THEN
        WRITE(NRES,100) IMAX
100     FORMAT(/,15X,' Object  comprised  of ',I6,'  particles, '  
     >  ,'distributed  on  ellipses,  as  follows :',/)
 
        WRITE(NRES,123) (CENTRE(J),J=2,MXJ1),CENTRE(1)
 123    FORMAT(15X,' Ellipse centres (m-rad): '
     >  ,/,11X,' HORIZONTAL    ( Yo, To ):',T50,1P,2G12.4
     >  ,/,11X,' VERTICAL      ( Zo, Po ):',T50,   2G12.4
     >  ,/,11X,' LONGITUDINAL  ( Xo, Do ):',T50,   2G12.4,/)
 
        WRITE(NRES,109) (ALP(J),BET(J),EPS(J),J=2,MXJ1,2)
 109    FORMAT(15X,
     >  ' Alpha(rad), Beta(m/rad), E/pi(m.rad) :'
     >  ,/,11X,' HORIZONTAL    :',T35,1P,3G12.4
     >  ,/,11X,' VERTICAL      :',T35,   3G12.4
     >  ,/,11X,' LONGITUDINAL  :',T35,   3G12.4,/)

        WRITE(NRES,FMT='(15X,'' Sampling on an invariant, IX/IY/IZ : '',
     >             3(I5,''/''))') (ISAM(I),I=1,3)

        WRITE(NRES,FMT='(15X,'' Number of invariants in X/Y/Z : '',
     >             3(I5,''/''))') (NVRNT(I),I=1,3)
      ENDIF

C----- CONSTITUTION DU FAISCEAU
      DO 1 J=2,MXJ1,2
        J1=J+1
        J2=J/2
        NVRNT2 = NVRNT(J2)
        IF(J1.EQ.MXJ) J1=1 
        IF(EPS(J).EQ.0.D0) THEN
          DO 11 I=IMI,IMA
            FO(J ,I)=0.D0
            FO(J1,I)=0.D0
 11       CONTINUE
        ELSE
          IF(ISAM(J2) .GT. 0) THEN 
            REB=SQRT(EPS(J)*BET(J))
            DA = 2.D0*PI/DBLE(ISAM(J2))
          ELSE
            REB=0.D0
            DA = 0.D0
          ENDIF
          DREB = REB / DBLE(NVRNT2)
          JIM = IMAX / NVRNT2
          IMI1 = IMI -  JIM
          IMA1 = IMI1 + JIM - 1
          REBIV = 0.D0
          DO IVRNT = 1, NVRNT2
            IMI1 = IMI1 + JIM            
            IMA1 = IMA1 + JIM            
            REBIV = REBIV + DREB
            ANG = 0.D0
            DO 12 I=IMI1,IMA1
              X = REBIV*COS(ANG)
              IF(I.GT.MXT) 
     >          CALL ENDJOB('ERROR, SBR OBJ8 : max  MXT is ',MXT)
              FO(J ,I) = X/UNIT(J-1)
              FO(J1,I) = (REBIV*SIN(ANG)-ALP(J)*X)/BET(J)/UNIT(J)
              ANG = ANG + DA
 12         CONTINUE         
          ENDDO
        ENDIF
 1    CONTINUE      
 
      IKAR = 1
      DO 5 I=IMI,IMA
        FO(1,I)=FO(1,I) + CENTRE(1)/UNIT(6)
        F(1,I)=FO(1,I)
        DO 4 J=2,MXJ1
          FO(J,I)=FO(J,I) + CENTRE(J)/UNIT(J-1)
          F(J,I)=FO(J,I)
 4      CONTINUE
        IREP(I)=I
        IEX(I)=1
        LET(I)= KAR(IKAR)
        IKAR=IKAR+1
        IF(IKAR.GT.41) IKAR=1
 5    CONTINUE
 
      RETURN

 98   TXT = '    Too  many  particles'
      CALL OBJERR(ABS(NRES),2,MXT,TXT)
      CALL ENDJOB('OBJERR,  Too  many  particles ',-99)

      RETURN
      END
