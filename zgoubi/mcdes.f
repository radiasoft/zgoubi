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
      SUBROUTINE MCDES(
     >   DS,KEX,Y,T,Z,P,DA,QBR,SAR,TAR,IT,AMT,Q,BORO,XAR,KART)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST/ CL9,CEL,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE "MAXTRA.H"
      COMMON/DESIN/ FDES(7,MXT),IFDES,KINFO,IRSAR,IRTET,IRPHI,NDES
     >,AMS,AMP,AM3,TDVM,TETPHI(2,MXT)
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      CHARACTER LET
      COMMON/FAISCT/ LET(MXT)
      COMMON/UNITS/ UNIT(MXJ)

      CHARACTER LTD(MXT)
      SAVE LTD

      IF(KEX .LT. -1) GOTO 99

      DP=QBR/(Q*BORO)
      DL = DS * COS(T)
      IF( LET(IT) .NE. 'S' .AND. DL .GT. 0D0 ) THEN
         DSAR = FDES(6,IT)-SAR
         IF( DSAR .LE. 0D0 ) THEN
C---------- Primary particle has reached or passed  its TOF limit. Therefore: 
C             * If field free section (CHANGREF, ESL... - B=0) then go back 
C             linearly to decay point (XAR=1.).
C             Test on time of flight of primary particles is performed at end of 
C             sextion, then back to decay point (by DSAR<0) for calculation of 
C             transport of secondary particle. 
C             * If B non zero, test on time of flight is performed after each 
C             integration step, and no corrective interpolation is 
C             performed (XAR=0) : hence an error on position of decay point < step. 
C             Save  coordinates of decay point.

           DSAR=DSAR*XAR
           FDES(6,IT) = SAR + DSAR
C           IF(IEX(IT).EQ.-4) THEN
CC Primary particle was stopped by COLLIMA or CHAMBR
C             IF(FDES(6,IT).LT.SORT(IT)) IEX(IT)=1
C           ENDIF
           Y = Y + DSAR*SIN(T)*COS(P)
           Z = Z + DSAR*SIN(P)
           FDES(3,IT) = T / UNIT(2)
           FDES(5,IT) = P / UNIT(4)
           FDES(1,IT) = DP
           LTD(IT) = LET(IT)

C------------------------------------------------------------------------
C DO NOT SCRATCH !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C---------- Genererate angles & parametres of new  particle
C           Angles wrt incoming direction,  COM frame :
C            TET = ACOS(1.D0 - 2.D0*RNDM(IRTET))
C            PHI = 2.D0*PI* RNDM(IRPHI)
C--------------------------------------------------------------------------
C--- It is necessary to catch TET, PHI from fixed sets if is FIT to be used
           TET = TETPHI(1,IT)
           PHI = TETPHI(2,IT)
C--------------------------------------------------------------------------
 
C---------- Ref. lab :
           BETAP = 1.D0/SQRT( 1.D0 +  (AMP/(CL9*QBR))**2  )
           GAMAP = 1.D0/SQRT(1.D0-BETAP*BETAP)

           AMP2 = AMP*AMP
           AMS2 = AMS*AMS
           AM32 = AM3*AM3
           ENSTAR = .5D0*(AMP2+AMS2-AM32)/AMP
           BSTAR  = SQRT(1.D0-AMS2/ENSTAR/ENSTAR)

           DP00=DP
C ERROR, corrctn FM Feb. 2002
C           DP=ENSTAR*GAMAP*(1.D0+BSTAR*BETAP*COS(TET))/(BORO*CL9)
           DP=SQRT((ENSTAR*GAMAP*(1.D0+BSTAR*BETAP*COS(TET)))**2 - 
     >       AMS2)       /(BORO*CL9)

CC---------- Test densities ("ligne verte")//
C           DP = .4d0
C           DP0=DP
C           costet = ( SQRT((DP*BORO*CL9)**2 + AMS2)/(ENSTAR*GAMAP) 
C     >             -1.D0) / (BSTAR*BETAP)
C           TET =  ACOS(costet)
CC-----------------------------------------//

C          *** TETA LAB PAR RAPPORT A L'INCIDENCE :
           TET = ATAN( SIN(TET) / GAMAP / (COS(TET) + BETAP/BSTAR) )

           IF(KINFO .EQ. 1) THEN
               WRITE(ABS(NRES),102) IT,T*1.D3, P*1.D3, DP00 ,BETAP
     >           , GAMAP, IT, FDES(6,IT), SAR
102            FORMAT(1P
     >         ,/,10X,'Decay  of  parent  particle  M1,  # ',I6
     >         ,',      with  following  parameters : '
     >         ,/,15X,'- current  direction     T = ',G12.4
     >         ,' mrd ,     P = ',G12.4,5X,' mrd '
     >         ,/,15X,'- relative  momentum  1.+dp/p = '
     >         ,G12.4,',     beta = ',G12.4,',    gamma = ',G12.4
     >         ,/,15X,'- native  life  flight  distance  of  #',I6
     >         ,'  is ',G12.4
     >         ,' cm,     current  path  length  is  ',G12.4,' cm')
 
               WRITE(ABS(NRES),100) TET*1.D3, PHI*1.D3, DP
100            FORMAT(1P
     >         ,/,10X,'Native  parameters  of  daugther  particle  M2 :'
     >         ,/,15X,'Lab  angles :  theta = ',G12.4
     >         ,' mrad ,     phi = ',G12.4,' mrad ,'
     >         ,/,15X,'Relative momentum  p/p_ref = ',F11.5)
           ENDIF
 
           CTET = COS(TET)
           STET = SIN(TET)
           CPHI = COS(PHI)
           SPHI = SIN(PHI)
           CT = COS(T)
           ST = SIN(T)
           CP = COS(P)
           SP = SIN(P)
 
           XP = CTET*CP*CT - STET*(CPHI*ST + SPHI*SP*CT)
           YP = STET*CPHI*CT + (CTET*CP - STET*SPHI*SP)*ST
           RP = SQRT(XP*XP + YP*YP)
 
C---------- ANGLES LAB APRES DESINTEGRATION
           T = ATAN2(YP,XP)
           P = ATAN2((CTET*SP + STET*SPHI*CP),RP)
 
C---------- COTES EN FIN DE DL
           Y = Y - DSAR*SIN(T)*COS(P)
           Z = Z - DSAR*SIN(P)
           SAR= FDES(6,IT)-DSAR
 
C----------- Calculate M1 TOF upstream of decay point
           QBRO = QBR*CL9
           DTAR = DSAR / (QBRO/SQRT(QBRO*QBRO+AMT*AMT)*CL9) *1.D-5
           TAR= TAR + DTAR
C---------- Assigns  mass 
           AMT = AMS
C----------- Calculate M2 TOF downstream of decay point
           IF(AMT .NE. 0.D0) THEN 
             QBRO = BORO*DP*CL9*Q
             DTAR = DSAR / (QBRO/SQRT(QBRO*QBRO+AMT*AMT)*CL9) *1.D-5
             TAR= TAR - DTAR
           ENDIF

C---------- Assign  lettre 'S'(-econdary) to decay particle
           LET(IT)='S'
C---------- FDES(7,IT) contains com life time (s) of S particle
           TOS = FDES(7,IT)
           IF(TOS.NE.0.D0) THEN
             U=RNDM()
             IF(U .EQ. 0.D0) U=1.D0
             FDES(6,IT) = FDES(6,IT)
     >                 -1.D-7*CEL*CEL * BORO*DP * TOS / AMS * LOG(U)
           ENDIF 
C---------- Decay counter
           NDES=NDES+1

           IF(KINFO .EQ. 1) THEN
             WRITE(ABS(NRES),FMT='(1P,/,10X
     >       ,''Parameters of daugther particle M2 in Zgoubi frame :''
     >       ,/,5X,''#, nDes, Y, T, Z, P, SAR, decay absissa :''
     >       ,2I6,2X,6G12.4,/)')
     >       IT, NDES, Y, T, Z, P, SAR, FDES(6,IT)
           ENDIF
 
C-----------------------
 
         ENDIF
 
      ELSEIF( LET(IT) .EQ. 'S') THEN
 
C------- Change  reference  at each  step
        IF(KART .EQ. 2) FDES(3,IT) = FDES(3,IT) - DA*1.D3
 
        IF( DL .LT. 0.D0) THEN
C--------- The particle has undergone linear backward motion (e.g., in CHANGREF, negative ESL)

C          DSAR = FDES(6,IT) - SAR 
C          IF( DSAR .GT. 0. ) THEN
CC----------- The parent particle decay point does stand somewhere along this backward step
C 
C            NDES = NDES - 1
CC----------- Go forth to decay point
C            CP = COS(P)
C            YD = Y + DSAR*SIN(T)*CP
C            ZD = Z + DSAR*SIN(P)
C            DX = DSAR*COS(T)*CP
CC----------- Get old  angles, old relative rigidity = those of parent particle
C            T = FDES(3,IT)*.001
C            P = FDES(5,IT)*.001
C            DP= FDES(1,IT)
CC-----------Compute new value of DSAR
C            CP = COS(P)
C            DSAR = DX /COS(T)/CP
CC----------- return to actual position with parent particle angles
C           Y = YD - DSAR*SIN(T)*CP
C            Z = ZD - DSAR*SIN(P)
C            SAR = FDES(6,IT) - DSAR
CC----------- Deduces TOF of M2
C            QBRO = BR*CL9*QT
C            DTAR = DSAR / (QBRO/SQRT(QBRO*QBRO+AMT*AMT)*CL9) *1.D-5
C            TAR= TAR-DTAR
CC----------- Resets  mass to parent, rigidity 
C            AMT = AMP
C            BR = BORO*DP
CCCCCCCCCC il manque la correction de dsar, en utilisant dx ou da
CC------------ Adds M1 TOF
C            QBRO = BR*CL9*QT   
C            DTAR = DSAR / (QBRO/SQRT(QBRO*QBRO+AMT*AMT)*CL9) *1.D-5
C            TAR= TAR + DTAR
C
CC----------- REAFFECTE SA LETTRE INITIALE A TOUTE PARTICULE REINTEGREE
C            LET(IT)=LTD(IT)
C 
C          ENDIF

        ELSE
C--------- FDES(7,IT) contains com life time (s) of S(econdary) particle
          IF(FDES(7,IT).GE.0.D0) THEN
            DSAR = FDES(6,IT)-SAR
            IF( DSAR .LE. 0D0 ) CALL KSTOP(10,IT,KEX,*99)
          ENDIF

        ENDIF
 
      ENDIF
      QBR=Q*BORO*DP
      BRI = QT/QBR
C-------- Test NuFact
CCCCCCCCCc         IF(BR.LT.666) DP=-999
 99   RETURN
      END
