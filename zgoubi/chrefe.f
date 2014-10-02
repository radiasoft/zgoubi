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
C  Upton, NY, 11973
C  -------
      SUBROUTINE CHREFE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     -------------------------------------------------
C     CHANGEMENT DE REFERENCE DE L'ENSEMBLE DU FAISCEAU
C     -------------------------------------------------
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      COMMON/CHAMBR/ LIMIT,IFORM,YL2,ZL2,SORT(MXT),FMAG,BMAX
     > ,YCH,ZCH
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXSTEP.H'
      INCLUDE 'CSR.H'
      COMMON/DESIN/ FDES(7,MXT),IFDES,KINFO,IRSAR,IRTET,IRPHI,NDES
     >,AMS,AMP,AM3,TDVM,TETPHI(2,MXT)
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER(80) TA
      PARAMETER (MXTA=45)
      COMMON/DONT/ TA(MXL,MXTA)
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      COMMON/GASC/ AI, DEN, KGA
      LOGICAL ZSYM
      COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)
      COMMON/SYNRA/ KSYN
      INCLUDE 'MXFS.H'
      COMMON/SCALP/ VPA(MXF,MXP),JPA(MXF,MXP)
      INCLUDE 'MXSCL.H'
      DIMENSION KFM(MXSCL)

      LOGICAL EVNT
      PARAMETER (MSR=10)
      CHARACTER(2) QSHRO(MSR)
      DIMENSION VSHRO(MSR)
      PARAMETER (I1=1)
 
      DUM = SCALER(I1,NOEL,
     >                     DUM)
      CALL SCALE9(
     >            KFM)
      DO IFM = 1, MXSCL
         IF(KFM(IFM) .LE. 0) GOTO 20
         DO I= 1 , JPA(KFM(IFM),MXP)
            A(NOEL,JPA(KFM(IFM),I)) = VPA(KFM(IFM) ,I)
         ENDDO
      ENDDO
 20   CONTINUE

      NSR = NINT(A(NOEL,9))
      IF(NSR .GE. MSR)
     >CALL ENDJOB('SBR CHREFE. Max. nmbr of transforms must be < ',MSR)
      VSHRO(MSR) = NSR
      DO I = 1, NSR
        QSHRO(I) = TA(NOEL,I)(1:2)
      ENDDO
C To allow for old style
      QSHRO(4) = TA(NOEL,4)(1:2)

      EVNT = KSPN.EQ.1 .OR. IFDES.EQ.1 .OR. KGA.EQ.1 .OR. 
     >  LIMIT.EQ.1 .OR. KSYN.GE.1 .OR. KCSR.EQ.1 

      IF( QSHRO(4) .EQ. 'OL') THEN
C Old style. TA(NOEL,4)(1:2) was set to 4 in rchang.f
        XC  = A(NOEL,1)
        YC  = A(NOEL,2)
        AA  = A(NOEL,3)*RAD
        VSHRO(1) = XC
        VSHRO(2) = YC
        VSHRO(3) = AA

        IF(NRES.GT.0) THEN
          WRITE(NRES,100) XC,YC,AA*DEG,AA
 100      FORMAT(/,' CHANGE  OF  REFERENCE  FRAME,   XC ='
     >    ,F10.3,' cm , YC ='
     >    ,F10.3,'  cm ,   A =',F12.5,' deg  (',F10.6,' rad)',/)
        ENDIF

        DO IT=1,IMAX
C--------- IEX<-1<=> PARTICULE STOPPEE
          IF( IEX(IT) .GE. -1) THEN
 
            IF(IT .EQ. IREP(IT) .OR. .NOT.ZSYM) THEN
              CALL INITRA(IT)
              CALL CHAREF(EVNT,XC,YC,AA)
              CALL MAJTRA(IT)
            ELSE
              CALL DEJACA(IT)
            ENDIF
          ENDIF
        ENDDO

      ELSE      
C New style

        DO I=1, NSR
          VSHRO(I) = A(NOEL,I)
          IF(QSHRO(I)(2:2).EQ.'R') VSHRO(I) = VSHRO(I)*RAD
        ENDDO

        IF(NRES.GT.0) THEN
          WRITE(NRES,FMT='(/,5X,''Change  of  reference, '',
     >    I2,''  transformations :'',/)') NSR
          DO I=1, NSR
            IF(QSHRO(I)(2:2).EQ.'S') 
     >        WRITE(NRES,110) QSHRO(I)(1:1),VSHRO(I)
 110          FORMAT(10X,
     >         'type : ',A1,'-shift,        value : ',1P,E14.6,' cm')
            IF(QSHRO(I)(2:2).EQ.'R')
     >        WRITE(NRES,111) QSHRO(I)(1:1),VSHRO(I)/RAD,VSHRO(I)
 111        FORMAT(10X,
     >         'type : ',A1,'-rotation,     value : ',1P,E14.6,' deg  ('
     >         ,E14.6,' rad)')
          ENDDO
        ENDIF

        DO IT=1,IMAX
C--------- IEX<-1<=> PARTICULE STOPPEE
          IF( IEX(IT) .GE. -1) THEN
 
            IF(IT .EQ. IREP(IT) .OR. .NOT.ZSYM) THEN
              CALL INITRA(IT)
C              CALL CHANRF(NSR,EVNT,QSHRO,VSHRO)
              CALL CHANRF(EVNT,QSHRO,VSHRO)
              CALL MAJTRA(IT)
            ELSE
              CALL DEJACA(IT)
            ENDIF
          ENDIF
        ENDDO

      ENDIF

      IF(NRES.GT.0) WRITE(NRES,101) IEX(1),(F(J,1),J=1,7)
 101  FORMAT(/,' Traj #1,  IEX,D,Y,T,Z,P,S,time :',I3,1P,5G12.4,2G17.5)

      RETURN
      END
