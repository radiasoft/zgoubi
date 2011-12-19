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
      SUBROUTINE ELCYLD(SCAL,
     >                       ORBL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON/AIM/ AE,AT,AS,RM,XI,XF,EN,EB1,EB2,EG1,EG2
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      PARAMETER(MCOEF=6)
      COMMON/CHAFUI/ XE,XS,CE(MCOEF),CS(MCOEF),QCE(MCOEF),QCS(MCOEF)
C      COMMON/CHAFUI/ XE,XS,CE(6),CS(6),QCE(6),QCS(6)
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QEL,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      PARAMETER (MDR=9)
      COMMON/DROITE/ CA(MDR),SA(MDR),CM(MDR),IDRT
      COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      PARAMETER(MPOL=10)
      COMMON/MULTPE/ EM(MPOL),QLE(MPOL),QLS(MPOL)
     >,QE(MPOL,MCOEF),QS(MPOL,MCOEF),RTQ(MPOL)
      LOGICAL ZSYM
      COMMON/OPTION/ KFLD,MG,LC,ML,ZSYM
      COMMON/PTICUL/ AM,Q,G,TO
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
 
      EQUIVALENCE (EB1,XLE), (EB2,XLS), (EG1,V0)
 
      DEV =A(NOEL,10)
      RM =A(NOEL,11)*1.D2
      EM(1)  = A(NOEL,12)*SCAL

      ORBL = RM * DEV
      IF (EM(1) .EQ. 0.D0) KFLD=0

      P0 = BORO*CL9 *Q
      EK = SQRT(P0*P0+AM*AM) - AM
      ETA = EK / 2.D0 / AM
      V0 = EM(1) * RM * (1D0+2.D0*ETA) / (1+ETA)

      XLE = A(NOEL,20)
      QLE(1) = A(NOEL,21)
      IF(QLE(1) .EQ. 0D0) XLE = 0.D0
      XLS = A(NOEL,40)
      QLS(1) = A(NOEL,41)
      IF(QLS(1) .EQ. 0D0) XLS = 0.D0
      DO 8 I=1,MCOEF
        QCE(I) = A(NOEL,30+I-1)
        QCS(I) = A(NOEL,50+I-1)
 8    CONTINUE
 
      AE = ATAN(XLE/RM)
      AS = ATAN(XLS/RM)
      AT = AE + DEV + AS

      XI = 0.D0
      XF = AT
      
      IF(NRES.GT.0) THEN
        WRITE(NRES,100) ' Electrical Cyl Deflector',
     >    DEV, RM,EM(1),P0/SQRT(P0*P0+AM*AM), 
     >    EM(1)*AM*RM/1.D2/P0/P0*1.D-6, ETA
 100    FORMAT(/,5X,' +++++  ', A30
     >        ,//,15X,' Deviation                 = ',F10.5 ,' rad'
     >        ,/ ,15X,' Layout  radius  R         = ',F10.5 ,' cm'
     >        ,/ ,15X,' E-field                   = ',G12.4,' V/m'
     >        ,//,15X,' v/c                       = ',G12.4
     >        ,/ ,15X,' Balance (E.Mass.R/p^2, normally 1)  = ',F12.6
     >        ,/ ,15X,' eta = Ekin / 2Mass                  = ',G12.4)
 
        WRITE(NRES,104) 'Entrance'
 104    FORMAT(/,15X,A8,'  face')
        WRITE(NRES,101) XLE,QLE(1)
 101    FORMAT(15X,' DX = ',F10.3,'    LAMBDA = ',F10.3)
        IF(QLE(1).NE.0D0) WRITE(NRES,132) (QCE(I),I=1,MCOEF)
 132    FORMAT(15X,' COEFFICIENTS :',6F9.5)
 
        WRITE(NRES,104) 'Exit    '
        WRITE(NRES,101) XLS,QLS(1)
        IF(QLS(1).NE.0D0) WRITE(NRES,132) (QCS(I),I=1,MCOEF)
 
      ENDIF
 
C----- E:MeV/cm
      EM(1)=EM(1)*1.D-8
      V0=V0*1.D-8

      IF(QLE(1) .NE. 0D0) THEN
C        QE(1,1)=-EM(1)/QLE(1)
        QE(1,1)=-V0/QLE(1)
        DO 10 I=2,4
          QE(1,I)=-QE(1,I-1)/QLE(1)
 10     CONTINUE
      ENDIF

      IF(QLS(1) .NE. 0D0) THEN
C        QS(1,1)=-EM(1)/QLS(1)
        QS(1,1)=-V0/QLS(1)
        DO 11 I=2,4
          QS(1,I)=-QS(1,I-1)/QLS(1)
 11     CONTINUE
      ENDIF
 
      RETURN
      END
