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
C  Brookhaven National Laboratory                                                               és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
C      SUBROUTINE HMTCGR(FNCT,NLOG,NV,V,NC,CX,CY,XM,SIG,MXV,MXC)
      SUBROUTINE HMTCGR(FNCT,NLOG,NV,V,NC,CX,CY,XM,SIG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(*), CX(*), CY(*)

C----- Plot the data and the matching function

      COMMON/CDF/ IES,IORDRE,LCHA,LIST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN    
      LOGICAL OKECH, OKVAR, OKBIN
      COMMON/ECHL/ OKECH, OKVAR, OKBIN
      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB

      PARAMETER(MXC=400)
      DIMENSION FCT(MXC)
      CHARACTER*80 TXT
C      CHARACTER TX1*80, TFM*80
      EXTERNAL FNCT

C      CALL HOMCLR

      KX = 2
      KY = 28
      OKVAR = .TRUE.

        IF(XMI .LT. XMA .AND. YMI .LT. YMA) THEN
          OKECH = .TRUE.
          CALL TRAXES(XMI,XMA,YMI,YMA*1.05,-1) 
        ELSE
          OKECH = .FALSE.
          WRITE(*,*) ' Scale problems ; no plot possible !!'
          GOTO 70
        ENDIF          

      CALL LINTYP(1)

      CALL FUNK(FNCT,V,CX(MXC-1),CX(MXC),NC,NV,CX,FCT)
      I = 1
      X = CX(I)
      Y = FCT(I)
      CALL VECTPL(X,Y,4)
      IF(LIS .EQ. 2) CALL IMPV(NLOG,I,X,Y,DUM,DUM,IDUM)
      DO 1 I = 2, NC
        X = CX(I)
        Y = FCT(I)
        CALL VECTPL(X,Y,2)
        IF(LIS .EQ. 2) CALL IMPV(NLOG,I,X,Y,DUM,DUM,IDUM)
 1    CONTINUE

      CALL LINTYP(12)
      I = 1
      X = CX(I)
      Y = CY(I)
      CALL VECTPL(X,Y,4)
      IF(LIS .EQ.  2) CALL IMPV(NLOG,I,X,Y,DUM,DUM,IDUM)
      DO 2 I = 2, NC
        X = CX(I)
        Y = CY(I)
        CALL VECTPL(X,Y,2)
        IF(LIS .EQ. 2) CALL IMPV(NLOG,I,X,Y,DUM,DUM,IDUM)
 2    CONTINUE
 
C------ Magnetic face
        CALL LINTYP(4)
        CALL VECTPL(XM,YMI,4)
        CALL VECTPL(XM,YMA,2)


C Hamel's graphic libraries,
C positionning of TXT is in screen units (x384/y256 Pix)
      CALL LOGO

      WRITE(TXT,FMT=
     >  '(''Profile  v.s.  distance (data & model)'')')
      CALL TRTXT(120.D0,245.D0,TXT,0)
      WRITE(TXT,FMT='(''XM, Sigma ='',2F11.6,'' cm'')') XM, SIG
      CALL TRTXT(10.D0,18.D0,TXT,0)

      WRITE(TXT,FMT='(''C0-C'',I1,'' ='',11F8.4)') NV-1,(V(I), I=1,NV)
      CALL TRTXT(10.D0,10.D0,TXT,0)

 70   CALL FBGTXT
      WRITE(*,FMT=
     >  '('' Profile (data & model)'')')
      WRITE(*,FMT='('' XM, Sigma ='',6F11.6,'' cm'')') XM, SIG
      WRITE(*,FMT='(''C0-C'',I1,'' ='',11F8.4)') NV-1,(V(I), I=1,NV)

      RETURN
      END
