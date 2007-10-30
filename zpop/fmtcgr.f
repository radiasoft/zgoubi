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
      SUBROUTINE FMTCGR(FNCT,NLOG,NV,V,NC,CX,CY,XFB,XLAMB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(*), CX(*), CY(*)

C----- Plot the data and the matched field fall-off curve.

      COMMON/CDF/ IES,IORDRE,LCHA,LIST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN    
      LOGICAL OKECH, OKVAR, OKBIN
      COMMON/ECHL/ OKECH, OKVAR, OKBIN
      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB

      PARAMETER(MXC=400)
      DIMENSION FCT(MXC)
      CHARACTER*80 TXT
C      CHARACTER TX1*80, TFM*80
      EXTERNAL FNCT

      KX = 7
      KY = 32
      OKVAR = .TRUE.
      YN = 1.D0

 11   CONTINUE

        IF(XMI .LT. XMA .AND. YMI .LT. YMA) THEN
          OKECH = .TRUE.
          CALL TRAXES(XMI,XMA,YMI,YMA*1.05,-1) 
        ELSE
          OKECH = .FALSE.
          WRITE(*,*) ' Scale problems ; no plot possible !!'
          RETURN
        ENDIF          

 12   CONTINUE
      CALL LINTYP(1)

      CALL FUNK(FNCT,V,CX(MXC-1),CX(MXC),NC,NV,CX,FCT)
      I = 1
      X = CX(I)
      Y = YMA * FCT(I)
      CALL VECTPL(X,Y,4)
      IF(LIS .EQ. 2) CALL IMPV(NLOG,I,X,Y,DUM,DUM,IDUM,IDUM)
      DO 1 I = 2, NC
        X = CX(I)
        Y = YMA * FCT(I)
        CALL VECTPL(X,Y,2)
        IF(LIS .EQ. 2) CALL IMPV(NLOG,I,X,Y,DUM,DUM,IDUM,IDUM)
 1    CONTINUE

      CALL LINTYP(12)
      I = 1
      X = CX(I)
      Y = CY(I) * YMA
      CALL VECTPL(X,Y,4)
      IF(LIS .EQ.  2) CALL IMPV(NLOG,I,X,Y,DUM,DUM,IDUM,IDUM)
      DO 2 I = 2, NC
        X = CX(I)
        Y = CY(I) * YMA
        CALL VECTPL(X,Y,2)
        IF(LIS .EQ. 2) CALL IMPV(NLOG,I,X,Y,DUM,DUM,IDUM,IDUM)
 2    CONTINUE
 
C------ Magnetic face
        CALL LINTYP(4)
        CALL VECTPL(XFB,YMI,4)
        CALL VECTPL(XFB,YMA,2)

C Hamel's graphic libraries,
C positionning of TXT is in screen units (x384/y256 Pix)
      CALL LOGO

      WRITE(TXT,FMT=
     >  '(''Field fall-off  v.s.  distance (m)  (data & model)'')')
      CALL TRTXT(120.D0,245.D0,TXT,80,0)
      WRITE(TXT,FMT='(''XFB, Lambda ='',2F11.6,'' m'')') XFB, XLAMB
      CALL TRTXT(10.D0,18.D0,TXT,80,0)
C      WRITE(TXT,FMT='(''C0-C'',I1,'' ='',6F9.5)') NV-1,(V(I), I=1,NV)
      WRITE(TXT,FMT='(''C0-C5 = '',6F9.5)') (V(I), I=1,6)
      CALL TRTXT(10.D0,10.D0,TXT,80,0)
      WRITE(TXT,FMT='(''I1.Gap = FINT.Gap ='',F11.6,'' m'')')
     >  GK1(FNCT,CX,NC,NV,V)
      CALL TRTXT(10.D0,2.D0,TXT,80,0)

      CALL FBGTXT
      WRITE(*,FMT=
     >  '('' Magnetic field  v.s.  distance (m)  (data & model)'')')
      WRITE(*,FMT='('' XFB, Lambda ='',6F11.6,'' m'')') XFB, XLAMB
C      WRITE(*,FMT='(''C0-C'',I1,'' ='',6F9.5)') NV-1,(V(I), I=1,NV)
      WRITE(*,FMT='(''C0-C5 = '',6F9.5)') (V(I), I=1,6)
      WRITE(*,FMT='('' I1.Gap = FINT.Gap ='',F11.6,'' m'')')
     >  GK1(FNCT,CX,NC,NV,V)

      RETURN
      END
