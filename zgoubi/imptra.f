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
      SUBROUTINE IMPTRA(IMAX1,IMAX2,NRES)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      INCLUDE "C.DESIN.H"     ! COMMON/DESIN/ FDES(7,MXT),IFDES,KINFO,IRSAR,IRTET,IRPHI,NDES
C     >,AMS,AMP,AM3,TDVM,TETPHI(2,MXT)
C     1,AMS ,AMP,ENSTAR,BSTAR,TDVM,TETPHI(2,MXT)
      INCLUDE "MAXCOO.H"
      INCLUDE "C.OBJET.H"     ! COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      CHARACTER(1) LET
      INCLUDE "C.FAISCT.H"     ! COMMON/FAISCT/ LET(MXT)
C      LOGICAL ZSYM
      INCLUDE "C.TYPFLD.H"     ! COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      INCLUDE "C.PTICUL.H"     ! COMMON/PTICUL/ AM,Q,G,TO
      INCLUDE "C.REBELO.H"   ! COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      INCLUDE "C.UNITS.H"     ! COMMON/UNITS/ UNIT(MXJ) 

      DIMENSION SIG(4,4)
      CHARACTER(5) TXT(3)
      CHARACTER(10) UU(3)

      DATA TXT / '(Y,T)', '(Z,P)', '(t,K)' /
      DATA UU / '(cm,rd)', '(cm,rd)', '(mu_s,MeV)' /

      WRITE(NRES,100) IMAX2-IMAX1+1
 100  FORMAT('0',45X,'TRACE DU FAISCEAU',//,45X,I6,' TRAJECTOIRES',//
     >,35X,'OBJET',50X,'FAISCEAU',//,2(10X,'D',7X,'Y(CM)',4X,'T(MR)'
     >,4X,'Z(CM)',4X,'P(MR)',6X,'S(CM)'),/)
C     1,35X,'OBJET',50X,'FAISCEAU',//,2(10X,'D',6X,'Y(CM)',4X,'T(MR)'
C     2,4X,'Z(CM)',4X,'P(MR)',3X,'S(CM)'),/)
 
      DO 1 I=IMAX1,IMAX2
        WRITE(NRES,101) LET(I),IEX(I),(FO(J,I),J=1,6)
     >  ,(F(J,I),J=1,5),F(6,I),I
C     >  ,(F(J,I),J=1,6),I
C 101    FORMAT(' ',A1,1X,I2,1X,F8.4,5F10.3,6X,F8.4,4F9.3,1X,F10.3,1X,I6)
C 101    FORMAT(A1,1X,I3,1X,F8.4,5F10.3,5X,F8.4,4F9.3,1X,F12.3,1X,I5)
C Aug 2013
C 101    FORMAT(A1,1X,I2,1X,F8.4,5F10.3,5X,F8.4,4F9.3,1X,F12.4,1X,I6)
 101    FORMAT(A1,1X,I2,1X,F8.4,4F10.3,1X,F12.4,
     >                     2X,F8.4,4F9.3,1X,(1P,E14.6),1X,I5)
        IF(AM .NE. 0D0) THEN
          IF(IFDES.EQ.1) THEN
            WRITE(NRES,FMT='(15X,''Time of flight (mus) :'',
     >      1P,G16.8,'' mass (MeV/c2) :'',G14.6, ''    decay at (m) :'',
     >           G14.6)') F(7,I),AMQ(1,I),FDES(6,I)*UNIT(5)
          ELSE
            WRITE(NRES,FMT='(15X,''Time of flight (mus) :'',
     >      1P,G16.8,'' mass (MeV/c2) :'', G14.6)') F(7,I),AMQ(1,I)
          ENDIF
        ENDIF
 1    CONTINUE

Compute rms ellipse
      WRITE(NRES,FMT='(//,A)') '------' 
      WRITE(NRES,FMT='(''  Characteristics of concentration ellipse '',
     >''(Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls'', 
     >'' inside ellips, ratio, space, pass#) : '',/)')
      PI = 4.D0 * ATAN(1.D0)
      DO JJ = 1, 3
        CALL LPSFIT(JJ, 
     >                 EMIT,ALP,BET,XM,XPM)
Compute number of particles alive and numberinside ellipse
        CALL CNTINL(JJ,PI*EMIT,ALP,BET,XM,XPM,
     >                                        NLIV,NINL)
        RATIN = DBLE(NINL)/DBLE(IMAX)
        WRITE(NRES,110)
     >  PI*EMIT,ALP,BET,XM,XPM,NLIV,NINL,RATIN,TXT(JJ),IPASS
 110    FORMAT(1P,3(1X,E12.4),2(1X,E14.6),2(1X,I8),1X,G12.4,2X,A,2X,I8)
      ENDDO

      DO JJ = 1, 3
        CALL LPSFIT(JJ, 
     >                 EMIT,ALP,BET,XM,XPM)
        WRITE(NRES,fmt='(1P,/,A,2(/,5X,A,E14.6))') 
     >  TXT(JJ)//'  space (units : '//UU(JJ)//') :  ',
     >  ' sigma_'//TXT(JJ)(2:2)//' = sqrt(Surface/pi * BET) = ',
     >  sqrt(emit * BET) ,
     >  ' sigma_'//TXT(JJ)(4:4)//
     >                   ' = sqrt(Surface/pi * (1+ALP^2)/BET) = ',
     >  sqrt(emit * (1.d0+alp**2)/BET) 
      ENDDO

Compute 4-D sigma matrix
C      WRITE(NRES,FMT='(//,''  Beam  characteristics '', 1X
C     >,'' SIGMA(4,4) : '',/)')
      WRITE(NRES,FMT='(//,''  Beam  sigma  matrix : '',/)')
      CALL LPSFI4( 
     >             SQX,SQZ,SIG)
      WRITE(NRES,120) ((SIG(I,J),J=1,4),I=1,4)
 120  FORMAT(1P,4(1X,E14.6))
      WRITE(NRES,fmt='(/,5X,1P,A,2(2X,E14.6),3X,A)') 
     >' sqrt(det_Y), sqrt(det_Z) : ', SQX, SQZ,
     >' (Note :  sqrt(determinant) = ellipse surface / pi)'

      RETURN
      END
