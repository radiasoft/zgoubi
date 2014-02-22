C  ZGOUBI, a program for computing the trajectories of charged particles
C  in electric and magnetic fields
C  Copyright (C) 1988-2007  Fran�ois M�ot
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
C  Fran�ois M�ot <fmeot@bnl.gov>
C  Brookhaven National Laboratory                    
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  -------
      SUBROUTINE IMPTRA(IMAX1,IMAX2,NRES)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      COMMON/DESIN/ FDES(7,MXT),IFDES,KINFO,IRSAR,IRTET,IRPHI,NDES
     >,AMS,AMP,AM3,TDVM,TETPHI(2,MXT)
C     1,AMS ,AMP,ENSTAR,BSTAR,TDVM,TETPHI(2,MXT)
      INCLUDE "MAXCOO.H"
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      CHARACTER LET
      COMMON/FAISCT/ LET(MXT)
      LOGICAL ZSYM
      COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      COMMON/PTICUL/ AM,Q,G,TO
      COMMON/REBELO/ NPASS,IPASS,KWRT,NNDES,STDVM
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      COMMON/UNITS/ UNIT(MXJ) 

      DIMENSION SIG(4,4)
      CHARACTER TXT*10, TXT2*2 

      WRITE(NRES,100) IMAX2-IMAX1+1
 100  FORMAT('0',45X,'TRACE DU FAISCEAU',//,45X,I6,' TRAJECTOIRES',//
     >,35X,'OBJET',50X,'FAISCEAU',//,2(10X,'D',7X,'Y(CM)',4X,'T(MR)'
     >,4X,'Z(CM)',4X,'P(MR)',6X,'S(CM)'),/)
C     1,35X,'OBJET',50X,'FAISCEAU',//,2(10X,'D',6X,'Y(CM)',4X,'T(MR)'
C     2,4X,'Z(CM)',4X,'P(MR)',3X,'S(CM)'),/)
 
      DO 1 I=IMAX1,IMAX2
        WRITE(NRES,101) LET(I),IEX(I),(FO(J,I),J=1,6)
     >  ,(F(J,I),J=1,6),I
C 101    FORMAT(' ',A1,1X,I2,1X,F8.4,5F10.3,6X,F8.4,4F9.3,1X,F10.3,1X,I6)
C 101    FORMAT(A1,1X,I3,1X,F8.4,5F10.3,5X,F8.4,4F9.3,1X,F12.3,1X,I5)
C Aug 2013
C 101    FORMAT(A1,1X,I2,1X,F8.4,5F10.3,5X,F8.4,4F9.3,1X,F12.4,1X,I6)
 101    FORMAT(A1,1X,I2,1X,F8.4,4F10.3,1X,F12.4,
     >                     2X,F8.4,4F9.3,1X,F14.4,1X,I5)
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
      WRITE(NRES,FMT='(//,''  Beam  characteristics '', 1X
     >,'' (EMIT,ALP,BET,XM,XPM,NLIV,NINL,RATIN) : '',/)')
      TXT = 'B-Dim '
      PI4 = 4.D0 *      4.D0 * ATAN(1.D0)
      DO 10 JJ = 1, 3
        CALL LPSFIT(JJ, 
     >                           EMIT,ALP,BET,XM,XPM)
Compute number of particles alive and numberinside ellipse
        CALL CNTINL(JJ,PI4*EMIT,ALP,BET,XM,XPM,
     >                                        NLIV,NINL)
        RATIN = DBLE(NINL)/DBLE(IMAX)
        WRITE(TXT2,FMT='(I1)') JJ
        WRITE(NRES,110)
     >       PI4*EMIT,ALP,BET,XM,XPM,NLIV,NINL,RATIN,TXT//TXT2,IPASS
 110    FORMAT(1P,5(1X,G12.4),2I8,1X,G12.4,A,I6)
 10   CONTINUE

Compute 4-D sigma matrix
      WRITE(NRES,FMT='(//,''  Beam  characteristics '', 1X
     >,'' SIGMA(4,4) : '',/)')
      CALL LPSFI4( 
     >             sqx,sqz,SIG)

      WRITE(NRES,fmt='(10X,1P,A,2E14.6)') ' Ex, Ez =', sqx,sqz
      WRITE(NRES,fmt='(10X,1P,A,2E14.6)') ' AlpX, BetX =', 
     >                  sig(1,2)/sqx, sig(1,1)/sqx
      WRITE(NRES,fmt='(10X,1P,A,2E14.6)') ' AlpZ, BetZ =', 
     >                  sig(3,4)/sqZ, sig(3,3)/sqx
      WRITE(NRES,120) ((SIG(I,J),J=1,4),I=1,4)
 120  FORMAT(/,1P,4E14.6)

      RETURN
      END
