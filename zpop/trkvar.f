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
C  Brookhaven National Laboratory                                               és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE TRKVAR(NOC,KVY,KDY,KVX,KDX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER  KVY*(*),KDY*(*),KVX*(*),KDX*(*)

      COMMON/REBELO/ NRBLT,IPASS,KWRI,NNDES,STDVM
      PARAMETER (MXJ=7)
      COMMON/UNITS/ UNIT(MXJ-1) 
      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB

      CHARACTER*130 TXT, TEMP
      CHARACTER*80 TITL
      SAVE TITL
      CHARACTER*14 TXTL
      CHARACTER*1 KLET

      CALL LOGO
      WRITE(TXT,FMT='(A54)') '* '//TITL(1:50)//' *'
      CALL TRTXT(60.D0,28.D0,TXT,0)

      CALL READC3(
     >            KL1,KL2)
      IF(KL1.EQ.-1) THEN
        WRITE(TXTL,FMT='(A5)') '* all'
      ELSE
        WRITE(TXTL,FMT='(I5,'' to '',I5)') KL1,KL2
      ENDIF
C          IF(LM.EQ.-1) THEN
C             WRITE(TXTL,FMT='(A5)') '* all'
C          ELSE
C             WRITE(TXTL,FMT='(I5)') LM
C          ENDIF

CCern, with HIGZ, 1994 - Positionning of TXT is in cm 
C      WRITE(TXT,101) TXTP,TXTL,NOC,IPASS
C      CALL TRTXT(.2D0,.1D0,TXT,50,0)
C      WRITE(TXT,100) KVY,KDY,KVX,KDX
C      CALL TRTXT(.2D0,.6D0,TXT,50,0)
C      WRITE(TXT,*) '   '
C      CALL TRTXT(.2D0,1.1D0,TXT,50,0)


C Hamel's graphic libraries,
C positionning of TXT is in screen units (x384/y256 Pix)
      WRITE(TXT,100) KVY,KDY,KVX,KDX
      CALL TRTXT(110.D0,248.D0,TXT,0)
      WRITE(TXT,107) XMI,XMA,YMI,YMA
C      CALL TRTXT(10.D0,18.D0,TXT,80,0)
      CALL TRTXT(10.D0,10.D0,TXT,0)
      CALL READC1(
     >            KP1,KP2,KP3)
      CALL READC5(
     >            KT1,KT2)
      CALL READC9(
     >            KEX,KLET)
      WRITE(TXT,101) KT1,KT2,KLET,TXTL,KP1,IPASS,KP3,NOC
      CALL TRTXT(10.D0,2.D0,TXT,0)

      CALL FBGTXT
C      WRITE(6,*)
CC----- Truc pour sortir de la fenetre graphique apres READCO, sur Alpha...
C      I=IDLG('('' Plot ended, press RETURN for more :'')','    ',1)

      CALL UNITR(KX,KY,
     >                 UX,UY)
      WRITE(6,*) 
      WRITE(6,100) KVY,KDY,KVX,KDX
 100  FORMAT(2A,'vs. ',2A) 
      WRITE(6,197) XMI,XMA,YMI,YMA,(XMI+XMA)/2.D0,(YMI+YMA)/2.D0,
     > KDX,KDY,(XMI+XMA)/2.D0/UX,(YMI+YMA)/2.D0/UY
 197  FORMAT('min.-max. ',T14,'Hor.:',1P,2G16.8,/,T14,' Ver.:',2G16.8,
     >   /,' (min+max)/2 Hor. & ver.:',2G16.8,'  ',A,'/',A,
     >   /,'                         ',2G16.8,'  (zgoubi units)')
 107  FORMAT('Mi-ma H/V: ',1P,2(1X,G10.3),'/',2(1X,G10.3))
      WRITE(6,*) ' ' 
      WRITE(6,101) KT1,KT2,KLET,TXTL,KP1,IPASS,KP3,NOC
 101  FORMAT('Part#',I3,'-',I6,' (',A1,'); Lmnt# ',A5,'; pass# ', 
     >    I6,'-',I6,', [',I4,'];',I7,' pnts') 

      RETURN

      ENTRY TRKVA2(TEMP)
      TITL = TEMP(1:50)
      RETURN
      END
