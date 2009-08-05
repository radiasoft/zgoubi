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
      SUBROUTINE PRDIP(NMAG,AT,RM,ACN,OP,XIE,OM,XIS,RE,TE,RS,TS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL IDLUNI
      character fnam*132, txt4*4
      save fnam
      data fnam /'zpop.out_PRDIP'/

      IF (IDLUNI(IUN)) THEN
        OPEN(UNIT=IUN,FILE=FNAM,ERR=997)
        CLOSE(unit=IUN,status='DELETE')
        OPEN(UNIT=IUN,FILE=FNAM,STATUS='NEW',ERR=997)
      ELSE
        WRITE(6,*) ' *** Problem in PRDIP : No idle unit number ! '
        GOTO 997
      ENDIF

      WRITE(IUN,FMT='(a2,i3,1p,11e14.6)')
     >  '% ',NMAG,AT,RM,ACN,OP,XIE,OM,XIS,RE,TE,RS,TS

      pi2 = 2.d0 * atan(1.d0)
      rd2deg = 90.d0/pi2
      deg2rd = 1.d0/rd2deg
      cm2m = 0.01d0
     
      rcm = rm * cm2m
 
c reference line (acn direction)
      x1seg = 0.d0
      y1seg = 0.d0
      ttseg = 0.d0
      sseg = 1.5d0 * rcm
      x2seg = sseg * cos(ttseg)
      y2seg = sseg * sin(ttseg)
      write(iun,*) 'LTYP 4'
      write(iun,fmt='(1p,2e14.6)') x1seg, y1seg
      write(iun,fmt='(1p,2e14.6)') x2seg, y2seg
      txt4 = 'REF.'
C      CALL TRTXT(x2seg,y2seg,TXT4,4,1)
c upstream limit
      x1seg = 0.d0
      y1seg = 0.d0
      ttseg = acn * deg2rd
      sseg = 1.5d0 * rcm
      x2seg = sseg * cos(ttseg)
      y2seg = sseg * sin(ttseg)
      write(iun,*) 'LTYP 1'
      write(iun,fmt='(1p,2e14.6)') x1seg, y1seg
      write(iun,fmt='(1p,2e14.6)') x2seg, y2seg
c acn arc
      xa = 0.d0
      ya = 0.d0
      tt1 = 0.
      dtt = acn * deg2rd
      call prArc(xa,ya,rcm/3.d0,tt1,dtt,iun) 
      txt4 = 'ACN'
C      CALL TRTXT(x2seg,y2seg,TXT4,4,1)
c downstream limit
      x1seg = 0.d0
      y1seg = 0.d0
      ttseg = (acn - at) * deg2rd
      sseg = 1.5d0 * rcm
      x2seg = sseg * cos(ttseg)
      y2seg = sseg * sin(ttseg)
      write(iun,*) 'LTYP 1'
      write(iun,fmt='(1p,2e14.6)') x1seg, y1seg
      write(iun,fmt='(1p,2e14.6)') x2seg, y2seg
c AT arc
      xa = 0.d0
      ya = 0.d0
      tt1 = ttseg
      dtt = at * deg2rd
      call prArc(xa,ya,rcm/2.d0,tt1,dtt,iun)
      txt4 = 'AT'
C      CALL TRTXT(x2seg,y2seg,TXT4,4,1)
c entrance sector limit
      x1seg = 0.d0
      y1seg = 0.d0
      ttseg = op * deg2rd
      sseg = 1.5d0 * rcm
      x2seg = sseg * cos(ttseg)
      y2seg = sseg * sin(ttseg)
      write(iun,*) 'LTYP 3'
      write(iun,fmt='(1p,2e14.6)') x1seg, y1seg
      write(iun,fmt='(1p,2e14.6)') x2seg, y2seg
      txt4 = 'EFBi'
c      CALL TRTXT(x2seg,y2seg,TXT4,4,1)
      xa = 0.d0
      ya = 0.d0
      tt1 = 0.d0
      dtt = op * deg2rd
      call prArc(xa,ya,rcm*1.3d0,tt1,dtt,iun)
c exit sector limit
      x1seg = 0.d0
      y1seg = 0.d0
      ttseg = om * deg2rd
      sseg = 1.5d0 * rcm
      x2seg = sseg * cos(ttseg)
      y2seg = sseg * sin(ttseg)
      write(iun,*) 'LTYP 3'
      write(iun,fmt='(1p,2e14.6)') x1seg, y1seg
      write(iun,fmt='(1p,2e14.6)') x2seg, y2seg
      txt4 = 'EFBo'
c      CALL TRTXT(x2seg,y2seg,TXT4,4,1)
      xa = 0.d0
      ya = 0.d0
      tt1 = 0.d0
      dtt = om * deg2rd
      call prArc(xa,ya,rcm*1.4d0,tt1,dtt,iun)
c entrance efb
      ttseg = (op+xie) * deg2rd
      sseg = 0.2d0 * rcm
      x2seg = rcm * cos(op*deg2rd) 
      y2seg = rcm * sin(op*deg2rd) 
      x1seg = x2seg + sseg * cos(ttseg)
      y1seg = y2seg + sseg * sin(ttseg)
      write(iun,*) 'LTYP  2'
      write(iun,fmt='(1p,2e14.6)') x1seg, y1seg
      write(iun,fmt='(1p,2e14.6)') x2seg, y2seg
      x1seg = x2seg - sseg * cos(ttseg)
      y1seg = y2seg - sseg * sin(ttseg)
      write(iun,fmt='(1p,2e14.6)') x1seg, y1seg
c exit efb
      ttseg = (om+xie) * deg2rd
      sseg = 0.2d0 * rcm
      x2seg = rcm * cos(om*deg2rd) 
      y2seg = rcm * sin(om*deg2rd) 
      x1seg = x2seg + sseg * cos(ttseg)
      y1seg = y2seg + sseg * sin(ttseg)
      write(iun,*) 'LTYP  2'
      write(iun,fmt='(1p,2e14.6)') x1seg, y1seg
      write(iun,fmt='(1p,2e14.6)') x2seg, y2seg
      x1seg = x2seg - sseg * cos(ttseg)
      y1seg = y2seg - sseg * sin(ttseg)
      write(iun,fmt='(1p,2e14.6)') x1seg, y1seg

      CLOSE(IUN)
      GOTO 99

 997  WRITE(6,*) ' Error upon OPEN statement '

 99   RETURN
      END 
