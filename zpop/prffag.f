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
      SUBROUTINE PRFFAG(NMAG,AT,RM,ACN,OP,XIE,OM,XIS,RE,TE,RS,TS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL IDLUNI
      character fnam*132, txt4*4
      logical first
      save fnam, iun, first, noc
      data fnam /'zpop.out_PRFFAG'/
      data first /.true./

      if(first) then
        noc = 1
        IF (IDLUNI(IUN)) THEN
          OPEN(UNIT=IUN,FILE=FNAM,ERR=997)
          CLOSE(unit=IUN,status='DELETE')
          OPEN(UNIT=IUN,FILE=FNAM,STATUS='NEW',ERR=997)
        ELSE
          WRITE(6,*) ' *** Problem in PRFFAG : No idle unit number ! '
          GOTO 997
        ENDIF
        first = .false.
      else
        noc = noc + 1
      ENDIF

      WRITE(IUN,FMT='(a2,2i3,1p,11e14.6)')
     >  '% ',noc,NMAG,AT,RM,ACN,OP,XIE,OM,XIS,RE,TE,RS,TS

      pi = 4.d0 * atan(1.d0)
      pi2 = pi/2.d0
      rd2deg = 90.d0/pi2
      deg2rd = 1.d0/rd2deg
      cm2m = 0.01d0

      rme = rm * cm2m
          DEV =  at*deg2rd + te - ts
      ttrf = -(2.d0 * ACN -AT/2.d0) * deg2rd 
     >      + dev * DBLE(noc-1)

      sign = 1.d0
c reference line (acn direction)
      x1seg = 0.d0
      y1seg = 0.d0
      ttseg = (2.d0 * ACN -AT/2.d0) * deg2rd + ttrf
      sseg = 1.5d0 * rme
      x2seg = sseg * cos(ttseg)
      y2seg = sseg * sin(ttseg)
      write(iun,*) 'LTYP 4    reference line (ACN dir.)'
      write(iun,fmt='(1p,2e14.6)') x1seg, y1seg
      write(iun,fmt='(1p,4e14.6)') x2seg, y2seg, sseg, ttseg
      txt4 = 'REF.'
c      CALL TRTXT(x2seg,y2seg,TXT4,4,1)
c upstream limit
      x1seg = 0.d0
      y1seg = 0.d0
      ttseg = sign * acn * deg2rd  + ttrf
      sseg = 1.5d0 * rme
      x2seg = sseg * cos(ttseg)
      y2seg = sseg * sin(ttseg)
      write(iun,*) 'LTYP 1       upstr. limit'
      write(iun,fmt='(1p,2e14.6)') x1seg, y1seg
      write(iun,fmt='(1p,4e14.6)') x2seg, y2seg, sseg, ttseg
c acn arc
      xa = 0.d0
      ya = 0.d0
      tt1 = 0.   + ttrf
      dtt = acn * deg2rd
      call prArc(xa,ya,rme/3.d0,tt1,dtt,iun) 
      txt4 = 'ACN'
c      CALL TRTXT(x2seg,y2seg,TXT4,4,1)
c downstream limit
      x1seg = 0.d0
      y1seg = 0.d0
      ttseg = sign * (acn - at) * deg2rd  + ttrf
      sseg = 1.5d0 * rme
      x2seg = sseg * cos(ttseg)
      y2seg = sseg * sin(ttseg)
      write(iun,*) 'LTYP 1       downstr. limit'
      write(iun,fmt='(1p,2e14.6)') x1seg, y1seg
      write(iun,fmt='(1p,4e14.6)') x2seg, y2seg, sseg, ttseg
c AT arc
      xa = 0.d0
      ya = 0.d0
      tt1 = ttseg 
      dtt = at * deg2rd
      call prArc(xa,ya,rme/2.d0,tt1,dtt,iun)
      txt4 = 'AT'
c      CALL TRTXT(x2seg,y2seg,TXT4,4,1)
c entrance sector limit
      x1seg = 0.d0
      y1seg = 0.d0
      ttseg = sign * op * deg2rd  +(2.d0 * ACN -AT/2.d0) * deg2rd + ttrf
      sseg = 1.5d0 * rme
      x2seg = sseg * cos(ttseg)
      y2seg = sseg * sin(ttseg)
      write(iun,*) 'LTYP 3     entr. sector'
      write(iun,fmt='(1p,2e14.6)') x1seg, y1seg
      write(iun,fmt='(1p,4e14.6)') x2seg, y2seg, sseg, ttseg
      txt4 = 'EFBi'
c      CALL TRTXT(x2seg,y2seg,TXT4,4,1)
      xa = 0.d0
      ya = 0.d0
      tt1 = (2.d0 * ACN -AT/2.d0) * deg2rd + ttrf
      dtt = op * deg2rd
      call prArc(xa,ya,rme*1.3d0,tt1,dtt,iun)
c      marker at rm
      xa = rme * cos(ttseg)
      ya = rme * sin(ttseg)
      call prArc(xa,ya,rme/50.d0,0.d0,4.*pi2,iun)
c exit sector limit
      x1seg = 0.d0
      y1seg = 0.d0
      ttseg = sign * om * deg2rd+(2.d0 * ACN -AT/2.d0) * deg2rd + ttrf
      sseg = 1.5d0 * rme
      x2seg = sseg * cos(ttseg)
      y2seg = sseg * sin(ttseg)
      write(iun,*) 'LTYP 3       exit sector'
      write(iun,fmt='(1p,2e14.6)') x1seg, y1seg
      write(iun,fmt='(1p,4e14.6)') x2seg, y2seg, sseg, ttseg
      txt4 = 'EFBo'
c      CALL TRTXT(x2seg,y2seg,TXT4,4,1)
      xa = 0.d0
      ya = 0.d0
      tt1 = (2.d0 * ACN -AT/2.d0) * deg2rd + ttrf
      dtt = om * deg2rd 
      call prArc(xa,ya,rme*1.4d0,tt1,dtt,iun)
c      marker at rm
      xa = rme * cos(ttseg)
      ya = rme * sin(ttseg)
      call prArc(xa,ya,rme/50.d0,0.d0,4.*pi2,iun)
c entrance efb
      ttseg = sign * (op+xie) * deg2rd +
     >   (2.d0 * ACN -AT/2.d0) * deg2rd + ttrf
      sseg = 0.2d0 * rme
      x2seg = rme * cos(op*deg2rd +(2.d0 * ACN -AT/2.d0) *deg2rd + ttrf)
      y2seg = rme * sin(op*deg2rd+(2.d0 * ACN -AT/2.d0) * deg2rd + ttrf)
      x1seg = x2seg + sseg * cos(ttseg)
      y1seg = y2seg + sseg * sin(ttseg)
      write(iun,*) 'LTYP  2       entr. EFB'
      write(iun,fmt='(1p,2e14.6)') x1seg, y1seg
      write(iun,fmt='(1p,4e14.6)') x2seg, y2seg, sseg, ttseg
      x1seg = x2seg - sseg * cos(ttseg)
      y1seg = y2seg - sseg * sin(ttseg)
      write(iun,fmt='(1p,2e14.6)') x1seg, y1seg
      tta0 = sign * op * deg2rd +(2.d0 * ACN -AT/2.d0) * deg2rd + ttrf
      write(iun,*) 'LTYP  1'
      call drawspi(tta0,rme,rme*.6d0,rme*1.3d0,sign * xie,iun)
c exit efb
      ttseg = sign * (om+xis) * deg2rd+(2.d0*ACN-AT/2.d0)*deg2rd+ttrf
      sseg = 0.2d0 * rme
      x2seg = rme * cos(om*deg2rd+(2.d0 * ACN -AT/2.d0) * deg2rd + ttrf)
      y2seg = rme * sin(om*deg2rd+(2.d0 * ACN -AT/2.d0) * deg2rd + ttrf)
      x1seg = x2seg + sseg * cos(ttseg)
      y1seg = y2seg + sseg * sin(ttseg)
      write(iun,*) 'LTYP  2        exit EFB'
      write(iun,fmt='(1p,2e14.6)') x1seg, y1seg
      write(iun,fmt='(1p,4e14.6)') x2seg, y2seg, sseg, ttseg
      x1seg = x2seg - sseg * cos(ttseg)
      y1seg = y2seg - sseg * sin(ttseg)
      write(iun,fmt='(1p,2e14.6)') x1seg, y1seg
      tta0 = sign * om * deg2rd+(2.d0 * ACN -AT/2.d0) * deg2rd + ttrf
      write(iun,*) 'LTYP  1'
      call drawspi(tta0,rme,rme*.6d0,rme*1.3d0,sign * xis,iun)

c      CLOSE(IUN)
      GOTO 99

 997  WRITE(6,*) ' Error upon OPEN statement '

 99   RETURN
      END 
