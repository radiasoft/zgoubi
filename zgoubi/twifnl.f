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
      SUBROUTINE TWIFNL(LUNW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER(800) TXT400
      CHARACTER(62) TXT62
      INTEGER DEBSTR, FINSTR
      LOGICAL STRCON

      CLOSE(LUNW)
      CALL SYSTEM('\cp zgoubi.TWISS.out zgoubi.TWISS.out_temp ')

      OK = IDLUNI(
     >            LUNW)
      OPEN(UNIT=LUNW,FILE='zgoubi.TWISS.out')
      OK = IDLUNI(
     >            LUNR)
      OPEN(UNIT=LUNR,FILE='zgoubi.TWISS.out_temp')

      CALL OPTIC1(
     >            DXMA, DYMA, DXMI, DYMI, 
     >            XMA, YMA, XMI, YMI, 
     >            BTXMA, BTYMA, BTXMI, BTYMI, 
     >            XM, YM, DXM, DYM, 
     >            XM2, YM2, DXM2, DYM2, DLTP, NC)
      XRMS = SQRT(XM2/DBLE(NC) - (XM/DBLE(NC))**2)
      YRMS = SQRT(YM2/DBLE(NC) - (YM/DBLE(NC))**2)
      DXRMS = SQRT(DXM2/DBLE(NC) - (DXM/DBLE(NC))**2)
      DYRMS = SQRT(DYM2/DBLE(NC) - (DYM/DBLE(NC))**2)

 1    CONTINUE
        READ(LUNR,FMT='(A)',ERR=77,END=88) TXT400
        IF    (STRCON(TXT400,'@ DXMAX            %le'
     >                    ,IS)) THEN
          WRITE(TXT62,FMT='(1P,E16.8,5X,A,3X,E16.8)') 
     >    DXMA, '@ DYMIN            %le',DXMI
          TXT400 = TXT400(DEBSTR(TXT400):22)//'   '
     >    //TXT62(debstr(TXT62):finstr(TXT62))
        ELSEIF(STRCON(TXT400,'@ DYMAX            %le'
     >                    ,IS)) THEN
          WRITE(TXT62,FMT='(1P,E16.8,5X,A,3X,E16.8)') 
     >    DYMA, '@ DYMIN            %le',DYMI
          TXT400 = TXT400(DEBSTR(TXT400):22)//'   '
     >    //TXT62(debstr(TXT62):finstr(TXT62))
        ELSEIF(STRCON(TXT400,'@ XCOMAX           %le'
     >                    ,IS)) THEN
          WRITE(TXT62,FMT='(1P,E16.8,5X,A,3X,E16.8)') 
     >    XMA,  '@ XMIN             %le',XMI
          TXT400 = TXT400(DEBSTR(TXT400):22)//'   '
     >    //TXT62(debstr(TXT62):finstr(TXT62))
        ELSEIF(STRCON(TXT400,'@ YCOMAX           %le'
     >                    ,IS)) THEN
          WRITE(TXT62,FMT='(1P,E16.8,5X,A,3X,E16.8)')
     >    YMA,  '@ YMIN             %le',YMI
          TXT400 = TXT400(DEBSTR(TXT400):22)//'   '
     >    //TXT62(debstr(TXT62):finstr(TXT62))
        ELSEIF(STRCON(TXT400,'@ BETXMAX          %le'
     >                    ,IS)) THEN
          WRITE(TXT62,FMT='(1P,E16.8,5X,A,3X,E16.8)') 
     >    BTXMA,'@ BTXMI            %le',BTXMI
          TXT400 = TXT400(DEBSTR(TXT400):22)//'   '
     >    //TXT62(debstr(TXT62):finstr(TXT62))
        ELSEIF(STRCON(TXT400,'@ BETYMAX          %le'
     >                    ,IS)) THEN
          WRITE(TXT62,FMT='(1P,E16.8,5X,A,3X,E16.8)') BT
     >    YMA,  '@ BTYMI            %le',BTYMI
          TXT400 = TXT400(DEBSTR(TXT400):22)//'   '
     >    //TXT62(debstr(TXT62):finstr(TXT62))
        ELSEIF(STRCON(TXT400,'@ XCORMS           %le'
     >                    ,IS)) THEN
          WRITE(TXT62,FMT='(1P,E16.8)') XRMS
          TXT400 = TXT400(DEBSTR(TXT400):22)//'   '
     >    //TXT62(debstr(TXT62):finstr(TXT62))
        ELSEIF(STRCON(TXT400,'@ YCORMS           %l
     >    e'
     >                    ,IS)) THEN
          WRITE(TXT62,FMT='(1P,E16.8)') YRMS
          TXT400 = TXT400(DEBSTR(TXT400):22)//'   '
     >    //TXT62(debstr(TXT62):finstr(TXT62))
        ELSEIF(STRCON(TXT400,'@ DXRMS            %le'
     >           
     >             ,IS)) THEN
          WRITE(TXT62,FMT='(1P,E16.8)') DXRMS
          TXT400 = TXT400(DEBSTR(TXT400):22)//'   '
     >    //TXT62(debstr(TXT62):finstr(TXT62))
        ELSEIF(STRCON(TXT400,'@ DYRMS            %le'
     >                    ,IS)) THEN
          WRITE(TXT62,FMT='(1P,E16
     >    .8)') DYRMS
          TXT400 = TXT400(DEBSTR(TXT400):22)//'   '
     >    //TXT62(debstr(TXT62):finstr(TXT62))
        ELSEIF(STRCON(TXT400,'@ DELTAP           %le'
     >                    ,IS)) THEN
          WRITE(TXT62,FMT='(1P,E16.8)') DLTP
          TXT400 = TXT400(
     >    DEBSTR(TXT400):22)//'   '
     >    //TXT62(debstr(TXT62):finstr(TXT62))
        ELSEIF(STRCON(TXT400,'XXXXXXXXX'
     >                    ,IS)) THEN
        ENDIF
        WRITE(LUNW,FMT='(A)') TXT400(DEBSTR(TXT400):FINSTR(TXT400))        
      GOTO 1

 77   CONTINUE
 88   CONTINUE
      CLOSE(LUNR)
      CLOSE(LUNW)
      RETURN
      END
