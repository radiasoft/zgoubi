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
C      SUBROUTINE LEVEL
C      PARAMETER (NLEV=5,NLEV2=2*NLEV)
C      CHARACTER KLEV*NLEV2, TXT*2
C      SAVE KLEV
C
C      ENTRY LEVIN(IOPT)
C      IF(IOPT.LE.9) THEN
C        IZERO = 0
C        WRITE(TXT,FMT='(2I1)') IZERO,IOPT
C      ELSE
C        WRITE(TXT,FMT='(I2)') IOPT
C      ENDIF
C      KLEV = KLEV(3:NLEV2)//OPT
C      RETURN
C
C      ENTRY LEVOUT(KLEV)
C      RETURN
C      END
      SUBROUTINE HELP(SHELP)
      CHARACTER*(*) SHELP

      WRITE(6,*) '      '
      WRITE(6,*) '      '
      WRITE(6,*) '------------------- HELP ------------------------'
      WRITE(6,*) '      '
      WRITE(6,*) '            MENU LEVEL : ',SHELP
      WRITE(6,*) '      '
      WRITE(6,*) '------------------- HELP ------------------------'
      WRITE(6,*) '      '

      IF(SHELP .EQ. '1/07') THEN
        WRITE(6,FMT='(A,/)') 
     > ' In order to get a plot rightaway it is sufficient to :',
     >'      - first, select option 1 and open a zgoubi output file,',
     >'      - next,  select option 7.'

C23456789012345678901234567890123456789012345678901234567890123456789012

        I=IDLG('('' Press RETURN for more :'')','    ',1)
        WRITE(6,FMT='(/,A)') 
     >'Option 1: ',
     >'   to chose the file data to be read and plotted from. ',
     >'   1/1 : zgoubi.plt. Contains particle coordinates, fields ',
     >'   along trajectories, etc., inside optical elements.  ', 
     >'   zgoubi.plt is filled upon IL=2 in optical elements. ',
     >'   1/2 :  zgoubi.fai. Contains particle coordinates at fixed',
     >'   positions along the beam line.', 
     >'   zgoubi.fai is filled upon use of FAISCNL. ',
     >'   1/3 zgoubi.spn. Contains particle and spin coordinates',
     >'   along trajectories, etc., inside optical elements.  ', 
     >'   zgoubi.spn is filled upon use of SPNTRK. ',
     >'   1/5 other... To specify any data file name different from ',
     >'   the above. It is then requested to provide the data file ',
     >'   type, in terms of the suffix .plt, .fai, .spn, as defined', 
     >'   above - even if this is not the suffix YOU use.' 

C23456789012345678901234567890123456789012345678901234567890123456789012

        I=IDLG('('' Press RETURN for more :'')','    ',1)
        WRITE(6,FMT='(/,A)') 
     >'Option 2: to select which variable is to be plotted against ',
     >'   which other. Attention ! available variables depend on the ',
     >'   content of the data file currently opened (presumably ',
     >'   by means of Option 1)'

        I=IDLG('('' Press RETURN for more :'')','    ',1)
        WRITE(6,FMT='(/,A)') 
     >'Option 3: to select various options ',
     >'   2/1 : in case REBELOTE is used, several passages through ',
     >'   the structure are performed, numbered IPASS=1,NPASS+1, ',
     >'   with NPASS being the argument under REBELOTE. 2/1 allows',
     >'   skiping data according to "IPASS=1,NPASS+1,Interval".',
     >'   2/3 : will list - or alternately stop listing of - the ',
     >'   plotted coordinates ',
     >'   2/4 : will scale H and V axis wrt. one another. Usefull when',
     >'   axis both are in identical units : a circle will appear as ',
     >'   circle.',
     >'   2/5 : to change the number of bins in the histograms ',
     >'   (histograms are plotted when variable #28 is entered ',
     >'   upon Option 3).',
     >'   2/6 : to get a linear extension of particle tracks inside',
     >'   field maps (upon use of Option 11), although these normally',
     >'   stop at the map limit.',
     >'   2/7 : when using .plt or .spn type data files (see Option 1)',
     >'   only those data generated within the specified optical ',
     >'   element will be plotted ; when using .fai type data file, ',
     >'   only those data generated right next to the specified ',
     >'   optical element will be plotted (as long as FAISCNL does ',
     >'   appear right there !). The optical element number',
     >'   of concern is the number appearing next to the optical ',
     >'   element keyword in zgoubi.res',
     >'   2/8 : to select the trajectory to be plotted. Trajectory',
     >'   numbers can be found out clearly in zgoubi.res by use of ',
     >'   keyword FAISCEAU',
     >'   2/10 : to get a tag attached with each plotted curve. Also ',
     >'   allows positionning of tags on the graph.',
     >'   2/11 : the path length S is continuously increased when ',
     >'   performing multiturn tracking (with REBELOTE). This option', 
     >'   allows reseting the plotted (image of) S to zero, turn after',
     >'   turn, which permits such plots as multiturn beam envelop ',
     >'   generation.'

      ELSE
        WRITE(6,FMT='(A,/)') ' To be provisionned...'
      ENDIF
      RETURN
      END
