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
C  Brookhaven National Laboratory            és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      PARAMETER (MXKLE=107)
      CHARACTER*(KSIZ)  KLE(MXKLE)
C           1           2            3            4            5
      DATA (KLE(I),I=1,50)/
     >     'ESL'  ,   'AIMANT' , 'QUADRUPO' , 'SEXTUPOL' ,   'IMAGE'  ,
     1   'IMAGES' , 'FAISCEAU' , 'FAISCNL'  ,   'CIBLE'  ,  'FOCALE'  ,
     > 'REBELOTE' , 'QUADISEX' , 'CHANGREF' ,  'SEXQUAD' ,  'TOROID'  ,
     2  'TOSCA3D' ,  'OBJETA'  ,  'MATRIX'  ,  'CHAMBR'  , 'QUADSPS2' ,
     > 'CARTEMES' ,    'YMY'   , 'OCTUPOLE' ,   'OBJET'  ,  'MCOBJET' ,
     3 'MCDESINT' ,   'HISTO'  , 'TRANSMAT' ,   'VENUS'  ,   'PS170'  ,
     >   'TOSCA'  ,  'AUTOREF' , 'COLLIMA'  , 'MULTIPOL' ,  'SEPARA'  ,
     4   'RESET'  ,  'ORDRE'   ,   'DRIFT'  , 'DIPOLE-M' ,  'TARGET'  ,
     >   'IMAGEZ' , 'IMAGESZ'  ,  'FOCALEZ' , 'PLOTDATA' ,  'BINARY'  ,
     5     'FIT'  ,  'SPES3'   ,  'CHALUT'  , 'SPNTRK'   ,  'SPNPRT'  /
      DATA (KLE(I),I=51,100)/
     >     'BEND' , 'SOLENOID' ,  'PLTLPS'  , 'SPNPRNL'  ,  'CAVITE'  ,
     6 'PARTICUL' ,  'SCALING' ,  'BREVOL'  , 'POISSON'  , 'DIPOLE'  ,
     >  'DANFGSI' ,   'MAP2D'  , 'DECAPOLE' , 'DODECAPO' , 'WIENFILT' ,
     7 'FAISTORE' , 'SPNSTORE' ,  'ELREVOL' ,  'EL2TUB'  ,  'UNIPOT'  ,
     >   'ELMULT' ,  'EBMULT'  ,  'ESLVIR'  , 'TROBJET'  , 'TRANSLAT' ,
     8 'POLARMES' ,  'TRAROT'  ,  'SRLOSS'  , 'PICKUPS'  ,  'OPTICS'  ,
     >   'GASCAT' ,'UNDULATOR' , 'ELCYLDEF' ,  'ELMIR'   ,  'ELCMIR'  ,
     9  'MAP2D_E' ,  'SRPRNT'  , 'BETATRON' ,  'TWISS'   ,   'END'    , 
     >    'FFAG'  ,  'HELIX'   ,   'CSR'    ,   'PATH'   ,   'COILS'  ,
     A 'GETFITVAL', 'SUPERPOSE',  'MARKER'  , 'DIPOLES'  , 'TRACKING' /
      DATA (KLE(I),I=101,MXKLE)/
     > 'FFAG-SPI' ,   'FIT2'   ,   'EMMA'   , 'DIPOLEC'  , 'REVERSE'  ,
     B  'SYSTEM'  ,  'SPINR'/

      GOTO 1000

 1001 GOTO(1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,31,17,18,19,20
     >,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37, 1,39, 9,41
     >,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62
     >,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83
     >,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103
     >,104,105,106,107
     > ),IKLE
      GOTO 999

 1000 CONTINUE
