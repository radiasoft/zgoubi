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
C  Upton, NY, 11973
C  -------
      SUBROUTINE ADPOL(ID,
     >             IK,B,DB,DDB,D3BX,D3BY,D4BX,D4BY,BT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C--------------------------------------------------------------------
C     Add up pole components of multipole lenses. Can be E or B field
C--------------------------------------------------------------------
      DIMENSION B(5,3),DB(3,3),DDB(3,3,3)
      DIMENSION D3BX(3,3,3), D3BY(3,3,3)
      DIMENSION D4BX(3,3,3,3) ,D4BY(3,3,3,3)
      DIMENSION BT(5,15) 
 
      GOTO(1,2,3) IK
 
 1    CONTINUE
        IF(ID.GE.0) THEN
          BT(1,1)= B(1,1)
          BT(1,2)= B(1,2)
          BT(1,3)= B(1,3)
 
          IF(ID.GE.1) THEN
            BT(2,1)= DB(1,1)
            BT(2,2)= DB(2,1)
            BT(2,3)= DB(3,1)
            BT(2,4)= DB(2,2)
            BT(2,5)= DB(3,2)          
 
            IF(ID.GE.2) THEN
              BT(3,1)= DDB(1,1,1)
              BT(3,2)= DDB(2,1,1)
              BT(3,3)= DDB(3,1,1)
              BT(3,4)= DDB(2,2,1)
              BT(3,5)= DDB(3,2,1)
              BT(3,7)= DDB(2,2,2)
              BT(3,8)= DDB(3,2,2)
 
              IF(ID.GE.3) THEN
                BT(4,1)= D3BX(1,1,1)
                BT(4,2)= D3BX(2,1,1)
                BT(4,3)= D3BX(3,1,1)
                BT(4,4)= D3BX(2,2,1)
                BT(4,5)= D3BX(3,2,1)
                BT(4,7)= D3BX(2,2,2)
                BT(4,8)= D3BX(3,2,2)
 
                BT(4,11)= D3BY(2,2,2)
                BT(4,12)= D3BY(3,2,2)
 
                IF(ID.GE.4) THEN
                  BT(5,1)= D4BX(3,1,1,1)
                  BT(5,2)= D4BX(3,2,1,1)
                  BT(5,3)= D4BX(3,2,2,1)
                  BT(5,4)= D4BX(3,3,3,1)
                  BT(5,5)= D4BX(3,2,2,2)
                  BT(5,6)= D4BX(3,3,3,2)
 
                  BT(5,7)= D4BY(3,2,2,2)
                  BT(5,8)= D4BY(3,3,3,2)
                ENDIF
C                  ID GE 4
              ENDIF
C                ID GE 3
            ENDIF
C              ID GE 2
          ENDIF
C            ID GE 1
        ENDIF
C          ID GE 0

C FM june 2015. Already in chamc        CALL RAZDRV(KFL)

        IK=2
 
      GOTO 99
 
 2    CONTINUE
        IF(ID.GE.0) THEN
          BT(1,1)=BT(1,1) + B(1,1)
          BT(1,2)=BT(1,2) + B(1,2)
          BT(1,3)=BT(1,3) + B(1,3)
 
          IF(ID.GE.1) THEN
            BT(2,1)=BT(2,1) + DB(1,1)
            BT(2,2)=BT(2,2) + DB(2,1)
            BT(2,3)=BT(2,3) + DB(3,1)
            BT(2,4)=BT(2,4) + DB(2,2)
            BT(2,5)=BT(2,5) + DB(3,2)
 
            IF(ID.GE.2) THEN
              BT(3,1)=BT(3,1) + DDB(1,1,1)
              BT(3,2)=BT(3,2) + DDB(2,1,1)
              BT(3,3)=BT(3,3) + DDB(3,1,1)
              BT(3,4)=BT(3,4) + DDB(2,2,1)
              BT(3,5)=BT(3,5) + DDB(3,2,1)
              BT(3,7)=BT(3,7) + DDB(2,2,2)
              BT(3,8)=BT(3,8) + DDB(3,2,2)
 
              IF(ID.GE.3) THEN
                BT(4,1)=BT(4,1) + D3BX(1,1,1)
                BT(4,2)=BT(4,2) + D3BX(2,1,1)
                BT(4,3)=BT(4,3) + D3BX(3,1,1)
                BT(4,4)=BT(4,4) + D3BX(2,2,1)
                BT(4,5)=BT(4,5) + D3BX(3,2,1)
                BT(4,7)=BT(4,7) + D3BX(2,2,2)
                BT(4,8)=BT(4,8) + D3BX(3,2,2)
 
                BT(4,11)=BT(4,11) + D3BY(2,2,2)
                BT(4,12)=BT(4,12) + D3BY(3,2,2)
 
                IF(ID.GE.4) THEN
                  BT(5,1)=BT(5,1) + D4BX(3,1,1,1)
                  BT(5,2)=BT(5,2) + D4BX(3,2,1,1)
                  BT(5,3)=BT(5,3) + D4BX(3,2,2,1)
                  BT(5,4)=BT(5,4) + D4BX(3,3,3,1)
                  BT(5,5)=BT(5,5) + D4BX(3,2,2,2)
                  BT(5,6)=BT(5,6) + D4BX(3,3,3,2)
 
                  BT(5,7)=BT(5,7) + D4BY(3,2,2,2)
                  BT(5,8)=BT(5,8) + D4BY(3,3,3,2)
 
                ENDIF
C                  ID GE 4
              ENDIF
C                ID GE 3
            ENDIF
C              ID GE 2
          ENDIF
C            ID GE 1
        ENDIF
C          ID GE 0

C FM june 2015. Already in chamc        CALL RAZDRV(KFL)

      GOTO 99
 
 3    CONTINUE
        IF(ID.GE.0) THEN
          B(1,1)=BT(1,1)
          B(1,2)=BT(1,2)
          B(1,3)=BT(1,3)
 
          IF(ID.GE.1) THEN
            DB(1,1)=BT(2,1) 
            DB(2,1)=BT(2,2) 
            DB(3,1)=BT(2,3) 
            DB(2,2)=BT(2,4) 
            DB(3,2)=BT(2,5) 
 
            IF(ID.GE.2) THEN
              DDB(1,1,1)=BT(3,1) 
              DDB(2,1,1)=BT(3,2) 
              DDB(3,1,1)=BT(3,3) 
              DDB(2,2,1)=BT(3,4) 
              DDB(3,2,1)=BT(3,5) 
              DDB(2,2,2)=BT(3,7) 
              DDB(3,2,2)=BT(3,8) 
c       write(89,*) b(1,3),DB(3,1),DB(3,2),DDB(3,2,2),DDB(3,1,1),
c     >         ' sbr adpol '
c      write(89,*)  DB(3,1),  DB(3,2), DB(1,1), DB(2,1),' sbr adpol '
 
              IF(ID.GE.3) THEN
                D3BX(1,1,1)=BT(4,1) 
                D3BX(2,1,1)=BT(4,2) 
                D3BX(3,1,1)=BT(4,3) 
                D3BX(2,2,1)=BT(4,4) 
                D3BX(3,2,1)=BT(4,5) 
                D3BX(2,2,2)=BT(4,7) 
                D3BX(3,2,2)=BT(4,8) 

                D3BY(2,2,2)=BT(4,11) 
                D3BY(3,2,2)=BT(4,12) 

                IF(ID.GE.4) THEN
                  D4BX(3,1,1,1)=BT(5,1) 
                  D4BX(3,2,1,1)=BT(5,2) 
                  D4BX(3,2,2,1)=BT(5,3) 
                  D4BX(3,3,3,1)=BT(5,4) 
                  D4BX(3,2,2,2)=BT(5,5) 
                  D4BX(3,3,3,2)=BT(5,6) 
 
                  D4BY(3,2,2,2)=BT(5,7) 
                  D4BY(3,3,3,2)=BT(5,8) 
 
                ENDIF
C                  ID GE 4
              ENDIF
C                ID GE 3
            ENDIF
C              ID GE 2
          ENDIF
C            ID GE 1
        ENDIF
C          ID GE 0
 
 99   CONTINUE

      RETURN
      END
