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
      SUBROUTINE DIAGNU(OKECH,NLOG)
      implicit double precision (a-h,o-z)
C----------------------------------------------------------------------------
C     DIAGRAMME NOMBRE D'ONDES
C----------------------------------------------------------------------------
      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB

      logical okech
      logical TRDROI, TRTOUR,xtrace
      logical SYSTEM, ALEAT

      save ires, irot, MM, nharm 

C nmail = number of superperiods 
C iqx, iqy = integer part of tune 
C AGS
C      data nmail, iqx, iqy /  12, 8, 8  /
C      DATA QXMI,XXMA,QYMI,YYMA/0.68d0,1.d0,0.68d0,1.d0/
C      DATA QXMI,XXMA,QYMI,YYMA/0.5d0,1.d0,0.5d0,1.d0/
C Luciano's cyclo
c      data nmail, iqx, iqy /  8, 0, 0  /
c      DATA QXMI,XXMA,QYMI,YYMA/1.0499d0,2.05001d0,0.299d0,1.03001d0/
C Fred's ETparam
      data nmail, iqx, iqy /  1, 0, 0  /
C      DATA QXMI,XXMA,QYMI,YYMA/ 0. , .5, 0., .5 /
C      DATA QXMI,XXMA,QYMI,YYMA/ 0.65 , 1., 0.65, 1. /
C      DATA QXMI,XXMA,QYMI,YYMA/ 0.65d0 , .85d0, 0.85d0, 1.d0 /
C      DATA QXMI,XXMA,QYMI,YYMA/ 0.65d0 , .9d0, 0.75d0, 1.d0 /
      DATA QXMI,XXMA,QYMI,YYMA/ 0.d0 , 1.8d0, 0.d0, 1.8d0 /
C      DATA QXMI,XXMA,QYMI,YYMA/ 3.5d0 , 4d0, 3.5d0, 4d0 /

      data ires, irot, MM, nharm / 3, 3, -4, 500 /

      write(*,*) ' ++ TUNE DIAGRAM MQx+NQy=P'
      write(*,*) ' '

      ires0 = ires
      irot0 = irot 
      MM0 = MM
      nharm0 = nharm

c      qxmin = qxmi + dble(iqx) 
c      xxmax = xxma + dble(iqx )
c      qymin = qymi + dble(iqy )
c      yymax = yyma + dble(iqy )

      WRITE(*,fmt='(a,1p,4e14.6)') 
     >' Scales (QX_min/max, Qy_min/max)       : '
     > , qxmi + dble(iqx), xxma + dble(iqx )
     > , qymi + dble(iqy ), yyma + dble(iqy ) 
C      read(*,fmt='(4e18.2)',err=5,end=5) QXMIi,XXMAi,QYMIi,YYMAi
C      goto 51
C 5    continue
      QXMIi = QXMI
      XXMAi = XXMA
      QYMIi = QYMI
      YYMAi = YYMA
      qxmin = qxmi + dble(iqx) 
      xxmax = xxma + dble(iqx )
      qymin = qymi + dble(iqy )
      yymax = yyma + dble(iqy )
C51    continue
      QXMI = QXMIi
      XXMA = XXMAi
      QYMI = QYMIi
      YYMA = YYMAi
      QXMIN = QXMI
      XXMAX = XXMA
      QYMIN = QYMI
      YYMAX = YYMA
c       write(*,*) ' diagnu ',okech,qXMIn,xXMAx,qYMIn,yYMAx
c                   read(*,*)

      DY=YYMAX-QYMIN
      DX=XXMAX-QXMIN
c      IF((DX-DY).LE.0.) THEN
c         QXMAX=QXMIN+DY
c         QYMAX=QYMIN+DY
c      ELSE
c         QXMAX=QXMIN+DX
c         QYMAX=QYMIN+DX
c      ENDIF
          QXMAX=QXMIN+DX
         QYMAX=QYMIN+DY
      XXMAX=QXMAX
      YYMAX=QYMAX

c      if(.not. okech) then 
        xmi = QXMIN
        xma = QXMAX
        ymi = QYMIN
        yma = QYMAX
        CALL TRAXES(XMI,XMA,YMI,YMA,1)
        okech = .true.
c      else
c        QXMIN = xmi
c        QXMAX = xma
c        QYMIN = ymi
c        QYMAX = yma
c        CALL TRAXES(XMI,XMA,YMI,YMA,1)
c      endif
c       write(*,*) ' diagnu ',okech,XMI,XMA,YMI,YMA
c                   read(*,*)

      CALL FBGTXT

      write(*,*) ' TYPE OF RESONANCE (SY,RA,AL) (1/2/3) : ',ires
      read(*,fmt='(i6)',err=1,end=1) ires
      if(ires.eq.0) ires = ires0
      goto 11
 1    continue
      ires = ires0
 11   continue
      ires0 = ires

      SYSTEM=(IRES.NE.2)
      ALEAT =(IRES.NE.1)
C
      write(*,*) ' REGULAR,SKEW,ALL (RE,SK,AL)  (1/2/3) : ',irot
      read(*,fmt='(i6)',err=2,end=2) irot
      if(irot.eq.0)  irot = irot0
      goto 21
 2    continue
      irot = irot0
 21   continue
      irot0 = irot

      TRDROI=(IROT.EQ.1).OR.(IROT.EQ.3)
      TRTOUR=(IROT.EQ.2).OR.(IROT.EQ.3)
 
      WRITE(*,*) ' Which order M+N (negative for 1 -> M+N) : ',MM
      read(*,fmt='(i6)',err=3,end=3)  MM
      if(MM.eq.0)  MM = MM0
      goto 31
 3    continue
      MM = MM0
 31   continue
      MM0 = MM

      if (MM.lt.0) then
        MM1 = 1
        MM2 = -MM
      else
        MM1 = MM
        MM2 = MM
      endif

      WRITE(*,*) ' MAX HARMONIC           : ',NHARM
      read(*,fmt='(i6)',err=4,end=4)  NHARM
      if(nharm.eq.0)  nharm = nharm0
      goto 41
 4    continue
      nharm=nharm0
 41   continue
      nharm0=nharm

c          write(*,*) ires,irot,MM1,MM2,nharm
c              stop
      DO 10 IH=0,NHARM
         IF(MOD(IH,NMAIL) .EQ. 0) THEN
            XTRACE=SYSTEM
         ELSE
            XTRACE=ALEAT
         ENDIF
         IF(XTRACE) THEN
 
            DO 20 LM = MM1, MM2
            DO 20 M=-LM,LM
              N=LM-IABS(M)
c  22           WRITE(6,*)
c               WRITE(6,*) ' Continue (Y/N) ?' 
c               READ(*,FMT='(A1)',ERR=22) REP        
c               IF(REP .NE. 'N' .AND. REP .NE. 'n') REP = 'y'
c               IF(REP.NE. 'y') goto 30
C               IF (INTRPT().NE.0) GO TO 30
 
              IF(MOD(N,2) .EQ. 0) THEN
                 XTRACE=TRDROI
              ELSE
                 XTRACE=TRTOUR
              ENDIF
              IF(XTRACE) THEN
                CALL FBGTXT
                CALL DIAGQ(nlog,lis,M,N,IH,
     >       sngl(QXMIN),sngl(QXMAX),sngl(QYMIN),sngl(QYMAX))
                IF(N.NE.0) CALL DIAGQ(nlog,lis,M,-N,IH,
     >       sngl(QXMIN),sngl(QXMAX),sngl(QYMIN),sngl(QYMAX))
              ENDIF
 
20          CONTINUE
         ENDIF
10    CONTINUE

      RETURN
      END   
