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
C  François Meot <fmeot@bnl.gov>
C  Brookhaven National Laboratory 
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  -------
      SUBROUTINE MATRIC(JORD,JFOC,KWR,OKCPLD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ------------------------------------
C     Compute transfer matrix coefficients
C     ------------------------------------
      LOGICAL OKCPLD
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
C      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
C      COMMON/DON/ A(09876,99),IQ(09876),IP(09876),NB,NOEL
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
 
C------         R_ref    +dp/p     -dp/p
      DIMENSION R(6,6), RPD(6,6), RMD(6,6) 
      DIMENSION T(6,6,6)
      DIMENSION T3(5,6) , T4(5,6)
      SAVE R,T, T3,       T4

C------        Beam_ref    +dp/p     -dp/p
      DIMENSION F0(6,6), F0PD(6,6), F0MD(6,6) 

      LOGICAL KWRMAT

      LOGICAL PRDIC

      DIMENSION RO(6,6)
C      CHARACTER(140) BUFFER
      INTEGER DEBSTR, FINSTR
      LOGICAL IDLUNI, OK

      CHARACTER(300) CMMND

      DATA KWRMAT / .FALSE. /

C      IF(NRES.LE.0) RETURN

      IF(.NOT. (KOBJ.EQ.5 .OR. KOBJ.EQ.6)) THEN
        IF(NRES.GT.0)
     >  WRITE(NRES,FMT='('' Matrix  cannot  be  computed :  need "OBJET" 
     >  with  KOBJ=5 or 6'')')
        RETURN
      ENDIF

C      IORD = A(NOEL,1)
      IORD = JORD
      IF(IORD .EQ. 0) THEN
        IF(NRES.GT.0) THEN
          WRITE(NRES,FMT='(/,9X,'' Matrix not computed : IORD = 0'',/)')
          CALL IMPTRA(1,IMAX,NRES)
        ENDIF
        RETURN
      ENDIF
      IF(KOBJ .EQ. 5) THEN
        IORD=1
      ELSEIF(KOBJ .EQ. 6) THEN
        IORD=2
      ENDIF 

C      IFOC = A(NOEL,2) 
      IFOC = JFOC
      PRDIC = IFOC .GT. 10
      NMAIL = IFOC-10

C      KWRMAT = NINT(A(NOEL,3)) .EQ. 1
      KWRMAT = KWR .EQ. 1
C      IF(KWRMAT) CALL MATIM6(KWRMAT)
      CALL MATIM6(KWRMAT)
      IF(KWRMAT .AND. NRES .GT. 0) THEN
        WRITE(NRES,*) 
        WRITE(NRES,*)  ' Matrix coefficients are printed in '
     >  // ' zgoubi.MATRIX.out.'
      ENDIF

      IF    (IORD .EQ. 1) THEN
        CALL OBJ51(
     >             NBREF)
        IREF = 0
 1      CONTINUE
          IREF = IREF + 1

          IT1 = 1 + 11 * (IREF-1)
          IT2 = IT1+3
          IT3 = IT1+4

C FM, Nov. 2008
C          CALL REFER(1,IORD,IFOC,IT1,IT2,IT3)
          IFC = IFOC
          IF(PRDIC) IFC = 0
          CALL REFER(1,IORD,IFC,IT1,IT2,IT3)
          CALL MAT1(IT1, 
     >                  R,T)
          CALL MKSA(IORD,R,T,T3,T4)
          IF(PRDIC) CALL TUNES(R,F0,NMAIL,IERY,IERZ,.TRUE.,
     >                                                YNU,ZNU,CMUY,CMUZ)
          CALL MATIMP(R,F0,YNU,ZNU,CMUY,CMUZ,NMAIL,PRDIC,IT1)

C FM, Nov. 2008
          CALL REFER(2,IORD,IFC,IT1,IT2,IT3)

          IF(IREF.LT.NBREF) GOTO 1
          
      ELSEIF(IORD .EQ. 2) THEN

C FM, Nov. 2008
C        CALL REFER(1,IORD,IFOC,1,6,7)
        IREF = 1
        IFC = IFOC
        IF(PRDIC) IFC = 0
        CALL REFER(1,IORD,IFC,1,6,7)
        CALL MAT2(R,T,T3,T4)
        CALL MKSA(IORD,R,T,T3,T4)
        IF(PRDIC) CALL TUNES(R,F0,NMAIL,IERY,IERZ,.TRUE.,
     >                                                YNU,ZNU,CMUY,CMUZ)
        CALL MATIMP(R,F0,YNU,ZNU,CMUY,CMUZ,NMAIL,PRDIC,iref)
        CALL MATIM2(R,T,T3)
        IF(PRDIC) THEN 
          CALL MAT2P(RPD,DP)
          CALL MKSA(IORD,RPD,T,T3,T4)
          CALL TUNES(RPD,F0PD,NMAIL,IERY,IERZ,.TRUE.,
     >                                              YNUP,ZNUP,CMUY,CMUZ)
          CALL MATIMP(RPD,F0PD,YNUP,ZNUP,CMUY,CMUZ,NMAIL,PRDIC,iref)
          CALL MAT2M(RMD,DP)
          CALL MKSA(IORD,RMD,T,T3,T4)
          CALL TUNES(RMD,F0MD,NMAIL,IERY,IERZ,.TRUE.,
     >                                              YNUM,ZNUM,CMUY,CMUZ)
          CALL MATIMP(RMD,F0MD,YNUM,ZNUM,CMUY,CMUZ,NMAIL,PRDIC,iref)
C Momentum detuning
          NUML = 1
C          DNUYDP = (YNUP-YNUM)/2.D0/A(NUML,25)
C          DNUZDP = (ZNUP-ZNUM)/2.D0/A(NUML,25)
          DNUYDP = (YNUP-YNUM)/2.D0/DP
          DNUZDP = (ZNUP-ZNUM)/2.D0/DP
          IF(NRES .GT. 0) WRITE(NRES,FMT='(/,34X,'' Chromaticities : '',
     >      //,30X,''dNu_y / dp/p = '',G15.8,/, 
     >         30X,''dNu_z / dp/p = '',G15.8)') DNUYDP, DNUZDP
C             write(nres,*) dp, a(numl,25)
        ENDIF
C        CALL REFER(2,IORD,IFOC,1,6,7)
        CALL REFER(2,IORD,IFC,1,6,7)
      ENDIF

c---------------------------------------------------------------------------------
c      write(*,*)
c     >'Exportation of the matrix coefficients (coupled formalism, Fred)'
c            read(*,*)
c---------------------------------------------------------------------------------


C       RETURN



      IF(.NOT. PRDIC) RETURN
      IF(.NOT. OKCPLD) RETURN

c      OK = IDLUNI(
c     >            LUNR)
c      IF(.NOT. OK) CALL ENDJOB(
c     >'SBR MATRIC. Problem open idle unit for READ. ',-99)
cC Just to spare former computation results
c      cmmnd = 'cp transfertM.dat transfertM_save.dat'
c      write(6,*) ' Pgm matric. Now doing ' 
c     > // cmmnd(debstr(cmmnd):finstr(cmmnd))
c      CALL SYSTEM(cmmnd)
c      OPEN(lunR,FILE='transfertM_save.dat',STATUS='UNKNOWN',IOSTAT=IOS1)
      OK = IDLUNI(
     >            LUNW)
      IF(.NOT. OK) CALL ENDJOB(
     >'SBR MATRIC. Problem open idle unit for WRITE. ',-99)
      OPEN(lunW,FILE='transfertM.dat',STATUS='UNKNOWN',IOSTAT=IOS2)
c----------------------------------------------

c 3    CONTINUE           
c        Reading the 1-turnM.dat
c        READ(lunR,FMT='(a)',END=39,ERR=38) BUFFER
c        Writing the new transfertM.dat
c        WRITE(lunW,FMT='(a)') BUFFER(DEBSTR(BUFFER):FINSTR(BUFFER))
c      GOTO 3

c 38   CONTINUE     
c      WRITE(6,FMT='(//,''SBR matric. 
c     >             ERROR while reading transfertm.dat.'',/,/)')
c 39   CONTINUE

      WRITE(lunW,FMT='(//)')
      WRITE(lunW,FMT='(
     >''TRANSPORT MATRIX (written by Zgoubi, for use by ETparam ):'')')
      DO I=1,4
         WRITE(lunW,FMT='(4(F15.8,1X))') (R(I,J),J=1,4)
      ENDDO

C      CLOSE(lunR,IOSTAT=IOS1)
      CLOSE(lunW,IOSTAT=IOS2)
c----------------------------------------------

c      cmmnd = '/home/owl/fmeot/zgoubi/current/coupling/ETparam'
c      CALL SYSTEM(cmmnd)

      call tunesc

c            read(*,*)

      IF(NRES .GT. 0) WRITE(NRES,FMT='(/,'' Pgm matric, '',
     >'' now calling et2res.  '')')

      call et2res(nres)
      call et2re1(
     >             F011,f012,f033,f034,phy,phz,Cstrn)
      IF(NRES .GT. 0) then
        WRITE(NRES,*)
        WRITE(NRES,*) '--------------------------------------'
        WRITE(NRES,*) ' Pgm matric. '
        WRITE(NRES,*)
        WRITE(NRES,*) ' Coupled modes : '
        WRITE(NRES,*) ' bet1, alf1 : ',    F011,-f012
        WRITE(NRES,*) ' bet2, alf2 : ',    f033,-f034
        WRITE(NRES,*) ' Q1, Q2 :     ',    phy,phz  
        WRITE(NRES,*) ' Coupling strength :     ',    Cstrn
        WRITE(NRES,*)
        WRITE(NRES,*) '--------------------------------------'
        WRITE(NRES,*)

      ENDIF

      RETURN
 
      ENTRY MATRI1(
     >             RO)
      DO IB = 1, 6
        DO IA = 1, 6
          RO(IA,IB)  = R(IA,IB)
        ENDDO
      ENDDO
      RETURN
      END
