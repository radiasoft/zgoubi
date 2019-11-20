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
      SUBROUTINE IMPTRA(IMAX1,IMAX2,LUN)
      use data_partition_ixfc, only : data_partition
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "MAXTRA.H"
      INCLUDE "C.DESIN.H"     ! COMMON/DESIN/ FDES(7,MXT),IFDES,KINFO,IRSAR,IRTET,IRPHI,NDES
C     >,AMS,AMP,AM3,TDVM,TETPHI(2,MXT)
C     1,AMS ,AMP,ENSTAR,BSTAR,TDVM,TETPHI(2,MXT)
      INCLUDE "MAXCOO.H"
      INCLUDE "C.OBJET.H"     ! COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT,KZOB
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      CHARACTER(1) LET
      INCLUDE "C.FAISCT.H"     ! COMMON/FAISCT/ LET(MXT)
C      LOGICAL ZSYM
      INCLUDE "C.TYPFLD.H"     ! COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      INCLUDE "C.PTICUL.H"     ! COMMON/PTICUL/ AMASS,Q,G,TO
      INCLUDE "C.REBELO.H"   ! COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,HDPRF,DP,QBR,BRI
      INCLUDE "C.UNITS.H"     ! COMMON/UNITS/ UNIT(MXJ)

      DIMENSION SIG(4,4)
      CHARACTER(5) TXT(3)
      CHARACTER(10) UU(3)
      type(data_partition) particle_set

      DATA TXT / '(Y,T)', '(Z,P)', '(t,K)' /
      DATA UU / '(cm,rd)', '(cm,rd)', '(mu_s,MeV)' /

      CALL ZGNOEL(
     >            NOEL)

      write(6,*) 'imptra(): calling gather(F) on image', this_image()
      flush(6)
      call particle_set%gather(F, result_image=1)
      write(6,*) 'imptra(): calling gather(FO) on image', this_image()
      flush(6)
      call particle_set%gather(FO, result_image=1)
      ! ^^^^^ Are we sure we need this?
      !call gather(FDES, result_image=1)
      ! ^^^^^ Do we also need this?

      IF(LUN .EQ. 0) RETURN

      block
      logical, parameter :: all_images_write=.true.
      associate( me => this_image() )
      image_1_writes: if (me==1 .or. all_images_write) then

      WRITE(LUN,100) NOEL-1,IMAX2-IMAX1+1
 100  FORMAT('0',45X,'TRACE DU FAISCEAU',
     > /,43X,'(follows element # ',I6,')',
     > /,45X,I6,' TRAJECTOIRES',//
     >,35X,'OBJET',50X,'FAISCEAU',//
     >,10X,'D',7X,'Y(cm)',5X,'T(mr)',5X,'Z(cm)',5X,'P(mr)',7X,'S(cm)'
     >,6X,'D-1',5X,'Y(cm)',4X,'T(mr)',4X,'Z(cm)',4X,'P(mr)',6X,'S(cm)'
     >,/)

      DO 1 I=IMAX1,IMAX2
        WRITE(LUN,101) LET(I),IEX(I),(FO(J,I),J=1,6)
     >  ,F(1,I)-1.D0,(F(J,I),J=2,5),F(6,I),I
 101    FORMAT(A1,1X,I2,1X,F8.4,4(1X,F9.3),1X,F11.4,
     >                     2X,F8.4,4(1X,F8.3),1X,(1P,E14.6),1X,I5)
        IF(AMASS .NE. 0D0) THEN
          IF(IFDES.EQ.1) THEN
            WRITE(LUN,FMT='(15X,''Time of flight (mus) :'',
     >      1P,G16.8,'' mass (MeV/c2) :'',G14.6, ''    decay at (m) :'',
     >           G14.6)') F(7,I),AMQ(1,I),FDES(6,I)*UNIT(5)
          ELSE
            WRITE(LUN,FMT='(15X,''Time of flight (mus) :'',
     >      1P,G16.8,'' mass (MeV/c2) :'', G14.6)') F(7,I),AMQ(1,I)
          ENDIF
        ENDIF
 1    CONTINUE

Compute rms ellipse
      WRITE(LUN,FMT='(//,A)') '---------------'//
     >'  Concentration ellipses : '
      WRITE(LUN,FMT='(T4,''surface'',T17,''alpha'',T30,''beta'',
     >T43,''<X>'',T58,''<XP>'',T73,''numb. of prtcls'',
     >''   ratio      space      pass# '')')
      WRITE(LUN,FMT='(T73,''in ellips,  out '')')
      PI = 4.D0 * ATAN(1.D0)
      DO JJ = 1, 3
        CALL LPSFIT(JJ,
     >                 EMIT,ALF,BET,XM,XPM)
Compute number of particles alive and numberinside ellipse
        CALL CNTINL(JJ,PI*EMIT,ALF,BET,XM,XPM,
     >                                        NLIV,NINL)
        RATIN = DBLE(NINL)/DBLE(IMAX)
        WRITE(LUN,110)
     >  PI*EMIT,ALF,BET,XM,XPM,NLIV,NINL,RATIN,TXT(JJ),IPASS
 110    FORMAT(1P,3(1X,E12.4),2(1X,E14.6),2(1X,I8),1X,G12.4,2X,A,2X,I8)
      ENDDO

      DO JJ = 1, 3
        CALL LPSFIT(JJ,
     >                 EMIT,ALF,BET,XM,XPM)
        WRITE(LUN,fmt='(1P,/,A,2(/,5X,A,E14.6))')
     >  TXT(JJ)//'  space (units : '//UU(JJ)//') :  ',
     >  ' sigma_'//TXT(JJ)(2:2)//' = sqrt(Surface/pi * BET) = ',
     >  sqrt(emit * BET) ,
     >  ' sigma_'//TXT(JJ)(4:4)//
     >                   ' = sqrt(Surface/pi * (1+ALF^2)/BET) = ',
     >  SQRT(EMIT * (1.D0+ALF**2)/BET)
      ENDDO

Compute 4-D sigma matrix
      WRITE(LUN,FMT='(//,''  Beam  sigma  matrix : '',/)')
      CALL LPSFI4(
     >             SQX,SQZ,SIG)
      WRITE(LUN,120) ((SIG(I,J),J=1,4),I=1,4)
 120  FORMAT(1P,4(1X,E14.6))
      WRITE(LUN,fmt='(/,5X,1P,A,2(2X,E14.6),3X,A)')
     >' sqrt(det_Y), sqrt(det_Z) : ', SQX, SQZ,
     >' (Note :  sqrt(determinant) = ellipse surface / pi)'

      end if image_1_writes
      end associate
      end block

      RETURN
      END
