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
      SUBROUTINE RSCAL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ----------------------------------------------
C     READS DATA FOR POWER SUPPLIES OF LMNT FAMILIES
C     ----------------------------------------------
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER(80) TA
      COMMON/DONT/ TA(MXL,40)
      INCLUDE 'MXFS.H'
      COMMON/SCAL/SCL(MXF,MXS),TIM(MXF,MXS),NTIM(MXF),KSCL
      PARAMETER (LBLSIZ=10)
      PARAMETER (KSIZ=10)
      CHARACTER FAM*(KSIZ),LBF*(LBLSIZ)
      COMMON/SCALT/ FAM(MXF),LBF(MXF,MLF),JPA(MXF,MXP)
 
      PARAMETER(MSTR=MLF+1)
      CHARACTER*10 STRA(MSTR)

      INTEGER DEBSTR, FINSTR

      LOGICAL STRCON, IDLUNI
      CHARACTER(132) TXT132
      DIMENSION MODSCL(MXF)
      SAVE MODSCL

      DIMENSION NTIM2(MXF),SCL2(MXF,MXS),TIM2(MXF,MXS)

      DATA MODSCL / MXF*0  /
      DATA FAC / 1.d0  /

C----- IOPT; NB OF DIFFRNT FAMILIES TO BE SCALED (<= MXF)
      READ(NDAT,*) A(NOEL,1),A(NOEL,2) 

      NFAM = A(NOEL,2)
      IF(NFAM .GT. MXF) 
     >  CALL ENDJOB('SBR RSCAL - Too many families, maximum allowed is '
     >  ,MXF)

      DO 1 IFM=1,NFAM

C------- Store name of family and label(s)
        READ(NDAT,FMT='(A)')  TXT132

       I =   0

        CALL STRGET(TXT132,MSTR,
     >                          NSTR,STRA)

        IF(NSTR .GT. MSTR) 
     >     CALL ENDJOB('SBR RSCAL - Too many labels per family, max is '
     >     ,MLF)

        FAM(IFM) = STRA(1)(1:KSIZ)

        IF(NSTR .GE. 2) THEN
          DO  KL=1,NSTR-1
            LBF(IFM,KL) =  STRA(KL+1)(1:LBLSIZ)
          ENDDO
        ENDIF

        DO KL=NSTR, MLF
          LBF(IFM,KL) = ' '
        ENDDO

C JPA tells which parameter in the element (IFM) the scaling is to be applied to
C Case ERR accounts for older versions of zgoubi, w/o JPA
C        READ(NDAT,*,ERR=121) NTIM(IFM)  !!, JPA(IFM)
C        GOTO 122
C 121    CONTINUE
C        JPA(IFM) = 0.D0
c 122    CONTINUE
        READ(NDAT,FMT='(A)') TXT132

        IF( STRCON(TXT132,'.', 
     >                        IS)) THEN

          READ(TXT132(DEBSTR(TXT132):IS-1),*) NTIM(IFM)
          READ(TXT132(IS+1:FINSTR(TXT132)),*) MODSCL(IFM)
c Possible additional scaling factor : 
          IF(MODSCL(IFM) .EQ. 10)
     >      READ(TXT132(IS+3:FINSTR(TXT132)),*,ERR=44,END=44) FAC
          GOTO 45

 44       CONTINUE
          FAC = 1.D0
 45       CONTINUE

          IF(NTIM(IFM) .GT. MXS-2) 
     >       CALL ENDJOB('SBR RSCAL - Too many timings, max is ',MXS-2)
        ELSE
          READ(TXT132,*) NTIM(IFM)
          IF(NTIM(IFM) .GT. MXD-2) 
     >    CALL ENDJOB('SBR RSCAL - Too many timings, max is ',MXD-2)
          DO J = 1, MXP-1
            READ(TXT132,*,ERR=46,END=46) DUM, (JPA(IFM,K), K=1, J)
          ENDDO
 46       CONTINUE
          JPA(IFM,MXP) = J-1
        ENDIF        

        IF(NTIM(IFM) .GE. 0) THEN

          NDSCL=NTIM(IFM)
          NDTIM=NTIM(IFM)
          MAX=NTIM(IFM)

        ELSEIF(NTIM(IFM) .LT. 0) THEN

          IF    (NTIM(IFM) .EQ. -1) THEN
            NDSCL=1
            NDTIM=1
C               max = max(NDSCL,NDTIM)
            MAX=NDSCL  

          ELSEIF(NTIM(IFM) .EQ. -2) THEN
C--------- Field law for scaling FFAG, LPSC, Sept. 2007
            NDSCL=1
            NDTIM=1
C               max = max(NDSCL,NDTIM)
            MAX=NDSCL  
          ELSEIF(NTIM(IFM) .EQ. -60) THEN
C--------- K1, K2 laws for AGS dipoles, FM, BNL, Jan. 2011
            NDSCL=1
            NDTIM=1
C               max = max(NDSCL,NDTIM)
            MAX=NDSCL  
          ELSEIF(NTIM(IFM) .EQ. -77) THEN
C---------- Field law protn driver, FNAL, Nov.2000 :
            NDSCL=4
            NDTIM=2
C               max = max(NDSCL,NDTIM)
            MAX=NDSCL  

          ELSEIF(NTIM(IFM) .EQ. -88) THEN
C--------- AC dipole at  BNL
            NDSCL=4
            NDTIM=3
C               max = max(NDSCL,NDTIM)
            MAX=NDSCL  

          ELSEIF(NTIM(IFM) .EQ. -87) THEN
C--------- AGS, Q-jump quads with snales
            NDSCL=1
            NDTIM=3
C               max = max(NDSCL,NDTIM)
            MAX=NDTIM

          ENDIF
        ENDIF

        IF    ( MODSCL(IFM) .LT. 10) then             

          IF(IFM .GT. MXD/10 -1) CALL ENDJOB(
     >    'SBR RSCAL: Too many families to scale. Max is',MXD/10-1)

C--------- SCL(IFM,IT)
          IT = NDSCL
          IF(IT.GT.10) CALL ENDJOB(
     >    'SBR RSCAL: Number of scaling factors cannot exceed',10)
          READ(NDAT,*) (A(NOEL,10*IFM+IT-1),IT=1,NDSCL)
C--------- TIM(IFM,IT)
          IT = NDTIM
          IF(IT.GT.10) CALL ENDJOB(
     >    'SBR RSCAL: Number of timings cannot exceed',10)
          READ(NDAT,*) (A(NOEL,10*IFM+NDSCL+IT-1),IT=1,NDTIM)
          IF(10*IFM+2*MAX.GT.MXD) CALL ENDJOB(
     >    'SBR RSCAL: Too maniy families or too many timings',-99)
          A(NOEL,10*IFM+2*MAX) = NSTR

        ELSEIF( MODSCL(IFM) .GE. 10) then             
C          Name of the storage file in the next line

          READ(NDAT,FMT='(A)') TA(NOEL,IFM)
          READ(NDAT,*) NTIMCOL, NSCALCOL
  
          IF(IDLUNI(
     >              LUN)) THEN
            OPEN(UNIT=LUN,FILE=TA(NOEL,IFM)(DEBSTR(TA(NOEL,IFM)):
     >      FINSTR(TA(NOEL,IFM))),STATUS='OLD',ERR=96)
          ELSE
            GOTO 97
          ENDIF

          NTIM(IFM) = 1 !number of timing is determined during lecture
 55       CONTINUE
          READ(LUN,*,ERR=55,END=56)
     >         (ANONE,ICOL=1,NTIMCOL-1),
     >         TIM(IFM,NTIM(IFM)),
     >         (ANONE,ICOL=1,NSCALCOL-NTIMCOL-1),
     >         SCL(IFM,NTIM(IFM))
          NTIM(IFM) = NTIM(IFM)+1
          IF(NTIM(IFM) .GT. MXS-2) 
     >    CALL ENDJOB('SBR RSCAL - Too many timings, max is ',MXS-2)
          GOTO 55

 56       CONTINUE
          CLOSE(LUN)
          NTIM(IFM) = NTIM(IFM)-1
          SCL(IFM,MXS) = fac

          IF( MODSCL(IFM) .EQ. 11) THEN             
            READ(NDAT,* ) NTIM2(IFM)
            IIT = NTIM2(IFM)
          IF(IIT .GT. MXS-2) 
     >    CALL ENDJOB('SBR RSCAL - Too many timings, max is ',MXS-2)
C            IF(IIT.GT.10) CALL ENDJOB(
C     >      'SBR RSCAL: Number of scaling factors cannot exceed',10)
C            IF(10*IFM+2*IIT-1.GT.MXD-1) CALL ENDJOB(
C     >      'SBR RSCAL: Too maniy families or too many timings',-99)
C----------- SCL2(IFM,IT)
            READ(NDAT,*) (SCL2(IFM,IT),IT=1,IIT)
C            READ(NDAT,*) (A(NOEL,10*(2*IFM-1)+IT-1),IT=1,IIT)
C----------- TIM2(IFM,IT)
            READ(NDAT,*) (TIM2(IFM,IT),IT=1,IIT)
C            READ(NDAT,*) (A(NOEL,10*(2*IFM  )+IT-1),IT=1,IIT)

            do it = 1, 10
              A(NOEL,10*(2*IFM-1)+IT-1) = SCL2(IFM,IT)
              A(NOEL,10*(2*IFM  )+IT-1) = TIM2(IFM,IT)
            enddo
            A(NOEL,10*(2*IFM  )+IIT) = NSTR

cC               if(IFM.eq.1 .or. IFM.eq.2)   then
c               if(IFM.eq.1 .or. IFM.eq.2)   then
c                  write(*,*) ' rscal '
c                  write(*,*) ' A(10*IFM / scl2) /  IFM :  ',IFM,noel
c                  write(*,*)   10*(2*IFM-1),'-',10*(2*IFM-1)+IIT-1
c                  write(*,*)  (A(NOEL,10*(2*IFM-1)+IT-1),IT=1,IIT)
c                  write(*,*) ' A(10*IFM... / tim2) : '
c                  write(*,*)   10*(2*IFM  ),'-',10*(2*IFM  )+IIT-1
c                  write(*,*)  (A(NOEL,10*(2*IFM  )+IT-1),IT=1,IIT)
c                   write(*,*) ' '
c                endif

            CALL SCALI4(
     >                  SCL2,TIM2,NTIM2,IFM)

          ENDIF

        ELSE

          CALL ENDJOB('SBR RSCAL: No such option MODSCL = '
     >    ,MODSCL(IFM))

        ENDIF

 1    CONTINUE
 
      close(lun)
      CALL SCALI6(MODSCL)

      RETURN

 96   CONTINUE
      WRITE(ABS(NRES),*) 'ERROR  OPEN  FILE ', 
     > TA(NOEL,IFM)(DEBSTR(TA(NOEL,IFM)):FINSTR(TA(NOEL,IFM)))
      WRITE(ABS(NRES),*) ' NOEL, IFM : ',NOEL,IFM
      CALL ENDJOB('SBR rscal. End job upon ERROR  OPEN  FILE.',-99)
      RETURN

 97   CONTINUE
      WRITE(ABS(NRES),*) 'SBR rscal. No idle unit  '
      CALL ENDJOB('SBR rscal. End job upon NO IDLE UNIT.',-99)
      RETURN
      END
