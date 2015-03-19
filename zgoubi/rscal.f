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
C  USA
C  -------
      SUBROUTINE RSCAL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ----------------------------------------------
C     READS DATA FOR POWER SUPPLIES OF LMNT FAMILIES
C     ----------------------------------------------
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER(80) TA
      PARAMETER (MXTA=45)
      INCLUDE "C.DONT.H"     ! COMMON/DONT/ TA(MXL,MXTA)
      INCLUDE 'MXFS.H'
      INCLUDE 'MXSCL.H'
      INCLUDE "C.SCAL.H"     ! COMMON/SCAL/ SCL(MXF,MXS,MXSCL),TIM(MXF,MXS),NTIM(MXF),KSCL
      INCLUDE "C.SCALP.H"     ! COMMON/SCALP/ VPA(MXF,MXP),JPA(MXF,MXP)
C      COMMON/SCAL/SCL(MXF,MXS),TIM(MXF,MXS),NTIM(MXF),JPA(MXF,MXP),KSCL
      PARAMETER (LBLSIZ=10)
      PARAMETER (KSIZ=10)
      CHARACTER FAM*(KSIZ),LBF*(LBLSIZ)
      INCLUDE "C.SCALT.H"     ! COMMON/SCALT/ FAM(MXF),LBF(MXF,MLF)
 
      PARAMETER (MSTR=MLF+1)
      CHARACTER(LBLSIZ) STRA(MSTR)
 
      PARAMETER (MPA = MXP-1, MSTRD=2*MXP, KSTRA=40)
      CHARACTER(KSTRA) STRAD(MSTRD)
 
      INTEGER DEBSTR, FINSTR
 
      LOGICAL STRCON, IDLUNI, OK
      CHARACTER(132) TXT132, TXT, TXTF
      DIMENSION MODSCL(MXF)
      SAVE MODSCL
 
      DIMENSION NTIM2(MXF),SCL2(MXF,MXS),TIM2(MXF,MXS),NSCALCOL(MXSCL)
 
      DATA MODSCL / MXF*0  /
      DATA FAC / 1.D0  /
      DATA LUN / 0 /
 
      IF(MXTA.LT.MXF)
     >CALL ENDJOB('SBR RSCAL. Change MXTA to same value as MXF',-99)

C----- IOPT; NB OF DIFFRNT FAMILIES TO BE SCALED (<= MXF)
      NP = 1
      READ(NDAT,*) A(NOEL,NP),NFAM
      NP = NP + 1
      A(NOEL,NP) = NFAM
 
      IF(NFAM .GT. MXF)
     >CALL ENDJOB('SBR RSCAL - Too many families. Max allowed is ',MXF)
      IF(NFAM .GT. MXTA)
     >CALL ENDJOB('SBR RSCAL - Too many TA. Max allowed is ',MXTA)
 
      DO 1 IFM=1,NFAM
 
        IF(NP .GT. MXD-2)
     >  CALL ENDJOB('SBR RSCAL - Too many data. Max allowed is ',MXD-2)
 
C------- Store name of family and label(s)
        READ(NDAT,FMT='(A)')  TXT132
C Remove possible comment trailer
        IF( STRCON(TXT132,'!',
     >                        IS)) TXT132 = TXT132(DEBSTR(TXT132):IS-1)
 
        I =   0
 
        CALL RAZS(STRA,MSTR)
        CALL STRGET(TXT132,MSTR,
     >                          NSTR,STRA)
 
        IF(NSTR .GT. MSTR)
     >  CALL ENDJOB
     >  ('SBR RSCAL - Too many labels in family. Max allowed is ',MLF)
 
        IF(KSIZ .GT. LBLSIZ) STOP ' Pgm rscal, ERR : KSIZ > LBLSIZ.'
        FAM(IFM) = STRA(1)(1:KSIZ)
        TA(NOEL,IFM) = FAM(IFM)
 
        IF(NSTR .GE. 2) THEN
          DO  KL=1,NSTR-1
            IF(KL+1 .GT. MSTR) STOP ' Pgm rscal, ERR : KL+1 > MSTR.'
            LBF(IFM,KL) =  STRA(KL+1)(1:LBLSIZ)
          ENDDO
        ENDIF
 
        DO KL=NSTR, MLF
          LBF(IFM,KL) = ' '
        ENDDO
 
C For the current family, get the number of timings or working mode and possible parameters
C (input data is of the form NT[.MODSCL].
        READ(NDAT,FMT='(A)') TXT132
C Remove possible comment trailer
        IF( STRCON(TXT132,'!',
     >                        IS)) TXT132 = TXT132(DEBSTR(TXT132):IS-1)
 
        READ(TXT132,*) TXT
 
        IF( STRCON(TXT,'.',
     >                     IS)) THEN
 
          OK = STRCON(TXT132,'.',
     >                           IS)
          READ(TXT132(DEBSTR(TXT132):IS-1),*) NTIM(IFM)
          READ(TXT132(IS+1:FINSTR(TXT132)),*) MODSCL(IFM)
          NP = NP + 1
          A(NOEL,NP) = NTIM(IFM)
          I_NTIM = NP
C Possible additional scaling factor :
          IF(MODSCL(IFM) .EQ. 10)
     >      READ(TXT132(IS+3:FINSTR(TXT132)),*,ERR=44,END=44) FAC
          GOTO 45
 
 44       CONTINUE
          FAC = 1.D0
 45       CONTINUE
 
          IF(NTIM(IFM) .GT. MXS-2)
     >       CALL ENDJOB('SBR RSCAL - Too many timings. Max is ',MXS-2)
        ELSE
 
          CALL RAZS(STRAD,MSTRD)
          CALL STRGET(TXT132,MSTRD,
     >                            KSTR,STRAD)
 
          IF(KSTR.GE.1) READ(STRAD(1),*) NTIM(IFM)
          NP = NP + 1
          A(NOEL,NP) = NTIM(IFM)
          IF(NTIM(IFM) .GT. MXD-2)
     >        CALL ENDJOB('SBR RSCAL - Too many timings. Max is ',MXD-2)
 
          IF(KSTR.GE.2) THEN
            READ(STRAD(2),*) NPA
            IF(NPA.GT.MPA) CALL
     >         ENDJOB('SBR RSCAL - Too many parameterss. Max is ',MPA)
            IF(NPA .GT. MXSCL)
     >       CALL ENDJOB('SBR rscal. Exceded size of SCL tab.',-99) 
          ELSE
            NPA = 0
          ENDIF
          JPA(IFM,MXP) = NPA
 
          DO J = 1, NPA
C            IF(2*J+1 .GT. KSTRA) STOP ' SBR rscal, ERR : 2J+1 > KSTRA.'
            IF(2*J+1 .GT. MSTRD) CALL ENDJOB(
     >      ' SBR rscal, ERR : 2J+1 > MSTRD.',-99)
            READ(STRAD(2*J+1),*) JPA(IFM,J)
            READ(STRAD(2*J+2),*) VPA(IFM,J)
            IF(NP.GT.MXD-2) CALL ENDJOB('SBR rscal. Too many data.',-99)
            NP = NP + 1
            IF(NP.GT.MXD) CALL ENDJOB('SBR rscal. NP >',MXD)
            A(NOEL,NP) = VPA(IFM,J)
 
          ENDDO
 
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
            IF(NP.GT.MXD-2) CALL ENDJOB('SBR RSCAL. Too many data.',-99)
            A(NOEL,MXD-IFM) = NTIM(IFM)
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
 
C--------- SCL(IFM,IT)
          NP = NP + 1
          IF(NP.GT.MXD-2) CALL ENDJOB('SBR RSCAL. Too many data.',-99)
          READ(NDAT,*) (A(NOEL,NP+IT-1),IT=1,NDSCL)
 
C             write(*,*) ' A(NOEL,NP+IT-1),IT=1,NDSCL : ',
C     >          (NP+IT-1,A(NOEL,NP+IT-1),IT=1,NDSCL)
 
          NP = NP + NDSCL -1
 
C--------- TIM(IFM,IT)
          NP = NP + 1
          IF(NP.GT.MXD-2) CALL ENDJOB('SBR RSCAL. Too many data.',-99)
          READ(NDAT,*) (A(NOEL,NP+IT-1),IT=1,NDTIM)
 
C             WRITE(*,*) ' A(NOEL,NP+IT-1),IT=1,NDtim : ',
C     >          (NP+IT-1,A(NOEL,NP+IT-1),IT=1,NDtim)
 
 
          NP = NP + NDTIM - 1
 
          NP = NP + 1
          IF(NP.GT.MXD-2) CALL ENDJOB('SBR RSCAL. Too many data.',-99)
          A(NOEL,NP) = NSTR
 
        ELSEIF( MODSCL(IFM) .GE. 10) then
C--------- Name of the storage file in the next line
c     yann : modif to setup the read of the cols if the external file
c     is used together with the new scaling method that point directly to the A table
          READ(NDAT,FMT='(A)') TXTF
          CALL SCALI8(TXTF, IFM)
          READ(NDAT,fmt='(a)') TXT132
          CALL RAZS(STRAD,MSTRD)
          CALL STRGET(TXT132,MSTRD,
     >         KSTR,STRAD)
 
          IF( KSTR .LE. 2 ) THEN !yann : this case is "as usual"
             READ(TXT132,*) NTIMCOL, NSCALCOL(1)
             NPA = 1
             JPA(IFM,MXP) = NPA
          ELSE                   !yann : this is the new case
             READ(TXT132,*) NPA, NTIMCOL
             IF( NPA.GT.10 ) CALL ENDJOB('SBR RSCAL - ' //
     >            'Too many parameter, max is ', 10)
             JPA(IFM,MXP) = NPA
             DO I=1, NPA         !yann : loop over the declared number of parameter
                READ(STRAD(2*I+1),*,err=95,end=95) JPA(IFM,I)     ! index in the table A
                READ(STRAD(2*I+2),*,err=95,end=95) NSCALCOL(I)    ! column number in the data file
             ENDDO
          ENDIF
c     yann : End of modif
 
          IF(IDLUNI(
     >              LUN)) THEN
            OPEN(UNIT=LUN,FILE=TXTF(DEBSTR(TXTF):FINSTR(TXTF))
     >            ,STATUS='OLD',ERR=96)
          ELSE
            GOTO 97
          ENDIF
 
          DO I=1, NPA ! yann : loop over the number of parameter, it is 1 for MOD .10 or .11
             NTIM(IFM) = 1      !number of timing is determined during lecture
             REWIND(LUN)
 55          CONTINUE
             READ(LUN,*,ERR=55,END=56)
     >            (ANONE,ICOL=1,NTIMCOL-1),
     >            TIM(IFM,NTIM(IFM)),
     >            (ANONE,ICOL=1,NSCALCOL(I)-NTIMCOL-1),
     >            SCL(IFM,NTIM(IFM),I)
             NTIM(IFM) = NTIM(IFM)+1
             IF(NTIM(IFM) .GT. MXS-2)
     >            CALL ENDJOB('SBR RSCAL - ' //
     >            'Too many timings, max is ',MXS-2)
             GOTO 55
 56          CONTINUE
          ENDDO       ! yann : end of loop
 
          CLOSE(LUN)
          NTIM(IFM) = NTIM(IFM)-1
          SCL(IFM,MXS,1) = fac
          A(NOEL,I_NTIM) = NTIM(IFM)
 
          IF( MODSCL(IFM) .EQ. 11) THEN
            READ(NDAT,* ) NTIM2(IFM)
            IIT = NTIM2(IFM)
            IF(IIT .GT. MXS-2)
     >      CALL ENDJOB('SBR RSCAL - Too many timings, max is ',MXS-2)
C            IF(IIT.GT.10) CALL ENDJOB(
C     >      'SBR RSCAL: Number of scaling factors cannot exceed',10)
C            IF(10*IFM+2*IIT-1.GT.MXD-1) CALL ENDJOB(
C     >      'SBR RSCAL: Too maniy families or too many timings',-99)
C----------- SCL2(IFM,IT)
            READ(NDAT,*) (SCL2(IFM,IT),IT=1,IIT)
            NP = NP + IIT
            IF(NP.GT.MXD-2)
     >        CALL ENDJOB('SBR RSCAL. Too many data.',-99)
            DO IT = 1, IIT
              A(NOEL,NP+IT-1) =     SCL2(IFM,IT)
            ENDDO
C----------- TIM2(IFM,IT)
            READ(NDAT,*) (TIM2(IFM,IT),IT=1,IIT)
            NP = NP + IIT
            IF(NP.GT.MXD-2)
     >        CALL ENDJOB('SBR RSCAL. Too many data.',-99)
            DO IT = 1, IIT
              A(NOEL,NP+IT-1) =  TIM2(IFM,IT)
            ENDDO
 
 
            CALL SCALI4(
     >                  SCL2,TIM2,NTIM2,IFM)
 
          ENDIF
 
        ELSE
 
          CALL ENDJOB
     >    ('SBR RSCAL: No such option MODSCL = ',MODSCL(IFM))
 
        ENDIF
 
 1    CONTINUE
 
      CLOSE(lun)
      CALL SCALI6(MODSCL)
 
      RETURN
 
 95   CONTINUE
      WRITE(ABS(NRES),*) 'ERROR READING SCALING'
      WRITE(ABS(NRES),*) '   Family number', IFM
      WRITE(ABS(NRES),*) '   at element : ', FAM(IFM)
      CALL ENDJOB('SBR rscal. End job upon reading scaling',-99)
      RETURN
 
 96   CONTINUE
      WRITE(ABS(NRES),*) 'ERROR  OPEN  FILE ',
     > TXTF(DEBSTR(TXTF):FINSTR(TXTF))
      WRITE(ABS(NRES),*) ' NOEL, IFM : ',NOEL,IFM
      CALL ENDJOB('SBR rscal. End job upon ERROR  OPEN  FILE.',-99)
      RETURN
 
 97   CONTINUE
      WRITE(ABS(NRES),*) 'SBR rscal. No idle unit  '
      CALL ENDJOB('SBR rscal. End job upon NO IDLE UNIT.',-99)
      RETURN
      END
