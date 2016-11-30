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
      SUBROUTINE RSCALE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ----------------------------------------------
C     READS DATA FOR POWER SUPPLIES OF LMNT FAMILIES
C     ----------------------------------------------
      COMMON/CDF/ IES,IORDRE,LCHA,LIST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER(80) TA
      COMMON/DONT/ TA(MXL,40)
      INCLUDE 'MXFS.H'
      COMMON/SCAL/SCL(MXF,MXS),TIM(MXF,MXS),KTI(MXF),KSCL
      PARAMETER (LBLSIZ=20)
      CHARACTER(LBLSIZ) LABEL
      COMMON /LABEL/ LABEL(MXL,2)
      PARAMETER (KSIZ=10)
      CHARACTER FAM*(KSIZ),LBF*(LBLSIZ)
      COMMON/SCALT/ FAM(MXF),LBF(MXF,MLF),JPA(MXF,MXP)
 
      PARAMETER(MSTR=MLF+1)
      CHARACTER(10) STRA(MSTR)
 
      INTEGER DEBSTR, FINSTR

C----- IOPT; NB OF DIFFRNT FAMILIES TO BE SCALED (<= MXF)
      LINE = 1
      READ(NDAT,*) A(NOEL,1),A(NOEL,2)

      NFAM = NINT(A(NOEL,2))
      IF(NFAM .GT. MXF) THEN
        WRITE(NRES,FMT=
     $ '('' To many families (maximum number allowed is '',I2,'')'')'
     $        )MXF
        STOP
      ENDIF

      DO 1 IF=1,NFAM

C------- Store name of family and label(s)
        LINE = LINE + 1
        READ(NDAT,100)  TA(NOEL,IF)
 100    FORMAT(A80)
        TA(NOEL,IF) = 
     >    TA(NOEL,IF)(DEBSTR(TA(NOEL,IF)):FINSTR(TA(NOEL,IF)))
        CALL STRGET(TA(NOEL,IF),MSTR,
     >                                 NST,STRA)
        IF(NST-1 .GT. MXLF) THEN
          WRITE(NRES,FMT='
     >  ('' To many labels per family (maximum number allowed is '',
     >         I2,'')'')') MXLF
          STOP
        ENDIF

        FAM(IF) = STRA(1)(1:KSIZ)
            
        IF(NST .GE. 2) THEN
          DO 11 KL=2,NST
 11         LBF(IF,KL-1) =  STRA(KL)(1:LBLSIZ)
        ENDIF
        DO 12 KL=NST+1, MSTR
 12       LBF(IF,KL-1) = ' '

        LINE = LINE + 1
        READ(NDAT,*) KTI(IF)
        IF(KTI(IF) .GT. MXS) THEN
          WRITE(NRES,FMT='
     >      ('' To many timings  (maximum number allowed is '',
     >                  I2,'')'')') MXS
          STOP
        ENDIF

C------- TIM(IF,IT)
        LINE = LINE + 1
        READ(NDAT,*) (A(NOEL,10*IF+IT-1),IT=1,KTI(IF))
C------- SCL(IF,IT)
        LINE = LINE + 1
        READ(NDAT,*) (A(NOEL,10*IF+KTI(IF)+IT-1),IT=1,KTI(IF))
        A(NOEL,10*IF+2*KTI(IF)) = NST

 1    CONTINUE
 
      RETURN
      END
