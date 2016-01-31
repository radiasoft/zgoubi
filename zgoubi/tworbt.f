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
      SUBROUTINE TWorbt
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      PARAMETER (LBLSIZ=10)
      CHARACTER(LBLSIZ) LABEL
      INCLUDE "C.LABEL.H"     ! COMMON/LABEL/ LABEL(MXL,2)
      INCLUDE "C.OBJET.H"     ! COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      INCLUDE "C.PTICUL.H"     ! COMMON/PTICUL/ AM,Q,G,TO
      INCLUDE "C.REBELO.H"   ! COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      
      DIMENSION RREF(6,6), T(6,6,6)
C      DIMENSION T3(5,6) , T4(5,6)

      LOGICAL READAT

      DATA PREC /  1E-4 /


      XAME = 1.D99
      DOWHILE(XAME .GT. PREC) 

          READAT = .FALSE.  ! OBJET has to be skipped (there may be another way).
                            ! Here, FO(j,1) has to be set to the new final coordinates
                            ! and the sampling of OBJET/KOBJ=5 has to be updated.
          CALL ZGOUBI(1,MXL,READAT,
     >                             NBEL)

        IF(KOBJ .EQ. 5) THEN
          IORD=1
        ELSEIF(KOBJ .EQ. 6) THEN
          IORD=2
        ELSE
          call endjob('SBR tworbt. Cannot run if kobj .ne. 5 or 6',-99)
        ENDIF 

        IF    (IORD .EQ. 1) THEN
          CALL REFER(1,1,0,1,4,5)
          CALL MAT1(1,F,IMAX,
     >                RREF,T)
          CALL REFER(2,1,0,1,4,5)
        ELSEIF(IORD .EQ. 2) THEN
          CALL REFER(1,2,0,1,6,7)
          CALL MAT2(
     >              RREF,T,TX3,TX4)
          CALL REFER(2,2,0,1,6,7)
        ENDIF

        XAME = (
     >  ABS(F(1,2) - FO(1,2)) + 
     >  ABS(F(1,3) - FO(1,3)) + 
     >  ABS(F(1,4) - FO(1,4)) + 
     >  ABS(F(1,5) - FO(1,5)) 
     >  )

      ENDDO
      
      IPASS = 1

      RETURN

      END
