C  ZGOUBI, a program for computing the trajectories of charged particles
C  in electric and magnetic fields
C  Copyright (C) 1988-2007  Fran�ois M�ot
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
C  Fran�ois M�ot <fmeot@bnl.gov>
C  Brookhaven National Laboratory    
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  -------
      SUBROUTINE TWorbt
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      PARAMETER (LBLSIZ=10)
      CHARACTER(LBLSIZ) LABEL
      COMMON /LABEL/ LABEL(MXL,2)
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      COMMON/PTICUL/ AM,Q,G,TO
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      
      DIMENSION RREF(6,6), T(6,6,6)
      DIMENSION T3(5,6) , T4(5,6)

      logical readat, endfit

      data prec /  1e-4 /


      xame = 1.d99
      dowhile(xame .gt. prec) 

          READAT = .FALSE.  ! OBJET has to be skipped (there may be another way).
                            ! Here, FO(j,1) has to be set to the new final coordinates
                            ! and the sampling of OBJET/KOBJ=5 has to be updated.
          ENDFIT = .FALSE.
          CALL ZGOUBI(1,MXL,READAT,
     >                             NBEL,ENDFIT)

        IF(KOBJ .EQ. 5) THEN
          IORD=1
        ELSEIF(KOBJ .EQ. 6) THEN
          IORD=2
        ELSE
          call endjob('SBR tworbt. Cannot run if kobj .ne. 5 or 6',-99)
        ENDIF 

        IF    (IORD .EQ. 1) THEN
          CALL REFER(1,1,0,1,4,5)
          CALL MAT1(1,
     >                RREF,T)
          CALL REFER(2,1,0,1,4,5)
        ELSEIF(IORD .EQ. 2) THEN
          CALL REFER(1,2,0,1,6,7)
          CALL MAT2(
     >              RREF,T,TX3,TX4)
          CALL REFER(2,2,0,1,6,7)
        ENDIF

        xame = (
     >  abs(f(1,2) - fo(1,2)) + 
     >  abs(f(1,3) - fo(1,3)) + 
     >  abs(f(1,4) - fo(1,4)) + 
     >  abs(f(1,5) - fo(1,5)) 
     >  )

      enddo
      
      ipass = 1

      RETURN

      END