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
C  USA
C  -------
      SUBROUTINE AVORPR(LUN,IOPT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXPUD=9,MXPU=5000)
      COMMON/CO/ FPU(MXPUD,MXPU),KCO,NPU,NFPU,IPU
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM

      DIMENSION CO1(7), CO2(7), COMA(7), COMI(7)
      
      IF    (IOPT .EQ. 1) THEN

        IF(IPASS .EQ. 1) WRITE(LUN,100) 
 100    FORMAT(/,1X,'Brute record at pick-ups :',/,
     >  2X,'PU#',T15,'Pos',T30,'Yco',T45,'Tco',T60,'Zco',T75,
     >  'Pco',T90,'L',T105,'D',T122,'t',
     >  T132,' #part ',T142,'Pass#'
     >  /,T15,'(cm)',T29,'(cm)',T44,'(mrad)',T59,'(cm)',T73,'(mrad)',
     >  T88,'(cm)',T103,'(dp/p)',T120,'(mu_s)'
     >  )

        CALL PCKUP1

      ELSEIF(IOPT .EQ. 2) THEN

        WRITE(LUN,100) 
        DO 251 J = 1, 7
          COMA(J) = -1.D10
          COMI(J) = 1.D10
          CO1(J) = 0.D0
          CO2(J) = 0.D0
          DO 261 I = 1, IPU
            DU = FPU(J,I) / IPASS
            CO1(J) = CO1(J) + DU
            DU2 = DU * DU
            CO2(J) = CO2(J) + DU2
C            IF( DU2 .GT. COMA(J)*COMA(J) ) COMA(J) = DU
            IF( DU .GT. COMA(J) ) COMA(J) = DU
            IF( DU .LT. COMI(J) ) COMI(J) = DU
 261      CONTINUE
          CO1(J) = CO1(J) / IPU
          IF(IPU.GT.1) CO2(J) = SQRT( CO2(J) / IPU - CO1(J) * CO1(J) )
 251    CONTINUE
        DO I = 1, IPU
          NT = NINT(FPU(8,I))
          WRITE(LUN,FMT= '(1P,2X,I4,6(1X,E14.6),2(1X,E17.9),2(1X,I6))')
     >    I,FPU(9,I),(FPU(J,I)/NT,J=2,6),FPU(1,I)/NT,FPU(7,I)/NT
     >    ,NT,IPASS
        ENDDO

      ENDIF

C Careful before changing output format : 
C - some post-processor use it
C - they may need high prec.!
        WRITE(LUN,FMT= '(/,1X,''PU_average (over partcl and pass) : '',
     >  /,T12,''Y'',T28,''T'',T45,''Z'',T60,''P'',
     >  T74,''path-L'',T90,''D'',T105,''time'', 
     >  /,2X,1P,7E16.8,/)' ) (CO1(J),J=2,6 ), CO1(1), CO1(7)
C     >  /,7X,1P,2E16.8,4E12.4,E16.8)' ) (CO1(J),J=2,6 ), CO1(1), CO1(7)
        IF(IPU.GT.1) WRITE(LUN,FMT= '(  1X,''Sigma  '',
     >  /,2X,1P,2E16.8,4E14.6,E16.8)' )
     >  (CO2(J),J=2,6 ), CO2(1), CO2(7)
        WRITE(LUN,FMT= '(  1X,''Max-PUSignal  '',
     >  /,2X,1P,2E16.8,4E14.6,E16.8)' )
     >  (COMA(J),J=2,6 ), COMA(1), COMA(7)
        WRITE(LUN,FMT= '(  1X,''Min-PUSignal  '',
     >  /,2X,1P,2E16.8,4E14.6,E16.8)' )
     >  (COMI(J),J=2,6 ), COMI(1), COMI(7)

      RETURN
      END
