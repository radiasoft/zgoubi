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
C  François Méot <meot@lpsc.in2p3.fr>
C  Service Accélerateurs
C  LPSC Grenoble
C  53 Avenue des Martyrs
C  38026 Grenoble Cedex
C  France
      SUBROUTINE SPNPRT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),IMAX,IEX(MXT),IREP(MXT),AMQLU
      CHARACTER LET
      COMMON/FAISCT/ LET(MXT)
      COMMON/PTICUL/ AM,Q,G,TO
      COMMON/RIGID/ BORO,DPREF,DP,BR
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)
 
      JDMAX=IDMAX
      JMAXT=IMAXT
      IF(JDMAX .GT. 1) WRITE(NRES,121)JDMAX
 121  FORMAT(/,25X,' .... ',I3
     >,'  GROUPS  OF  MOMENTA  FOLLOW    ....')
 
      DO 3 ID=1,JDMAX
        IMAX1=1+(ID-1)*JMAXT
        IMAX2=IMAX1+JMAXT-1
        SX = 0D0
        SY = 0D0
        SZ = 0D0
        SXF = 0D0
        SYF = 0D0
        SZF = 0D0
        II=0
        DO 1 I=IMAX1,IMAX2
          IF( IEX(I) .LT. -1 ) GOTO 1
          II=II+1
          SX = SX + SI(1,I)
          SY = SY + SI(2,I)
          SZ = SZ + SI(3,I)
          SXF = SXF + SF(1,I)
          SYF = SYF + SF(2,I)
          SZF = SZF + SF(3,I)
 1      CONTINUE
 
        IF(NRES.GT.0) THEN
          SM = SQRT(SX*SX+SY*SY+SZ*SZ)/II
          SMF = SQRT(SXF*SXF+SYF*SYF+SZF*SZF)/II
          WRITE(NRES,120) II,SX/II,SY/II,SZ/II,SM
     >    ,SXF/II,SYF/II,SZF/II,SMF
 120      FORMAT(//,25X,' POLARISATION  MOYENNE  DU'
     >    ,2X,'FAISCEAU  DE  ',I3,'  PARTICULES :'
     >    ,//,T20,'INITIALE',T70,'FINALE'
     >    ,//,T12,'<SX>',T22,'<SY>',T32,'<SZ>',T42,'<S>'
     >    ,T61,'<SX>',T71,'<SY>',T81,'<SZ>',T91,'<S>'
     >    ,/,5X,4F10.4,10X,4F10.4)
 
          WRITE(NRES,110) JMAXT
 110      FORMAT(///,15X,' SPIN  COMPONENTS  OF  EACH  OF  THE '
     >    ,I5,'  PARTICLES :',//,T20,'INITIAL',T70,'FINAL'
     >    ,//,T12,'SX',T22,'SY',T32,'SZ',T42,'S'
     >    ,T61,'SX',T71,'SY',T81,'SZ',T91,'S',T99,'GAMMA',/)
 
          DO 10 I=IMAX1,IMAX2
            IF( IEX(I) .LT. -1 ) GOTO 10
C            P = BORO*CL*1D-9 *F(1,I) *(Q/QE)
C            P = BORO*CL*1D-9 *F(1,I) *Q
            P = BORO*CL9 *F(1,I) *Q
            GAMA = SQRT(P*P + AM*AM)/AM
            WRITE(NRES,101) LET(I),IEX(I),(SI(J,I),J=1,4)
     X      ,(SF(J,I),J=1,4),GAMA,I
 101        FORMAT(1X,A1,1X,I2,4F10.4,10X,4F10.4,F10.6,I4)
 10       CONTINUE
 
        ENDIF
 3    CONTINUE
 
      RETURN
      END
