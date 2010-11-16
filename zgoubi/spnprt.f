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
C  Brookhaven National Laboratory               és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE SPNPRT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      CHARACTER LET
      COMMON/FAISCT/ LET(MXT)
      COMMON/PTICUL/ AM,Q,G,TO
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)

C      DIMENSION SMI(4,MXT), SMA(4,MXT)
      LOGICAL IDLUNI

      DIMENSION SPMI(4,MXT), SPMA(4,MXT)
      PARAMETER (ICMXT=4*MXT)
      DATA SPMI, SPMA / ICMXT*1D10, ICMXT* -1D10 /

      JDMAX=IDMAX
      JMAXT=IMAX/IDMAX
      IF(JDMAX .GT. 1) WRITE(NRES,121) JDMAX
 121  FORMAT(/,25X,' .... ',I3
     >,'  Groups  of  momenta  follow    ....')
 
      DO 3 ID=1,JDMAX
        IMAX1=1+(ID-1)*JMAXT
        IMAX2=IMAX1+JMAXT-1
C             write(abs(nres),*) ' IMAX2,IMAX1,JMAXT  ',IMAX2,IMAX1,JMAXT

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

          IF(SF(1,I).LT.SPMI(1,I)) SPMI(1,I) = SF(1,I)          
          IF(SF(2,I).LT.SPMI(2,I)) SPMI(2,I) = SF(2,I)          
          IF(SF(3,I).LT.SPMI(3,I)) SPMI(3,I) = SF(3,I)          
          IF(SF(4,I).LT.SPMI(4,I)) SPMI(4,I) = SF(4,I)          
          IF(SF(1,I).GT.SPMA(1,I)) SPMA(1,I) = SF(1,I)          
          IF(SF(2,I).GT.SPMA(2,I)) SPMA(2,I) = SF(2,I)          
          IF(SF(3,I).GT.SPMA(3,I)) SPMA(3,I) = SF(3,I)          
          IF(SF(4,I).GT.SPMA(4,I)) SPMA(4,I) = SF(4,I)          
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
     >    ,//,T15,'SX',T25,'SY',T35,'SZ',T45,'S'
     >    ,T67,'SX',T77,'SY',T87,'SZ',T97,'S',T106,'GAMMA',/)
          DO I=IMAX1,IMAX2
            IF( IEX(I) .GE. -1 ) THEN
              P = BORO*CL9 *F(1,I) *Q
              GAMA = SQRT(P*P + AM*AM)/AM
              WRITE(NRES,101) LET(I),IEX(I),(SI(J,I),J=1,4)
     X        ,(SF(J,I),J=1,4),GAMA,I
 101          FORMAT(1X,A1,1X,I2,4(1X,F10.6),9X,4(1X,F10.6),
     >                                         1X,F12.6,1X,I4)
C              WRITE(NRES,*)'ATN(sy/sx)=',ATAN(SF(2,I)/SF(1,I))*DEG,'deg'
            ENDIF
          ENDDO
 
          CALL SPNTR4(IMAX,SPMI,SPMA)
C          CALL SPNTR3(IMAX,
C     >                     SMI, SMA)
          WRITE(NRES,130) JMAXT
 130      FORMAT(///,15X,' Min/Max  components  of  each  of  the '
     >    ,I5,'  particles :'
     >    ,//,T3,'SX_mi',T15,'SX_ma',T27,'SY_mi',T39,'SY_ma'
     >    ,T51,'SZ_mi',T63,'SZ_ma',T75,'|S|_mi',T87,'|S|_ma'
     >    ,T99,'p/p_0',T111,'GAMMA',T125,'I  IEX',/)
          DO I=IMAX1,IMAX2
            IF( IEX(I) .GE. -1 ) THEN
              P = BORO*CL9 *F(1,I) *Q
              GAMA = SQRT(P*P + AM*AM)/AM
              WRITE(NRES,131) (SPMI(J,I),SPMA(J,I),J=1,4),F(1,I)
     >           ,GAMA,I,IEX(I)
 131          FORMAT(1P,8E12.4,2E13.5,I5,I3)
            ENDIF
          ENDDO
 
          IF(IPASS.EQ.NRBLT+1) THEN 

            IF(IDLUNI(
     >                LUN)) THEN
              OPEN(UNIT=LUN,FILE='zgoubi.SPNPRT.out',ERR=96)
            ELSE
              GOTO 96
            ENDIF

            WRITE(LUN,130) JMAXT
            DO I=IMAX1,IMAX2
              IF( IEX(I) .GE. -1 ) THEN
                P = BORO*CL9 *F(1,I) *Q
                GAMA = SQRT(P*P + AM*AM)/AM
                WRITE(LUN,131) (SPMI(J,I),SPMA(J,I),J=1,4),F(1,I)
     >             ,GAMA,I,IEX(I)
              ENDIF
            ENDDO
            CLOSE(LUN)

          ENDIF 
        ENDIF
 3    CONTINUE
 
      RETURN

 96       CONTINUE
          WRITE(ABS(NRES),FMT='(/,''SBR SPNPRT : '',
     >               ''Error open file zgoubi.SPNPRT.out'')')
          WRITE(*        ,FMT='(/,''SBR SPNPRT : '',
     >               ''Error open file zgoubi.SPNPRT.out'')')

      RETURN
      END
