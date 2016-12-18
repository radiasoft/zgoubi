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
      SUBROUTINE SPNPRT(LBL1, LBL2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER(*) LBL1, LBL2
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      CHARACTER(1) LET
      INCLUDE "C.FAISCT.H"     ! COMMON/FAISCT/ LET(MXT)
      INCLUDE "C.OBJET.H"     ! COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      INCLUDE "C.PTICUL.H"     ! COMMON/PTICUL/ AM,Q,G,TO
      INCLUDE "C.REBELO.H"   ! COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      INCLUDE "C.SPIN.H"     ! COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)

C      DIMENSION SMI(4,MXT), SMA(4,MXT)
      LOGICAL IDLUNI

      DIMENSION SPMI(4,MXT), SPMA(4,MXT)
      PARAMETER (ICMXT=4*MXT)

      DIMENSION AA(3),BB(3)  !,XX(3)
      DIMENSION PHI(MXT), PHIX(MXT)
      SAVE SXMF, SYMF, SZMF

      INTEGER DEBSTR, FINSTR
      LOGICAL FIRST
      SAVE LUNPRT, FIRST

      DIMENSION SMAT(3,3), DLT(3,3), PROD(3,3)

      DATA SXMF, SYMF, SZMF /  3 * 0.D0 /
      DATA SPMI, SPMA / ICMXT*1D10, ICMXT* -1D10 /
      DATA FIRST / .TRUE. /

      JDMAX=IDMAX
      JMAXT=IMAX/IDMAX

      IF(NRES.GT.0) THEN
        IF(JDMAX .GT. 1) WRITE(NRES,121) JDMAX
 121  FORMAT(/,25X,' .... ',I3,'  Groups  of  momenta  follow   ....')
      ENDIF
 
      DO 3 ID=1,JDMAX
        IMAX1=1+(ID-1)*JMAXT
        IMAX2=IMAX1+JMAXT-1
        SX =   0D0;  SY =   0D0;  SZ =   0D0
        SXM =  0D0;  SYM =  0D0;  SZM =  0D0
        SXF =  0D0;  SYF =  0D0;  SZF =  0D0 
        SXMF = 0D0;  SYMF = 0D0;  SZMF = 0D0
        PHIM = 0.D0 ; GGM = 0.D0

        II=0
        DO I=IMAX1,IMAX2
         IF( IEX(I) .GT. 0 ) THEN
          II=II+1
          SX = SX + SI(1,I)
          SY = SY + SI(2,I)
          SZ = SZ + SI(3,I)
          SXF = SXF + SF(1,I)
          SYF = SYF + SF(2,I)
          SZF = SZF + SF(3,I)
          SXMF = SXMF + SXF
          SYMF = SYMF + SYF
          SZMF = SZMF + SZF
          P = BORO*CL9 *F(1,I) *Q
          GAMA = SQRT(P*P + AM*AM)/AM
          GGM = GGM + G * GAMA 

          IF(SF(1,I).LT.SPMI(1,I)) SPMI(1,I) = SF(1,I)          
          IF(SF(2,I).LT.SPMI(2,I)) SPMI(2,I) = SF(2,I)          
          IF(SF(3,I).LT.SPMI(3,I)) SPMI(3,I) = SF(3,I)          
          IF(SF(4,I).LT.SPMI(4,I)) SPMI(4,I) = SF(4,I)          
          IF(SF(1,I).GT.SPMA(1,I)) SPMA(1,I) = SF(1,I)          
          IF(SF(2,I).GT.SPMA(2,I)) SPMA(2,I) = SF(2,I)          
          IF(SF(3,I).GT.SPMA(3,I)) SPMA(3,I) = SF(3,I)          
          IF(SF(4,I).GT.SPMA(4,I)) SPMA(4,I) = SF(4,I)          

              AA(1) = SI(1,I)
              AA(2) = SI(2,I)
              AA(3) = SI(3,I)
              BB(1) = SF(1,I)
              BB(2) = SF(2,I)
              BB(3) = SF(3,I)
              APSCAL = ACOS(VSCAL(AA,BB,3))

              PHIM = PHIM + APSCAL

         ENDIF
        ENDDO

        PHIM = PHIM / DBLE(II)
        GGM = GGM / DBLE(II)
 
        PHIM2 = 0.D0
        II=0
        DO I=IMAX1,IMAX2
         IF( IEX(I) .GT. 0 ) THEN
          II=II+1
              AA(1) = SI(1,I)
              AA(2) = SI(2,I)
              AA(3) = SI(3,I)
              BB(1) = SF(1,I)
              BB(2) = SF(2,I)
              BB(3) = SF(3,I)
              APSCAL = ACOS(VSCAL(AA,BB,3))

              PHIM2 = PHIM2 + (APSCAL -PHIM)**2

         ENDIF
        ENDDO
 
        PHIM2 = PHIM2 / DBLE(II)
        SIGPHI = SQRT(PHIM2) 
        

        IF(NRES.GT.0) THEN
          SM = SQRT(SX*SX+SY*SY+SZ*SZ)/DBLE(II)
          SMF = SQRT(SXF*SXF+SYF*SYF+SZF*SZF)/DBLE(II)
          PHIM = PHIM * DEG
          SIGPHI = SIGPHI * DEG
          WRITE(NRES,120) II,SX/DBLE(II),SY/DBLE(II),SZ/DBLE(II),SM
     >    ,SXF/DBLE(II),SYF/DBLE(II),SZF/DBLE(II),SMF,GGM,phim,sigphi
 120      FORMAT(//,25X,' Average  over  particles at this pass ; '
     >    ,2X,'beam with  ',I6,'  particles :'
     >    ,//,T20,'INITIAL',T70,'FINAL'
     >    ,//,T7,'<SX>',T18,'<SY>',T29,'<SZ>',T40,'<|S|>'
     >    ,T57,'<SX>',T68,'<SY>',T79,'<SZ>',T89,'<|S|>',T98,'<G.gma>'
     >    ,T109,'<(SI,SF)>',T120,'sigma_(SI,SF)'
     >    ,/,T110,'  (deg)',T121,'   (deg)'
     >    ,/,4(1x,F10.6),6X,7(1x,F10.6))
 
          WRITE(NRES,110) JMAXT
 110      FORMAT(//,15X,' Spin  components  of  each  of  the '
     >    ,I6,'  particles,  and  rotation  angle :'
     >    ,//,T20,'INITIAL',T70,'FINAL'
     >    ,//,T12,'SX',T22,'SY',T32,'SZ',T42,'|S|'
     >    ,T60,'SX',T70,'SY',T80,'SZ',T90,'|S|',T101,'GAMMA'
     >    ,T110,'(Si,Sf)',T120,'(Si,Sf_x)')
          WRITE(NRES,FMT='(
     >    T110,'' (deg.)'',T120,''  (deg.)'')')
          WRITE(NRES,fmt='(T92,a,/)') 
     >           '(Sf_x : projection of Sf on plane x=0)'

          DO I=IMAX1,IMAX2
            IF( IEX(I) .GE. -1 ) THEN
              P = BORO*CL9 *F(1,I) *Q
              GAMA = SQRT(P*P + AM*AM)/AM
              AA(1) = SI(1,I)
              AA(2) = SI(2,I)
              AA(3) = SI(3,I)
              BB(1) = SF(1,I)
              BB(2) = SF(2,I)
              BB(3) = SF(3,I)

              CPHI = VSCAL(AA,BB,3)
              PHI(I) = ACOS(CPHI) * DEG

C If initial spin is // Z
              PHIZF = ATAN(SQRT(SF(1,I)**2+SF(2,I)**2)/SF(3,I)) * DEG
              BB(1)= 0.D0

              CPHIX = VSCAL(AA,BB,3)/XNORM(BB,3)
              PHIX(I) = ACOS(CPHIX) * DEG 

              WRITE(NRES,101) LET(I),IEX(I),(SI(J,I),J=1,4)
     >        ,(SF(J,I),J=1,4),GAMA,PHI(I),PHIX(I),I
 101          FORMAT(1X,A1,1X,I2,4(1X,F9.6),9X,4(1X,F9.6),1X,F11.4,
     >        2(1X,F9.4),1X,I4)
C              WRITE(NRES,*)'ATN(sy/sx)=',ATAN(SF(2,I)/SF(1,I))*DEG,'deg'
            ENDIF
          ENDDO
 
          CALL SPNTR4(IMAX,SPMI,SPMA)

          WRITE(NRES,130) JMAXT
 130      FORMAT(///,15X,' Min/Max  components  of  each  of  the '
     >    ,I6,'  particles :'
     >    ,//,T3,'SX_mi',T15,'SX_ma',T27,'SY_mi',T39,'SY_ma'
     >    ,T51,'SZ_mi',T63,'SZ_ma',T75,'|S|_mi',T87,'|S|_ma'
     >    ,T99,'p/p_0',T112,'GAMMA',T127,'I  IEX',/)
          DO I=IMAX1,IMAX2
            IF( IEX(I) .GE. -1 ) THEN
              P = BORO*CL9 *F(1,I) *Q
              GAMA = SQRT(P*P + AM*AM)/AM
              WRITE(NRES,131) (SPMI(J,I),SPMA(J,I),J=1,4),F(1,I)
     >           ,GAMA,I,IEX(I)
 131          FORMAT(1P,8E12.4,2E13.5,1X,I6,1X,I3)
            ENDIF
          ENDDO
 
        ENDIF
 3    CONTINUE
 
C Matrix computation
C Pattern in OBJET has to be 
C  3 identical particles on orbit with their spin components resp. 0,0,1, 0,1,0, 1,0,0
      IF  (LBL1(DEBSTR(LBL1):FINSTR(LBL1)) .EQ. 'MATRIX'
     >.OR. LBL2(DEBSTR(LBL2):FINSTR(LBL2)) .EQ. 'MATRIX') THEN

        IT = 1
        SMAT(1,IT) = SF(1,IT) 
        SMAT(2,IT) = SF(2,IT) 
        SMAT(3,IT) = SF(3,IT) 
        IT = 2
        SMAT(1,IT) = SF(1,IT) 
        SMAT(2,IT) = SF(2,IT) 
        SMAT(3,IT) = SF(3,IT) 
        IT = 3
        SMAT(1,IT) = SF(1,IT) 
        SMAT(2,IT) = SF(2,IT) 
        SMAT(3,IT) = SF(3,IT) 

        TRM = SMAT(1,1) + SMAT(2,2) + SMAT(3,3) 
        SROT = ACOS((TRM-1.D0)/2.D0)

        CALL RAZ(DLT,3*3)
        DLT(2,3) = -1.D0 
        DLT(3,2) = +1.D0 
        PROD = MATMUL(DLT,SMAT)
        TR1 = (PROD(1,1) + PROD(2,3) + PROD(3,3))/(2.D0*SROT)
        DLT(2,3) = 0.D0 
        DLT(3,2) = 0.D0 
        DLT(1,3) = +1.D0 
        DLT(3,1) = -1.D0 
        PROD = MATMUL(DLT,SMAT)
        TR2 = (PROD(1,1) + PROD(2,3) + PROD(3,3))/(2.D0*SROT)
        DLT(1,3) = 0.D0 
        DLT(3,1) = 0.D0 
        DLT(1,2) = -1.D0 
        DLT(2,1) = +1.D0 
        PROD = MATMUL(DLT,SMAT)
        TR2 = (PROD(1,1) + PROD(2,3) + PROD(3,3))/(2.D0*SROT)
        REN = SQRT(TR1*TR1+TR2*TR2+TR3*TR3)
        TR1 = TR1 / REN ; TR2 = TR2 / REN ; TR3 = TR3 / REN 

        IF(NRES.GT.0) THEN
          WRITE(NRES,103) 
 103      FORMAT(//,18X,'Spin  transfer  matrix  ',/)
          WRITE(NRES,104) (( SMAT(IA,IB) , IB=1,3) , IA=1,3)
 104      FORMAT(6X,1P,3G16.6)
          WRITE(NRES,112) TRM, SROT*180.D0/(4.D0*ATAN(1.D0)),TR1,TR2,TR3
112       FORMAT(/,10X,'Trace = ',F18.10,',',4X,
     >    ';   spin precession acos((trace-1)/2) = ',F18.10,' deg',
     >    /,10X,'Rotation axis :   (',F7.4,', ',F7.4,', ',F7.4,')')
        ENDIF
      ENDIF

C Print to zgoubi.SPNPRT.Out
      IF  (LBL1(DEBSTR(LBL1):FINSTR(LBL1)) .EQ. 'PRINT'
     >.OR. LBL2(DEBSTR(LBL2):FINSTR(LBL2)) .EQ. 'PRINT') THEN

        IF(FIRST) THEN
          FIRST = .FALSE.
          IF(IDLUNI(
     >              LUNPRT)) THEN
            OPEN(UNIT=LUNPRT,FILE='zgoubi.SPNPRT.Out',ERR=96)
          ELSE
            GOTO 96
          ENDIF
          WRITE(LUNPRT,FMT='(3(A,/),A)') 
     >    '# spin data. PRINT by spnprt.f ',
     >    '#  1  2  3  4  5  6    7    8     9 10 11 12   '//
     >    '    13 14 15 16      17       18   19 '//
     >    '   20     21     22  23  24  25  26  27     '//
     >    ' etc.                                               '//
     >    '           ',
     >    '# Y, T, Z, P, S, D, TAG, IEX, (SI(J,I),J=1,4), '//
     >    '(SF(J,I),J=1,4), gamma, G.gamma, PHI, '//
     >    'PHIX, ITRAJ, IPASS, Yo, To, Zo, Po, So, Do '//
     >    ' !spnprt.f ',
     >    '# cm mr cm mr cm  -   -     -    - - - -       '//
     >    '    - - - -            -        -     '//
     >    '  etc.     '
        ENDIF
        DO I=IMAX1,IMAX2
          IF( IEX(I) .GE. -1 ) THEN
            P = BORO*CL9 *F(1,I) *Q
            GAMA = SQRT(P*P + AM*AM)/AM
            WRITE(LUNPRT,111) 
     >      (F(J,I),J=2,6),F(1,I)
     >      ,'''',LET(I),'''',IEX(I),(SI(J,I),J=1,4)
     >      ,(SF(J,I),J=1,4),GAMA,G*GAMA,PHI(I),PHIX(I),I,IPASS,NOEL
     >      ,(FO(J,I),J=2,6),FO(1,I)
 111        FORMAT(1X,1P,6(1X,E14.6),1X,3A1,1X,I2,12(1X,E14.6),3(1X,I6)
     >      ,6(1X,E14.6))
          ENDIF
        ENDDO
C Leaving unclosed allows stacking when combined use of FIT and REBELOTE
C        CLOSE(LUNPRT)
      ENDIF 

      CALL FLUSH2(LUNPRT,.FALSE.)

      RETURN

 96   CONTINUE
      WRITE(*        ,FMT='(/,''SBR SPNPRT : '',
     >           ''Error open file zgoubi.SPNPRT.Out'')')
      CALL ENDJOB('Pgm spnprt. Error open file zgoubi.SPNPRT.Out',-99)

      RETURN
      END
