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
      SUBROUTINE SPNPRT(LBL1, LBL2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (LBLSIZ=20)
      CHARACTER(LBLSIZ) LBL1, LBL2
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL,PI,RAD,DEG,QEL,AMPROT,CM2M
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      CHARACTER(1) LET
      INCLUDE "C.FAISCT.H"     ! COMMON/FAISCT/ LET(MXT)
      INCLUDE "C.OBJET.H"     ! COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT,KZOB
      INCLUDE "C.PTICUL.H"     ! COMMON/PTICUL/ AMASS,Q,G,TO
      INCLUDE "C.REBELO.H"   ! COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,HDPRF,DP,QBR,BRI
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

      DIMENSION SMAT(3,3)
      LOGICAL FITIN

      DATA SXMF, SYMF, SZMF /  3 * 0.D0 /
      DATA SPMI, SPMA / ICMXT*1D10, ICMXT* -1D10 /
      DATA FIRST / .TRUE. /
      DATA TR1, TR2, TR3 / 3 * 0.D0 /

      JDMAX=IDMAX
      JMAXT=IMAX/IDMAX

      IF(NRES.GT.0) THEN
        IF(JDMAX .GT. 1) WRITE(NRES,121) JDMAX
 121  FORMAT(/,25X,' -- ',I0,'  GROUPS  OF  MOMENTA  FOLLOW   --')
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
          GAMA = SQRT(P*P + AMASS*AMASS)/AMASS
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
          WRITE(NRES,120) ID,II,SX/DBLE(II),SY/DBLE(II),SZ/DBLE(II),SM
     >    ,SXF/DBLE(II),SYF/DBLE(II),SZF/DBLE(II),SMF,GGM,phim,sigphi
 120      FORMAT(//,25X,
     >    'Momentum  group  #',I0,' ; '
     >    ,'average  over ',I0,' particles at this pass : '
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
     >    ,T54,'SX',T64,'SY',T74,'SZ',T84,'|S|',T94,'GAMMA'
     >    ,T102,'(Si,Sf)',T112,'(Z,Sf_X)',T122,'(Z,Sf)')
          WRITE(NRES,FMT='(
     >    T102,'' (deg.)'',T112,'' (deg.)'',T122,'' (deg.)'')')
          WRITE(NRES,fmt='(T92,a,/)')
     >           '(Sf_X : projection of Sf on YZ plane)'

          DO I=IMAX1,IMAX2
            IF( IEX(I) .GE. -1 ) THEN
              P = BORO*CL9 *F(1,I) *Q
              GAMA = SQRT(P*P + AMASS*AMASS)/AMASS
              AA(1) = SI(1,I)
              AA(2) = SI(2,I)
              AA(3) = SI(3,I)
              BB(1) = SF(1,I)
              BB(2) = SF(2,I)
              BB(3) = SF(3,I)

C Angle between SI and SF
              CPHI = VSCAL(AA,BB,3) /XNORM(AA,3) /XNORM(BB,3)
              IF(ABS(CPHI) .GT.  1.D0) THEN
                PHI(I) = 0.D0
              ELSE
                PHI(I) = SIGN(ACOS(CPHI) * DEG,CPHI)
              ENDIF
C Angle between SF and Z axis
              PHIZF = ATAN2(SQRT(BB(1)**2+BB(2)**2),BB(3)) * DEG

C Angle between 'projection of SF on YZ plane' and Z axis
              PHIX(I) = ATAN2(BB(2), BB(3)) * DEG

              WRITE(NRES,101) LET(I),IEX(I),(SI(J,I),J=1,4)
C     >        ,(SF(J,I),J=1,4),CPHI,PHI(I),PHIX(I),PHIZF,I
     >        ,(SF(J,I),J=1,4),GAMA,PHI(I),PHIX(I),PHIZF,I
 101          FORMAT(1X,A1,1X,I2,4(1X,F9.6),3X,4(1X,F9.6),1X,F11.4,
     >        3(1X,F8.3),1X,I4)

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
              GAMA = SQRT(P*P + AMASS*AMASS)/AMASS
              WRITE(NRES,131) (SPMI(J,I),SPMA(J,I),J=1,4),F(1,I)
     >           ,GAMA,I,IEX(I)
 131          FORMAT(1P,8E12.4,2E13.5,1X,I6,1X,I3)
            ENDIF
          ENDDO

        ENDIF

C Matrix computation
C Pattern in OBJET has to be
C  3 identical particles on orbit with their spin components resp. 0,0,1, 0,1,0, 1,0,0
        IF  (LBL1(DEBSTR(LBL1):FINSTR(LBL1)) .EQ. 'MATRIX'
     >  .OR. LBL2(DEBSTR(LBL2):FINSTR(LBL2)) .EQ. 'MATRIX') THEN

          IF(JMAXT.NE.3) THEN
            IF(NRES .NE. 0) WRITE(ABS(NRES),FMT='(/,5X,
     >      ''Pgm spnprt. For computation of spin matrices, '',
     >      ''OBJET must create groups of 3 same-momentum particles, '',
     >      /,5X,''with respective spins in direction X, Y, Z.'',
     >      /,5X,''Hint :  use  OBJET with KOBJ = 1 or 2'')')

          ELSE

            CALL SPNMAT(ID,
     >                   SMAT,TRM, SROT,TR1,TR2,TR3,QS)

            IF(NRES.GT.0) THEN
              A2XY = ATAN2(TR3,SQRT(TR1*TR1+TR2*TR2))*DEG
              A2X = ATAN2(TR2,TR1)*DEG
              WRITE(NRES,103) ID
 103          FORMAT(//,18X,'Spin transfer matrix, momentum group # '
     >        ,I0,' :',/)
              WRITE(NRES,104) (( SMAT(IA,IB) , IB=1,3) , IA=1,3)
 104          FORMAT(6X,1P,3G16.6)
              WRITE(NRES,112)TRM,SROT*DEG,TR1,TR2,TR3,A2XY,A2X,QS
112           FORMAT(/,5X,'Trace = ',F18.10,',',4X,
     >        ';   spin precession acos((trace-1)/2) = ',F18.10,' deg',
     >        /,5X,'Rotation axis :   (',F7.4,', ',F7.4,', ',F7.4,')',
     >        '   ->   angle to (X,Y) plane,  angle to X axis : ',
     >         F10.4,', ',F10.4,'  degree',
     >        /,5X,'Spin  tune  Qs =   ',1P,E12.4)
            ENDIF
          ENDIF

      ENDIF
 3    CONTINUE

C Print to zgoubi.SPNPRT.Out
      IF  (LBL1(DEBSTR(LBL1):FINSTR(LBL1)) .EQ. 'PRINT'
     >.OR. LBL2(DEBSTR(LBL2):FINSTR(LBL2)) .EQ. 'PRINT') THEN
        CALL FITSTA(5,
     >                FITIN)

        IF(.NOT. FITIN) THEN    ! This can be past 'FIT' keyword
          IF(FIRST) THEN
            FIRST = .FALSE.
            IF(IDLUNI(
     >                LUNPRT)) THEN
              OPEN(UNIT=LUNPRT,FILE='zgoubi.SPNPRT.Out',ERR=96)
            ELSE
              GOTO 96
            ENDIF
            WRITE(LUNPRT,FMT='(3(A,/),A)')
     >      '# spin data. PRINT by spnprt.f ',
     >      '# 1  2  3  4  5  6    7    8     9 10 11 12   '//
     >      '    13 14 15 16      17       18   19 '//
     >      '   20     21     22   23   24  25  26  27  28  29 '//
     >      ' 30-32     etc.                               '//
     >      '           ',
     >      '# Y, T, Z, P, S, D, TAG, IEX, (SI(J,I),J=1,4), '//
     >      '(SF(J,I),J=1,4), gamma, G.gamma, PHI, '//
     >      'PHIX, ITRAJ, IPASS, NOEL, Yo, To, Zo, Po, So, Do, '//
     >      'AXE(1,3) Qs !spnprt.f ',
     >      '# cm mr cm mr cm  -   -     -    - - - -       '//
     >      '    - - - -            -         Spin angles in deg. '//
     >      '                 -         -           -          '//
     >      '     Qs !spnprt.f  -     '
          ENDIF
          DO I=IMAX1,IMAX2
            IF( IEX(I) .GE. -1 ) THEN
              P = BORO*CL9 *F(1,I) *Q
              GAMA = SQRT(P*P + AMASS*AMASS)/AMASS
              WRITE(LUNPRT,111)
     >        (F(J,I),J=2,6),F(1,I)
     >        ,'''',LET(I),'''',IEX(I),(SI(J,I),J=1,4)
     >        ,(SF(J,I),J=1,4),GAMA,G*GAMA,PHI(I),PHIX(I),I,IPASS,NOEL
     >        ,(FO(J,I),J=2,6),FO(1,I),TR1,TR2,TR3,QS
     >        ,'!spnprt.f',LBL1,LBL2
 111          FORMAT(1X,1P,6(1X,E14.6),1X,3A1,1X,I2,12(1X,E16.8),
     >        3(1X,I6),9(1X,E16.8),1X,0P,F9.6,3(1X,A))
            ENDIF
          ENDDO
C Leaving unclosed allows stacking when combined use of FIT and REBELOTE
C        CLOSE(LUNPRT)
        ENDIF
      ENDIF

      CALL FLUSH2(LUNPRT,.FALSE.)

      RETURN

 96   CONTINUE
      IF(NRES .NE. 0)
     >WRITE(ABS(NRES),FMT='(/,''SBR SPNPRT : '',
     >           ''Error open file zgoubi.SPNPRT.Out'')')
      WRITE(*        ,FMT='(/,''SBR SPNPRT : '',
     >           ''Error open file zgoubi.SPNPRT.Out'')')
      CALL ENDJOB('Pgm spnprt. Error open file zgoubi.SPNPRT.Out',-99)

      RETURN
      END
