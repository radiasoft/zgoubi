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
      SUBROUTINE SPNPRT(LABEL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER(*) LABEL
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      CHARACTER(1) LET
      COMMON/FAISCT/ LET(MXT)
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      COMMON/PTICUL/ AM,Q,G,TO
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)

C      DIMENSION SMI(4,MXT), SMA(4,MXT)
      LOGICAL IDLUNI

      DIMENSION SPMI(4,MXT), SPMA(4,MXT)
      PARAMETER (ICMXT=4*MXT)

      DIMENSION AA(3),BB(3)  !,XX(3)
      DIMENSION PHI(MXT), PHIX(MXT)
      SAVE SXMF, SYMF, SZMF

      INTEGER DEBSTR, FINSTR
      LOGICAL FIRST
      SAVE LUN 

      DATA SXMF, SYMF, SZMF /  3 * 0.D0 /
      DATA SPMI, SPMA / ICMXT*1D10, ICMXT* -1D10 /
      DATA FIRST / .TRUE. /

      JDMAX=IDMAX
      JMAXT=IMAX/IDMAX
      IF(JDMAX .GT. 1) WRITE(NRES,121) JDMAX
 121  FORMAT(/,25X,' .... ',I3
     >,'  Groups  of  momenta  follow    ....')
 
      DO 3 ID=1,JDMAX
        IMAX1=1+(ID-1)*JMAXT
        IMAX2=IMAX1+JMAXT-1
        SX =   0D0;  SY =   0D0;  SZ =   0D0
        SXM =  0D0;  SYM =  0D0;  SZM =  0D0
        SXF =  0D0;  SYF =  0D0;  SZF =  0D0 
        SXMF = 0D0;  SYMF = 0D0;  SZMF = 0D0
        phim = 0.D0
C        IF(IPASS.EQ.1) THEN
C          SXMT = 0D0
C          SYMT = 0D0
C          SZMT = 0D0
C        ENDIF

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

          IF(SF(1,I).LT.SPMI(1,I)) SPMI(1,I) = SF(1,I)          
          IF(SF(2,I).LT.SPMI(2,I)) SPMI(2,I) = SF(2,I)          
          IF(SF(3,I).LT.SPMI(3,I)) SPMI(3,I) = SF(3,I)          
          IF(SF(4,I).LT.SPMI(4,I)) SPMI(4,I) = SF(4,I)          
          IF(SF(1,I).GT.SPMA(1,I)) SPMA(1,I) = SF(1,I)          
          IF(SF(2,I).GT.SPMA(2,I)) SPMA(2,I) = SF(2,I)          
          IF(SF(3,I).GT.SPMA(3,I)) SPMA(3,I) = SF(3,I)          
          IF(SF(4,I).GT.SPMA(4,I)) SPMA(4,I) = SF(4,I)          

              aa(1) = si(1,i)
              aa(2) = si(2,i)
              aa(3) = si(3,i)
              bb(1) = sf(1,i)
              bb(2) = sf(2,i)
              bb(3) = sf(3,i)
              apscal = acos(vscal(aa,bb,3))

              phim = phim + apscal

         ENDIF
        ENDDO

        PHIM = PHIM / DBLE(II)
 
        phim2 = 0.D0
        II=0
        DO I=IMAX1,IMAX2
         IF( IEX(I) .GT. 0 ) THEN
          II=II+1
              aa(1) = si(1,i)
              aa(2) = si(2,i)
              aa(3) = si(3,i)
              bb(1) = sf(1,i)
              bb(2) = sf(2,i)
              bb(3) = sf(3,i)
              apscal = acos(vscal(aa,bb,3))

              phim2 = phim2 + (apscal -phim)**2

         ENDIF
        ENDDO
 
        PHIM2 = PHIM2 / DBLE(II)
        SIGPHI = SQRT(PHIM2) 
        

        IF(NRES.GT.0) THEN
          SM = SQRT(SX*SX+SY*SY+SZ*SZ)/DBLE(II)
          SMF = SQRT(SXF*SXF+SYF*SYF+SZF*SZF)/DBLE(II)
          PHIM = PHIM * deg
          sigphi = sigphi * deg
          WRITE(NRES,120) II,SX/DBLE(II),SY/DBLE(II),SZ/DBLE(II),SM
     >    ,SXF/DBLE(II),SYF/DBLE(II),SZF/DBLE(II),SMF,phim,sigphi
 120      FORMAT(//,25X,' Average  over  particles at this pass ; '
     >    ,2X,'beam with  ',I6,'  particles :'
     >    ,//,T20,'INITIAL',T70,'FINAL'
     >    ,//,T9,'<SX>',T21,'<SY>',T33,'<SZ>',T45,'<S>'
     >    ,T67,'<SX>',T78,'<SY>',T91,'<SZ>',T104,'<S>'
     >    ,t109,'<(SI,SF)>',t120,'sigma_(SI,SF)'
     >    ,/,t110,'  (deg)',t121,'   (deg)'
     >    ,/,4(2x,F10.6),10X,6(2x,F10.6))
 
cc          WRITE(NRES,140) II,SXF/DBLE(II),SYF/DBLE(II),SZF/DBLE(II),SMF
c          WRITE(NRES,140) II,SXmt/DBLE(II),SYmt/DBLE(II),SZmt/DBLE(II),SMF
c 140      FORMAT(//,25X,' Average  over  particles and pass, '
c     >    ,'at this pass ;  beam with  ',I6,'  particles :'
c     >    ,//,T20,'FINAL'
c     >    ,//,T13,'<SX>',T24,'<SY>',T35,'<SZ>'
c     >    ,/,5X,4(1X,F10.4))
 
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
              aa(1) = si(1,i)
              aa(2) = si(2,i)
              aa(3) = si(3,i)
              bb(1) = sf(1,i)
              bb(2) = sf(2,i)
              bb(3) = sf(3,i)

              cphi = vscal(aa,bb,3)
              PHI(I) = acos(cphi) * deg

C If initial spin is // Z
              phizf = atan(sqrt(sf(1,i)**2+sf(2,i)**2)/sf(3,i)) * deg
c                  write(*,*) ' '
c                  write(*,*) ' spnprt it, phizf : ',i,phizf
C              call vvect(aa,bb,xx)
C              sphi = xnorm(xx)
C Sfx=(0,sfy,sfz) = projection de Sf sur le plan (y,z)
              bb(1)= 0.d0

              cphix = vscal(aa,bb,3)/xnorm(bb,3)
              PHIX(I) = acos(cphix) * deg 

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
 
      IF(LABEL(DEBSTR(LABEL):FINSTR(LABEL)) .EQ. 
     >                                    'PRINT') THEN
        IF(FIRST) THEN
          FIRST = .FALSE.
          IF(IDLUNI(
     >              LUN)) THEN
            OPEN(UNIT=LUN,FILE='zgoubi.SPNPRT.Out',ERR=96)
          ELSE
            GOTO 96
          ENDIF
          WRITE(LUN,fmt='(A)') 
     >    '# Y, T, Z, P, S, D, TAG, IEX, (SI(J,I),J=1,4)
     >    , (SF(J,I),J=1,4), gamma, G.gamma, PHI, 
     >    , PHIX, ITRAJ, IPASS, Yo, To, Zo, Po, So, Do'
        ENDIF
        DO I=IMAX1,IMAX2
          IF( IEX(I) .GE. -1 ) THEN
            P = BORO*CL9 *F(1,I) *Q
            GAMA = SQRT(P*P + AM*AM)/AM
            WRITE(LUN,111) 
     >      (F(J,I),J=2,6),F(1,I)
     >      ,'''',LET(I),'''',IEX(I),(SI(J,I),J=1,4)
     >      ,(SF(J,I),J=1,4),GAMA,G*GAMA,PHI(I),PHIX(I),I,ipass,noel
     >      ,(FO(J,I),J=2,6),FO(1,I)
 111        FORMAT(1X,1p,6(1X,E14.6),1X,3A1,1X,I2,12(1X,e14.6),3(1X,I4)
     >      ,6(1X,E14.6))
          ENDIF
        ENDDO
C Leaving unclosed allows stacking when combined use of FIT and REBELOTE
C        CLOSE(LUN)
      ENDIF 

      RETURN

 96       CONTINUE
          WRITE(ABS(NRES),FMT='(/,''SBR SPNPRT : '',
     >               ''Error open file zgoubi.SPNPRT.Out'')')
          WRITE(*        ,FMT='(/,''SBR SPNPRT : '',
     >               ''Error open file zgoubi.SPNPRT.Out'')')

      RETURN
      END
