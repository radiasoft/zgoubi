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
      SUBROUTINE AGSKS(NOEL,P,
     >                        AK1,AK2,AKS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AK1(*), AK2(*), AKS(*)
      INCLUDE 'MXFS.H'
      INCLUDE 'MXLD.H'
      PARAMETER (LBLSIZ=10)
      CHARACTER(LBLSIZ) LABEL
      COMMON /LABEL/ LABEL(MXL,2)
      PARAMETER (KSIZ=10)
      CHARACTER FAM*(KSIZ),LBF*(LBLSIZ)
      COMMON/SCALT/ FAM(MXF),LBF(MXF,MLF),JPA(MXF,MXP)

      DIMENSION AP(50), ADBD(50), ADBF(50)

      DIMENSION PNOW(6)
      SAVE PNOW

      CHARACTER(LBLSIZ) LBL1

      PARAMETER (BDOT=0.D0)
      PARAMETER(CONV = 0.0254D0)
      PARAMETER(ALA = 94D0*CONV)        !LENGTH IN METERS
      PARAMETER(ALB = 79D0*CONV)
      PARAMETER(ALC = 94D0*CONV)

      DIMENSION SAK1(6), SAK2(6)
      SAVE SAK1, SAK2

      PARAMETER(   UKAFM3 =  0.0D0)
      PARAMETER(   UKAFM2 = -9.9628200D-6)
      PARAMETER(   UKAFM1 = -2.7124400D-4)
      PARAMETER(   UKAF0 =   0.0487709D0)
      PARAMETER(   UKAF1 =   4.2519700D-5)
      PARAMETER(   UKAF2 =  -1.5027310D-5 )
      PARAMETER(   UKAF3 =   1.9900440D-6)
      PARAMETER(   UKAF4 =  -1.2474890D-7)
      PARAMETER(   UKAF5 =   3.7239700D-9)
      PARAMETER(   UKAF6 =  -4.3303320D-11)

      PARAMETER(    UKADM3 =  0.0D0)
      PARAMETER(   UKADM2 =  1.0494200D-5)
      PARAMETER(   UKADM1 =  2.7764700D-4)
      PARAMETER(   UKAD0 =  -0.0487187D0)
      PARAMETER(   UKAD1 =  -4.4223160D-5)
      PARAMETER(   UKAD2 =   1.5432060D-5)
      PARAMETER(   UKAD3 =  -2.0235660D-6)
      PARAMETER(   UKAD4 =   1.2581210D-7)
      PARAMETER(   UKAD5 =  -3.7265440D-9)
      PARAMETER(   UKAD6 =   4.3053950D-11)

      PARAMETER(    UKBFM3 =  0.0D0)
      PARAMETER(   UKBFM2 = -1.0370300D-5)
      PARAMETER(   UKBFM1 = -2.7908700D-4)
      PARAMETER(   UKBF0 =   0.0485774D0)
      PARAMETER(   UKBF1 =   4.1807340D-5)
      PARAMETER(   UKBF2 =  -1.4546960D-5)
      PARAMETER(   UKBF3 =   1.9047810D-6)
      PARAMETER(   UKBF4 =  -1.1846940D-7)
      PARAMETER(   UKBF5 =   3.5161530D-9)
      PARAMETER(   UKBF6 =  -4.0853520D-11)

      PARAMETER(    UKBDM3 =  0.0D0)
      PARAMETER(   UKBDM2 =  1.0855000D-5)
      PARAMETER(   UKBDM1 =  2.8391800D-4)
      PARAMETER(   UKBD0 =  -0.048532D0)
      PARAMETER(   UKBD1 =  -4.3289400D-5)
      PARAMETER(   UKBD2 =   1.4902690D-5)
      PARAMETER(   UKBD3 =  -1.9351350D-6)
      PARAMETER(   UKBD4 =   1.1951780D-7)
      PARAMETER(   UKBD5 =  -3.5235550D-9)
      PARAMETER(   UKBD6 =   4.0708920D-11)

      PARAMETER(    UKCFM3 =  0.0D0)
      PARAMETER(   UKCFM2 =  1.9389200D-5)
      PARAMETER(   UKCFM1 =  3.7263000D-4)
      PARAMETER(   UKCF0 =   0.0485356D0)
      PARAMETER(   UKCF1 =   1.5650290D-5)
      PARAMETER(   UKCF2 =  -7.0134400D-6)
      PARAMETER(   UKCF3 =   1.1275010D-6)
      PARAMETER(   UKCF4 =  -8.2173200D-8)
      PARAMETER(   UKCF5 =   2.7713320D-9)
      PARAMETER(   UKCF6 =  -3.5625590D-11)

      PARAMETER(    UKCDM3 =  0.0D0)
      PARAMETER(   UKCDM2 = -1.8159700D-5)
      PARAMETER(   UKCDM1 = -3.6821400D-4)
      PARAMETER(   UKCD0 =  -0.0484683D0)
      PARAMETER(   UKCD1 =  -1.4804670D-5)
      PARAMETER(   UKCD2 =   6.8025790D-6)
      PARAMETER(   UKCD3 =  -1.1017970D-6)
      PARAMETER(   UKCD4 =   8.0248420D-8)
      PARAMETER(   UKCD5 =  -2.6963500D-9)

      PARAMETER(    UKCD6 =   3.4569950D-11)

      PARAMETER(   TKADM3 = -4.54103D-5)
      PARAMETER(   TKADM2 =  3.86648D-4)
      PARAMETER(   TKADM1 =  -5.15221D-3)
      PARAMETER(   TKAD0  = -6.23676D-3)
      PARAMETER(   TKAD1  =  -8.21074D-5)
      PARAMETER(   TKAD2  =  2.94841D-5)
      PARAMETER(   TKAD3  =  -2.63597D-6  )
      PARAMETER(   TKAD4  =  2.17817D-9)
      PARAMETER(   TKAD5  =  6.02362D-9)
      PARAMETER(   TKAD6  = -1.60702D-10)
      PARAMETER(   TKAFM3 =  -2.11163D-5)
      PARAMETER(   TKAFM2 =  2.31252D-4)
      PARAMETER(   TKAFM1 =  -4.98909D-3)
      PARAMETER(   TKAF0  =  -6.35613D-3)
      PARAMETER(   TKAF1  =  -1.3545D-4)
      PARAMETER(   TKAF2  =  5.20196D-5)
      PARAMETER(   TKAF3  =  -5.93495D-6)
      PARAMETER(   TKAF4  =  2.12422D-7)
      PARAMETER(   TKAF5  =  6.36497D-11)
      PARAMETER(   TKAF6  =  -9.88187D-11)
      PARAMETER(   TKBDM3 =  -5.33308D-5  )
      PARAMETER(   TKBDM2 =  4.61419D-4)
      PARAMETER(   TKBDM1 =  -5.36794D-3)
      PARAMETER(   TKBD0  =  -7.54801D-3)
      PARAMETER(   TKBD1  =  -8.46257D-5)
      PARAMETER(   TKBD2  =  3.25779D-5)
      PARAMETER(   TKBD3  =  -3.46201D-6)
      PARAMETER(   TKBD4  =  7.60854D-8)
      PARAMETER(   TKBD5  =  3.2396D-9)
      PARAMETER(   TKBD6  =  -1.24659D-10)
      PARAMETER(   TKBFM3 =  -3.35980D-5)
      PARAMETER(   TKBFM2 =  3.42004D-4)
      PARAMETER(   TKBFM1 =  -5.26172D-3)
      PARAMETER(   TKBF0  =  -7.65876D-3)
      PARAMETER(   TKBF1  =  -1.20210D-4)
      PARAMETER(   TKBF2  =  4.83680D-5)
      PARAMETER(   TKBF3  =  -5.85176D-6)
      PARAMETER(   TKBF4  =  2.31758D-7)
      PARAMETER(   TKBF5  =  -1.23709D-9)
      PARAMETER(   TKBF6  =  -7.77011D-11)
      PARAMETER(   TKCDM3 =  5.46173D-5)
      PARAMETER(   TKCDM2 =  -3.1184D-4)
      PARAMETER(   TKCDM1 =  4.53022D-3)
      PARAMETER(   TKCD0  =  -1.03323D-2)
      PARAMETER(   TKCD1  =  -5.10398D-4)
      PARAMETER(   TKCD2  =  1.74831D-4)
      PARAMETER(   TKCD3  =  -1.90748D-5)
      PARAMETER(   TKCD4  =  9.02456D-7)
      PARAMETER(   TKCD5  =  -1.75112D-8)
      PARAMETER(   TKCD6  =  7.24016D-11)
      PARAMETER(   TKCFM3 =  3.24138D-5)
      PARAMETER(   TKCFM2 =  -1.49126D-4)
      PARAMETER(   TKCFM1 =  4.2301D-3)
      PARAMETER(   TKCF0  =  -1.03626D-2)
      PARAMETER(   TKCF1  =  -4.55906D-4)
      PARAMETER(   TKCF2  =  1.55926D-4)
      PARAMETER(   TKCF3  =  -1.70233D-5)
      PARAMETER(   TKCF4  =  7.91965D-7)
      PARAMETER(   TKCF5  =  -1.44494D-8)
      PARAMETER(   TKCF6  =  3.8317D-11)

C CORRECTION FOR B-DOT EFFECTS IS ADDED IN HERE.  K.BROWN (8/31/98)
      PARAMETER(   TKDBDCOEF = -0.0025D0)
      PARAMETER(   TKFBDCOEF = -0.0025D0)

      INTEGER FINSTR

      DATA  AP(1) /    2.1633D0 / !170 
      DATA  ADBD(1) /  -0.1D-3 /
      DATA  ADBF(1) /   1.0D-3 /
      DATA   AP(2) /    2.3083D0 / !190 
      DATA  ADBD(2) /  -1.2D-3 /
      DATA  ADBF(2) /   0.6D-3 /
      DATA   AP(3) /    2.6140D0 / !210 
      DATA  ADBD(3) /  -1.2D-3 /
      DATA  ADBF(3) /   0.9D-3 /
      DATA   AP(4) /    3.9661D0 / !250 
      DATA  ADBD(4) /  -1.0D-3 /
      DATA  ADBF(4) /   1.2D-3 /
      DATA   AP(5) /    4.4980D0 / !260,MATCHED WITH INTERPOLATED TUNES 
      DATA  ADBD(5) /  -0.9D-3 /
      DATA  ADBF(5) /   1.8D-3 /
      DATA   AP(6) /    5.1099D0 / !270 
      DATA  ADBD(6) /  -0.4D-3 /
      DATA  ADBF(6) /   1.5D-3 /
      DATA   AP(7) /    7.6973D0 / !310 
      DATA  ADBD(7) /  -0.7D-3 /
      DATA  ADBF(7) /   0.2D-3 /
      DATA   AP(8) /    10.231D0 / !350 
      DATA  ADBD(8) /   0.4D-3 /
      DATA  ADBF(8) /   1.6D-3 /
      DATA   AP(9) /    15.180D0 / !430 
      DATA  ADBD(9) /   0.5D-3 /
      DATA  ADBF(9) /   1.7D-3 /
      DATA   AP(10) /    23.732D0 / !590 
      DATA  ADBD(10) /   1.1D-3 /
      DATA  ADBF(10) /   2.7D-3 /
        
      DATA  ISCLMAX / 10 /

      DATA PNOW / 6 * -1D10 /
      DATA DB1D, DB1F / 2 * 0.D0 /

      LBL1 = LABEL(NOEL,1)
      LAST = FINSTR(LBL1)

      IF    (LBL1(LAST-1:LAST) .EQ. 'BF') THEN
        IMM = 1
        IF(PNOW(1) .EQ. P) THEN
          AK1(1)=SAK1(1)
          AK2(1)=SAK2(1)
        ELSE
          PNOW(1) = P
          IF( P .LT. AP(1) ) THEN
            DB1F = ADBF(1)
          ELSEIF( P .GT. AP(ISCLMAX)) THEN 
            DB1F = ADBF(ISCLMAX)
          ELSE
            DO I=1, ISCLMAX-1
               IF( P .GE. AP(I) .AND. P .LE. AP(I+1)) THEN
                 DB1F = ADBF(I) + 
     >           ((ADBF(I) - ADBF(I+1))/(AP(I) - AP(I+1)))*(P-AP(I))
                 GOTO 861
               ENDIF
            ENDDO
 861        CONTINUE
          ENDIF

          DB1F = 0.D0
          P2 = P*P
          P3 = P*P2
          P4 = P*P3
          P5 = P*P4
          P6 = P*P5
C  K1BF :
          AK1(1)= UKBFM3/(P3)+ UKBFM2/(P2)+ UKBFM1/P+ UKBF0+UKBF1*P 
     >    +UKBF2*P2+UKBF3*P3+UKBF4*P4+UKBF5*P5 
     >    +UKBF6*P6 
          AK1(1) = AK1(1)*(DB1F + 1.D0)
C  K2BF :
          AK2(1)= TKBFM3/(P3)+ TKBFM2/(P2)+ TKBFM1/P + TKBF0+TKBF1*P 
     >    +TKBF2*P2+TKBF3*P3+TKBF4*P4+TKBF5*P5  
     >    +TKBF6*P6 - TKFBDCOEF*BDOT*ALB/ALA
  
          SAK1(1)=AK1(1)
          SAK2(1)=AK2(1)

        ENDIF

      ELSEIF(LBL1(LAST-1:LAST) .EQ. 'CD') THEN
        IMM = 2
        IF(PNOW(2) .EQ. P) THEN
          AK1(2)=SAK1(2)
          AK2(2)=SAK2(2)
        ELSE
          PNOW(2) = P
          IF( P .LT. AP(1) ) THEN
            DB1D = ADBD(1)
          ELSEIF( P .GT. AP(ISCLMAX)) THEN 
            DB1D = ADBD(ISCLMAX)
          ELSE
            DO I=1, ISCLMAX-1
               IF( P .GE. AP(I) .AND. P .LE. AP(I+1)) THEN
                 DB1D = ADBD(I) + 
     >           ((ADBD(I) - ADBD(I+1))/(AP(I) - AP(I+1)))*(P-AP(I))
                 GOTO 862
               ENDIF
            ENDDO
 862        CONTINUE
          ENDIF

          DB1D = 0.D0
          P2 = P*P
          P3 = P*P2
          P4 = P*P3
          P5 = P*P4
          P6 = P*P5
C  K1CD :
          AK1(2)= UKCDM3/(P3)+ UKCDM2/(P2)+ UKCDM1/P+ UKCD0+UKCD1*P 
     >    +UKCD2*P2+UKCD3*P3+UKCD4*P4+UKCD5*P5 
     >    +UKCD6*P6 
          AK1(2) = AK1(2)*(DB1D + 1.D0)
C  K2CD :
          AK2(2)= TKCDM3/(P3)+ TKCDM2/(P2)+ TKCDM1/P + TKCD0+TKCD1*P 
     >    +TKCD2*P2+TKCD3*P3+TKCD4*P4+TKCD5*P5  
     >    +TKCD6*P6 - TKDBDCOEF*BDOT*ALC/ALA

          SAK1(2)=AK1(2)
          SAK2(2)=AK2(2)

        ENDIF
      ELSEIF(LBL1(LAST-1:LAST) .EQ. 'AF') THEN
        IMM = 3
        IF(PNOW(3) .EQ. P) THEN
          AK1(3)=SAK1(3)
          AK2(3)=SAK2(3)
        ELSE
          PNOW(3) = P
          IF( P .LT. AP(1) ) THEN
            DB1F = ADBF(1)
          ELSEIF( P .GT. AP(ISCLMAX)) THEN 
            DB1F = ADBF(ISCLMAX)
          ELSE
            DO I=1, ISCLMAX-1
               IF( P .GE. AP(I) .AND. P .LE. AP(I+1)) THEN
                 DB1F = ADBF(I) + 
     >           ((ADBF(I) - ADBF(I+1))/(AP(I) - AP(I+1)))*(P-AP(I))
                 GOTO 863
               ENDIF
            ENDDO
 863        CONTINUE
          ENDIF

          DB1F = 0.D0
          P2 = P*P
          P3 = P*P2
          P4 = P*P3
          P5 = P*P4
          P6 = P*P5
C  K1AF :
          AK1(3)= UKAFM3/(P3)+ UKAFM2/(P2)+ UKAFM1/P+ UKAF0+UKAF1*P 
     >    +UKAF2*P2+UKAF3*P3+UKAF4*P4+UKAF5*P5 
     >    +UKAF6*P6 
           AK1(3) = AK1(3)*(DB1F + 1.D0)
C  K2AF :
          AK2(3)= TKAFM3/(P3)+ TKAFM2/(P2)+ TKAFM1/P + TKAF0+TKAF1*P 
     >    +TKAF2*P2+TKAF3*P3+TKAF4*P4+TKAF5*P5  
     >    +TKAF6*P6 - TKFBDCOEF*BDOT*ALA/ALA

          SAK1(3)=AK1(3)
          SAK2(3)=AK2(3)

        ENDIF
      ELSEIF(LBL1(LAST-1:LAST) .EQ. 'BD') THEN
        IMM = 4
        IF(PNOW(4) .EQ. P) THEN
          AK1(4)=SAK1(4)
          AK2(4)=SAK2(4)
        ELSE
          PNOW(4) = P
          IF( P .LT. AP(1) ) THEN
            DB1D = ADBD(1)
          ELSEIF( P .GT. AP(ISCLMAX)) THEN 
            DB1D = ADBD(ISCLMAX)
          ELSE
            DO I=1, ISCLMAX-1
               IF( P .GE. AP(I) .AND. P .LE. AP(I+1)) THEN
                 DB1D = ADBD(I) + 
     >           ((ADBD(I) - ADBD(I+1))/(AP(I) - AP(I+1)))*(P-AP(I))
                 GOTO 864
               ENDIF
            ENDDO
 864        CONTINUE
          ENDIF

          DB1D = 0.D0
          P2 = P*P
          P3 = P*P2
          P4 = P*P3
          P5 = P*P4
          P6 = P*P5
C  K1BD :
          AK1(4)= UKBDM3/(P3)+ UKBDM2/(P2)+ UKBDM1/P+ UKBD0+UKBD1*P 
     >    +UKBD2*P2+UKBD3*P3+UKBD4*P4+UKBD5*P5 
     >    +UKBD6*P6 
          AK1(4) = AK1(4)*(DB1D + 1.D0)
C  K2BD :
          AK2(4)= TKBDM3/(P3)+ TKBDM2/(P2)+ TKBDM1/P + TKBD0+TKBD1*P 
     >    +TKBD2*P2+TKBD3*P3+TKBD4*P4+TKBD5*P5  
     >    +TKBD6*P6 - TKDBDCOEF*BDOT*ALB/ALA

          SAK1(4)=AK1(4)
          SAK2(4)=AK2(4)

        ENDIF
      ELSEIF(LBL1(LAST-1:LAST) .EQ. 'CF') THEN
        IMM = 5
        IF(PNOW(5) .EQ. P) THEN
          AK1(5)=SAK1(5)
          AK2(5)=SAK2(5)
        ELSE
          PNOW(5) = P
          IF( P .LT. AP(1) ) THEN
            DB1F = ADBF(1)
          ELSEIF( P .GT. AP(ISCLMAX)) THEN 
            DB1F = ADBF(ISCLMAX)
          ELSE
            DO I=1, ISCLMAX-1
               IF( P .GE. AP(I) .AND. P .LE. AP(I+1)) THEN
                 DB1F = ADBF(I) + 
     >           ((ADBF(I) - ADBF(I+1))/(AP(I) - AP(I+1)))*(P-AP(I))
                 GOTO 865
               ENDIF
            ENDDO
 865        CONTINUE
          ENDIF

          DB1F = 0.D0
          P2 = P*P
          P3 = P*P2
          P4 = P*P3
          P5 = P*P4
          P6 = P*P5
C  K1CF :
          AK1(5)= UKCFM3/(P3)+ UKCFM2/(P2)+ UKCFM1/P+ UKCF0+UKCF1*P 
     >    +UKCF2*P2+UKCF3*P3+UKCF4*P4+UKCF5*P5 
     >    +UKCF6*P6 
          AK1(5) = AK1(5)*(DB1F + 1.D0)
C  K2CF :
          AK2(5)= TKCFM3/(P3)+ TKCFM2/(P2)+ TKCFM1/P + TKCF0+TKCF1*P 
     >    +TKCF2*P2+TKCF3*P3+TKCF4*P4+TKCF5*P5  
     >    +TKCF6*P6 - TKFBDCOEF*BDOT*ALC/ALA

          SAK1(5)=AK1(5)
          SAK2(5)=AK2(5)

        ENDIF
      ELSEIF(LBL1(LAST-1:LAST) .EQ. 'AD') THEN
        IMM = 6
        IF(PNOW(6) .EQ. P) THEN
          AK1(6)=SAK1(6)
          AK2(6)=SAK2(6)
        ELSE
          PNOW(6) = P
          IF( P .LT. AP(1) ) THEN
            DB1D = ADBD(1)
          ELSEIF( P .GT. AP(ISCLMAX)) THEN 
            DB1D = ADBD(ISCLMAX)
          ELSE
            DO I=1, ISCLMAX-1
               IF( P .GE. AP(I) .AND. P .LE. AP(I+1)) THEN
                 DB1D = ADBD(I) + 
     >           ((ADBD(I) - ADBD(I+1))/(AP(I) - AP(I+1)))*(P-AP(I))
                 GOTO 866
               ENDIF
            ENDDO
 866        CONTINUE
          ENDIF

          DB1D = 0.D0
          P2 = P*P
          P3 = P*P2
          P4 = P*P3
          P5 = P*P4
          P6 = P*P5
C  K1AD :
          AK1(6)= UKADM3/(P3)+ UKADM2/(P2)+ UKADM1/P+ UKAD0+UKAD1*P 
     >    +UKAD2*P2+UKAD3*P3+UKAD4*P4+UKAD5*P5  
     >    +UKAD6*P6   
          AK1(6) = AK1(6)*(DB1D + 1.D0)
C  K2AD :
          AK2(6)= TKADM3/(P3)+ TKADM2/(P2)+ TKADM1/P + TKAD0+TKAD1*P 
     >    +TKAD2*P2+TKAD3*P3+TKAD4*P4+TKAD5*P5  
     >    +TKAD6*P6 - TKDBDCOEF*BDOT*ALA/ALA

          SAK1(6)=AK1(6)
          SAK2(6)=AK2(6)

        ENDIF
      ENDIF

      AKS(2) = AK1(IMM) 
      AKS(3) = AK2(IMM) 
      RETURN
      END
