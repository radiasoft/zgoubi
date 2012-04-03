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
      SUBROUTINE AGSKS(P,
     >                   AK1,AK2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension ak1(*), ak2(*)
      dimension ap(50), adbd(50), adbf(50)
      parameter (bdot=0.d0)
      parameter(conv = 0.0254d0)
      parameter(ALA = 94d0*CONV)        !length in meters
      parameter(ALB = 79d0*CONV)
      parameter(ALC = 94d0*CONV)
      SAVE sak11, sak12, sak13, sak14, sak15, sak16 
      SAVE sak21, sak22, sak23, sak24, sak25, sak26
      SAVE Pnow
      data Pnow / -1d10 /
      data db1d, db1f / 2*0.d0 /

      if(Pnow .EQ. P) then
         ak1(1)=sak11
         ak1(2)=sak12
         ak1(3)=sak13
         ak1(4)=sak14
         ak1(5)=sak15
         ak1(6)=sak16
         ak2(1)=sak21
         ak2(2)=sak22
         ak2(3)=sak23
         ak2(4)=sak24
         ak2(5)=sak25
         ak2(6)=sak26
         goto 695
      endif

      Pnow = P
c      write(*,*)'agsks.f', p
c      read(*,*)

      UKAFM3 =  0.0D0
      UKAFM2 = -9.9628200D-6
      UKAFM1 = -2.7124400D-4
      UKAF0 =   0.0487709D0
      UKAF1 =   4.2519700D-5
      UKAF2 =  -1.5027310D-5 
      UKAF3 =   1.9900440D-6
      UKAF4 =  -1.2474890D-7
      UKAF5 =   3.7239700D-9
      UKAF6 =  -4.3303320D-11

      UKADM3 =  0.0D0
      UKADM2 =  1.0494200D-5
      UKADM1 =  2.7764700D-4
      UKAD0 =  -0.0487187D0
      UKAD1 =  -4.4223160D-5
      UKAD2 =   1.5432060D-5
      UKAD3 =  -2.0235660D-6
      UKAD4 =   1.2581210D-7
      UKAD5 =  -3.7265440D-9
      UKAD6 =   4.3053950D-11

      UKBFM3 =  0.0D0
      UKBFM2 = -1.0370300D-5
      UKBFM1 = -2.7908700D-4
      UKBF0 =   0.0485774D0
      UKBF1 =   4.1807340D-5
      UKBF2 =  -1.4546960D-5
      UKBF3 =   1.9047810D-6
      UKBF4 =  -1.1846940D-7
      UKBF5 =   3.5161530D-9
      UKBF6 =  -4.0853520D-11

      UKBDM3 =  0.0D0
      UKBDM2 =  1.0855000D-5
      UKBDM1 =  2.8391800D-4
      UKBD0 =  -0.048532D0
      UKBD1 =  -4.3289400D-5
      UKBD2 =   1.4902690D-5
      UKBD3 =  -1.9351350D-6
      UKBD4 =   1.1951780D-7
      UKBD5 =  -3.5235550D-9
      UKBD6 =   4.0708920D-11

      UKCFM3 =  0.0D0
      UKCFM2 =  1.9389200D-5
      UKCFM1 =  3.7263000D-4
      UKCF0 =   0.0485356D0
      UKCF1 =   1.5650290D-5
      UKCF2 =  -7.0134400D-6
      UKCF3 =   1.1275010D-6
      UKCF4 =  -8.2173200D-8
      UKCF5 =   2.7713320D-9
      UKCF6 =  -3.5625590D-11

      UKCDM3 =  0.0D0
      UKCDM2 = -1.8159700D-5
      UKCDM1 = -3.6821400D-4
      UKCD0 =  -0.0484683D0
      UKCD1 =  -1.4804670D-5
      UKCD2 =   6.8025790D-6
      UKCD3 =  -1.1017970D-6
      UKCD4 =   8.0248420D-8
      UKCD5 =  -2.6963500D-9

      UKCD6 =   3.4569950D-11


      ap(1) =    2.1633d0 !170
      adbd(1) =  -0.1d-3
      adbf(1) =   1.0d-3

      ap(2) =    2.3083d0 !190
      adbd(2) =  -1.2d-3
      adbf(2) =   0.6d-3

      ap(3) =    2.6140d0 !210
      adbd(3) =  -1.2d-3
      adbf(3) =   0.9d-3

      ap(4) =    3.9661d0 !250
      adbd(4) =  -1.0d-3
      adbf(4) =   1.2d-3

      ap(5) =    4.4980d0 !260,matched with interpolated tunes
      adbd(5) =  -0.9d-3
      adbf(5) =   1.8d-3

      ap(6) =    5.1099d0 !270
      adbd(6) =  -0.4d-3
      adbf(6) =   1.5d-3

      ap(7) =    7.6973d0 !310
      adbd(7) =  -0.7d-3
      adbf(7) =   0.2d-3

      ap(8) =    10.231d0 !350
      adbd(8) =   0.4d-3
      adbf(8) =   1.6d-3

      ap(9) =    15.180d0 !430
      adbd(9) =   0.5d-3
      adbf(9) =   1.7d-3

      ap(10) =    23.732d0 !590
      adbd(10) =   1.1d-3
      adbf(10) =   2.7d-3
      


      isclMAX = 10

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if( p .LT. ap(1) ) then
         db1d = adbd(1)
         db1f = adbf(1)
      elseif( p .GT. ap(isclMAX)) then 
         db1d = adbd(isclMAX)
         db1f = adbf(isclMAX)
      else
         do i=1, isclMAX-1
            if( p .GE. ap(i) .AND. p .LE. ap(i+1)) then
               db1d = adbd(i) + 
     > ((adbd(i) - adbd(i+1))/(ap(i) - ap(i+1)))*(p-ap(i))
               db1f = adbf(i) + 
     > ((adbf(i) - adbf(i+1))/(ap(i) - ap(i+1)))*(p-ap(i))
               goto 865
            endif
         enddo
 865     continue
      endif
!!!!!!!!!
      db1d = 0.d0
      db1f = 0.d0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C  K1CD :
       AK1(2)= UKCDM3/(p*p*p)+ UKCDM2/(p*p)+ UKCDM1/p + UKCD0+UKCD1*p 
     >  +UKCD2*p*p+UKCD3*p*p*p+UKCD4*p*p*p*p+UKCD5*p*p*p*p*p 
     >  +UKCD6*p*p*p*p*p*p 


C  K1BD :
       AK1(4)= UKBDM3/(p*p*p)+ UKBDM2/(p*p)+ UKBDM1/p + UKBD0+UKBD1*p 
     >  +UKBD2*p*p+UKBD3*p*p*p+UKBD4*p*p*p*p+UKBD5*p*p*p*p*p 
     >  +UKBD6*p*p*p*p*p*p 


C  K1AD :
       AK1(6)= UKADM3/(p*p*p)+ UKADM2/(p*p)+ UKADM1/p + UKAD0+UKAD1*p 
     >  +UKAD2*p*p+UKAD3*p*p*p+UKAD4*p*p*p*p+UKAD5*p*p*p*p*p  
     >  +UKAD6*p*p*p*p*p*p   

       AK1(2) = AK1(2)*(db1d + 1.d0)
       AK1(4) = AK1(4)*(db1d + 1.d0)
       AK1(6) = AK1(6)*(db1d + 1.d0)




C  K1BF :
       AK1(1)= UKBFM3/(p*p*p)+ UKBFM2/(p*p)+ UKBFM1/p + UKBF0+UKBF1*p 
     >  +UKBF2*p*p+UKBF3*p*p*p+UKBF4*p*p*p*p+UKBF5*p*p*p*p*p 
     >  +UKBF6*p*p*p*p*p*p 


C  K1AF :
       AK1(3)= UKAFM3/(p*p*p)+ UKAFM2/(p*p)+ UKAFM1/p + UKAF0+UKAF1*p 
     >  +UKAF2*p*p+UKAF3*p*p*p+UKAF4*p*p*p*p+UKAF5*p*p*p*p*p 
     >  +UKAF6*p*p*p*p*p*p 


C  K1CF :
       AK1(5)= UKCFM3/(p*p*p)+ UKCFM2/(p*p)+ UKCFM1/p + UKCF0+UKCF1*p 
     >  +UKCF2*p*p+UKCF3*p*p*p+UKCF4*p*p*p*p+UKCF5*p*p*p*p*p 
     >  +UKCF6*p*p*p*p*p*p 

       AK1(1) = AK1(1)*(db1f + 1.d0)
       AK1(3) = AK1(3)*(db1f + 1.d0)
       AK1(5) = AK1(5)*(db1f + 1.d0)


 
c         print *, p, ak1
c         stop
!!4.2 Main Magnet Sextupole strength vs momentum
!===========================

      TKADM3 = -4.54103D-5
      TKADM2 =  3.86648D-4
      TKADM1 =  -5.15221D-3
      TKAD0  = -6.23676D-3
      TKAD1  =  -8.21074D-5
      TKAD2  =  2.94841D-5
      TKAD3  =  -2.63597D-6  
      TKAD4  =  2.17817D-9
      TKAD5  =  6.02362D-9
      TKAD6  = -1.60702D-10

      TKAFM3 =  -2.11163D-5
      TKAFM2 =  2.31252D-4
      TKAFM1 =  -4.98909D-3
      TKAF0  =  -6.35613D-3
      TKAF1  =  -1.3545D-4
      TKAF2  =  5.20196D-5
      TKAF3  =  -5.93495D-6
      TKAF4  =  2.12422D-7
      TKAF5  =  6.36497D-11
      TKAF6  =  -9.88187D-11

      TKBDM3 =  -5.33308D-5  
      TKBDM2 =  4.61419D-4
      TKBDM1 =  -5.36794D-3
      TKBD0  =  -7.54801D-3
      TKBD1  =  -8.46257D-5
      TKBD2  =  3.25779D-5
      TKBD3  =  -3.46201D-6
      TKBD4  =  7.60854D-8
      TKBD5  =  3.2396D-9
      TKBD6  =  -1.24659D-10

      TKBFM3 =  -3.35980D-5
      TKBFM2 =  3.42004D-4
      TKBFM1 =  -5.26172D-3
      TKBF0  =  -7.65876D-3
      TKBF1  =  -1.20210D-4
      TKBF2  =  4.83680D-5
      TKBF3  =  -5.85176D-6
      TKBF4  =  2.31758D-7
      TKBF5  =  -1.23709D-9
      TKBF6  =  -7.77011D-11

      TKCDM3 =  5.46173D-5
      TKCDM2 =  -3.1184D-4
      TKCDM1 =  4.53022D-3
      TKCD0  =  -1.03323D-2
      TKCD1  =  -5.10398D-4
      TKCD2  =  1.74831D-4
      TKCD3  =  -1.90748D-5
      TKCD4  =  9.02456D-7
      TKCD5  =  -1.75112D-8
      TKCD6  =  7.24016D-11

      TKCFM3 =  3.24138D-5
      TKCFM2 =  -1.49126D-4
      TKCFM1 =  4.2301D-3
      TKCF0  =  -1.03626D-2
      TKCF1  =  -4.55906D-4
      TKCF2  =  1.55926D-4
      TKCF3  =  -1.70233D-5
      TKCF4  =  7.91965D-7
      TKCF5  =  -1.44494D-8
      TKCF6  =  3.8317D-11

!! Correction for B-dot effects is added in here.  K.Brown (8/31/98)
      TKDBDCOEF = -0.0025d0
      TKFBDCOEF = -0.0025d0

C  K2AD :
      AK2(6)= TKADM3/(p*p*p)+ TKADM2/(p*p)+ TKADM1/p + TKAD0+TKAD1*p 
     >  +TKAD2*p*p+TKAD3*p*p*p+TKAD4*p*p*p*p+TKAD5*p*p*p*p*p  
     >  +TKAD6*p*p*p*p*p*p - TKDBDCOEF*BDOT*ALA/ALA
      
C  K2AF :
      AK2(3)= TKAFM3/(p*p*p)+ TKAFM2/(p*p)+ TKAFM1/p + TKAF0+TKAF1*p 
     >  +TKAF2*p*p+TKAF3*p*p*p+TKAF4*p*p*p*p+TKAF5*p*p*p*p*p  
     >  +TKAF6*p*p*p*p*p*p - TKFBDCOEF*BDOT*ALA/ALA

C  K2BD :
      AK2(4)= TKBDM3/(p*p*p)+ TKBDM2/(p*p)+ TKBDM1/p + TKBD0+TKBD1*p 
     >  +TKBD2*p*p+TKBD3*p*p*p+TKBD4*p*p*p*p+TKBD5*p*p*p*p*p  
     >  +TKBD6*p*p*p*p*p*p - TKDBDCOEF*BDOT*ALB/ALA

C  K2BF :
      AK2(1)= TKBFM3/(p*p*p)+ TKBFM2/(p*p)+ TKBFM1/p + TKBF0+TKBF1*p 
     >  +TKBF2*p*p+TKBF3*p*p*p+TKBF4*p*p*p*p+TKBF5*p*p*p*p*p  
     >  +TKBF6*p*p*p*p*p*p - TKFBDCOEF*BDOT*ALB/ALA

C  K2CD :
      AK2(2)= TKCDM3/(p*p*p)+ TKCDM2/(p*p)+ TKCDM1/p + TKCD0+TKCD1*p 
     >  +TKCD2*p*p+TKCD3*p*p*p+TKCD4*p*p*p*p+TKCD5*p*p*p*p*p  
     >  +TKCD6*p*p*p*p*p*p - TKDBDCOEF*BDOT*ALC/ALA

C  K2CF :
      AK2(5)= TKCFM3/(p*p*p)+ TKCFM2/(p*p)+ TKCFM1/p + TKCF0+TKCF1*p 
     >  +TKCF2*p*p+TKCF3*p*p*p+TKCF4*p*p*p*p+TKCF5*p*p*p*p*p  
     >  +TKCF6*p*p*p*p*p*p - TKFBDCOEF*BDOT*ALC/ALA


      sak11=ak1(1)
      sak12=ak1(2)
      sak13=ak1(3)
      sak14=ak1(4)
      sak15=ak1(5)
      sak16=ak1(6)
      sak21=ak2(1)
      sak22=ak2(2)
      sak23=ak2(3)
      sak24=ak2(4)
      sak25=ak2(5)
      sak26=ak2(6)
         

 695  continue

      RETURN
      END
