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
      SUBROUTINE AGSQKS(NOEL,A1,A2,A3,XL,
     >                                   BM)
      IMPLICIT DOUBLE PRECISION (A-Z)
      INTEGER NOEL
      INCLUDE 'MXFS.H'
      INCLUDE 'MXLD.H'
      INTEGER LBLSIZ, KSIZ 
      PARAMETER (LBLSIZ=10)
      PARAMETER (KSIZ=10)
      CHARACTER FAM*(KSIZ),LBF*(LBLSIZ),KLEY*(KSIZ),LABEL*(LBLSIZ)
      COMMON/SCALT/ FAM(MXF),LBF(MXF,2),KLEY,LABEL(MXL,2)

      CHARACTER(LBLSIZ) LBL1
      INTEGER INIL
      LOGICAL DEBSTR, FINSTR

C  In xags.madx :
      PARAMETER (PQH= 1.D0)
      PARAMETER (PQV=-1.D0)
      
      PARAMETER (EPS = 1D-30)   !  A TINY NUMBER TO PREVENT DIVISION BY ZERO WHEN CALC'ING POLARITY BITS.

      !  THE FOLLOWING ARE THE POLYNOMIAL COEFFICIENTS FOR A TUNE QUAD
      !       TRANSFER FUNCTION WITH NON-ZERO Y-OFFSET.
      PARAMETER (NKQHC6 =  -1.398198D-19  )
      PARAMETER (NKQHC5 =   2.580390D-16  )
      PARAMETER (NKQHC4 =  -1.111090D-13  )
      PARAMETER (NKQHC3 =  -6.402958D-11  )
      PARAMETER (NKQHC2 =   6.384814D-8   )
      PARAMETER (NKQHC1 =   1.719269D-3   )
      PARAMETER (NKQHC0 =   1.771798D-3   )
      
      !  THE FOLLOWING ARE THE POLYNOMIAL COEFFICIENTS FOR A TUNE QUAD
      !       TRANSFER FUNCTION WITH A ZERO Y-OFFSET       
      PARAMETER (KQHC6 = -3.641203D-18  )
      PARAMETER (KQHC5 =  8.233111D-15  )
      PARAMETER (KQHC4 = -7.140673D-12  )
      PARAMETER (KQHC3 =  2.941244D-9   )
      PARAMETER (KQHC2 = -5.725403D-7   )
      PARAMETER (KQHC1 =  1.779240D-3   )
      PARAMETER (KQHC0 =  0.D0)
      
      !  COLD SNAKE COMPENSATION QUADRUPOLES.      
      PARAMETER (A19C6 =  1.917404D-15 )
      PARAMETER (A19C5 = -1.800937D-12 )
      PARAMETER (A19C4 =  5.946546D-10 )
      PARAMETER (A19C3 = -8.898595D-8  )
      PARAMETER (A19C2 =  5.948189D-6 )
      PARAMETER (A19C1 =  2.984878D-3 )
      PARAMETER (A19C0 =  1.025763D-3 *0.D0)
      ! A19C4 = -8.938974D-12 
      ! A19C3 = -4.360751D-9  
      ! A19C2 =  1.787175D-6 
      ! A19C1 =  2.985794D-3 
      ! A19C0 =  2.195629D-3 

      PARAMETER (B1C6 =  1.917404D-15 )
      PARAMETER (B1C5 = -1.800937D-12 )
      PARAMETER (B1C4 =  5.946546D-10 )
      PARAMETER (B1C3 = -8.898595D-8  )
      PARAMETER (B1C2 =  5.948189D-6 )
      PARAMETER (B1C1 =  2.984878D-3 )
      PARAMETER (B1C0 =  1.025763D-3 *0D0)

      !  WARM SNAKE COMPENSATION
      PARAMETER (E19C6 =  1.917404D-15 )
      PARAMETER (E19C5 = -1.800937D-12 )
      PARAMETER (E19C4 =  5.946546D-10 )
      PARAMETER (E19C3 = -8.898595D-8  )
      PARAMETER (E19C2 =  5.948189D-6 )
      PARAMETER (E19C1 =  2.984878D-3 )
      PARAMETER (E19C0 =  1.025763D-3 *0D0)

      PARAMETER (F1C6 =  1.917404D-15 )
      PARAMETER (F1C5 = -1.800937D-12 )
      PARAMETER (F1C4 =  5.946546D-10 )
      PARAMETER (F1C3 = -8.898595D-8  )
      PARAMETER (F1C2 =  5.948189D-6 )
      PARAMETER (F1C1 =  2.984878D-3 )
      PARAMETER (F1C0 =  1.025763D-3 *0D0)
           
      LBL1 = LABEL(NOEL,1)
      INIL = DEBSTR(LBL1)

      ! ! 6.0 TUNE QUADRUPOLES
      ! =================
      IF    (LBL1(INIL:INIL+2) .EQ. 'QH_') THEN
        LENQ = XL
        XQH = A1
        XQHP2 = XQH*XQH 
        XQHP3 = XQH*XQHP2 
        XQHP4 = XQH*XQHP3 
        XQHP5 = XQH*XQHP4 
        XQHP6 = XQH*XQHP5 
        
      ! ==================================================================
      !  WE'RE GOING TO HAVE TO USE THE POLYNOMIAL COEFFICIENTS OVER AND OVER
      !    LATER IN THE FILE, SO I'M GOING TO DEFINE A SET THAT WE'LL USE 
      !    THROUGHOUT THE REST OF THE FILE.  IF ONE WANTS TO CHANGE WHICH POLYNOMIAL
      !    IS USED, ONE ONLY HAS TO CHANGE THE VALUES OF Q0-Q6
      ! ==================================================================
      
        Q6 = KQHC6 
        Q5 = KQHC5 
        Q4 = KQHC4 
        Q3 = KQHC3 
        Q2 = KQHC2 
        Q1 = KQHC1 
        Q0 = KQHC0 
      
      !  ===================================================================
      !  ALMOST ALL OF THE TUNE QUADS HAVE MULTIPLE SUPPLIES WIRED TO THEM.
      !    HERE WE COMBINE THE CURRENTS IN EACH OF THE QUADS, DETERMINE
      !    THE POLARITY BIT AND THEN APPLY THE TRANSFER FUNCTIONS.
      !  ===================================================================
      
      !  HORIZONTAL TUNE QUAD CURRENTS (MADE POSITIVE FOR USE IN THE TRANSFER FUNC. POLYNOMIALS)

        IF    (LBL1(INIL:INIL+5) .EQ. 'QH_A17') THEN
          A17_QUAD_I = A2
          XQH_A17 = ABS(XQH - A17_QUAD_I             ) 
        !  HORIZONTAL TUNE QUAD POLARITY BITS = CURRENT/ABS(CURRENT)
          PQH_A17 =    (XQH - A17_QUAD_I             ) / (XQH_A17+EPS) 
      !  NOW MAKE B STRENGTHS OUT OF THE CURRENTS USING THE TRANSFER FUNCTION (ABOVE)      
          BQH_A17 =   Q6*XQH_A17**6 + Q5*XQH_A17**5 + Q4*XQH_A17**4 
     >            + Q3*XQH_A17**3 
     >            + Q2*XQH_A17**2 + Q1*XQH_A17**1 + Q0*XQH_A17**0 
      !  AND NOW CALCULATE THE K VALUES FROM THE BQH VALUES....
          KQH_A17 =  PQH_A17*BQH_A17/(LENQ)
          BM =  KQH_A17 

c                write(*,*)' agsks. A17. xqh a2 '
c     >        ,xqh,A17_QUAD_I,PQH_A17,BQH_A17,KQH_A17
c                   read(*,*)
          
        ELSEIF(LBL1(INIL:INIL+5) .EQ. 'QH_B17') THEN
          XQH_B17 = ABS(XQH        - 0.D0 * XSB_FHLB17) 
          PQH_B17 =    (XQH        - 0.D0 * XSB_FHLB17) / (XQH_B17+EPS) 
          BQH_B17 =   Q6*XQH_B17**6 + Q5*XQH_B17**5 + Q4*XQH_B17**4 
     >            + Q3*XQH_B17**3 
     >            + Q2*XQH_B17**2 + Q1*XQH_B17**1 + Q0*XQH_B17**0 
          KQH_B17 =  PQH_B17*BQH_B17/(LENQ) 
          BM =  KQH_B17 

c                write(*,*)' agsks. B17. xqh a2 '
c     >        ,xqh,B17_QUAD_I,PQH_B17,BQH_B17,KQH_B17
c                   read(*,*)
          
        ELSEIF(LBL1(INIL:INIL+5) .EQ. 'QH_C17') THEN
          XQH_C17 = ABS(XQH        + 0.D0 * XSB_CEIK17) 
          PQH_C17 =    (XQH        + 0.D0 * XSB_CEIK17) / (XQH_C17+EPS) 
          BQH_C17 =   Q6*XQH_C17**6 + Q5*XQH_C17**5 + Q4*XQH_C17**4 
     >            + Q3*XQH_C17**3 
     >            + Q2*XQH_C17**2 + Q1*XQH_C17**1 + Q0*XQH_C17**0 
          KQH_C17 =  PQH_C17*BQH_C17/(LENQ) 
          BM =  KQH_C17
      
        ELSEIF(LBL1(INIL:INIL+5) .EQ. 'QH_D17') THEN
          XQH_D17 = ABS(XQH                          ) 
          PQH_D17 =    (XQH                          ) / (XQH_D17+EPS) 
          BQH_D17 =   Q6*XQH_D17**6 + Q5*XQH_D17**5 + Q4*XQH_D17**4 
     >            + Q3*XQH_D17**3 
     >            + Q2*XQH_D17**2 + Q1*XQH_D17**1 + Q0*XQH_D17**0 
          KQH_D17 =  PQH_D17*BQH_D17/(LENQ) 
          BM =  KQH_D17
      
        ELSEIF(LBL1(INIL:INIL+5) .EQ. 'QH_E17') THEN
          E17_QUAD_I = A2      
          XQH_E17 = ABS(XQH - E17_QUAD_I + 0.D0 * XSB_CEIK17) 
          PQH_E17 =(XQH - E17_QUAD_I + 0.D0 * XSB_CEIK17) /(XQH_E17+EPS) 
          BQH_E17 =   Q6*XQH_E17**6 + Q5*XQH_E17**5 + Q4*XQH_E17**4 
     >            + Q3*XQH_E17**3 
     >            + Q2*XQH_E17**2 + Q1*XQH_E17**1 + Q0*XQH_E17**0 
          KQH_E17 =  PQH_E17*BQH_E17/(LENQ)
          BM =  KQH_E17
      
c       write(*,*)' agsks. xq a2 ',xqh,E17_QUAD_I,PQH_e17,BQH_e17,KQH_e17
c                   read(*,*)

        ELSEIF(LBL1(INIL:INIL+5) .EQ. 'QH_F17') THEN
          XQH_F17 = ABS(XQH        + 0.D0 * XSB_FHLB17) 
          PQH_F17 =    (XQH        + 0.D0 * XSB_FHLB17) / (XQH_F17+EPS) 
          BQH_F17 =   Q6*XQH_F17**6 + Q5*XQH_F17**5 + Q4*XQH_F17**4 
     >            + Q3*XQH_F17**3 
     >            + Q2*XQH_F17**2 + Q1*XQH_F17**1 + Q0*XQH_F17**0 
          KQH_F17 =  PQH_F17*BQH_F17/(LENQ) 
          BM =  KQH_F17
      
        ELSEIF(LBL1(INIL:INIL+5) .EQ. 'QH_G17') THEN
          XQH_G17 = ABS(XQH                          ) 
          PQH_G17 =    (XQH                          ) / (XQH_G17+EPS) 
          BQH_G17 =   Q6*XQH_G17**6 + Q5*XQH_G17**5 + Q4*XQH_G17**4 
     >            + Q3*XQH_G17**3 
     >            + Q2*XQH_G17**2 + Q1*XQH_G17**1 + Q0*XQH_G17**0 
          KQH_G17 =  PQH_G17*BQH_G17/(LENQ) 
          BM =  KQH_G17
      
        ELSEIF(LBL1(INIL:INIL+5) .EQ. 'QH_H17') THEN
          XQH_H17 = ABS(XQH        + 0.D0 * XSB_FHLB17) 
          PQH_H17 =    (XQH        + 0.D0 * XSB_FHLB17) / (XQH_H17+EPS) 
          BQH_H17 =   Q6*XQH_H17**6 + Q5*XQH_H17**5 + Q4*XQH_H17**4 
     >            + Q3*XQH_H17**3 
     >            + Q2*XQH_H17**2 + Q1*XQH_H17**1 + Q0*XQH_H17**0 
          KQH_H17 =  PQH_H17*BQH_H17/(LENQ) 
          BM =  KQH_H17
      
        ELSEIF(LBL1(INIL:INIL+5) .EQ. 'QH_I17') THEN
          XQH_I17 = ABS(XQH        - 0.D0 * XSB_CEIK17) 
          PQH_I17 =    (XQH        - 0.D0 * XSB_CEIK17) / (XQH_I17+EPS) 
          BQH_I17 =   Q6*XQH_I17**6 + Q5*XQH_I17**5 + Q4*XQH_I17**4 
     >            + Q3*XQH_I17**3 
     >            + Q2*XQH_I17**2 + Q1*XQH_I17**1 + Q0*XQH_I17**0 
          KQH_I17 =  PQH_I17*BQH_I17/(LENQ) 
          BM =  KQH_I17
      
        ELSEIF(LBL1(INIL:INIL+5) .EQ. 'QH_J17') THEN
          XQH_J17 = ABS(XQH                          ) 
          PQH_J17 =    (XQH                          ) / (XQH_J17+EPS) 
          BQH_J17 =   Q6*XQH_J17**6 + Q5*XQH_J17**5 + Q4*XQH_J17**4 
     >            + Q3*XQH_J17**3 
     >            + Q2*XQH_J17**2 + Q1*XQH_J17**1 + Q0*XQH_J17**0 
          KQH_J17 =  PQH_J17*BQH_J17/(LENQ) 
          BM =  KQH_J17
      
        ELSEIF(LBL1(INIL:INIL+5) .EQ. 'QH_K17') THEN
          XQH_K17 = ABS(XQH        - 0.D0 * XSB_CEIK17) 
          PQH_K17 =    (XQH        - 0.D0 * XSB_CEIK17) / (XQH_K17+EPS) 
          BQH_K17 =   Q6*XQH_K17**6 + Q5*XQH_K17**5 + Q4*XQH_K17**4 
     >            + Q3*XQH_K17**3 
     >            + Q2*XQH_K17**2 + Q1*XQH_K17**1 + Q0*XQH_K17**0 
          KQH_K17 =  PQH_K17*BQH_K17/(LENQ) 
          BM =  KQH_K17
      
        ELSEIF(LBL1(INIL:INIL+5) .EQ. 'QH_L17') THEN
          XQH_L17 = ABS(XQH        - 0.D0 * XSB_FHLB17)              
          PQH_L17 =    (XQH        - 0.D0 * XSB_FHLB17) / (XQH_L17+EPS) 
          BQH_L17 =   Q6*XQH_L17**6 + Q5*XQH_L17**5 + Q4*XQH_L17**4 
     >            + Q3*XQH_L17**3 
     >            + Q2*XQH_L17**2 + Q1*XQH_L17**1 + Q0*XQH_L17**0 
          KQH_L17 =  PQH_L17*BQH_L17/(LENQ) 
          BM =  KQH_L17
            
        ENDIF        

      ELSEIF(LBL1(INIL:INIL+2) .EQ. 'QV_') THEN
        LENQ = XL
        XQV = A1      

        XQVP2 = XQV*XQV 
        XQVP3 = XQV*XQVP2 
        XQVP4 = XQV*XQVP3 
        XQVP5 = XQV*XQVP4 
        XQVP6 = XQV*XQVP5 
      
      ! ==================================================================
      !  WE'RE GOING TO HAVE TO USE THE POLYNOMIAL COEFFICIENTS OVER AND OVER
      !    LATER IN THE FILE, SO I'M GOING TO DEFINE A SET THAT WE'LL USE 
      !    THROUGHOUT THE REST OF THE FILE.  IF ONE WANTS TO CHANGE WHICH POLYNOMIAL
      !    IS USED, ONE ONLY HAS TO CHANGE THE VALUES OF Q0-Q6
      ! ==================================================================
      
        Q6 = KQHC6 
        Q5 = KQHC5 
        Q4 = KQHC4 
        Q3 = KQHC3 
        Q2 = KQHC2 
        Q1 = KQHC1 
        Q0 = KQHC0 
      
        IF    (LBL1(INIL:INIL+5) .EQ. 'QV_A03') THEN
      !  VERTICAL TUNE QUAD CURRENTS (MADE POSITIVE FOR USE IN THE TRANSFER FUNC. POLYNOMIALS)
          XQV_A03 = ABS(XQV                          ) 
        !  VERTICAL TUNE QUAD POLARITY BITS = CURRENT/ABS(CURRENT)
          PQV_A03 =    (XQV                          ) / (XQV_A03+EPS) 
      !  NOW MAKE B STRENGTHS OUT OF THE CURRENTS USING THE TRANSFER FUNCTION (ABOVE)
          BQV_A03 =   Q6*XQV_A03**6 + Q5*XQV_A03**5 + Q4*XQV_A03**4 
     >            + Q3*XQV_A03**3 
     >            + Q2*XQV_A03**2 + Q1*XQV_A03**1 + Q0*XQV_A03**0 
        !  NOW CALCULATE K VALUES FROM THE BQV      
          KQV_A03 =  PQV_A03*BQV_A03/(LENQ) 
          BM =  KQV_A03
      
        ELSEIF(LBL1(INIL:INIL+5) .EQ. 'QV_B03') THEN
          B3_QUAD_I = A2    
          XQV_B03 = ABS(XQV        - 0.D0 * XSB_FHLB03) 
          PQV_B03 =    (XQV        - 0.D0 * XSB_FHLB03) / (XQV_B03+EPS) 
          BQV_B03 =   Q6*XQV_B03**6 + Q5*XQV_B03**5 + Q4*XQV_B03**4 
     >            + Q3*XQV_B03**3 
     >            + Q2*XQV_B03**2 + Q1*XQV_B03**1 + Q0*XQV_B03**0 
          KQV_B03 =  PQV_B03*BQV_B03/(LENQ) 
          BM =  KQV_B03
      
        ELSEIF(LBL1(INIL:INIL+5) .EQ. 'QV_C03') THEN
          XQV_C03 = ABS(XQV        + 0.D0 * XSB_CEIK03) 
          PQV_C03 =    (XQV        + 0.D0 * XSB_CEIK03) / (XQV_C03+EPS) 
          BQV_C03 =   Q6*XQV_C03**6 + Q5*XQV_C03**5 + Q4*XQV_C03**4 
     >            + Q3*XQV_C03**3 
     >            + Q2*XQV_C03**2 + Q1*XQV_C03**1 + Q0*XQV_C03**0 
          KQV_C03 =  PQV_C03*BQV_C03/(LENQ) 
          BM =  KQV_C03
      
        ELSEIF(LBL1(INIL:INIL+5) .EQ. 'QV_D03') THEN
          XQV_D03 = ABS(XQV                          ) 
          PQV_D03 =    (XQV                          ) / (XQV_D03+EPS) 
          BQV_D03 =   Q6*XQV_D03**6 + Q5*XQV_D03**5 + Q4*XQV_D03**4 
     >            + Q3*XQV_D03**3 
     >            + Q2*XQV_D03**2 + Q1*XQV_D03**1 + Q0*XQV_D03**0 
          KQV_D03 =  PQV_D03*BQV_D03/(LENQ) 
          BM =  KQV_D03
      
        ELSEIF(LBL1(INIL:INIL+5) .EQ. 'QV_E03') THEN
          XQV_E03 = ABS(XQV        + 0.D0 * XSB_CEIK03) 
          PQV_E03 =    (XQV        + 0.D0 * XSB_CEIK03) / (XQV_E03+EPS) 
          BQV_E03 =   Q6*XQV_E03**6 + Q5*XQV_E03**5 + Q4*XQV_E03**4 
     >            + Q3*XQV_E03**3 
     >            + Q2*XQV_E03**2 + Q1*XQV_E03**1 + Q0*XQV_E03**0 
          KQV_E03 =  PQV_E03*BQV_E03/(LENQ) 
          BM =  KQV_E03
      
        ELSEIF(LBL1(INIL:INIL+5) .EQ. 'QV_F03') THEN
          F3_QUAD_I = A2
          XQV_F03 = ABS(XQV        + 0.D0 * XSB_FHLB03) 
          PQV_F03 =    (XQV        + 0.D0 * XSB_FHLB03) / (XQV_F03+EPS) 
          BQV_F03 =   Q6*XQV_F03**6 + Q5*XQV_F03**5 + Q4*XQV_F03**4 
     >            + Q3*XQV_F03**3 
     >            + Q2*XQV_F03**2 + Q1*XQV_F03**1 + Q0*XQV_F03**0 
          KQV_F03 =  PQV_F03*BQV_F03/(LENQ) 
          BM =  KQV_F03
      
        ELSEIF(LBL1(INIL:INIL+5) .EQ. 'QV_G03') THEN
          XQV_G03 = ABS(XQV                          ) 
          PQV_G03 =    (XQV                          ) / (XQV_G03+EPS) 
          BQV_G03 =   Q6*XQV_G03**6 + Q5*XQV_G03**5 + Q4*XQV_G03**4 
     >            + Q3*XQV_G03**3 
     >            + Q2*XQV_G03**2 + Q1*XQV_G03**1 + Q0*XQV_G03**0 
          KQV_G03 =  PQV_G03*BQV_G03/(LENQ) 
          BM =  KQV_G03
      
        ELSEIF(LBL1(INIL:INIL+5) .EQ. 'QV_H03') THEN
          XQV_H03 = ABS(XQV        + 0.D0 * XSB_FHLB03) 
          PQV_H03 =    (XQV        + 0.D0 * XSB_FHLB03) / (XQV_H03+EPS) 
          BQV_H03 =   Q6*XQV_H03**6 + Q5*XQV_H03**5 + Q4*XQV_H03**4 
     >            + Q3*XQV_H03**3 
     >            + Q2*XQV_H03**2 + Q1*XQV_H03**1 + Q0*XQV_H03**0 
          KQV_H03 =  PQV_H03*BQV_H03/(LENQ) 
          BM =  KQV_H03
      
        ELSEIF(LBL1(INIL:INIL+5) .EQ. 'QV_I03') THEN
          XQV_I03 = ABS(XQV        - 0.D0 * XSB_CEIK03) 
          PQV_I03 =    (XQV        - 0.D0 * XSB_CEIK03) / (XQV_I03+EPS) 
          BQV_I03 =   Q6*XQV_I03**6 + Q5*XQV_I03**5 + Q4*XQV_I03**4 
     >            + Q3*XQV_I03**3 
     >            + Q2*XQV_I03**2 + Q1*XQV_I03**1 + Q0*XQV_I03**0 
          KQV_I03 =  PQV_I03*BQV_I03/(LENQ) 
          BM =  KQV_I03
      
        ELSEIF(LBL1(INIL:INIL+5) .EQ. 'QV_J03') THEN
          XQV_J03 = ABS(XQV                          ) 
          PQV_J03 =    (XQV                          ) / (XQV_J03+EPS) 
          BQV_J03 =   Q6*XQV_J03**6 + Q5*XQV_J03**5 + Q4*XQV_J03**4 
     >            + Q3*XQV_J03**3 
     >            + Q2*XQV_J03**2 + Q1*XQV_J03**1 + Q0*XQV_J03**0 
          KQV_J03 =  PQV_J03*BQV_J03/(LENQ) 
          BM =  KQV_J03
      
        ELSEIF(LBL1(INIL:INIL+5) .EQ. 'QV_K03') THEN
          XQV_K03 = ABS(XQV        - 0.D0 * XSB_CEIK03) 
          PQV_K03 =    (XQV        - 0.D0 * XSB_CEIK03) / (XQV_K03+EPS) 
          BQV_K03 =   Q6*XQV_K03**6 + Q5*XQV_K03**5 + Q4*XQV_K03**4 
     >            + Q3*XQV_K03**3 
     >            + Q2*XQV_K03**2 + Q1*XQV_K03**1 + Q0*XQV_K03**0 
          KQV_K03 =  PQV_K03*BQV_K03/(LENQ) 
          BM =  KQV_K03
      
        ELSEIF(LBL1(INIL:INIL+5) .EQ. 'QV_L03') THEN
          XQV_L03 = ABS(XQV        - 0.D0 * XSB_FHLB03)    
          PQV_L03 =    (XQV        - 0.D0 * XSB_FHLB03) / (XQV_L03+EPS)    
          BQV_L03 =   Q6*XQV_L03**6 + Q5*XQV_L03**5 + Q4*XQV_L03**4 
     >            + Q3*XQV_L03**3 
     >            + Q2*XQV_L03**2 + Q1*XQV_L03**1 + Q0*XQV_L03**0 
          KQV_L03 =  PQV_L03*BQV_L03/(LENQ) 
          BM =  KQV_L03
        
        ENDIF
      
            
      ELSEIF(LBL1(INIL:INIL+2) .EQ. 'QP_') THEN
      ! =====================================================================
      !  POLARIZATION QUADS
      !   THESE WERE WIRED TO THE VERTICAL TUNE QUAD STRING AS OF THE FY07 RUN
      !    THEY ARE IDENTICAL TO THE VERTICAL TUNE QUADS EXCEPT THAT THEY HAVE
      !    HALF THE NUMBER OF WINDINGS (HENCE THE FACTOR OF /2 IN CURRENT BELOW) 
      !    THIS MEANS THAT F03 AND B03 EACH HAVE TWO CURRENTS (SNAKE COMPENSATION
      !    AND VERTICAL TUNE CURRENT.
      ! =====================================================================
        LENQ = XL
        XQV = A1      
        IF    (LBL1(INIL:INIL+5) .EQ. 'QP_B03') THEN
          B3_QUAD_I = A2
        ELSEIF(LBL1(INIL:INIL+5) .EQ. 'QP_F03') THEN
          F3_QUAD_I = A2
        ENDIF

        IF    (LBL1(INIL:INIL+5) .EQ. 'QP_A03') THEN
          IPOL_A03 = ABS(XQV             ) /2.D0 
        !  POLARITY BITS (FACTOR OF TWO HERE ALSO...):
          PPOL_A03 = (XQV             ) / (IPOL_A03+EPS) /2.D0 
      !  SAME TRANSFER FUNC AS TUNE QUADS.      
          BPOL_A03 =   Q6*IPOL_A03**6 + Q5*IPOL_A03**5 + Q4*IPOL_A03**4 
     >            + Q3*IPOL_A03**3 
     >             + Q2*IPOL_A03**2 + Q1*IPOL_A03**1 + Q0*IPOL_A03**0 
        !   NOW GET THE K'S FROM THE B'S
          KPOL_A03 = PPOL_A03*BPOL_A03/(LENQ) 
          BM =  KPOL_A03 
      
        ELSEIF(LBL1(INIL:INIL+5) .EQ. 'QP_B03') THEN
          IPOL_B03 = ABS(XQV - B3_QUAD_I ) /2.D0  
          PPOL_B03 = (XQV - B3_QUAD_I ) / (IPOL_B03+EPS) /2.D0 
          BPOL_B03 =   Q6*IPOL_B03**6 + Q5*IPOL_B03**5 + Q4*IPOL_B03**4 
     >            + Q3*IPOL_B03**3 
     >             + Q2*IPOL_B03**2 + Q1*IPOL_B03**1 + Q0*IPOL_B03**0 
          KPOL_B03 = PPOL_B03*BPOL_B03/(LENQ) 
          BM =  KPOL_B03 
        
        ELSEIF(LBL1(INIL:INIL+5) .EQ. 'QP_C03') THEN
          IPOL_C03 = ABS(XQV             ) /2.D0 
          PPOL_C03 = (XQV             ) / (IPOL_C03+EPS) /2.D0 
          BPOL_C03 =   Q6*IPOL_C03**6 + Q5*IPOL_C03**5 + Q4*IPOL_C03**4 
     >            + Q3*IPOL_C03**3 
     >             + Q2*IPOL_C03**2 + Q1*IPOL_C03**1 + Q0*IPOL_C03**0 
          KPOL_C03 = PPOL_C03*BPOL_C03/(LENQ) 
          BM =  KPOL_C03 
        
        ELSEIF(LBL1(INIL:INIL+5) .EQ. 'QP_D03') THEN
          IPOL_D03 = ABS(XQV             ) /2.D0 
          PPOL_D03 = (XQV             ) / (IPOL_D03+EPS) /2.D0 
          BPOL_D03 =   Q6*IPOL_D03**6 + Q5*IPOL_D03**5 + Q4*IPOL_D03**4 
     >            + Q3*IPOL_D03**3 
     >             + Q2*IPOL_D03**2 + Q1*IPOL_D03**1 + Q0*IPOL_D03**0 
          KPOL_D03 = PPOL_D03*BPOL_D03/(LENQ) 
          BM =  KPOL_D03 
        
        ELSEIF(LBL1(INIL:INIL+5) .EQ. 'QP_E03') THEN
          IPOL_E03 = ABS(XQV             ) /2.D0 
          PPOL_E03 = (XQV             ) / (IPOL_E03+EPS) /2.D0 
          BPOL_E03 =   Q6*IPOL_E03**6 + Q5*IPOL_E03**5 + Q4*IPOL_E03**4 
     >            + Q3*IPOL_E03**3 
     >             + Q2*IPOL_E03**2 + Q1*IPOL_E03**1 + Q0*IPOL_E03**0 
          KPOL_E03 = PPOL_E03*BPOL_E03/(LENQ) 
          BM =  KPOL_E03 
        
        ELSEIF(LBL1(INIL:INIL+5) .EQ. 'QP_F03') THEN
          IPOL_F03 = ABS(XQV - F3_QUAD_I ) /2.D0 
          PPOL_F03 = (XQV - F3_QUAD_I ) / (IPOL_F03+EPS) /2.D0 
          BPOL_F03 =   Q6*IPOL_F03**6 + Q5*IPOL_F03**5 + Q4*IPOL_F03**4 
     >            + Q3*IPOL_F03**3 
     >             + Q2*IPOL_F03**2 + Q1*IPOL_F03**1 + Q0*IPOL_F03**0 
          KPOL_F03 = PPOL_F03*BPOL_F03/(LENQ) 
          BM =  KPOL_F03 
        
        ELSEIF(LBL1(INIL:INIL+5) .EQ. 'QP_G03') THEN
          IPOL_G03 = ABS(XQV             ) /2.D0 
          PPOL_G03 = (XQV             ) / (IPOL_G03+EPS) /2.D0 
          BPOL_G03 =   Q6*IPOL_G03**6 + Q5*IPOL_G03**5 + Q4*IPOL_G03**4 
     >            + Q3*IPOL_G03**3 
     >             + Q2*IPOL_G03**2 + Q1*IPOL_G03**1 + Q0*IPOL_G03**0 
          KPOL_G03 = PPOL_G03*BPOL_G03/(LENQ) 
          BM =  KPOL_G03 
        
        ELSEIF(LBL1(INIL:INIL+5) .EQ. 'QP_H03') THEN
          IPOL_H03 = ABS(XQV             ) /2.D0 
          PPOL_H03 = (XQV             ) / (IPOL_H03+EPS) /2.D0 
          BPOL_H03 =   Q6*IPOL_H03**6 + Q5*IPOL_H03**5 + Q4*IPOL_H03**4 
     >            + Q3*IPOL_H03**3 
     >             + Q2*IPOL_H03**2 + Q1*IPOL_H03**1 + Q0*IPOL_H03**0 
          KPOL_H03 = PPOL_H03*BPOL_H03/(LENQ) 
          BM =  KPOL_H03 
        
        ELSEIF(LBL1(INIL:INIL+5) .EQ. 'QP_I03') THEN
          IPOL_I03 = ABS(XQV             ) /2.D0 
          PPOL_I03 = (XQV             ) / (IPOL_I03+EPS) /2.D0 
          BPOL_I03 =   Q6*IPOL_I03**6 + Q5*IPOL_I03**5 + Q4*IPOL_I03**4 
     >            + Q3*IPOL_I03**3 
     >             + Q2*IPOL_I03**2 + Q1*IPOL_I03**1 + Q0*IPOL_I03**0 
          KPOL_I03 = PPOL_I03*BPOL_I03/(LENQ) 
          BM =  KPOL_I03 
        
        ELSEIF(LBL1(INIL:INIL+5) .EQ. 'QP_J03') THEN
          IPOL_J03 = ABS(XQV             ) /2.D0 
          PPOL_J03 = (XQV             ) / (IPOL_J03+EPS) /2.D0 
          BPOL_J03 =   Q6*IPOL_J03**6 + Q5*IPOL_J03**5 + Q4*IPOL_J03**4 
     >            + Q3*IPOL_J03**3 
     >             + Q2*IPOL_J03**2 + Q1*IPOL_J03**1 + Q0*IPOL_J03**0 
          KPOL_J03 = PPOL_J03*BPOL_J03/(LENQ) 
          BM =  KPOL_J03 
        
        ELSEIF(LBL1(INIL:INIL+5) .EQ. 'QP_K03') THEN
          IPOL_K03 = ABS(XQV             ) /2.D0 
          PPOL_K03 = (XQV             ) / (IPOL_K03+EPS) /2.D0 
          BPOL_K03 =   Q6*IPOL_K03**6 + Q5*IPOL_K03**5 + Q4*IPOL_K03**4 
     >            + Q3*IPOL_K03**3 
     >             + Q2*IPOL_K03**2 + Q1*IPOL_K03**1 + Q0*IPOL_K03**0 
          KPOL_K03 = PPOL_K03*BPOL_K03/(LENQ) 
          BM =  KPOL_K03 
        
        ELSEIF(LBL1(INIL:INIL+5) .EQ. 'QP_L03') THEN
          IPOL_L03 = ABS(XQV             ) /2.D0 
          PPOL_L03 = (XQV             ) / (IPOL_L03+EPS) /2.D0 
          BPOL_L03 =   Q6*IPOL_L03**6 + Q5*IPOL_L03**5 + Q4*IPOL_L03**4 
     >            + Q3*IPOL_L03**3 
     >             + Q2*IPOL_L03**2 + Q1*IPOL_L03**1 + Q0*IPOL_L03**0 
          KPOL_L03 = PPOL_L03*BPOL_L03/(LENQ) 
          BM =  KPOL_L03 
        
        ENDIF

        

      ELSEIF(LBL1(INIL:INIL+5) .EQ. 'QTHIN_') THEN
        LQTHIN = XL

      ! ===================================================================
      !  COLD SNAKE COMPENSATION QUADRUPOLES.
      ! 
      !  THE A17 TUNE QUAD HAS A FLOATING POWER SUPPLY WIRED UP TO ADD COMPENSATION
      !    FOR THE COLD SNAKE. HAVE TO CALCULATE STRENGTH SEPARATE FROM OTHER TUNE
      !      QUADS
      ! 
      !  THE B03 POLARIZATION QUAD IS POWERED BY TWO SEPARATE SUPPLIES.
      !       SINCE 09/06 THE POLARIZATION QUADS HAVE BEEN WIRED IN SERIES WITH
      !       THE VERTICAL TUNE QUADS.  SINCE 09/05, B03 HAS HAD IT'S OWN INDIVIDUAL SUPPLY.
      !  A19 AND B01 THIN QUADS HAVE THEIR OWN SUPPLIES AND TRANSFER FUNCTIONS 
      ! ==============================================================
        IF    (LBL1(INIL:INIL+8) .EQ. 'QTHIN_A19') THEN
          A19_QUAD_I = A1      
          XA19 = A19_QUAD_I

          KTHIN_A19_RAW =   A19C6*(XA19*XA19*XA19*XA19*XA19*XA19)   +
     >                   A19C5*(XA19*XA19*XA19*XA19*XA19)      +
     >                   A19C4*(XA19*XA19*XA19*XA19)             +
     >                   A19C3*(XA19*XA19*XA19)                   +
     >                   A19C2*(XA19*XA19)                   +
     >                   A19C1*(XA19)                         +
     >                   A19C0                                
      
          KTHIN_A19 = -1*KTHIN_A19_RAW/(LQTHIN)    !  DEFOCUSING.
          BM =  KTHIN_A19 

        ELSEIF(LBL1(INIL:INIL+8) .EQ. 'QTHIN_B01') THEN
          B1_QUAD_I = A1      
          XB1 = B1_QUAD_I 
      
          KTHIN_B01_RAW =    B1C6*(XB1*XB1*XB1*XB1*XB1*XB1)    +
     >                   B1C5*(XB1*XB1*XB1*XB1*XB1)      +
     >                   B1C4*(XB1*XB1*XB1*XB1)             +
     >                   B1C3*(XB1*XB1*XB1)             +
     >                   B1C2*(XB1*XB1)                   +
     >                   B1C1*(XB1)                   +
     >                   B1C0                         

          KTHIN_B01 = -1*KTHIN_B01_RAW/(LQTHIN)    !  DEFOCUSING.
          BM =  KTHIN_B01 


      ! ===============================================================
      !  WARM SNAKE COMPENSATION
      ! ===============================================================
      
      !  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      !  NOTE:  E19 AND F1 THIN QUADS ARE WIRED IN SERIES BUT WITH 
      !              OPPOSITE POLARITIES. E19 IS FOCUSING,F01 IS DEFOCUSING.
      
        ELSEIF(LBL1(INIL:INIL+8) .EQ. 'QTHIN_E19') THEN
          E19_F1_QUAD_I = A1      
          XE19 = E19_F1_QUAD_I 
      
          KTHIN_E19_RAW = E19C6*(XE19*XE19*XE19*XE19*XE19*XE19)  +
     >                E19C5*(XE19*XE19*XE19*XE19*XE19)      +
     >                E19C4*(XE19*XE19*XE19*XE19)             +
     >                E19C3*(XE19*XE19*XE19)               +
     >                E19C2*(XE19*XE19)                   +
     >                E19C1*(XE19)                         +
     >                E19C0                                
      
          KTHIN_E19 = 1.D0* KTHIN_E19_RAW / (LQTHIN) 
          BM =  KTHIN_E19 
            
        ELSEIF(LBL1(INIL:INIL+8) .EQ. 'QTHIN_F01') THEN
          E19_F1_QUAD_I = A1       ! E19 AND F1 THIN QUADS ARE WIRED IN SERIES 
          XF1 = E19_F1_QUAD_I 
      
          KTHIN_F01_RAW =       F1C6*(XF1*XF1*XF1*XF1*XF1*XF1)       +
     >             F1C5*(XF1*XF1*XF1*XF1*XF1)      +
     >             F1C4*(XF1*XF1*XF1*XF1)             +
     >             F1C3*(XF1*XF1*XF1)             +
     >             F1C2*(XF1*XF1)                   +
     >             F1C1*(XF1)                   +
     >            F1C0
      
          KTHIN_F01 = -1.D0* KTHIN_F01_RAW/ (LQTHIN) 
          BM =  KTHIN_F01 

        ENDIF

      ENDIF


      RETURN
      END

