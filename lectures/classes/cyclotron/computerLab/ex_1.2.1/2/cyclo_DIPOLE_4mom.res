Cyclotron, classical.
 'OBJET'                                                                                                      1
64.62444403717985                                   ! 200keV
2
4 1
12.9248888074 0.  0.   0.   0.  1. 'm'              ! 200keV. R=Brho/B=*/.5
28.9070891209 0.  0.   0.   0.  2.23654451125 'm'   ! 1 MeV. R=Brho/B=*/.5
50.           0.  0.   0.   0.  3.86850523397 'o'   ! at RM  (B*rho=0.5*0.5=0.25T.m, 2.9885 MeV)
64.7070336799 0.  0.   0.   0.  5.0063899693  'M'   ! 5 MeV. R=Brho/B=*/.5
1 1 1 1
 
 'DIPOLE'                                                                                                     2
0
60. 50.
30.  5. 0. 0. 0.
0.  0.                                        ! EFB 1  hard-edge
4  .1455   2.2670  -.6395  1.1558  0. 0.  0.
30. 0.  1.E6  -1.E6  1.E6  1.E6
0.   0.                                       ! EFB 2
4  .1455   2.2670  -.6395  1.1558  0. 0.  0.
-30. 0.  1.E6  -1.E6  1.E6  1.E6
0. 0.                                         ! EFB 3
0  0.      0.      0.      0.      0. 0.  0.
0. 0.  1.E6  -1.E6  1.E6  1.E6 0.
4   10.
0.03                 ! 0.03 rsults in Y==Y0 at all E. If .ne.0.03, IT=1 or IT=2 do not close.
2  0. 0. 0. 0.       ! could also be, e.g., 2 50. 0. 50. 0. with Y0 amended accordingly in OBJET
 'FAISCEAU'                                                                                                   3
 'FAISTORE'                                                                                                   4
zgoubi.fai
1
 'END'                                                                                                        5

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 

                          MAGNETIC  RIGIDITY =         64.624 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       4 POINTS 



************************************************************************************************************************************
      2  Keyword, label(s) :  DIPOLE                                                

                    Dipole  magnet

           ANGLES : A.TOTAL =   60.00     degrees     A.CENTRAL =   30.00     degrees
           RM =   50.00     cm
           HNORM =   5.000     kGauss     COEF.N =   0.000         COEF.B =   0.000         COEF.G=   0.000    

     Entrance  EFB 
          Fringe  field  : LAMBDA =   0.00 CM     QSI=   0.00
           COEFFICIENTS :  4   0.14550   2.26700  -0.63950   1.15580   0.00000   0.00000
           Shift  of  EFB  =    0.000     CM

          OMEGA =  30.00 deg.     Wedge  angle  =   0.00 deg.
           Radius 1 =  1.00E+06 CM
           Straight  segment 1 = -1.00E+06 CM
           Straight  segment 2 =  1.00E+06 CM
           Radius 2 =  1.00E+06 CM

     Exit  EFB     
          Fringe  field  : LAMBDA =   0.00 CM     QSI=   0.00
           COEFFICIENTS :  4   0.14550   2.26700  -0.63950   1.15580   0.00000   0.00000
           Shift  of  EFB  =    0.000     CM

          OMEGA = -30.00 deg.     Wedge  angle  =   0.00 deg.
           Radius 1 =  1.00E+06 CM
           Straight  segment 1 = -1.00E+06 CM
           Straight  segment 2 =  1.00E+06 CM
           Radius 2 =  1.00E+06 CM

     Lateral  EFB  
          Fringe  field  : LAMBDA =   0.00 CM     QSI=   0.00
           COEFFICIENTS :  0   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000
           Shift  of  EFB  =    0.000     CM

                     Face centred on direction ACENT+OMEGA, A    0.000     CM

          OMEGA =   0.00 deg.     Wedge  angle  =   0.00 deg.
           Radius 1 =  1.00E+06 CM
           Straight  segment 1 = -1.00E+06 CM
           Straight  segment 2 =  1.00E+06 CM
           Radius 2 =  1.00E+06 CM

                     Interpolation  option : 4
                    25-point  interpolation, size of flying mesh :   STEP /   10.0

                    Integration step :  3.0000E-02 cm   (i.e.,   6.0000E-04 rad  at mean radius RM =    50.00    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad

  A    1  1.0000    12.925     0.000     0.000     0.000            1.047    12.925     0.000     0.000     0.000            1
  A    1  2.2365    28.907     0.000     0.000     0.000            1.047    28.907    -0.000     0.000     0.000            2
  A    1  3.8685    50.000     0.000     0.000     0.000            1.047    50.000     0.000     0.000     0.000            3
  A    1  5.0064    64.707     0.000     0.000     0.000            1.047    64.707    -0.000     0.000     0.000            4

 Cumulative length of optical axis =   0.523598776     m ;  Time  (for ref. rigidity & particle) =   8.460221E-08 s 

************************************************************************************************************************************
      3  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      2)
                                                  4 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    1.0000   12.925    0.000    0.000    0.000   1.353491E+01     1
m  1   2.2365    28.907     0.000     0.000     0.000       0.0000    2.2365   28.907   -0.000    0.000    0.000   3.027143E+01     2
o  1   3.8685    50.000     0.000     0.000     0.000       0.0000    3.8685   50.000    0.000    0.000    0.000   5.235988E+01     3
M  1   5.0064    64.707     0.000     0.000     0.000       0.0000    5.0064   64.707   -0.000    0.000    0.000   6.776105E+01     4


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   1.1409E-12   4.2309E-02   1.0764E+11   3.913475E-01  -4.526658E-13        4        4    1.000      (Y,T)         1
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        4        4    1.000      (Z,P)         1
          NaN   0.0000E+00   1.0000E+00   1.367006E-03   5.866152E+01        4        0    0.000      (t,K)         1

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   1.977092E-01
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   1.838416E-12

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =            NaN
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =            NaN


  Beam  sigma  matrix : 

   3.908893E-02  -1.536432E-14   0.000000E+00   0.000000E+00
  -1.536432E-14   3.379774E-24   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     3.631469E-13    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      4  Keyword, label(s) :  FAISTORE                                              


                OPEN FILE zgoubi.fai                                                                      
                FOR PRINTING COORDINATES 

               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      5  Keyword, label(s) :  END                                                   


                             4 particles have been launched
                     Made  it  to  the  end :      4

************************************************************************************************************************************

          Pgm zgoubi : Execution ended upon key  END       

************************************************************************************************************************************
   
            ZGOUBI RUN COMPLETED. 

  Zgoubi, author's dvlpmnt version.
  Job  started  on  19-05-0017,  at  06:18:05 
  JOB  ENDED  ON    19-05-0017,  AT  06:18:05 

   CPU time, total :    4.800200000000000E-002
