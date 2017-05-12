Cyclotron, classical.
 'OBJET'                                                                                                      1
64.62444403717985                                   ! 200keV
2
1 1
28.9070891209 0.  0.   0.   0.  2.23654451125 'm'   ! 1 MeV. R=Brho/B=*/.5
1
 'FAISCEAU'                                                                                                   2
 
 'DIPOLE'                                                                                                     3
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
 'FAISCEAU'                                                                                                   4
 
 'FIT'                                                                                                        5
1   nofinal
1 30 0 [28.,65.]
2
3.1 1 2 #End 0. 1. 0
3.1 1 3 #End 0. 1. 0
 
 'FAISTORE'                                                                                                   6
zgoubi.fai
1
 
 'REBELOTE'                                                                                                   7
20 0.2 99 1
1
OBJET 35 2.23654451125:5.0063899693
 
 'END'                                                                                                        8

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 

                          MAGNETIC  RIGIDITY =         64.624 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       1 POINTS 



************************************************************************************************************************************
      2  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      1)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

m  1   2.2365    28.907     0.000     0.000     0.000       0.0000    2.2365   28.907    0.000    0.000    0.000   0.000000E+00     1


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   2.890709E-01   0.000000E+00        1        1    1.000      (Y,T)         1
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)         1
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   4.333064E+01        1        1    1.000      (t,K)         1

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      3  Keyword, label(s) :  DIPOLE                                                

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

  A    1  2.2365    28.907     0.000     0.000     0.000            1.047    28.907    -0.000     0.000     0.000            1

 Cumulative length of optical axis =   0.523598776     m ;  Time  (for ref. rigidity & particle) =   8.460221E-08 s 

************************************************************************************************************************************
      4  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      3)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

m  1   2.2365    28.907     0.000     0.000     0.000       0.0000    2.2365   28.907   -0.000    0.000    0.000   3.027143E+01     1


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   2.890709E-01  -3.406018E-12        1        1    1.000      (Y,T)         1
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)         1
   0.0000E+00   0.0000E+00   1.0000E+00   1.009746E-03   4.333064E+01        1        1    1.000      (t,K)         1

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      5  Keyword, label(s) :  FIT                                                   

     FIT procedure launched. Method is 1

           variable #            1       IR =            1 ,   ok.
           variable #            1       IP =           30 ,   ok.
           constraint #            1       IR =            4 ,   ok.
           constraint #            1       I  =            1 ,   ok.
           constraint #            2       IR =            4 ,   ok.
           constraint #            2       I  =            1 ,   ok.

                    FIT  variables  and  constraints  in  good  order,  FIT  will proceed. 

                    Final FIT status will NOT be saved. For so, use the 'save [FileName]' command

 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    28.0        28.9       28.907089       65.0      0.123      OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    5.730172E-11    2.83E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    3.406018E-09    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   1.1604E-17

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


                OPEN FILE zgoubi.fai                                                                      
                FOR PRINTING COORDINATES 

               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        1 through the optical structure 

                     Total of          1 particles have been launched

     Multiple pass, 
          from element #     1 : OBJET     /label1=                    /label2=                    
                             to  REBELOTE  /label1=                    /label2=                    
     ending at pass #      21 at element #     7 : REBELOTE  /label1=                    /label2=                    


     Parameter #   35 in element #    1 will be modified at each pass. 
     list of requested parameter values :
                   1   2.23654451E+00
                   2   2.38232585E+00
                   3   2.52810719E+00
                   4   2.67388853E+00
                   5   2.81966987E+00
                   6   2.96545121E+00
                   7   3.11123255E+00
                   8   3.25701389E+00
                   9   3.40279523E+00
                  10   3.54857657E+00
                  11   3.69435791E+00
                  12   3.84013925E+00
                  13   3.98592059E+00
                  14   4.13170193E+00
                  15   4.27748327E+00
                  16   4.42326461E+00
                  17   4.56904595E+00
                  18   4.71482729E+00
                  19   4.86060863E+00
                  20   5.00638997E+00

 Pgm rebel. At pass #    1/  21.  In element #    1,  parameter # 35  changed to    2.23654451E+00   (was    2.23654451E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.


  *** ERROR in object definition : 
    momentum value 0 found.
    Only  non-zero  dp/p  allowed...

 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    28.0        28.9       28.907089       65.0      0.123      OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
 Fit reached penalty value   0.0000E+00

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        2 through the optical structure 

                     Total of          1 particles have been launched

 Pgm rebel. At pass #    2/  21.  In element #    1,  parameter # 35  changed to    2.38232585E+00   (was    2.23654451E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    28.0        28.9       28.907089       65.0      0.123      OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
 Fit reached penalty value   0.0000E+00

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        3 through the optical structure 

                     Total of          1 particles have been launched

 Pgm rebel. At pass #    3/  21.  In element #    1,  parameter # 35  changed to    2.52810719E+00   (was    2.38232585E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    28.0        28.9       28.907089       65.0      0.123      OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
 Fit reached penalty value   0.0000E+00

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        4 through the optical structure 

                     Total of          1 particles have been launched

 Pgm rebel. At pass #    4/  21.  In element #    1,  parameter # 35  changed to    2.67388853E+00   (was    2.52810719E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    28.0        28.9       28.907089       65.0      0.123      OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
 Fit reached penalty value   0.0000E+00

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        5 through the optical structure 

                     Total of          1 particles have been launched

 Pgm rebel. At pass #    5/  21.  In element #    1,  parameter # 35  changed to    2.81966987E+00   (was    2.67388853E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    28.0        28.9       28.907089       65.0      0.123      OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
 Fit reached penalty value   0.0000E+00

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        6 through the optical structure 

                     Total of          1 particles have been launched

 Pgm rebel. At pass #    6/  21.  In element #    1,  parameter # 35  changed to    2.96545121E+00   (was    2.81966987E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    28.0        28.9       28.907089       65.0      0.123      OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
 Fit reached penalty value   0.0000E+00

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        7 through the optical structure 

                     Total of          1 particles have been launched

 Pgm rebel. At pass #    7/  21.  In element #    1,  parameter # 35  changed to    3.11123255E+00   (was    2.96545121E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    28.0        28.9       28.907089       65.0      0.123      OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
 Fit reached penalty value   0.0000E+00

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        8 through the optical structure 

                     Total of          1 particles have been launched

 Pgm rebel. At pass #    8/  21.  In element #    1,  parameter # 35  changed to    3.25701389E+00   (was    3.11123255E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    28.0        28.9       28.907089       65.0      0.123      OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
 Fit reached penalty value   0.0000E+00

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        9 through the optical structure 

                     Total of          1 particles have been launched

 Pgm rebel. At pass #    9/  21.  In element #    1,  parameter # 35  changed to    3.40279523E+00   (was    3.25701389E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    28.0        28.9       28.907089       65.0      0.123      OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
 Fit reached penalty value   0.0000E+00

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       10 through the optical structure 

                     Total of          1 particles have been launched

 Pgm rebel. At pass #   10/  21.  In element #    1,  parameter # 35  changed to    3.54857657E+00   (was    3.40279523E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    28.0        28.9       28.907089       65.0      0.123      OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
 Fit reached penalty value   0.0000E+00

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       11 through the optical structure 

                     Total of          1 particles have been launched

 Pgm rebel. At pass #   11/  21.  In element #    1,  parameter # 35  changed to    3.69435791E+00   (was    3.54857657E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    28.0        28.9       28.907089       65.0      0.123      OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
 Fit reached penalty value   0.0000E+00

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       12 through the optical structure 

                     Total of          1 particles have been launched

 Pgm rebel. At pass #   12/  21.  In element #    1,  parameter # 35  changed to    3.84013925E+00   (was    3.69435791E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    28.0        28.9       28.907089       65.0      0.123      OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
 Fit reached penalty value   0.0000E+00

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       13 through the optical structure 

                     Total of          1 particles have been launched

 Pgm rebel. At pass #   13/  21.  In element #    1,  parameter # 35  changed to    3.98592059E+00   (was    3.84013925E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    28.0        28.9       28.907089       65.0      0.123      OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
 Fit reached penalty value   0.0000E+00

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       14 through the optical structure 

                     Total of          1 particles have been launched

 Pgm rebel. At pass #   14/  21.  In element #    1,  parameter # 35  changed to    4.13170193E+00   (was    3.98592059E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    28.0        28.9       28.907089       65.0      0.123      OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
 Fit reached penalty value   0.0000E+00

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       15 through the optical structure 

                     Total of          1 particles have been launched

 Pgm rebel. At pass #   15/  21.  In element #    1,  parameter # 35  changed to    4.27748327E+00   (was    4.13170193E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    28.0        28.9       28.907089       65.0      0.123      OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
 Fit reached penalty value   0.0000E+00

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       16 through the optical structure 

                     Total of          1 particles have been launched

 Pgm rebel. At pass #   16/  21.  In element #    1,  parameter # 35  changed to    4.42326461E+00   (was    4.27748327E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    28.0        28.9       28.907089       65.0      0.123      OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
 Fit reached penalty value   0.0000E+00

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       17 through the optical structure 

                     Total of          1 particles have been launched

 Pgm rebel. At pass #   17/  21.  In element #    1,  parameter # 35  changed to    4.56904595E+00   (was    4.42326461E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    28.0        28.9       28.907089       65.0      0.123      OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
 Fit reached penalty value   0.0000E+00

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       18 through the optical structure 

                     Total of          1 particles have been launched

 Pgm rebel. At pass #   18/  21.  In element #    1,  parameter # 35  changed to    4.71482729E+00   (was    4.56904595E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    28.0        28.9       28.907089       65.0      0.123      OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
 Fit reached penalty value   0.0000E+00

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       19 through the optical structure 

                     Total of          1 particles have been launched

 Pgm rebel. At pass #   19/  21.  In element #    1,  parameter # 35  changed to    4.86060863E+00   (was    4.71482729E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    28.0        28.9       28.907089       65.0      0.123      OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
 Fit reached penalty value   0.0000E+00

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       20 through the optical structure 

                     Total of          1 particles have been launched


      Next  pass  is  #    21 and  last  pass  through  the  optical  structure


 Pgm rebel. At pass #   20/  21.  In element #    1,  parameter # 35  changed to    5.00638997E+00   (was    4.86060863E+00)

     WRITE statements to zgoubi.res are re-established from now on.

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      1)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

m -6   0.0000     0.000     0.000     0.000     0.000       0.0000    0.0000    0.000    0.000    0.000    0.000   0.000000E+00     1


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

          NaN   0.0000E+00   1.0000E+00            NaN            NaN        0        0    0.000      (Y,T)        21
          NaN   0.0000E+00   1.0000E+00            NaN            NaN        0        0    0.000      (Z,P)        21
          NaN   0.0000E+00   1.0000E+00            NaN            NaN        0        0    0.000      (t,K)        21

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =            NaN
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =            NaN

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =            NaN
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =            NaN

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =            NaN
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =            NaN


  Beam  sigma  matrix : 

            NaN            NaN            NaN            NaN
            NaN            NaN            NaN            NaN
            NaN            NaN            NaN            NaN
            NaN            NaN            NaN            NaN

      sqrt(det_Y), sqrt(det_Z) :              NaN             NaN    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      3  Keyword, label(s) :  DIPOLE                                                

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


 Cumulative length of optical axis =   0.523598776     m ;  Time  (for ref. rigidity & particle) =   1.692044E-07 s 

************************************************************************************************************************************
      4  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      3)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

m -6   0.0000     0.000     0.000     0.000     0.000       0.0000    0.0000    0.000    0.000    0.000    0.000   0.000000E+00     1


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

          NaN   0.0000E+00   1.0000E+00            NaN            NaN        0        0    0.000      (Y,T)        21
          NaN   0.0000E+00   1.0000E+00            NaN            NaN        0        0    0.000      (Z,P)        21
          NaN   0.0000E+00   1.0000E+00            NaN            NaN        0        0    0.000      (t,K)        21

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =            NaN
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =            NaN

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =            NaN
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =            NaN

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =            NaN
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =            NaN


  Beam  sigma  matrix : 

            NaN            NaN            NaN            NaN
            NaN            NaN            NaN            NaN
            NaN            NaN            NaN            NaN
            NaN            NaN            NaN            NaN

      sqrt(det_Y), sqrt(det_Z) :              NaN             NaN    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      5  Keyword, label(s) :  FIT                                                   

     FIT procedure launched. Method is 1


                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      1)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

m -6   0.0000     0.000     0.000     0.000     0.000       0.0000    0.0000    0.000    0.000    0.000    0.000   0.000000E+00     1


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

          NaN   0.0000E+00   1.0000E+00            NaN            NaN        0        0    0.000      (Y,T)        21
          NaN   0.0000E+00   1.0000E+00            NaN            NaN        0        0    0.000      (Z,P)        21
          NaN   0.0000E+00   1.0000E+00            NaN            NaN        0        0    0.000      (t,K)        21

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =            NaN
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =            NaN

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =            NaN
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =            NaN

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =            NaN
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =            NaN


  Beam  sigma  matrix : 

            NaN            NaN            NaN            NaN
            NaN            NaN            NaN            NaN
            NaN            NaN            NaN            NaN
            NaN            NaN            NaN            NaN

      sqrt(det_Y), sqrt(det_Z) :              NaN             NaN    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      3  Keyword, label(s) :  DIPOLE                                                

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


 Cumulative length of optical axis =   0.523598776     m ;  Time  (for ref. rigidity & particle) =   8.460221E-08 s 

************************************************************************************************************************************
      4  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      3)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

m -6   0.0000     0.000     0.000     0.000     0.000       0.0000    0.0000    0.000    0.000    0.000    0.000   0.000000E+00     1


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

          NaN   0.0000E+00   1.0000E+00            NaN            NaN        0        0    0.000      (Y,T)        21
          NaN   0.0000E+00   1.0000E+00            NaN            NaN        0        0    0.000      (Z,P)        21
          NaN   0.0000E+00   1.0000E+00            NaN            NaN        0        0    0.000      (t,K)        21

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =            NaN
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =            NaN

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =            NaN
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =            NaN

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =            NaN
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =            NaN


  Beam  sigma  matrix : 

            NaN            NaN            NaN            NaN
            NaN            NaN            NaN            NaN
            NaN            NaN            NaN            NaN
            NaN            NaN            NaN            NaN

      sqrt(det_Y), sqrt(det_Z) :              NaN             NaN    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      5  Keyword, label(s) :  FIT                                                   

     FIT procedure launched. Method is 1

           variable #            1       IR =            1 ,   ok.
           variable #            1       IP =           30 ,   ok.
           constraint #            1       IR =            4 ,   ok.
           constraint #            1       I  =            1 ,   ok.
           constraint #            2       IR =            4 ,   ok.
           constraint #            2       I  =            1 ,   ok.

                    FIT  variables  and  constraints  in  good  order,  FIT  will proceed. 

                    Final FIT status will NOT be saved. For so, use the 'save [FileName]' command

 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    28.0        28.9       28.907089       65.0      0.123      OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    0.000000E+00         NaN FAISCEAU   -                    -                    0
 Fit reached penalty value   0.0000E+00

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                         ****  End  of  'REBELOTE'  procedure  ****

      There  has  been         21  passes  through  the  optical  structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      8  Keyword, label(s) :  END                                                   


************************************************************************************************************************************

          Pgm zgoubi : Execution ended upon key  END       

************************************************************************************************************************************
   
            ZGOUBI RUN COMPLETED. 

  Zgoubi, author's dvlpmnt version.
  Job  started  on  11-05-0017,  at  21:17:20 
  JOB  ENDED  ON    11-05-0017,  AT  21:17:31 

   CPU time, total :    0.228014000000000     
