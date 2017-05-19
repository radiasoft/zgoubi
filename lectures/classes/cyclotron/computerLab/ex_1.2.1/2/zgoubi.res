Cyclotron, classical.
 'OBJET'                                                                                                      1
64.62444403717985                                   ! 200keV
2
1 1
12.9248888074 0.  0.   0.   0.  1. 'm'              ! 200keV. R=Brho/B=*/.5
1
 
 'PARTICUL'   ! This is required to get time-of-flight. Otherwise not: zgoubi only needs rigidity.            2
PROTON
 
 'FAISCEAU'                                                                                                   3
 
 'DIPOLE'                                                                                                     4
0
60. 50.                                       ! AT, RM
30.  5. 0. 0. 0.                              ! ACN, B0, N, B, G
0.  0.                                        ! EFB 1  hard-edge (can only work with sector magnet))
4  .1455   2.2670  -.6395  1.1558  0. 0.  0.
30. 0.  1.E6  -1.E6  1.E6  1.E6
0.   0.                                       ! EFB 2
4  .1455   2.2670  -.6395  1.1558  0. 0.  0.
-30. 0.  1.E6  -1.E6  1.E6  1.E6
0. 0.                                         ! EFB 3
0  0.      0.      0.      0.      0. 0.  0.
0. 0.  1.E6  -1.E6  1.E6  1.E6 0.
4   10.
0.01                 ! 0.03 rsults in Y==Y0 at all E. If .ne.0.03, IT=1 or IT=2 do not close.
2  0. 0. 0. 0.       ! could also be, e.g., 2 50. 0. 50. 0. with Y0 amended accordingly in OBJET
 'FAISCEAU'                                                                                                   5
 
 'FIT'                                                                                                        6
2   nofinal
1 30 0 [12.,65.]
4 60 0 2.
2
3.1 1 2 #End 0. 1. 0
3.1 1 3 #End 0. 1. 0
 
 'FAISTORE'                                                                                                   7
zgoubi.fai
1
 
 'REBELOTE'                                                                                                   8
20 0.2  0 1        ! 20 different rigidities ;  ; coordinates as found in OBJET ; change parameter(s)
1
OBJET 35     1:5.0063899693    ! 0.2 MeV to 5 MeV
 
 'END'                                                                                                        9

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 

                          MAGNETIC  RIGIDITY =         64.624 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       1 POINTS 



************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


     Particle  properties :
     PROTON
                     Mass          =    938.272        MeV/c2
                     Charge        =   1.602176E-19    C     
                     G  factor     =    1.79285              
                     COM life-time =   1.000000E+99    s     

              Reference  data :
                    mag. rigidity (kG.cm)   :   64.624444      =p/q, such that dev.=B*L/rigidity
                    mass (MeV/c2)           :   938.27203    
                    momentum (MeV/c)        :   19.373921    
                    energy, total (MeV)     :   938.47203    
                    energy, kinetic (MeV)   :  0.20000000    
                    beta = v/c              :  2.0644111178E-02
                    gamma                   :   1.000213158    
                    beta*gamma              :  2.0648511631E-02
                    G*gamma                 :   1.793229509    
                    electric rigidity (MeV) :  0.3999573775      =T[eV]*(gamma+1)/gamma, such that dev.=E*L/rigidity
  
 I, AMQ(1,I), AMQ(2,I)/QE, P/Pref, v/c, time, s :
  
     1   9.38272030E+02  1.00000000E+00  1.00000000E+00  2.06441112E-02  0.00000000E+00  0.00000000E+00

************************************************************************************************************************************
      3  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      2)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    1.0000   12.925    0.000    0.000    0.000   0.000000E+00     1
               Time of flight (mus) :   0.0000000     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01   0.000000E+00        1        1    1.000      (Y,T)         1
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)         1
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   2.000000E-01        1        1    1.000      (t,K)         1

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
      4  Keyword, label(s) :  DIPOLE                                                

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

                    Integration step :  1.0000E-02 cm   (i.e.,   2.0000E-04 rad  at mean radius RM =    50.00    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad

  A    1  1.0000    12.925     0.000     0.000     0.000            1.047    12.925    -0.000     0.000     0.000            1

 Cumulative length of optical axis =   0.523598776     m ;  Time  (for ref. rigidity & particle) =   8.460221E-08 s 

************************************************************************************************************************************
      5  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      4)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    1.0000   12.925   -0.394    0.000    0.000   1.353491E+01     1
               Time of flight (mus) :  2.18694843E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -3.936662E-04        1        1    1.000      (Y,T)         1
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)         1
   0.0000E+00   0.0000E+00   1.0000E+00   2.186948E-02   2.000000E-01        1        1    1.000      (t,K)         1

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
      6  Keyword, label(s) :  FIT                                                   

     FIT procedure launched. Method is 1

           variable #            1       IR =            1 ,   ok.
           variable #            1       IP =           30 ,   ok.
           variable #            2       IR =            4 ,   ok.
           variable #            2       IP =           60 ,   ok.
           constraint #            1       IR =            5 ,   ok.
           constraint #            1       I  =            1 ,   ok.
           constraint #            2       IR =            5 ,   ok.
           constraint #            2       I  =            1 ,   ok.

                    FIT  variables  and  constraints  in  good  order,  FIT  will proceed. 

                    Final FIT status will NOT be saved. For so, use the 'save [FileName]' command

 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        12.9       12.924889       65.0      0.177      OBJET      -                    -                   
   4   2    60  -1.000E-02   9.200E-03  9.20000000E-03  3.000E-02  1.333E-04  DIPOLE     -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     5    0.000000E+00    1.000E+00    1.825917E-11    5.86E-05 FAISCEAU   -                    -                    0
  3   1   3     5    0.000000E+00    1.000E+00    2.384760E-09    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   5.6874E-18

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


                OPEN FILE zgoubi.fai                                                                      
                FOR PRINTING COORDINATES 

               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        1 through the optical structure 

                     Total of          1 particles have been launched

     Multiple pass, 
          from element #     1 : OBJET     /label1=                    /label2=                    
                             to  REBELOTE  /label1=                    /label2=                    
     ending at pass #      21 at element #     8 : REBELOTE  /label1=                    /label2=                    


     Parameter #   35 in element #    1 will be modified at each pass. 
     list of requested parameter values :
                   1   1.00000000E+00
                   2   1.21086263E+00
                   3   1.42172526E+00
                   4   1.63258789E+00
                   5   1.84345052E+00
                   6   2.05431315E+00
                   7   2.26517578E+00
                   8   2.47603841E+00
                   9   2.68690104E+00
                  10   2.89776367E+00
                  11   3.10862630E+00
                  12   3.31948893E+00
                  13   3.53035156E+00
                  14   3.74121419E+00
                  15   3.95207682E+00
                  16   4.16293945E+00
                  17   4.37380208E+00
                  18   4.58466471E+00
                  19   4.79552734E+00
                  20   5.00638997E+00

 Pgm rebel. At pass #    1/  21.  In element #    1,  parameter # 35  changed to    1.00000000E+00   (was    1.00000000E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        12.9       12.924889       65.0      0.177      OBJET      -                    -                   
   4   2    60  -1.000E-02   9.200E-03  9.20000000E-03  3.000E-02  1.333E-04  DIPOLE     -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     5    0.000000E+00    1.000E+00    1.825917E-11    5.86E-05 FAISCEAU   -                    -                    0
  3   1   3     5    0.000000E+00    1.000E+00    2.384760E-09    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   5.6874E-18

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        2 through the optical structure 

                     Total of          2 particles have been launched

 Pgm rebel. At pass #    2/  21.  In element #    1,  parameter # 35  changed to    1.21086263E+00   (was    1.00000000E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        15.7       15.650265       65.0      3.324E-07  OBJET      -                    -                   
   4   2    60  -1.000E-02   1.130E-02  1.13013057E-02  3.000E-02  2.509E-10  DIPOLE     -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     5    0.000000E+00    1.000E+00    8.244742E-09    8.17E-05 FAISCEAU   -                    -                    0
  3   1   3     5    0.000000E+00    1.000E+00    9.123843E-07    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   8.3251E-13

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        3 through the optical structure 

                     Total of          3 particles have been launched

 Pgm rebel. At pass #    3/  21.  In element #    1,  parameter # 35  changed to    1.42172526E+00   (was    1.21086263E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        18.4       18.375641       65.0      3.324E-07  OBJET      -                    -                   
   4   2    60  -1.000E-02   1.363E-02  1.36278488E-02  3.000E-02  2.509E-10  DIPOLE     -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     5    0.000000E+00    1.000E+00    1.650725E-08    1.13E-04 FAISCEAU   -                    -                    0
  3   1   3     5    0.000000E+00    1.000E+00    1.555819E-06    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   2.4208E-12

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        4 through the optical structure 

                     Total of          4 particles have been launched

 Pgm rebel. At pass #    4/  21.  In element #    1,  parameter # 35  changed to    1.63258789E+00   (was    1.42172526E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        21.1       21.097924       65.0      8.211E-99  OBJET      -                    -                   
   4   2    60  -1.000E-02   1.594E-02  1.59429633E-02  3.000E-02  6.197-102  DIPOLE     -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     5    0.000000E+00    1.000E+00    1.546247E-03    8.46E-01 FAISCEAU   -                    -                    0
  3   1   3     5    0.000000E+00    1.000E+00    6.601042E-04    1.54E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   2.8266E-06

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        5 through the optical structure 

                     Total of          5 particles have been launched

 Pgm rebel. At pass #    5/  21.  In element #    1,  parameter # 35  changed to    1.84345052E+00   (was    1.63258789E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        23.8       23.826393       65.0      3.324E-07  OBJET      -                    -                   
   4   2    60  -1.000E-02   1.640E-02  1.64005163E-02  3.000E-02  2.509E-10  DIPOLE     -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     5    0.000000E+00    1.000E+00    5.845306E-08    1.89E-04 FAISCEAU   -                    -                    0
  3   1   3     5    0.000000E+00    1.000E+00    4.249287E-06    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   1.8060E-11

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        6 through the optical structure 

                     Total of          6 particles have been launched

 Pgm rebel. At pass #    6/  21.  In element #    1,  parameter # 35  changed to    2.05431315E+00   (was    1.84345052E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        26.5       26.549303       65.0      4.172-103  OBJET      -                    -                   
   4   2    60  -1.000E-02   1.845E-02  1.84504131E-02  3.000E-02  3.149-106  DIPOLE     -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     5    0.000000E+00    1.000E+00    1.232902E-03    8.49E-01 FAISCEAU   -                    -                    0
  3   1   3     5    0.000000E+00    1.000E+00    5.209476E-04    1.51E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   1.7914E-06

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        7 through the optical structure 

                     Total of          7 particles have been launched

 Pgm rebel. At pass #    7/  21.  In element #    1,  parameter # 35  changed to    2.26517578E+00   (was    2.05431315E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        29.3       29.277145       65.0      1.108E-07  OBJET      -                    -                   
   4   2    60  -1.000E-02   1.891E-02  1.89133371E-02  3.000E-02  8.363E-11  DIPOLE     -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     5    0.000000E+00    1.000E+00    8.209812E-08    2.86E-04 FAISCEAU   -                    -                    0
  3   1   3     5    0.000000E+00    1.000E+00    4.856888E-06    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   2.3596E-11

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        8 through the optical structure 

                     Total of          8 particles have been launched

 Pgm rebel. At pass #    8/  21.  In element #    1,  parameter # 35  changed to    2.47603841E+00   (was    2.26517578E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        32.0       32.002521       65.0      3.324E-07  OBJET      -                    -                   
   4   2    60  -1.000E-02   2.264E-02  2.26437034E-02  3.000E-02  2.509E-10  DIPOLE     -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     5    0.000000E+00    1.000E+00    9.036153E-08    3.41E-04 FAISCEAU   -                    -                    0
  3   1   3     5    0.000000E+00    1.000E+00    4.890409E-06    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   2.3924E-11

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        9 through the optical structure 

                     Total of          9 particles have been launched

 Pgm rebel. At pass #    9/  21.  In element #    1,  parameter # 35  changed to    2.68690104E+00   (was    2.47603841E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        34.7       34.727897       65.0      3.324E-07  OBJET      -                    -                   
   4   2    60  -1.000E-02   2.310E-02  2.30988170E-02  3.000E-02  2.509E-10  DIPOLE     -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     5    0.000000E+00    1.000E+00    9.862366E-08    4.02E-04 FAISCEAU   -                    -                    0
  3   1   3     5    0.000000E+00    1.000E+00    4.918690E-06    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   2.4203E-11

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       10 through the optical structure 

                     Total of         10 particles have been launched

 Pgm rebel. At pass #   10/  21.  In element #    1,  parameter # 35  changed to    2.89776367E+00   (was    2.68690104E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        37.5       37.453273       65.0      3.324E-07  OBJET      -                    -                   
   4   2    60  -1.000E-02   2.523E-02  2.52334283E-02  3.000E-02  2.509E-10  DIPOLE     -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     5    0.000000E+00    1.000E+00    1.068852E-07    4.67E-04 FAISCEAU   -                    -                    0
  3   1   3     5    0.000000E+00    1.000E+00    4.942869E-06    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   2.4443E-11

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       11 through the optical structure 

                     Total of         11 particles have been launched

 Pgm rebel. At pass #   11/  21.  In element #    1,  parameter # 35  changed to    3.10862630E+00   (was    2.89776367E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        40.2       40.178649       65.0      3.324E-07  OBJET      -                    -                   
   4   2    60  -1.000E-02   2.726E-02  2.72577878E-02  3.000E-02  2.509E-10  DIPOLE     -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     5    0.000000E+00    1.000E+00    1.151469E-07    5.38E-04 FAISCEAU   -                    -                    0
  3   1   3     5    0.000000E+00    1.000E+00    4.963761E-06    9.99E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   2.4652E-11

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       12 through the optical structure 

                     Total of         12 particles have been launched

 Pgm rebel. At pass #   12/  21.  In element #    1,  parameter # 35  changed to    3.31948893E+00   (was    3.10862630E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        42.9       42.904026       65.0      3.324E-07  OBJET      -                    -                   
   4   2    60  -1.000E-02   2.945E-02  2.94497188E-02  3.000E-02  2.509E-10  DIPOLE     -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     5    0.000000E+00    1.000E+00    1.234087E-07    6.13E-04 FAISCEAU   -                    -                    0
  3   1   3     5    0.000000E+00    1.000E+00    4.981998E-06    9.99E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   2.4836E-11

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       13 through the optical structure 

                     Total of         13 particles have been launched

 Pgm rebel. At pass #   13/  21.  In element #    1,  parameter # 35  changed to    3.53035156E+00   (was    3.31948893E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        45.6       45.629402       65.0      3.324E-07  OBJET      -                    -                   
   4   2    60  -1.000E-02   3.000E-02  2.99977059E-02  3.000E-02  2.509E-10  DIPOLE     -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     5    0.000000E+00    1.000E+00    1.316726E-07    6.94E-04 FAISCEAU   -                    -                    0
  3   1   3     5    0.000000E+00    1.000E+00    4.998036E-06    9.99E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   2.4998E-11

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       14 through the optical structure 

                     Total of         14 particles have been launched

 Pgm rebel. At pass #   14/  21.  In element #    1,  parameter # 35  changed to    3.74121419E+00   (was    3.53035156E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        48.4       48.354778       65.0      3.324E-07  OBJET      -                    -                   
   4   2    60  -1.000E-02   3.000E-02  2.99990245E-02  3.000E-02  2.509E-10  DIPOLE     -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     5    0.000000E+00    1.000E+00    1.399354E-07    7.79E-04 FAISCEAU   -                    -                    0
  3   1   3     5    0.000000E+00    1.000E+00    5.012269E-06    9.99E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   2.5142E-11

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       15 through the optical structure 

                     Total of         15 particles have been launched

 Pgm rebel. At pass #   15/  21.  In element #    1,  parameter # 35  changed to    3.95207682E+00   (was    3.74121419E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        51.1       51.079675       65.0      1.108E-07  OBJET      -                    -                   
   4   2    60  -1.000E-02   2.967E-02  2.96650616E-02  3.000E-02  8.363E-11  DIPOLE     -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     5    0.000000E+00    1.000E+00    2.391882E-04    1.00E+00 FAISCEAU   -                    -                    0
  3   1   3     5    0.000000E+00    1.000E+00    3.828126E-07    2.56E-06 FAISCEAU   -                    -                    0
 Fit reached penalty value   5.7211E-08

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       16 through the optical structure 

                     Total of         16 particles have been launched

 Pgm rebel. At pass #   16/  21.  In element #    1,  parameter # 35  changed to    4.16293945E+00   (was    3.95207682E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        53.8       53.805530       65.0      3.324E-07  OBJET      -                    -                   
   4   2    60  -1.000E-02   2.978E-02  2.97825053E-02  3.000E-02  2.509E-10  DIPOLE     -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     5    0.000000E+00    1.000E+00    1.564579E-07    9.64E-04 FAISCEAU   -                    -                    0
  3   1   3     5    0.000000E+00    1.000E+00    5.036452E-06    9.99E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   2.5390E-11

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       17 through the optical structure 

                     Total of         17 particles have been launched

 Pgm rebel. At pass #   17/  21.  In element #    1,  parameter # 35  changed to    4.37380208E+00   (was    4.16293945E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        56.5       56.530906       65.0      3.324E-07  OBJET      -                    -                   
   4   2    60  -1.000E-02   2.905E-02  2.90541151E-02  3.000E-02  2.509E-10  DIPOLE     -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     5    0.000000E+00    1.000E+00    1.647208E-07    1.06E-03 FAISCEAU   -                    -                    0
  3   1   3     5    0.000000E+00    1.000E+00    5.046775E-06    9.99E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   2.5497E-11

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       18 through the optical structure 

                     Total of         18 particles have been launched

 Pgm rebel. At pass #   18/  21.  In element #    1,  parameter # 35  changed to    4.58466471E+00   (was    4.37380208E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        59.3       59.255571       65.0      1.108E-07  OBJET      -                    -                   
   4   2    60  -1.000E-02   2.999E-02  2.99886841E-02  3.000E-02  8.363E-11  DIPOLE     -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     5    0.000000E+00    1.000E+00    3.555050E-04    9.99E-01 FAISCEAU   -                    -                    0
  3   1   3     5    0.000000E+00    1.000E+00    8.840993E-06    6.18E-04 FAISCEAU   -                    -                    0
 Fit reached penalty value   1.2646E-07

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       19 through the optical structure 

                     Total of         19 particles have been launched

 Pgm rebel. At pass #   19/  21.  In element #    1,  parameter # 35  changed to    4.79552734E+00   (was    4.58466471E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        62.0       61.981658       65.0      3.324E-07  OBJET      -                    -                   
   4   2    60  -1.000E-02   3.000E-02  2.99999669E-02  3.000E-02  2.509E-10  DIPOLE     -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     5    0.000000E+00    1.000E+00    1.502789E-08    1.28E-03 FAISCEAU   -                    -                    0
  3   1   3     5    0.000000E+00    1.000E+00    4.199344E-07    9.99E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   1.7657E-13

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       20 through the optical structure 

                     Total of         20 particles have been launched


      Next  pass  is  #    21 and  last  pass  through  the  optical  structure


 Pgm rebel. At pass #   20/  21.  In element #    1,  parameter # 35  changed to    5.00638997E+00   (was    4.79552734E+00)

     WRITE statements to zgoubi.res are re-established from now on.

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 

                          MAGNETIC  RIGIDITY =         64.624 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       1 POINTS 



************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      2)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

m  1   5.0064    61.982     0.000     0.000     0.000       0.0000    5.0064   61.982    0.000    0.000    0.000   0.000000E+00     1
               Time of flight (mus) :   0.0000000     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   6.198166E-01   0.000000E+00        1        1    1.000      (Y,T)        21
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        21
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   5.000000E+00        1        1    1.000      (t,K)        21

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
      4  Keyword, label(s) :  DIPOLE                                                

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

  A    1  5.0064    61.982     0.000     0.000     0.000            1.047    63.301     0.036     0.000     0.000            1

 Cumulative length of optical axis =   0.523598776     m ;  Time  (for ref. rigidity & particle) =   1.692044E-07 s 

************************************************************************************************************************************
      5  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      4)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

m  1   5.0064    61.982     0.000     0.000     0.000       0.0000    5.0064   63.301   36.484    0.000    0.000   6.540028E+01     1
               Time of flight (mus) :  2.12155187E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   6.330129E-01   3.648396E-02        1        1    1.000      (Y,T)        21
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        21
   0.0000E+00   0.0000E+00   1.0000E+00   2.121552E-02   5.000000E+00        1        1    1.000      (t,K)        21

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
      6  Keyword, label(s) :  FIT                                                   

     FIT procedure launched. Method is 1


                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 

                          MAGNETIC  RIGIDITY =         64.624 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       1 POINTS 



************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      2)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

m  1   5.0064    61.982     0.000     0.000     0.000       0.0000    5.0064   61.982    0.000    0.000    0.000   0.000000E+00     1
               Time of flight (mus) :   0.0000000     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   6.198166E-01   0.000000E+00        1        1    1.000      (Y,T)        21
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        21
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   5.000000E+00        1        1    1.000      (t,K)        21

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
      4  Keyword, label(s) :  DIPOLE                                                

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

  A    1  5.0064    61.982     0.000     0.000     0.000            1.047    63.301     0.036     0.000     0.000            1

 Cumulative length of optical axis =   0.523598776     m ;  Time  (for ref. rigidity & particle) =   8.460221E-08 s 

************************************************************************************************************************************
      5  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      4)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

m  1   5.0064    61.982     0.000     0.000     0.000       0.0000    5.0064   63.301   36.484    0.000    0.000   6.540028E+01     1
               Time of flight (mus) :  2.12155187E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   6.330129E-01   3.648396E-02        1        1    1.000      (Y,T)        21
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        21
   0.0000E+00   0.0000E+00   1.0000E+00   2.121552E-02   5.000000E+00        1        1    1.000      (t,K)        21

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
      6  Keyword, label(s) :  FIT                                                   

     FIT procedure launched. Method is 1

           variable #            1       IR =            1 ,   ok.
           variable #            1       IP =           30 ,   ok.
           variable #            2       IR =            4 ,   ok.
           variable #            2       IP =           60 ,   ok.
           constraint #            1       IR =            5 ,   ok.
           constraint #            1       I  =            1 ,   ok.
           constraint #            2       IR =            5 ,   ok.
           constraint #            2       I  =            1 ,   ok.

                    FIT  variables  and  constraints  in  good  order,  FIT  will proceed. 

                    Final FIT status will NOT be saved. For so, use the 'save [FileName]' command

 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        64.7       64.707034       65.0      3.324E-07  OBJET      -                    -                   
   4   2    60  -1.000E-02   3.000E-02  2.99967146E-02  3.000E-02  2.509E-10  DIPOLE     -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     5    0.000000E+00    1.000E+00    2.329175E-08    1.39E-03 FAISCEAU   -                    -                    0
  3   1   3     5    0.000000E+00    1.000E+00    6.234029E-07    9.99E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   3.8917E-13

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                         ****  End  of  'REBELOTE'  procedure  ****

      There  has  been         21  passes  through  the  optical  structure 

                     Total of         21 particles have been launched

************************************************************************************************************************************
      9  Keyword, label(s) :  END                                                   


************************************************************************************************************************************

          Pgm zgoubi : Execution ended upon key  END       

************************************************************************************************************************************
   
            ZGOUBI RUN COMPLETED. 

  Zgoubi, author's dvlpmnt version.
  Job  started  on  19-05-0017,  at  07:08:56 
  JOB  ENDED  ON    19-05-0017,  AT  07:09:50 

   CPU time, total :     53.6993560000000     
